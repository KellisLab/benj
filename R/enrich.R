

#' Binomial CI estimation wilson score interval
#' as taken from Carles Boix
#'
#' @param ns number successes
#' @param n total number
#' @param z standard deviations
#' @param max whether lower or upper tail
#' @export
calc_binomial_CI <- function(ns, n, z, max=TRUE) {
    nf = n - ns
    zsq = z^2
    p = (ns + zsq/2) / (n + zsq)
    pm = (z / (n + zsq)) * sqrt((ns * nf) / n + zsq / 4)
    return(p + (2 * max - 1) * pm)
}

#' Summarize JASPAR motifs as taken from ArchR
#'
#' @param motifs JASPAR motifs
#' @return list of motifs and a summary dataframe
.summarizeJASPARMotifs <- function(motifs) {
    mnames = sapply(seq_along(motifs), function(x) {
        namex = gsub("::", "_", motifs[[x]]@name)
        namex = make.names(namex)
        return(paste0(namex, "_", x))
    })
    mdf = lapply(seq_along(motifs), function(x) {
        data.frame(
            row.names = mnames[x],
            name = motifs[[x]]@name[[1]],
            ID = motifs[[x]]@ID,
            strand = motifs[[x]]@strand,
            symbol = ifelse(!is.null(motifs[[x]]@tags$symbol[1]), motifs[[x]]@tags$symbol[1], NA) ,
            family = ifelse(!is.null(motifs[[x]]@tags$family[1]), motifs[[x]]@tags$family[1], NA),
            alias = ifelse(!is.null(motifs[[x]]@tags$alias[1]), motifs[[x]]@tags$alias[1], NA),
            stringsAsFactors = FALSE) })
    names(motifs) = mnames
    return(list(motifs = motifs, motifSummary = Reduce(rbind, mdf)))
}

#' Motif overlap of genomic ranges
#'
#' @param peakSet GenomicRanges
#' @param jaspar JASPAR2022::JASPAR2022 or similar object
#' @param species Species to select from JASPAR
#' @param collection Collection from JASPAR
#' @param width Width for motifmatchr. Default=7
#' @param BSgenome BSgenome object to match peaks against. NULL means hg38
#' @param cutoff P-value cutoff for motifmatchr. Default=5e-05
#' @return RangedSummarizedExperiment of each motif's enrichment
#' @export
motif_overlap <- function(peakSet, jaspar, species="Homo sapiens", collection="CORE",
                          BSgenome=NULL,
                          width=7, cutoff=5e-05, counts=FALSE) {
    if (is.null(BSgenome)) {
        warning("BSgenome unspecified. Using BSgenome.Hsapiens.UCSC.hg38")
        BSgenome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    }
    motifs = TFBSTools::getMatrixSet(jaspar, list(species=species, collection=collection))
    obj = .summarizeJASPARMotifs(motifs)
    ### get matches
    mpos = motifmatchr::matchMotifs(pwms=obj$motifs, subject=peakSet, genome=BSgenome, out="positions", p.cutoff=cutoff, w=width)
    ### filter motifs with no matches
    novp = sapply(mpos, length)
    no_match = names(which(novp == 0))
    mpos = mpos[novp > 0]
    obj$motifSummary = obj$motifSummary[names(mpos),]
    obj$motifs = obj$motifs[names(mpos)]
    ### find overlaps with peaks
    apos = unlist(mpos)
    ovp = GenomicRanges::findOverlaps(peakSet, apos, ignore.strand=TRUE)
    motifMat = Matrix::sparseMatrix(
                           i = S4Vectors::from(ovp),
                           j = match(names(apos), names(mpos))[S4Vectors::to(ovp)],
                           x = rep(ifelse(counts, 1, TRUE), length(ovp)),
                           dimnames = list(peaks=names(peakSet), motifs=names(mpos)),
                           dims = c(length(peakSet), length(mpos)))
    return(SummarizedExperiment::SummarizedExperiment(assays=list(matches=motifMat), rowRanges=peakSet, colData=obj$motifSummary[colnames(motifMat),]))
}


#' Enrich against a custom range
#'
#' @param peakSet The peaks to enrich
#' @param gr The custom enrichment ranges to enrich against
#' @param gr.col Column in gr to split upon
#' @param counts logical whether to accumulate counts or not
#' @return RangedSummarizedExperiment
#' @export
enrich_overlap_custom <- function(peakSet, gr, gr.col=NULL, counts=FALSE) {
    ovp = GenomicRanges::findOverlaps(peakSet, gr)
    if (is.null(gr.col)) {
        cols = as.factor(rep("CustomEnrichment", length(ovp)))
    } else {
        cols = S4Vectors::mcols(gr)[S4Vectors::to(ovp), gr.col]
        cols = as.factor(cols)
    }
    motifMat = Matrix::sparseMatrix(
                           i = S4Vectors::from(ovp),
                           j = as.integer(cols),
                           x=rep(ifelse(counts, 1, TRUE), length(ovp)),
                           dims=c(length(peakSet),
                                  length(levels(cols))),
                           dimnames=list(names(peakSet), levels(cols)))
    return(SummarizedExperiment::SummarizedExperiment(assays=list(matches=motifMat), rowRanges=peakSet))
}
#' Enrich SummarizedExperiment by groupby column
#'
#' @param se SummarizedExperiment to enrich
#' @param groupby Column in rowData(se) to split upon
#' @param feature_name Feature column in returned dataframe
#' @param group_name Group name in returned dataframe
#' @param cols.colData Columns to add from colData(se)
#' @return dataframe with group, log2FC, and mlog10p columns along with $feature_name
#' @export
enrich_se <- function(se, groupby, feature_name="features", group_name="group", cols.colData="name") {
    P = t(make_pseudobulk(SummarizedExperiment::rowData(se)[[groupby]]))
    M = as.matrix(P %*% SummarizedExperiment::assays(se)[[1]])
    mdf = reshape2::melt(M, value.name="pfg", varnames=c(group_name, feature_name))
    cd = SummarizedExperiment::colData(se)
    for (col in cols.colData[cols.colData %in% names(cd)]) {
        mdf[[col]] = cd[mdf[[feature_name]], col]
    }
    mdf$cfg = rowSums(P)[mdf[[group_name]]]
    mdf$pbg = colSums(M)[mdf[[feature_name]]]
    mdf$cbg = sum(P) - mdf$cfg
    ### Compute log2FC with wilson score interval
    mdf$log2FC = with(mdf, log2(pfg / pbg) - log2(cfg / cbg))
    pind = which(mdf$log2FC >= 0)
    nind = which(mdf$log2FC < 0)
    mdf$pfrac = 1
    mdf$cfrac = 1
    mdf$pfrac[pind] = with(mdf, sapply(pind, function(i) {
        calc_binomial_CI(pfg[i], pbg[i], 1.5, max=F)
    }))
    mdf$cfrac[pind] = with(mdf, sapply(pind, function(i) {
        calc_binomial_CI(cfg[i], cbg[i], 1.5, max=T)
    }))
    mdf$pfrac[nind] = with(mdf, sapply(nind, function(i) {
        calc_binomial_CI(pfg[i], pbg[i], 1.5, max=T)
    }))
    mdf$cfrac[nind] = with(mdf, sapply(nind, function(i) {
        calc_binomial_CI(cfg[i], cbg[i], 1.5, max=F)
    }))
    mdf$pfrac[is.na(mdf$pfrac)] = 0
    mdf$cfrac[is.na(mdf$cfrac)] = 0
    mdf$log2FC = log2(mdf$pfrac / mdf$cfrac)
    ### Calculate p-value
    enr_mlog10p = apply(mdf[,c("pfg", "pbg", "cfg", "cbg")], 1, function(y) {
        -phyper(q=y[1] - 1, m=y[3], n=y[4] - y[3], k=y[2], lower.tail=FALSE, log.p=TRUE) / log(10)
    })
    dep_mlog10p = apply(mdf[,c("pfg", "pbg", "cfg", "cbg")], 1, function(y) {
        -phyper(q=y[1], m=y[3], n=y[4] - y[3], k=y[2], lower.tail=TRUE, log.p=TRUE) / log(10)
    })
    mdf$mlog10p = apply(cbind(enr_mlog10p, dep_mlog10p), 1, max)
    return(mdf[order(mdf$mlog10p, decreasing=TRUE), ])
}
