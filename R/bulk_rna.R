

#' Aggregate bulk RNA-seq into SummarizedExperiment
#' @param dname Directory name containing bulk RNA-seq data e.g. 230210Kel
#' @param use_gene_names Use gene names instead of Ensembl ids for genes
#' @return SummarizedExperiment object
#' @export
bulk_aggregate <- function(dname, use_gene_names=TRUE) {
    rname = list.files(dname, pattern="^RNAseqQC-[A-Z0-9]+$", full.names=TRUE)
    suffix = gsub("^RNAseqQC-", "", basename(rname))
    pat = paste0("-", suffix, "_geneexp.txt$")
    FL = list.files(rname, full.names=TRUE, pattern=pat)
    df = do.call(rbind, lapply(FL, function(fname) {
        df = read.table(fname, sep="\t", col.names=c("gene_id", "count", "gene_name"))
        df$library_id = gsub(pat, "", basename(fname))
        df = df[!is.na(df$gene_name),]
        return(df)
    }))
    M = pivot(df, "gene_id", "library_id", "count")
    rd = data.frame(row.names=rownames(M), gene_id=rownames(M))
    gf = df[!duplicated(df$gene_id), c("gene_id", "gene_name")]
    rownames(gf) = gf$gene_id
    rd$gene_name = var_names_make_unique(gf[rownames(rd),"gene_name"])
    if (use_gene_names) {
        rownames(M) = rd[rownames(M),]$gene_name
        rownames(rd) = rd$gene_name
    }
    cd = data.frame(row.names=colnames(M), batch=rep(basename(dname), ncol(M)))
    return(SummarizedExperiment::SummarizedExperiment(list(counts=M), rowData=rd, colData=cd))
}

#' Calculate TMM assay in SummarizedExperiment
#'
#' @param se SummarizedExperiment
#' @param method Method used for normalization, e.g. TMM or RLE
#' @param log Whether to use log-CPM or CPM after normalization
#' @export
se_tmm <- function(se, method="TMM", log=FALSE) {
    dgel = edgeR::calcNormFactors(se, method)
    tmm = edgeR::cpm(dgel, log=log)
    SummarizedExperiment::assays(se)[[method]] = tmm
    return(se)
}

#' Add genomicranges of gene to
#' @export
se_gene_ranges <- function(se, gtf, by="gene_id") {
    rd = SummarizedExperiment::rowData(se)
    gf = as.data.frame(rtracklayer::readGFF(gtf))
    gf = gf[gf$type == "gene",]
    rownames(gf) = var_names_make_unique(gf[[by]])
    gr = with(gf, GenomicRanges::GRanges(seqid, ranges=IRanges::IRanges(start, end)))
    names(gr) = rownames(gf)
    gr = gr[rd[[by]]]
    for (name in colnames(rd)) {
        S4Vectors::mcols(gr)[[name]] = rd[[name]]
    }
    names(gr) = rownames(rd)
    SummarizedExperiment::rowRanges(se) = gr
    return(se)
}
#' Load STAR ReadsPerGene into SummarizedExperiment
#'
#' @param star.prefix Vector of directories with STAR ReadsPerGene.out.tab inside. Can be named
#' @param index STAR index used. Will use geneInfo.tab to load in gene names
#' @param strand By default, use unstranded. Can use either first or second strand as well.
#' @param featureCounts use featureCounts output file if exists. Will use instead of STAR ReadsPerGene.out.tab
#' @export
read_star <- function(star.prefix, index="/net/bmc-lab5/data/kellis/group/Benjamin/ref/STAR_gencode43/", featureCounts="featureCounts", strand=NULL, ...) {
    gf = NULL
    if (file.exists(paste0(index, "/geneInfo.tab"))) {
        gf = read.table(paste0(index, "/geneInfo.tab"), skip=1,sep="\t")
        if (ncol(gf) == 3) {
            colnames(gf) = c("gene_id", "gene_name", "gene_type")
        } else {
            colnames(gf) = "gene_id"
        }
        rownames(gf) = gf$gene_id
    }
    xf = NULL
    for (i in seq_along(star.prefix)) {
        if (all(file.exists(paste0(star.prefix, featureCounts)))) {
            print(paste0("Using ", star.prefix[i], "featureCounts"))
            ### Only use featureCounts if exists for all in prefix.
            df = read.table(paste0(star.prefix[i], featureCounts), sep="\t", skip="#", header=TRUE)
            if (ncol(df) == 7) {
                df = df[c("Geneid", colnames(df)[ncol(df)])]
                colnames(df) = c("gene_id", "unstranded")
            } else {
                stop(paste0("Incorrect number of columns in ", star.prefix[i], featureCounts, ", ", ncol(df)))
            }
        } else {
            print(paste0("Using ", star.prefix[i], "ReadsPerGene.out.tab"))
            df = read.table(paste0(star.prefix[i], "ReadsPerGene.out.tab"), sep="\t", col.names=c("gene_id", "unstranded", "first", "second"))
        }
        df$sample = ifelse(is.null(names(star.prefix)), basename(star.prefix[i]), names(star.prefix)[i])
        if (is.null(xf)) {
            xf = df
        } else {
            xf = rbind(xf, df)
        }
    }
    if (is.null(strand)) {
        xf$count = xf$unstranded
    } else {
        xf$count = xf[[strand]]
    }
    if (is.null(gf)) {
        M = pivot(xf, "gene_id", "sample", "count")
        se = SummarizedExperiment::SummarizedExperiment(list(counts=M), ...)
    } else {
        M = pivot(xf[xf$gene_id %in% gf$gene_id,], "gene_id", "sample", "count")
        se = SummarizedExperiment::SummarizedExperiment(list(counts=M), rowData=gf[rownames(M),], ...)
        if ("gene_name" %in% colnames(gf)) {
            rownames(se) = var_names_make_unique(gf[rownames(se),]$gene_name)
        }
    }
    return(se)
}

#' List BMC filenames in batch directory
#' @param df Dataframe with rownames for sample names, and colnames for colData columns
#' @param seq.dir Outer directory where sample names are located
#' @param extra Boolean indicating whether BMC style suffixes should be added
#' @export
bulk_aggregate_star <- function(df, seq.dir, extra=TRUE) {
    sample.dir.list = sapply(rownames(df), function(rn) {
        sample.dir = list.files(paste0(seq.dir, "/"),
                                pattern=paste0("^", rn,
                                               ifelse(extra, "-[A-Z0-9]+", ""),
                                               "$"),
                                full.names=TRUE)
        stopifnot(length(sample.dir) == 1)
        return(sample.dir)
    })
    stopifnot(length(sample.dir.list) == nrow(df))
    se = read_star(setNames(paste0(sample.dir.list, "/"), rownames(df)))
    for (i in seq_along(colnames(df))) {
        SummarizedExperiment::colData(se)[[ colnames(df)[i] ]] = df[[ i ]]
    }
    return(se)
}

#' @export
plotPCA <- function(se, title, method="TMM", top=500, cpm.frac=0.25, cpm.cutoff=100, gene.selection="common", correct=FALSE, use_label=TRUE, force=0.2, text.size=4, max.overlaps=10) {
    require(ggplot2)
    dgel = edgeR::calcNormFactors(se, method)
    to_keep = rowMeans(edgeR::cpm(dgel) > cpm.cutoff) >= cpm.frac
    dgel = dgel[to_keep,,keep.lib.sizes=FALSE]
    if (correct) {
        dgel = limma::removeBatchEffect(dgel, dgel$batch)
    }
    pl = limma::plotMDS(dgel, top=top, gene.selection=gene.selection, plot=FALSE)
    pf = as.data.frame(SummarizedExperiment::colData(se))
    pf$x = pl$x
    pf$y = pl$y
    if (length(unique(pf$batch))==1) {
        g = ggplot(pf, aes(x=x, y=y, label=title))
    } else {
        g = ggplot(pf, aes(x=x, y=y, color=batch, label=title)) + scale_color_brewer(palette="Set2")
    }
    g = g + geom_point()
    if (use_label) {
        g = g + ggrepel::geom_text_repel(force=force, size=text.size, max.overlaps=max.overlaps)
    }
    g = g + xlab(paste0("PC1: ", round(100 * pl$var.explained[1]), "% of variance"))
    g = g + ylab(paste0("PC2: ", round(100 * pl$var.explained[2]), "% of variance"))
    g = g + ggpubr::theme_pubr() + ggtitle(title) + theme(plot.title=element_text(size=10), legend.text=element_text(size=8))
    return(g)
}
