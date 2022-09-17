
.summarizeJASPARMotifs <- function(motifs) {
    mnames = sapply(seq_along(motifs), function(x) {
        namex = gsub("::", "_", motifs[[x]]@name)
        namex = make.names(namex)
        return(paste0(namex, "_", x))
    })
    mdf = lapply(seq_along(motifs), function(x){
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
