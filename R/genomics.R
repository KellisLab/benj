
#' Parse string range to seqnames, start, end
#'
#' @param string Input string in format seqnames:start-end
#' @param ranges Boolean whether to return as GRanges
#' @param name Column to keep original name as
#' @param sep Suffix to name duplicates. "-" works well as it makes idempotent
#' @return data.frame or GRanges
#' @export
parse.range <- function(string, ranges=TRUE, name=NULL, sep="-") {
    S = strsplit(string, "[:-]")
    df = data.frame(seqnames = sapply(S, "[[", 1),
                    start = as.integer(sapply(S, "[[", 2)),
                    end = as.integer(sapply(S, "[[", 3)),
                    row.names = make.unique(string, sep))
    if (ranges) {
        gr = with(df, GenomicRanges::GRanges(seqnames=seqnames,
                                             ranges=IRanges::IRanges(start, end)))
        names(gr) = rownames(df)
        if (!is.null(name)) {
            S4Vectors::mcols(gr)[[name]] = string
        }
        return(gr)
    } else {
        if (!is.null(name)) {
            df[[name]] = string
        }
        return(df)
    }
}

#' Slurp BED file into a GenomicRanges object
#'
#' @param fname Filename
#' @param row.names logical whether to set rownames from 4th column
#' @return GenomicRanges object
#' @export
bed.slurp <- function(fname, row.names=FALSE, name.col="name", score.col="score") {
    df = read.table(fname, sep="\t")
    if (ncol(df) > 6) {
        warning("BED7+ formats not supported, ignoring extra columns")
    }
    gr = with(df, GenomicRanges::GRanges(seqnames=V1, ranges=IRanges::IRanges(V2, V3),
                                         ifelse(ncol(df) >= 6, V6, "*")))
    names(gr) = with(df, paste0(V1, ":", V2, "-", V3))
    if (ncol(df) >= 4) {
        if (row.names) {
            names(gr) = make.unique(df$V4, "-")
        } else {
            S4Vectors::mcols(gr)[[name.col]] = df$V4
        }
    }
    if (ncol(df) >= 5) {
        S4Vectors::mcols(gr)[[score.col]] = df$V5
    }
    return(gr)
}


#' Annotate peaks based on GFF
#'
#' gff can be a dataframe for filtering, e.g. gff[gff$gene_type == "protein_coding",]
#'
#' @param gr Ranges to annotate
#' @param gff GFF filename or data.frame with seq info
#' @param upstream Base pairs upstream for TSS centric
#' @param downstream Base pairs downstream for TSS centric
#' @return Updated ranges
#' @export
peak_annotation <- function(gr, gff, upstream=2000, downstream=200) {
    if (!is.data.frame(gff)) {
        gff = rtracklayer::readGFF(gff)
    }
    exons = with(gff[gff$type == "exon", ],
                 GenomicRanges::GRanges(seqnames=seqid,
                                        ranges=IRanges::IRanges(start, end),
                                        strand=strand,
                                        gene_name=var_names_make_unique(gene_name)))
    genes = with(gff[gff$type == "gene", ],
                 GenomicRanges::GRanges(seqnames=seqid,
                                        ranges=IRanges::IRanges(start, end),
                                        strand=strand,
                                        gene_name=var_names_make_unique(gene_name)))
    prom = GenomicRanges::promoters(genes, 0, 0)
    gr$nearestGene = genes[GenomicRanges::nearest(gr, prom)]$gene_name
    gr$distToTSS = S4Vectors::mcols(GenomicRanges::distanceToNearest(gr, prom))$distance
    prom = GenomicRanges::promoters(genes, upstream=upstream, downstream=downstream)
    op = S4Vectors::from(GenomicRanges::findOverlaps(gr, prom))
    og = S4Vectors::from(GenomicRanges::findOverlaps(gr, genes))
    oe = S4Vectors::from(GenomicRanges::findOverlaps(gr, exons))
    peakType = rep("Distal", length(gr))
    peakType[intersect(og, oe)] = "Exonic"
    peakType[setdiff(og, oe)] = "Intronic"
    peakType[unique(op)] = "TSS-proximal"
    gr$peakType = peakType
    return(gr)
}

run.great <- function(gr, split=NULL, bg=TRUE) {
    grcol = colnames(S4Vectors::mcols(gr))
    if (!is.null(split) & all(split %in% grcol)) {

    }
}
