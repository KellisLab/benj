
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

#'
#' Dump genomic ranges to BED6+
#'
#' gr Ranges
#' file File or connection for output
#' name Name of column to put as 4th column in BED6+
#' score Score column to put as 5th column in BED6+
#' extra Extra columns
#' default.score Number from 0-1000 inclusive to put as default for score if none exists
#' @export
bed.dump <- function(gr, file="", name="name", score="score", extra=NULL, default.score=0) {
    df = as.data.frame(gr)
    if ((!is.null(name)) & (name %in% colnames(df))) {
        df$name = df[[name]]
    } else if (!is.null(names(gr))) {
        df$name = names(gr)
    } else {
        df$name = with(df, paste0(seqnames, ":", start, "-", end))
    }
    extra = intersect(extra, colnames(df))
    if ((!is.null(score)) & (score %in% colnames(df))) {
        df$score = df[[score]]
    } else {
        df$score = default.score
    }
    df = df[c("seqnames", "start", "end", "name", "score", "strand", extra)]
    write.table(df, file=file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
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
        gff = as.data.frame(rtracklayer::readGFF(gff))
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

#' Run GREAT via rGREAT
#'
#' @param gr Ranges to send
#' @param spt Column(s) in ranges to split upon. Optional
#' @param bg logical set to use rest all of ranges as a background when split
#' @param species Species to send to GREAT
#' @return A dataframe with aggregated GREAT output and correct columns set if split
#' @export
run.great <- function(gr, spt=NULL, bg=TRUE, species="hg38") {
    .rgt <- function(x, y=NULL) {
        if (!is.null(y)) {
            job = rGREAT::submitGreatJob(x, y, species=species)
        } else {
            job = rGREAT::submitGreatJob(x, species=species)
        }
        tbl = rGREAT::getEnrichmentTables(job, download_by="tsv")
        if (is.list(tbl)) {
            tbl = Reduce(rbind, tbl)
        }
        return(tbl)
    }
    grcol = colnames(S4Vectors::mcols(gr))
    if (!is.null(spt) & all(spt %in% grcol)) {
        res = dplyr::bind_rows(lapply(split(gr, S4Vectors::mcols(gr)[[spt]]), function(lr) {
            if (bg) {
                tbl = .rgt(lr, gr)
            } else {
                tbl = .rgt(lr)
            }
            print(str(tbl))
            if (is.null(tbl)) {
                tbl = data.frame()
            }
            if (nrow(tbl) > 0) {
                ### Add cols
                for (col in spt) {
                    tbl[[col]] = rep(unique(S4Vectors::mcols(lr)[col]), nrow(tbl))
                }
            }
            return(tbl)
        }))
    } else {
        res = .rgt(gr)
    }
    return(res)
}
