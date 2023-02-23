
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

#' Load STAR ReadsPerGene into SummarizedExperiment
#'
#' @param star.prefix Vector of directories with STAR ReadsPerGene.out.tab inside. Can be named
#' @param index STAR index used. Will use geneInfo.tab to load in gene names
#' @param strand By default, use unstranded. Can use either first or second strand as well.
#' @export
read_star <- function(star.prefix, index="/net/bmc-lab5/data/kellis/group/Benjamin/ref/STAR_gencode43/", strand=NULL, ...) {
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
        df = read.table(paste0(star.prefix[i], "ReadsPerGene.out.tab"), sep="\t", col.names=c("gene_id", "unstranded", "first", "second"))
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
