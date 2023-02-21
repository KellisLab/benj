
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
