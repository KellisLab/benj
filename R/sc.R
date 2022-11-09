
#' Pseudobulk a SummarizedExperiment
#'
#' Determines constant columns across dataframe
#' and summarizes counts
#'
#' @param se SummarizedExperiment
#' @param on Column in colData() to split upon
#' @export
se_make_pseudobulk <- function(se, on, missing.levels=FALSE, unlevel=TRUE) {
    cd = SummarizedExperiment::colData(se)
    rd = SummarizedExperiment::rowData(se)
    A = SummarizedExperiment::assays(se)
    to_P = Matrix::Matrix(make_pseudobulk(cd[[on]], unlevel=unlevel))
    P = lapply(A, function(mat) {
        mat %*% to_P
    })
    cd[[on]] = factor(as.character(cd[[on]]), levels=colnames(to_P))
    pcd = data.frame(row.names=levels(cd[[on]]))
    for (cn in colnames(cd)) {
### for each "on", ensure all equal
        u = lapply(split(cd[[cn]], cd[[on]]), unique)
        if (missing.levels) {
            good = all(sapply(u, function(x) { length(x) <= 1 }))
        } else {
            good = all(sapply(u, function(x) { length(x) == 1 }))
        }
        if (good) {
            pcd[[cn]] = sapply(u, identity)[rownames(pcd)]
            if (is.factor(cd[[cn]])) {
                pcd[[cn]] = as.factor(as.character(pcd[[cn]]))
            }
        }
    }
    return(SummarizedExperiment::SummarizedExperiment(assays=P,
                                                      rowData=rd,
                                                      colData=pcd))
}

#' Calculate qc metrics similar to Scanpy
#'
#' @param se SummarizedExperiment
#' @param assay Assay to calculate upon
#' @param log1p logical whether to use log statistics as well
#' @export
calculate_qc_metrics <- function(se, assay=NULL, log1p=TRUE) {
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(se)[[1]]
    } else {
        X = SummarizedExperiment::assays(se)[[assay]]
    }
    rs = Matrix::rowSums(X)
    SummarizedExperiment::rowData(se)$total_counts = rs
    SummarizedExperiment::rowData(se)$n_cells_by_counts = Matrix::rowSums(X > 0)
    cs = Matrix::colSums(X)
    csz = Matrix::colSums(X > 0)
    SummarizedExperiment::colData(se)$total_counts = cs
    SummarizedExperiment::colData(se)$n_genes_by_counts = csz
    SummarizedExperiment::colData(se)$n_genes = csz
    if (log1p) {
        SummarizedExperiment::rowData(se)$log1p_total_counts = base::log1p(rs)
        SummarizedExperiment::colData(se)$log1p_total_counts = base::log1p(cs)
        SummarizedExperiment::colData(se)$log1p_n_genes_by_counts = base::log1p(csz)
    }
    return(se)
}

#' Save SummarizedExperiment to CellRanger style H5
#' for transfer to ScanPy
#'
#' Format: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices
#'
#' @param se SummarizedExperiment
#' @param h5 HDF5 output filename
#' @param assay Assay to export
#' @param chunk Chunk size for HDF5
#' @param genome Genome to set
#' @param compression integer GZIP compression level used
#' @export
se_to_cellranger_h5 <- function(se, h5, assay=NULL, chunk=10000, genome="GRCh38", compression=6) {
    rhdf5::h5createFile(h5)
    rhdf5::h5createGroup(h5, "matrix")
    rhdf5::h5createGroup(h5, "matrix/features")
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(se)[[1]]
    } else {
        X = SummarizedExperiment::assays(se)[[assay]]
    }
    if (!all(c("x", "i", "p") %in% slotNames(X))) {
        X = Matrix::Matrix(X, nrow=nrow(X), ncol=ncol(X), byrow=FALSE, sparse=TRUE) ### CSC
    }
    ml = list("matrix/data"=slot(X, "x"),
              "matrix/indices"=slot(X, "i"),
              "matrix/indptr"=slot(X, "p"),
              "matrix/barcodes"=colnames(X),
              "matrix/features/name"=var_names_make_unique(rownames(X)),
              "matrix/features/id"=rownames(X),
              "matrix/features/feature_type"=ifelse(strtrim(rownames(X), 3) == "chr",
                                                    "Peaks",
                                                    "Gene Expression"),
              "matrix/features/genome"=rep(genome, nrow(X)),
              "matrix/shape"=t(dim(X)))
    for (dname in names(ml)) {
        rhdf5::h5createDataset(h5, dname,
                               dims=length(ml[[dname]]),
                               storage.mode=storage.mode(ml[[dname]]),
                               chunk=min(length(ml[[dname]]), chunk),
                               level=compression)
        rhdf5::h5write(ml[[dname]], file=h5, name=dname)
    }
}

#' Normalize data similar to ScanPy
#'
#' @param se SummarizedExperiment object
#' @param target_sum Per UMI sum to target
#' @param assay Assay to normalize
#' @param normalized New assay name
#' @export
se_normalize_total <- function(se, target_sum=10000, assay=NULL, normalized="normalized") {
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(se)[[1]]
    } else {
        X = SummarizedExperiment::assays(se)[[assay]]
    }
    SummarizedExperiment::assays(se)[[normalized]] = X %*% Matrix::Diagonal(x=target_sum / Matrix::colSums(X))
    return(se)
}

se_log1p <- function(se, assay="normalized") {
    X = SummarizedExperiment::assays(se)[[assay]]
    SummarizedExperiment::assays(se)[[paste0("log1p_", assay)]] = log1p(X)
    return(se)
}
