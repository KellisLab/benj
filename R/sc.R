
#' Pseudobulk a SummarizedExperiment
#'
#' Determines constant columns across dataframe
#' and summarizes counts
#'
#' @param se SummarizedExperiment
#' @param on Column in colData() to split upon
#' @export
se_make_pseudobulk <- function(se, on, missing.levels=FALSE, unlevel=TRUE, ncell_col="ncell") {
    cd = SummarizedExperiment::colData(se)
    rd = SummarizedExperiment::rowData(se)
    A = SummarizedExperiment::assays(se)
    to_P = Matrix::Matrix(make_pseudobulk(cd[[on]], unlevel=unlevel))
    P = lapply(A, function(mat) {
        mat %*% to_P
    })
    if ("counts" %in% names(A)) {
        Percent <- as.matrix((A$counts > 0) %*% to_P) %*% diag(100/Matrix::colSums(to_P))
        dimnames(Percent) <- dimnames(P$counts)
        P$Percent <- Percent
    }
    cd[[on]] = factor(as.character(cd[[on]]), levels=colnames(to_P))
    pcd = data.frame(row.names=levels(cd[[on]]))
    pcd[[ncell_col]] <- Matrix::colSums(to_P)
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
    if (attr(se, "class") == "SingleCellExperiment") {
        to_AP = make_average(cd[[on]], unlevel=unlevel)
        newred = lapply(SingleCellExperiment::reducedDims(se), function(x) {
            t(to_AP) %*% x
        })
        return(SingleCellExperiment::SingleCellExperiment(assays=P,
                                                          rowData=rd,
                                                          colData=pcd,
                                                          metadata=S4Vectors::metadata(se),
                                                          reducedDims=S4Vectors::SimpleList(newred)))
    } else {
        return(SummarizedExperiment::SummarizedExperiment(assays=P,
                                                          rowData=rd,
                                                          colData=pcd,
                                                          metadata=S4Vectors::metadata(se)))
    }
}


#' Runs PCA on the SummarizedExperiment object
#'
#' @param se SummarizedExperiment
#' @param n Number of PCs
#' @param assay Assay to use for PCA
#' @param retx Whether to return rotation
#' @param center Whether to center observations
#' @param scale. Whether to scale the data
#' @param ... Extra parameters passed to irlba::prcomp_irlba
#' @return Updated SummarizedExperiment object
#' @export
se_prcomp <- function(se, n=3, assay="TMM", retx=TRUE, center=TRUE, scale.=FALSE, ...) {
    X = SummarizedExperiment::assays(se)[[assay]]
    stopifnot(!is.null(X))
    n = min(c(n, dim(se)-1))
    pca = irlba::prcomp_irlba(X, n=n, retx=retx, center=center, scale.=scale., ...)
    if (retx) {
        for (cn in colnames(pca$retx)) {
            SummarizedExperiment::rowData(se)[[cn]] = pca$rotation[,cn]
        }
    }
    for (cn in colnames(pca$x)) {
        SummarizedExperiment::colData(se)[[cn]] = pca$x[,cn]
    }
    SummarizedExperiment::metadata(se)$pca = list(scale=pca$scale,
                                                  totalvar=pca$totalvar,
                                                  sdev=pca$sdev,
                                                  center=pca$center,
                                                  n=n)
    return(se)
}
#' Calculate qc metrics similar to Scanpy
#'
#' @param se SummarizedExperiment
#' @param assay Assay to calculate upon
#' @param log1p logical whether to use log statistics as well
#' @export
calculate_qc_metrics <- function(se, assay=NULL, log1p=TRUE, qc_vars=c()) {
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
    for (qcv in intersect(qc_vars, colnames(SummarizedExperiment::rowData(se)))) {
### pct_counts_${qcv}
### check if logical
        vals = SummarizedExperiment::rowData(se)[[qcv]]
        if (is.logical(vals)) {
            tc_qc = setNames(Matrix::colSums(X[vals,]), NULL)
            SummarizedExperiment::colData(se)[[paste0("total_counts_", qcv)]] = tc_qc
            SummarizedExperiment::colData(se)[[paste0("pct_counts_", qcv)]] = 100 * tc_qc / cs
            if (log1p) {
                SummarizedExperiment::colData(se)[[paste0("log1p_total_counts_", qcv)]] = log1p(tc_qc)
            }
        }
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
              "matrix/features/feature_type"=ifelse(1:nrow(X) %in% grep(".*[:][0-9]+-[0-9]+$", rownames(X)),
                                                    "Chromatin Accessibility",
                                                    "Gene Expression"),
              "matrix/features/genome"=rep(genome, nrow(X)),
              "matrix/shape"=t(dim(X)))
    if (!is.null(SummarizedExperiment::rowRanges(se))) {
        ml[["matrix/features/interval"]] = with(as.data.frame(SummarizedExperiment::rowRanges(se)), paste0(seqnames, ":", start, "-", end))
    }
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
    out = X %*% Matrix::Diagonal(x=target_sum / Matrix::colSums(X))
    dimnames(out) = dimnames(X)
    SummarizedExperiment::assays(se)[[normalized]] = out
    return(se)
}

#' Log1p assay similar to Scanpy
#'
#' @param se SummarizedExperiment object
#' @param assay Base assay
#' @export
se_log1p <- function(se, assay="normalized") {
    X = SummarizedExperiment::assays(se)[[assay]]
    SummarizedExperiment::assays(se)[[paste0("log1p_", assay)]] = log1p(X)
    return(se)
}

#' Un-normalize SummarizedExperiment using total_counts column
#' @param se SummarizedExperiment object
#' @param total_counts total_counts column in colData(se) to extract real counts from
se_unnormalize <- function(se, assay=NULL, unnormalized="raw", eps=0.01, total_counts="total_counts") {
    if (is.null(assay)) {
        E = SummarizedExperiment::assays(se)[[1]]
    } else {
        E = SummarizedExperiment::assays(se)[[assay]]
    }
    cs = Matrix::colSums(E)
    if (sd(cs) < eps) {
        ### first check if is just CPM normalized
        M = E %*% Matrix::Diagonal(x=colData(se)[[total_counts]]/mean(cs))
    } else {
        ### otherwise assume log1p-CPM normalization
        cs = Matrix::colSums(expm1(E))
        M = NULL
        if (sd(cs) < eps) {
            M = expm1(E) %*% Matrix::Diagonal(x=colData(se)[[total_counts]]/mean(cs))
        }
    }
    if (!is.null(M)) {
        diff = abs(colSums(E) - colData(se)$total_counts)
        print(paste0("max diff: ", mean(diff), " counts per cell"))
        print(paste0("Saving to '", unnormalized, "'"))
        dimnames(M) = dimnames(E)
        assays(se)[[unnormalized]] = round(M)
    }
    return(se)
}

#' Concatenate SummarizedExperiment objects, throwing out non-matching rowData columns
#'
#' @param se.list List of SummarizedExperiment objects
#' @export
se_concat <- function(se.list) {
    cat(paste0("Concatenating ", length(se.list), " objects\n"))
    se.list = Filter(function(x) { S4Vectors::ncol(x) > 0 }, se.list)
    cat(paste0("Concatenating ", length(se.list), " nonzero objects\n"))
    rowDatas = lapply(se.list, SummarizedExperiment::rowData)
    colDatas = lapply(se.list, SummarizedExperiment::colData)
    stopifnot(do.call(all.equal, lapply(rowDatas, rownames)))
    metaDatas = lapply(se.list, S4Vectors::metadata)
    good.row.cols = as.logical(sapply(setNames(colnames(rowDatas[[1]]), colnames(rowDatas[[1]])), function(cn) {
        all(sapply(rowDatas, function(rd) {
            identical(rd[[cn]], rowDatas[[1]][[cn]])
        }))
    }))
    good.col.cols = as.logical(Reduce(intersect, lapply(colDatas, colnames)))
    se = do.call(SummarizedExperiment::cbind, lapply(se.list, function(se) {
      rd = as.data.frame(SummarizedExperiment::rowData(se))
      if (sum(good.row.cols, na.rm=TRUE) > 0) { 
        SummarizedExperiment::rowData(se) = rd[,names(which(good.row.cols)),drop=FALSE]
      }
      cd = as.data.frame(SummarizedExperiment::colData(se))
      if (sum(good.col.cols, na.rm=TRUE) > 0) {
        SummarizedExperiment::colData(se) = S4Vectors::DataFrame(cd[,good.col.cols,drop=FALSE])
      }
      return(se)
    }))
    S4Vectors::metadata(se)$h5ad = sapply(metaDatas, function(md) {
        md$h5ad
    })
    return(se)
}

se_as_Seurat <- function(sce) {
    C = SummarizedExperiment::assays(sce)$counts
    dimnames(C) <- dimnames(sce)
    cd = SummarizedExperiment::colData(sce)
    rownames(cd) <- colnames(sce)
    s.obj = Seurat::CreateSeuratObject(counts=C, meta.data=as.data.frame(cd))
    s.obj[["UMAP"]] <- Seurat::CreateDimReducObject(
         embeddings = SingleCellExperiment::reducedDims(sce)$X_umap,
         key = "UMAP_", assay = "RNA")
    s.obj[["X_pca"]] <- Seurat::CreateDimReducObject(
         embeddings = SingleCellExperiment::reducedDims(sce)$X_pca,
         key = "PC_", assay = "RNA")
    s.obj = Seurat::NormalizeData(s.obj)
    return(s.obj)
}

ad_as_Seurat <- function(h5ad) {
  s.obj = anndataR::read_h5ad(h5ad, to="Seurat")
  obsm = rhdf5::h5read(h5ad, "obsm")
  colnames(obsm$X_umap) <- rhdf5::h5read(h5ad, "/obs/_index")
  s.obj[["UMAP"]] <- Seurat::CreateDimReducObject(embeddings=t(obsm$X_umap), key="UMAP_", assay="RNA")
}
