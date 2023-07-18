
#' Parse H5AD dataframe
#'
#' This function converts an object from rhdf5::h5read into a usable dataframe
#'
#' @param obj List of objects as read from rhdf5::h5read(h5, "obs") or rhdf5::h5read(h5, "var")
#' @return A properly formatted dataframe
.parse_h5ad_dataframe <- function(obj, index_name="_index") {
    if (!is.list(obj)) {
        stop("Not a list!")
    }
    if (index_name %in% names(obj)) {
        df = data.frame(row.names=var_names_make_unique(obj[[index_name]]))
        bad.cols = c("_index", "__categories")
    } else { ## old version of anndata (pbmc3k)
        df = data.frame(row.names=var_names_make_unique(obj[["index"]]))
        bad.cols = c("index", "__categories")
    }
    for (cname in setdiff(names(obj), bad.cols)) {
        item = obj[[cname]]
        if (is.list(item)) {
            categories = item$categories
            codes = item$codes + 1
            codes[codes <= 0] = NA
            df[[cname]] = factor(categories[codes], levels=categories)
        } else {
            df[[cname]] = c(item)
        }
    }
    if (("__categories" %in% names(obj)) & (is.list(obj[["__categories"]]))) {
        for (cname in names(obj[["__categories"]])) {
            categories = c(obj[["__categories"]][[cname]])
            codes = df[[cname]] + 1
            codes[codes <= 0] = NA
            df[[cname]] = factor(categories[codes], levels=categories)
        }
    }
    return(df)
}


#' Read H5AD obs
#'
#' This function takes an H5 filename and reads the .obs dataframe
#'
#' @param h5ad An H5AD filename path
#' @return A properly formatted dataframe
#' @export
read_h5ad_obs <- function(h5ad, base="/") {
    if (length(h5ad) == 1) {
        df = rhdf5::h5read(h5ad, paste0(base, "/obs"))
        return(.parse_h5ad_dataframe(df))
    } else {
        df.list = parallel::mclapply(h5ad, function(h5) {
            .parse_h5ad_dataframe(rhdf5::h5read(h5, paste0(base, "/obs")))
        })
### Now to combine obs
        df = do.call(rbind, lapply(df.list, function(dfl) {
            dfl$df
        }))
        stopifnot(all(do.call(c, lapply(df.list, rownames)) %in% rownames(df)))
        return(df)
    }
}

#' Read H5AD var
#'
#' This function takes an H5 filename and reads the .var dataframe
#'
#' @param h5ad An H5AD filename path
#' @return A properly formatted dataframe
#' @export
read_h5ad_var <- function(h5ad, base="/") {
    df = rhdf5::h5read(h5ad, paste0(base, "/var"))
    if (!("_index" %in% names(df)) & ("symbol" %in% names(df))) {
        return(.parse_h5ad_dataframe(df, "symbol"))
    } else {
        return(.parse_h5ad_dataframe(df))
    }
}

#' Parse observed barcode dataframe
#'
#' This function reads .obs, optionally
#' subsetting to barcodes passed from
#' a dataframe or a vector of valid barcodes
#'
#' @param h5ad An H5AD filename path
#' @param obs Either NULL (all data), a vector of valid UMIs, or a dataframe of .obs with UMI rownames
#' @param subset A list of column-value(s) mappings to subset OBS to
#' @return A list of the dataframe (df) and the integer index with df's rownames as names
.parse_h5ad_obs <- function(h5ad, obs=NULL, subset=list(), base="/", refactor=TRUE) {
    obs_df = read_h5ad_obs(h5ad, base=base)
    obs_index = 1:nrow(obs_df)
    if (!is.null(obs)) {
        if (is.data.frame(obs)) {
            obs_index = match(rownames(obs), rownames(obs_df))
            obs_df = obs_df[rownames(obs),]
            for (cn in colnames(obs)) {
                ### Add new columns as well
                obs_df[[cn]] = obs[[cn]]
            }
        } else {
            obs_index = match(obs, rownames(obs_df))
            obs_df = obs_df[obs,]
        }
    }
    if (length(subset) > 0) {
        for (i in seq_along(subset)) {
            if (names(subset)[i] %in% colnames(obs_df)) {
                cn = names(subset)[i]
                if (any(subset[[i]] %in% obs_df[[cn]])) {
                    flag = obs_df[[cn]] %in% subset[[i]]
                    obs_df = obs_df[flag, ]
                    obs_index = obs_index[flag]
                } else {
                    warning(paste0("No object(s) ", subset[[i]], " in column ", cn))
                }
            } else {
                warning(paste0("Column ", cn, " is not in .obs"))
            }
        }
    }
    for (cn in colnames(obs_df)) {
        if (refactor & is.factor(obs_df[[cn]]) & (min(table(obs_df[[cn]]))==0)) {
            obs_df[[cn]] = as.factor(as.character(obs_df[[cn]]))
        }
    }
    names(obs_index) = rownames(obs_df)
    return(list(df=obs_df, index=obs_index))
}


#' Parse variables dataframe
#'
#' This function reads .var, optionally
#' subsetting to features passed from
#' a dataframe or a vector of valid features
#'
#' @param h5ad An H5AD filename path
#' @param var Either NULL (all features), a vector of valid features, or a dataframe of .var with feature rownames
#' @return A list of the dataframe (df) and the integer index with df's rownames as names
.parse_h5ad_var <- function(h5ad, var=NULL, base="/") {
    if (is.null(var)) {
        var_df = read_h5ad_var(h5ad, base=base)
        var_index = 1:nrow(var_df)
    } else if (is.data.frame(var)) {
        varnames = rownames(read_h5ad_var(h5ad, base=base))
        var_df = var
        var_index = match(rownames(var_df), varnames)
    } else {
        var_df = read_h5ad_var(h5ad, base=base)
        var_index = match(var, rownames(var_df))
        var_df = var_df[var_index,]
    }
    names(var_index) = rownames(var_df)
    return(list(df=var_df, index=var_index))
}


#' Take indptr from H5AD and turn into a dense index
#' Take in indptr assuming it is integer64
#' But use regular integer in rep() as number of compressed indices is not likely going to overflow
.decompress_indptr <- function(indptr, obs_index=NULL) {
    return(do.call(c, lapply(seq(1, length(indptr)-1), function(i) {
        begin = indptr[i] + 1
        end = indptr[i + 1]
        rep(match(i, obs_index), as.integer(end - begin + 1))
    })))
}

#' Parse a sparse matrix given indices
#'
#' THis function takes a matrix and a list of indices to transform
#' into a Matrix::sparseMatrix
#'
#' @param h5 An H5AD filename path
#' @param matrix A matrix object in the H5 file
#' @param obs_index A named vector of observed (barcode) indices
#' @param var_index A named vector of variable (gene) indices
#' @return A sparse matrix
.get_sparse <- function(h5, matrix, obs_index, var_index) {
    indptr = rhdf5::h5read(h5, paste0(matrix, "/indptr"), bit64conversion="bit64")
    if (is.null(obs_index)) {
        obs_index = seq(1, length(indptr)-1)
    }
    if (is.null(var_index)) {
        stop("Not implemented")
    }
    VI = match(rhdf5::h5read(h5, paste0(matrix, "/indices")) + 1,
                    var_index)
    data = rhdf5::h5read(h5, paste0(matrix, "/data"), bit64conversion="bit64")
    OI = .decompress_indptr(indptr, obs_index)
    flag = (!is.na(OI)) & (!is.na(VI))
    return(Matrix::sparseMatrix(i=VI[flag],
                                j=OI[flag],
                                x=c(data[flag]),
                                dims=c(length(var_index),
                                       length(obs_index))))
}

#' Parse a matrix given indices
#'
#' THis function takes a matrix and a list of indices to transform
#' into a matrix or Matrix::sparseMatrix
#'
#' @param h5 An H5AD filename path
#' @param matrix A matrix object in the H5 file
#' @param obs_index A named vector of observed (barcode) indices
#' @param var_index A named vector of variable (gene) indices
#' @return A named matrix
.parse_h5ad_X <- function(h5, matrix, obs_index=NULL, var_index=NULL) {
    #matrix = paste0(dirname(matrix), "/", basename(matrix))
    if (matrix %in% rhdf5::h5ls(h5)$group) {
        M = .get_sparse(h5, matrix, obs_index, var_index)
    } else {
        M = rhdf5::h5read(h5, matrix, index=list(var_index, obs_index))
    }
    if (!is.null(obs_index)) {
        colnames(M) = names(obs_index)
    }
    if (!is.null(var_index)) {
        rownames(M) = names(var_index)
    }
    return(M)
}

.get_dimreduc <- function(h5, group, obs_index=NULL, var_index=NULL) {
    h5list = rhdf5::h5ls(h5)
    rd_names = h5list[(h5list$group == group) & (h5list$otype %in% c("H5I_GROUP", "H5I_DATASET")),]$name
    names(rd_names) = rd_names
    return(lapply(rd_names, function(name) {
#        print(paste0("DR name:", name))
        .parse_h5ad_X(h5, paste0(group, "/", name),
                      obs_index=obs_index,
                      var_index=var_index)
   }))
}

.sparse2selfhits <- function(M) {
    N = nrow(M)
    if (is.matrix(M)) {
        M = reshape2::melt(M, varnames=c("i", "j"), value.name="x")
    } else {
        M = Matrix::summarize(M)
    }
    hits = SelfHits(
        mdf$i,
        mdf$j,
        nnode=N
    )
    mcols(hits)$value = mdf$x
    return(hits)
}
#' Read H5AD matrix
#'
#' This function reads the matrix from an H5AD file
#' and optionally subsets to only certain genes.
#' If raw, will only take raw variables and then return,
#' since .raw can contain a different set of genes.
#'
#' @param h5ad An H5AD filename path
#' @param obs A dataframe of .obs or vector of obsnames to subset by
#' @param var A dataframe of .var or vector of varnames to subset by
#' @param subset A list of column-value(s) pairs to subset .obs by, e.g. list("CellType"=c("Exc","Inh"))
#' @param layer A vector of layer(s) to subset
#' @return A SingleCellExperiment object filled with assays
#' @export
read_h5ad <- function(h5ad, obs=NULL, var=NULL, layer=NULL, raw=FALSE, obsm=TRUE, obsp=TRUE, varp=TRUE, base="/", subset=list(), refactor=TRUE) {
    obs = .parse_h5ad_obs(h5ad, obs, subset=subset, base=base, refactor=refactor)
    var = .parse_h5ad_var(h5ad, var, base=paste0(base, ifelse(raw, "/raw/", "/")))
    if (raw) {
        assays = list(counts=.parse_h5ad_X(h5ad, paste0(base, "/raw/X"),
                                           obs_index=obs$index,
                                           var_index=var$index))
        sce = SingleCellExperiment::SingleCellExperiment(assays, colData=obs$df, rowData=var$df)
    } else {
        assays = list(counts=.parse_h5ad_X(h5ad, "/X",
                                           obs_index=obs$index,
                                           var_index=var$index))
        for (L in layer) {
            assays[[L]] = .parse_h5ad_X(h5ad, paste0("/layers/", L),
                                        obs_index=obs$index,
                                        var_index=var$index)
        }
        sce = SingleCellExperiment::SingleCellExperiment(assays, colData=obs$df, rowData=var$df)
    }
    if (obsm) {
        rd = .get_dimreduc(h5ad, group="/obsm", obs_index=obs$index)
        SingleCellExperiment::reducedDims(sce) = lapply(rd, Matrix::t)
    }
    if (obsp) {
        ## print("OBSP:")
        cp = .get_dimreduc(h5ad, group="/obsp", obs_index=obs$index, var_index=obs$index)
        ## print(str(cp))
        ## SingleCellExperiment::colPairs(sce, asSparse=TRUE) = lapply(cp, Matrix::t)
    }
    if (varp & !raw) {
        #print("reading varp")
    }
    S4Vectors::metadata(sce)$subset = subset
    S4Vectors::metadata(sce)$h5ad = h5ad
    S4Vectors::metadata(sce)$layer = layer
    return(sce)
}
