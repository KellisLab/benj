#!/usr/bin/env Rscript


## 1. get args
## 2. get obs
## 3. benj::read_h5ad
## 4. set metadata (obs=path)
## 5. benj::deg

parse_subset <- function(subset.args, sep=Sys.getenv("DEG_sep", ","), assign=Sys.getenv("DEG_assign", "=")) {
    subset_str = paste0(make.names(gsub("_", "-", subset.args)), collapse="_")
    names(subset.args) = sapply(strsplit(subset.args, assign), "[[", 1)
    return(list(subset_str=subset_str,
                subset=lapply(subset.args, function(item) {
                    strsplit(gsub(paste0("^[^", assign, "]+="), "", item), sep)[[1]]
    })))
}

get_args <- function(args) {
### options: aggr
    params = list()
    params$cpm.cutoff = 10
    params$sample.col = "Sample"
    params$NRUV = 0
    params$IQR.factor = 1.5
    params$outlier.covariates = c("log1p_total_counts", "n_genes_by_counts")
    params$min.total.counts.per.sample = 100
    params$ncores = as.integer(Sys.getenv("OMP_NUM_THREADS", getOption("mc.cores", 2)))
    i = 1
    help = function() {
        cat("deg.R [options]
options:
  -h, --help Show this message and exit.
  -i, --h5ad List of H5AD object(s) to read in.
  -a, --annotation Annotation if necessary, to add to colData to subset
  -s, --subset Subset the H5AD objects with key=value assignments. Multiple values could be comma separated, like key1=valueA,valueB,valueC key2=valueX
  -o, --output XLSX output file
  -j, --ncores Number of cores to use. By default, use OMP_NUM_THREADS
  --pathology Column in colData/.obs that separates pathology/diagnosis.
  --case  Value within --pathology column that specifices cells that are cases.
  --control Value within --pathology column that specifies cells that are controls.
  --method Method(s) (whitespace separated) to be used for DEG analysis. Supported are DESeq2, edgeR-lrt, edgeR-QL, nebula, mast, mast-re.
  --covariates Covariates (in colData/.obs) to be used as covariates for DEG calling. If pseudobulk, may be removed if not a subject-level variable
  --sample-col Column corresponding to sample for pseudobulk
  --ruv Number of subject-level RUV terms to include
  --cpm-cutoff Counts per million cutoff for filtering genes.
  --min-total-counts Number of total counts required per cell
  --iqr-factor Number of IQR above 75% away from outlier covariates to be a bad sample.
  --outlier-covariates Covariates (in colData/.obs) to be used as covariates for pseudobulk outlier detection. If not present, or not calculated by benj::calculate_qc_metrics, will be thrown out.
")
        stop("Help")
    }
    while (i <= length(args)) {
        arg = args[[i]]
        if (arg %in% c("-h", "--help")) {
            help()
        } else if (arg %in% c("-i", "--h5ad")) {
            i = i + 1
            while(i <= length(args) && !startsWith(args[[i]], "-")) {
                params$h5ad = c(params$h5ad, args[[i]])
                i = i + 1
            }
        } else if (arg %in% c("-a", "--annotation")) {
            i = i + 1
            params$annotation = args[[i]]
            i = i + 1
        } else if (arg %in% c("-s", "--subset")) {
            assign = Sys.getenv("DEG_assign", "=")
            i = i + 1
            subset.args = NULL
            while(i <= length(args) & !startsWith(args[[i]], "-")) {
                subset.args = c(subset.args, args[[i]])
                i = i + 1
            }
            names(subset.args) = sapply(strsplit(subset.args, assign), "[[", 1)
            params$subset = lapply(subset.args, function(item) {
                strsplit(gsub(paste0("^[^", assign, "]+", assign), "", item), Sys.getenv("DEG_sep", ","))[[1]]
            })
        } else if (arg %in% c("-o", "--output")) {
            i = i + 1
            params$output = args[[i]]
            i = i + 1
        } else if (arg %in% c("-j", "--cores")) {
            i = i + 1
            params$ncores = as.integer(args[[i]])
            i = i + 1
        } else if (arg %in% c("--path", "--pathology")) {
            i = i + 1
            params$pathology = args[[i]]
            i = i + 1
        } else if (arg == "--case") {
            i = i + 1
            params$case = args[[i]]
            i = i + 1
        } else if (arg == "--control") {
            i = i + 1
            params$control = args[[i]]
            i = i + 1
        } else if (arg == "--covariates") {
            i = i + 1
            while(i <= length(args) && !startsWith(args[[i]], "-")) {
                params$covariates = c(params$covariates, args[[i]])
                i = i + 1
            }
        } else if (arg %in% c("-m", "--method")) {
            i = i + 1
            while(i <= length(args) && !startsWith(args[[i]], "-")) {
                params$method = c(params$method, args[[i]])
                i = i + 1
            }
        } else if (arg == "--sample-col") {
            i = i + 1
            params$sample.col = args[[i]]
            i = i + 1
        } else if (arg %in% c("-r", "--ruv")) {
            i = i + 1
            params$NRUV = as.integer(args[[i]])
            i = i + 1
        } else if (arg == "--cpm-cutoff") {
            i = i + 1
            params$cpm.cutoff = as.numeric(args[[i]])
            i = i + 1
        } else if (arg == "--min-total-counts") {
            i = i + 1
            params$min.total.counts.per.sample = as.numeric(args[[i]])
            i = i + 1
        } else if (arg == "--iqr-factor") {
            i = i + 1
            params$IQR.factor = as.numeric(args[[i]])
            i = i + 1
        } else if (arg == "--outlier-covariates") {
            i = i + 1
            params$outlier.covariates = NULL
            while(i <= length(args) && !startsWith(args[[i]], "-")) {
                params$outlier.covariates = c(params$outlier.covariates, args[[i]])
                i = i + 1
            }
        } else {
            stop(paste0("Unknown option ", args[[i]]))
        }
    }
    if (is.null(params$subset)) {
        params$subset = list()
    }
    print(jsonlite::toJSON(params, pretty=TRUE, auto_unbox=TRUE))
    stopifnot(!is.null(params$pathology))
    stopifnot(!is.null(params$case))
    stopifnot(!is.null(params$control))
    stopifnot(!is.null(params$method))
    stopifnot(!is.null(params$output))
    stopifnot(!is.null(params$h5ad))
    return(params)
}

load.anndata <- function(params) {
    obs = NULL
    if (!is.null(params$annotation)) {
        stopifnot(file.exists(params$annotation))
        obs = read.table(params$annotation, sep="\t", header=TRUE, comment.char="", row.names=1)
    }

    if (length(params$h5ad) > 1) {
        adata = benj::read_h5ad_parallel(h5ad=params$h5ad, obs=obs, subset=params$subset, ncores=params$ncores)
    } else {
        adata = benj::read_h5ad(h5ad=params$h5ad, obs=obs, subset=params$subset)
    }
    if (!is.null(params$annotation)) {
        S4Vectors::metadata(adata)$annotation = params$annotation
    }
    return(adata)
}

### TODO filter h5ad in read_h5ad_parallel by removing bad samples
params = get_args(commandArgs(trailing=TRUE))
adata = load.anndata(params)
adata = benj::deg(adata,
                  pathology=params$pathology,
                  case=params$case,
                  control=params$control,
                  covariates=params$covariates,
                  method=params$method,
                  output=params$output,
                  sample.col=params$sample.col,
                  cpm.cutoff=params$cpm.cutoff,
                  NRUV=params$NRUV,
                  min.total.counts.per.sample=params$min.total.counts.per.sample,
                  IQR.factor=params$IQR.factor,
                  outlier.covariates=params$outlier.covariates)
