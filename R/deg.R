

#' Filtering of design matrix for DEG comparisons
#' This method modifies a design matrix such that combinations or unused comparisons are removed.
#' It also renames items with spaces because some DEG methods complain about spaces in the design matrix.
#' @param design Design matrix.
#' @param rename Whether to use make.names to rename columns
#' @export
deg.filter.design <- function(design, rename=TRUE) {
    nzv = caret::nearZeroVar(design, saveMetrics=TRUE)
    design = design[,!nzv$nzv | rownames(nzv) == "(Intercept)",drop=FALSE]
    tryCatch({
        linear_combos = caret::findLinearCombos(design)
        if (!is.null(linear_combos$remove)) {
            design = design[,-linear_combos$remove,drop=FALSE]
        }
    }, error=function(e) {
        message("Disregarded error: ", e)
        print(str(design))
    })
    if (rename & is.matrix(design)) {
        if (ncol(design) > 1) {
            colnames(design) = make.names(colnames(design)) ### need spaces to be removed
        }
    }
    return(design)
}


#' When a design matrix is filtered, we need to be able to extract the covariates used for MAST::zlm.
#' Here, we use the filtered design to find if any columns have been themselves filtered.
#' @param design The design matrix
#' @param covariates.used A vector of the column names used for constructing the design matrix
#' @return A filtered list of covariates
#' @export
deg.extract.covariates.from.design <- function(design, covariates.used) {
    contrasts = names(attr(design, "contrasts"))
    return(covariates.used[covariates.used %in% c(contrasts, colnames(design))])
}
#' Filter outlying samples from a pseudobulk SummarizedExperiment.
#' This method uses covariates only, and finds distance from the centroid in PCA space
#' @param se Pseudobulked SummarizedExperiment object
#' @param covariates Covariates that are scaled and PCA'd to compute outliers
#' @param IQR.factor Number of IQR above 75%ile to use for exclusion of samples. Set to Inf for no exclusion
#' @return Updated SummarizedExperiment object with bad samples removed
#' @export
deg.filter.outliers <- function(se, covariates=c("log1p_total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"), IQR.factor=1.5) {
    cd = as.data.frame(SummarizedExperiment::colData(se))
    covariates = covariates[covariates %in% colnames(deg.filter.design(cd[c(covariates)]))]
    design = model.matrix(as.formula(paste0(c("~0", covariates), collapse="+")),
                          data=cd)
    design = deg.filter.design(design)
    if (ncol(design)==0) {
        ### Remove in case no covariates are passed
        return(se)
    }
    pca = prcomp(design, scale.=TRUE)$x
    centroid = colMeans(pca)
    distances = sqrt(rowSums((pca - centroid)^2))
    Q1 = quantile(distances, 0.25)
    Q3 = quantile(distances, 0.75)
    upper_bound = Q3 + IQR.factor * (Q3 - Q1)
    outliers = which(distances > upper_bound)
    if (length(outliers) > 0) {
        outliers = colnames(se)[outliers]
        cat("Filtering", length(outliers), "outlier(s):", outliers, "\n")
        se = se[,setdiff(colnames(se), outliers)]
    }
    return(se)
}
#' Prepare SummarizedExperiment object for DEG calling.
#' @export
deg.prepare <- function(se, pathology, case, control, sample.col, filter_only_case_control=TRUE,
                        cpm.cutoff=10, min.total.counts.per.sample=100, IQR.factor=1.5,
                        outlier.covariates=c("log1p_total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"),
                        ensure.integer.counts=TRUE) {
    if (filter_only_case_control) {
        if (pathology %in% colnames(SummarizedExperiment::colData(se))) {
            se = se[,SummarizedExperiment::colData(se)[[pathology]] %in% c(case, control)]
        } else {
            stop(paste0("Pathology not present in colData"))
        }
    }
    stopifnot("counts" %in% names(SummarizedExperiment::assays(se)))
    pb = calculate_qc_metrics(se_make_pseudobulk(se, sample.col), assay="counts", qc_vars=c("mt", "ribo", "pc"))
    stopifnot(S4Vectors::ncol(pb) <= 1000) ### we can't have uber-large sample numbers
    if (any(SummarizedExperiment::colData(pb)$total_counts < min.total.counts.per.sample)) {
        cd = SummarizedExperiment::colData(pb)
        cat("Bad samples:", rownames(cd)[cd$total_counts < min.total.counts.per.sample], "\n")
        se = se[, SummarizedExperiment::colData(se)[[sample.col]] %in% rownames(cd)[cd$total_counts >= min.total.counts.per.sample]]
        pb = se_make_pseudobulk(se, sample.col)
    }
    ### Outlier detection
    pb = deg.filter.outliers(pb, covariates=outlier.covariates, IQR.factor=IQR.factor)
    cpm.flag = rownames(se)
    if (!is.na(cpm.cutoff) & !is.null(cpm.cutoff)) {
### Get minimum samples for CPM threshold to pass
        cpm.n = min(table(SummarizedExperiment::colData(pb)[[pathology]]))
    ### Filter genes above a cpm cutoff
        cpm.flag = edgeR::cpm(SummarizedExperiment::assays(pb)$counts) > cpm.cutoff
    ### Then extract across number of samples
        cpm.flag = rownames(cpm.flag)[Matrix::rowSums(cpm.flag) >= cpm.n]
    }
    pb = pb[cpm.flag, ]
    se = se[rownames(pb), SummarizedExperiment::colData(se)[[sample.col]] %in% colnames(pb)]

### Get logCPM info
    logCPM = edgeR::cpmByGroup(pb, SummarizedExperiment::colData(pb)[[pathology]], log=TRUE)
    colnames(logCPM) = paste0("logCPM_", colnames(logCPM))
    SummarizedExperiment::rowData(se) = cbind(SummarizedExperiment::rowData(se), as.data.frame(logCPM))
### Check matrix
    X = SummarizedExperiment::assays(se)$counts
    if (inherits(X, "Matrix")) {
        is_integer_matrix = all(X@x == round(X@x))
    } else {
        is_integer_matrix = all(apply(X, c(1, 2), function(x) { round(x) == x }))
    }
    if (ensure.integer.counts) {
        stopifnot(is_integer_matrix)
    }
    S4Vectors::metadata(se)$deg = list(pathology=pathology,
                                       case=case,
                                       control=control,
                                       sample_column=sample.col,
                                       filter_only_case_control=filter_only_case_control,
                                       cpm_cutoff=cpm.cutoff,
                                       min_total_counts_per_sample=min.total.counts.per.sample,
                                       IQR_factor=IQR.factor,
                                       outlier_covariates=paste0(outlier.covariates, collapse=" + "))

### Re-order to group by sample
    SummarizedExperiment::colData(se)[[sample.col]] = as.factor(as.character(SummarizedExperiment::colData(se)[[sample.col]]))
    se = se[,order(SummarizedExperiment::colData(se)[[sample.col]])]

    ### Remove levels with no comparisons.
    SummarizedExperiment::colData(se) = droplevels(SummarizedExperiment::colData(se))
    return(se)
}
#' Run DEGs, switching on the method
#'
#' Note: It is recommended to save subset information in metadata(se)$subset
#' since there is *technically* no way to recover that. However, from sample_metadata this should be inferrable.
#' @param se SummarizedExperiment object, with integer counts in the "counts" assay.
#' @param pathology Column in colData(se) corresponding to diagnosis or experimental condition
#' @param case Column in colData(se) corresponding to case values, which must be in "pathology"
#' @param control Column in colData(se) corresponding to control values in "pathology"
#' @param covariates Covariates to be used. Non-existent variables will be removed, and RUV variables will be added (if NRUV>0).
#' @param sample.col Sample column in colData(se) corresponding to the experiment. This is the column on which data is pseudobulked.
#' @param method Method(s) of DEGs to be used. Valid values are DESeq2, edgerR-GLM, edgeR-QL, nebula, mast-hurdle, mast-random-effect
#' @param output Output XLSX
#' @param cpm.cutoff Counts per Million cutoff used for filtering lowly-expressed genes. Genes must be expressed in at least min(case, control) samples at least at this log(CPM) level to be used. This n=min(case, control) is after outlier removal and low-count data.
#' @param NRUV Number of RUVSeq components to include
#' @param filter_only_case_control Filter se to only case+control only. Default TRUE
#' @param IQR.factor Number of IQR factors above 75% percentile to be removed, according to scaled PCA of outlier.covariates. Set to Inf to remove none.
#' @param outlier.covariates Sample-level covariates which are scaled and PCA transformed, to get distances from the centroid of the PCA embedding. Outliers above 75% + IQR.factor * IQR are removed.
#' @export
deg <- function(se, pathology, case, control, covariates,
                method, output=NULL,
                sample.col="Sample", cpm.cutoff=10,
                filter_only_case_control=TRUE, NRUV=0,
                min.total.counts.per.sample=100, IQR.factor=1.5,
                outlier.covariates=c("log1p_total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"),
                verbose=TRUE,
                ensure.integer.counts=TRUE) {
    method = match.arg(gsub(" ", "-", tolower(method)),
                       c("deseq2", "edger", "edger-lrt", "edger-ql", "nebula", "mast", "mast-re", "lmer"), several.ok=TRUE)
    se = deg.prepare(se, pathology=pathology,
                     case=case, control=control,
                     sample.col=sample.col,
                     filter_only_case_control=filter_only_case_control,
                     min.total.counts.per.sample=min.total.counts.per.sample,
                     cpm.cutoff=cpm.cutoff, outlier.covariates=outlier.covariates,
                     ensure.integer.counts=ensure.integer.counts)
    ### RUVSeq
    covariates = covariates[covariates %in% names(SummarizedExperiment::colData(se))]
    if (NRUV > 0) {
        se = deg.ruvseq(se, pathology=pathology, covariates=covariates, sample.col=sample.col, NRUV=NRUV, verbose=verbose)
        covariates = c(covariates, paste0("RUV_", 1:NRUV))
        covariates = covariates[covariates %in% names(SummarizedExperiment::colData(se))]
    }
    S4Vectors::metadata(se)$deg$covariates = paste0(covariates, collapse=" + ")
    ### Iterate through methods
    for (meth in method) {
        if (meth == "deseq2") {
            se = deg.deseq2(se, pathology=pathology,
                            case=case, control=control,
                            sample.col=sample.col,
                            covariates=covariates)
        } else if (grepl("^edger", meth)) {
            se = deg.edger(se, pathology=pathology,
                           case=case, control=control,
                           sample.col=sample.col,
                           covariates=covariates,
                           method=ifelse(grepl("lrt$",meth), "LRT", "QL"))
        } else if (meth == "nebula") {
            se = deg.nebula(se, pathology=pathology,
                            case=case, control=control,
                            sample.col=sample.col,
                            covariates=covariates,
                            offset="total_counts", model="NBGMM")
        } else if (grepl("^mast", meth)) {
            se = deg.mast(se, pathology=pathology,
                          case=case, control=control,
                          sample.col=sample.col,
                          covariates=covariates,
                          method=ifelse(grepl("re$", meth), "RE", "Hurdle"))
        } else if (meth == "lmer") {
            se = deg.lmer(se, pathology=pathology,
                          case=case, control=control,
                          sample.col=sample.col,
                          covariates=covariates)
        } else {
            stop(paste0("Invalid method ", meth))
        }
        if (!is.null(output)) {
            dump.excel(se_make_pseudobulk(se, sample.col), output)
        }
    }
    return(se)
}

#' Run edgeR on a SummarizedExperiment
#'
#' @param se SummarizedExperiment to run
#' @param pathology Column in colData(se) that contains case and control
#' @param case Case string within pathology column. Example "AD"
#' @param control Control string within pathology column. Example "control"
#' @param covariates Vector of covariates to add to design matrix
#' @param method Method used for edgeR. Either LRT or QL
#' @param sample.col Sample column to pseudobulk on
#' @param prefix Prefix to add in rowData(se) where edgeR info is added
#' @export
deg.edger <- function(se, pathology, case, control,
                      sample.col, covariates=NULL,
                      method=c("LRT", "QL"), prefix="edgeR") {
    method = match.arg(gsub("^EDGER[^A-Z]*","", toupper(method)), c("LRT", "QL"))
    pb = se_make_pseudobulk(se, sample.col)
    X = SummarizedExperiment::assays(pb)$counts
    cd = as.data.frame(SummarizedExperiment::colData(pb))
    cd[[pathology]] = relevel(as.factor(cd[[pathology]]), ref=control)
    covariates = covariates[covariates %in% colnames(cd)]
    dgel = edgeR::DGEList(X, group=cd[[pathology]], remove.zeros=TRUE)
    dgel = edgeR::calcNormFactors(dgel, method="TMM")
    covariates = covariates[covariates %in% colnames(deg.filter.design(cd[c(covariates)]))]
    design = model.matrix(as.formula(paste0(c("~0", pathology, covariates), collapse=" + ")),
                          data=cd)
    design = deg.filter.design(design)

    contrasts = limma::makeContrasts(contrasts=paste0(make.names(paste0(pathology, case)),
                                                      "-",
                                                      make.names(paste0(pathology, control))),
                                     levels=design)
    if (method == "QL") {
        dgel = edgeR::estimateDisp(dgel, design, robust=TRUE)
        fit = edgeR::glmQLFit(dgel, design, robust=TRUE)
        res = edgeR::glmQLFTest(fit, contrast=contrasts)
    } else if (method == "LRT") {
        dgel = edgeR::estimateGLMCommonDisp(dgel, design)
        dgel = edgeR::estimateGLMTagwiseDisp(dgel, design)
        fit = edgeR::glmFit(dgel, design, robust=TRUE)
        res = edgeR::glmLRT(fit, contrast=contrasts)
    }
    tbl = edgeR::topTags(res, n=Inf, sort.by="none")$table
    colnames(tbl) = gsub("^logFC$", "log2FC", colnames(tbl))
    rd = as.data.frame(SummarizedExperiment::rowData(se))
    for (suffix in intersect(c("log2FC", "FDR", "F", "LR"), colnames(tbl))) {
        rd[[paste0(prefix, "_", method, "_", suffix)]] = NA
        rd[rownames(tbl), paste0(prefix, "_", method, "_", suffix)] = tbl[[suffix]]
    }
    SummarizedExperiment::rowData(se) = rd
    return(se)
}

#' Add RUVseq columns to SummarizedExperiment
#' Based on Carles Boix's RUVSeq workflow
#'
#' @param sce SummarizedExperiment object
#' @param sample.col Sample on which to pseudobulk
#' @param pathology Main variable of interest
#' @param NRUV Number of RUV components to generate. Non-variable RUV components are removed.
#' @param norm edgeR norm method for cal`<cNormFactors
#' @export
deg.ruvseq <- function(sce, sample.col, pathology, covariates=NULL, NRUV=3, norm="TMM", verbose=TRUE) {
    pb = se_make_pseudobulk(sce, sample.col)
    cd = SummarizedExperiment::colData(pb)
    X = SummarizedExperiment::assays(pb)$counts
    dgel = edgeR::DGEList(X, group=cd[[pathology]], remove.zeros=TRUE)
    covariates = covariates[covariates %in% names(SummarizedExperiment::colData(pb))]
    if (verbose) {
        cat("Design matrix:\n")
        print(tibble::as_tibble(cd[c(pathology, covariates)]), n=nrow(cd))
        print(str(cd))
    }
    covariates = covariates[covariates %in% colnames(deg.filter.design(cd[c(covariates)]))]
    design = model.matrix(as.formula(paste0("~", paste0(c(pathology, covariates), collapse="+"))), data=cd)
    design = deg.filter.design(design)
### Use LRT workflow
    dgel = edgeR::calcNormFactors(dgel, method=norm)
    dgel = edgeR::estimateGLMCommonDisp(dgel, design)
    dgel = edgeR::estimateGLMTagwiseDisp(dgel, design)
    fit1 = edgeR::glmFit(dgel, design)
    res1 = RUVSeq:::residuals.DGEGLM(fit1, type="deviance")
    ruv = RUVSeq::RUVr(dgel$counts,
                       cIdx=as.character(rownames(dgel$counts)),
                       k=NRUV,
                       residuals=res1)
    ruvw = ruv$W
    ruvw = ruvw[,apply(ruvw, 2, sd) > 0]
    colnames(ruvw) = paste0("RUV", gsub("^W", "", colnames(ruvw)))
    rownames(ruvw) = colnames(ruv$normalizedCounts)
    S = as.character(SummarizedExperiment::colData(sce)[[sample.col]])
    for (cn in colnames(ruvw)) {
        SummarizedExperiment::colData(sce)[[cn]] = ruvw[S, cn]
    }
    return(sce)
}

#' Run Nebula on a SummarizedExperiment
#'
#' @param sce SummarizedExperiment object
#' @param sample Sample which to model dispersions
#' @param pathology Pathology column to study
#' @param case Case to look at within pathology column
#' @param control Control value to look at within pathology column
#' @param covariates Vector of covariate columns within colData(sce)
#' @param offset Offset from which to weight cells in model
#' @param assay Assay within sce. Default is "counts"
#' @param model Model to pass to nebula::nebula
#' @param filter_only_case_control logical telling whether to filter sce to only case and control values within pathology
#' @param NRUV Number of RUVSeq::RUVr components. 0 means no RUVr components
#' @param cpm.cutoff Cutoff of CPM to utilize for RUV
#' @param cpm.count Number of observations passing CPM cutoff for filter for RUV
#' @export
deg.nebula <- function(sce, pathology, case, control, sample.col, covariates=NULL,
                       offset="total_counts", model="NBGMM",
                       ncore=getOption("mc.cores", 2),
                       cpc=0, reml=1,
                       prefix="nebula") {
    cd = SummarizedExperiment::colData(sce)
    if (!is.null(offset)) {
        stopifnot(offset %in% colnames(cd))
        offset = cd[[offset]]
    }
    X = SummarizedExperiment::assays(sce)$counts
    cd[[pathology]] = as.factor(as.character(cd[[pathology]]))
    cd[[pathology]] = relevel(cd[[pathology]], control)
    covariates = covariates[covariates %in% colnames(deg.filter.design(cd[c(covariates)]))]
    design = model.matrix(as.formula(paste0("~", paste0(c(pathology, covariates), collapse="+"))), data=cd)
    design = deg.filter.design(design)
    neb = nebula::nebula(X,
                         cd[[sample.col]],
                         design,
                         offset=offset,
                         cpc=cpc,
                         model=model, reml=ifelse(model=="NBLMM", reml, 0),
                         ncore=ncore)
    ### Now to parse output
    name.case = paste0(pathology, case)
    name.control = paste0(pathology, control)
    if (paste0("p_", name.case) %in% colnames(neb$summary)) {
        neb$summary$FDR = p.adjust(neb$summary[[paste0("p_", name.case)]], "fdr")
        neb$summary$log2FC = neb$summary[[paste0("logFC_", name.case)]] / log(2)
    } else if (paste0("p_", name.control) %in% colnames(neb$summary)) {
        neb$summary$FDR = p.adjust(neb$summary[[paste0("p_", name.control)]], "fdr")
        neb$summary$log2FC = -1 * neb$summary[[paste0("logFC_", name.control)]] / log(2)
    }
    neb_df = neb$summary
    ovr = neb$overdispersion
    colnames(ovr) = paste0("overdispersion_", colnames(ovr))
    print(str(ovr))
    neb_df = cbind(neb_df, ovr)
    neb_df$convergence = neb$convergence
    neb_df$algorithm = neb$algorithm
    rownames(neb_df) = neb_df$gene
    rd = as.data.frame(SummarizedExperiment::rowData(sce))
    for (suffix in intersect(c("log2FC", "FDR", "convergence", "algorithm", colnames(ovr)), colnames(neb_df))) {
        rd[[paste0(prefix, "_", suffix)]] = NA
        rd[rownames(neb_df), paste0(prefix, "_", suffix)] = neb_df[[suffix]]
    }
    SummarizedExperiment::rowData(sce) = rd
    return(sce)
}

#' Run DESeq2 on pseudobulked data
#'
#' @param se SummarizedExperiment data, expected pseudobulk
#' @param pathology Column in colData(se) from which to study differential changes
#' @param case Case value within pathology column
#' @param control Control value within pathology column
#' @param sample.col Sample column to pseudobulk on
#' @param covariates Covariates to pass to DESeq2
#' @param prefix Prefix to add in rowData(se) where edgeR info is added
#' @export
deg.deseq2 <- function(se,
                       pathology,
                       case,
                       control,
                       sample.col,
                       covariates=NULL,
                       prefix="DESeq2") {
    pb = se_make_pseudobulk(se, sample.col)
    X = SummarizedExperiment::assays(pb)$counts
    cd = as.data.frame(SummarizedExperiment::colData(pb))
    cd[[pathology]] = relevel(as.factor(cd[[pathology]]), ref=control)
    covariates = covariates[covariates %in% colnames(cd)]
    covariates = covariates[covariates %in% colnames(deg.filter.design(cd[c(covariates)]))]
    design = model.matrix(as.formula(paste0("~", c(pathology, covariates), collapse=" + ")),
                          data=cd)
    design = deg.filter.design(design)
    dds = DESeq2::DESeqDataSetFromMatrix(
                      X,
                      colData=cd,
                      design=design)
    out = DESeq2::DESeq(dds)
    df = DESeq2::results(out, list(make.names(paste0(pathology, case))), independentFiltering=FALSE)
    rd = as.data.frame(SummarizedExperiment::rowData(se))
    rd[[paste0(prefix, "_log2FC")]] = NA
    rd[rownames(df), paste0(prefix, "_log2FC")] = df$log2FoldChange
    rd[[paste0(prefix, "_FDR")]] = NA
    rd[rownames(df), paste0(prefix, "_FDR")] = df$padj
    rd[[paste0(prefix, "_stat")]] = NA
    rd[rownames(df), paste0(prefix, "_stat")] = df$stat
    SummarizedExperiment::rowData(se) = rd
    return(se)
}

deg.lmer <- function(se, sample.col, pathology, case, control, covariates=NULL,
                     weights="log1p_total_counts", mc.cores=2,
                     prefix="lmer") {
    stop("not implemented yet")
    M = SummarizedExperiment::assays(se)$counts
    cd = as.data.frame(SummarizedExperiment::colData(se))
    cd[[pathology]] = as.factor(as.character(cd[[pathology]]))
    cd[[pathology]] = relevel(cd[[pathology]], control)
    covariates = covariates[covariates %in% colnames(deg.filter.design(cd[c(covariates)]))]
    design = model.matrix(as.formula(paste0(c("~0", pathology, covariates), collapse="+")), data=cd)
    design = deg.filter.design(design, rename=FALSE)
    covariates = deg.extract.covariates.from.design(design, covariates)
    if (is.null(cd[[weights]])) {
        weights = rep(1, nrow(cd))
    } else {
        weights = cd[[weights]]
    }
    result = lapply(rownames(se), function(gene) {
        cat(gene,"\n")
        print(table(se@colData$Pathology))
        lme4::glmer(as.formula(paste0(pathology, "~ gene + (1|", sample.col, ") + ",
                                              paste0(covariates, collapse="+"))),
                   data=cbind(cd, data.frame(gene=M[gene,])),
                   weights=weights,
                   family=stats::binomial(link="logit"))

    })
    coefs = t(sapply(setNames(result, rownames(se)), function(model) {
        coef(summary(model))["gene",]
    }))
    return(data.frame(coef=coefs[,"Estimate"],
                      pvalue=coefs[,"Pr(>|z|)"],
                      row.names=rownames(coefs)))

}

#' Run MAST on a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object
#' @param sample Sample which to model dispersions
#' @param pathology Pathology column to study
#' @param case Case to look at within pathology column
#' @param control Control value to look at within pathology column
#' @param covariates Vector of covariate columns within colData(sce)
#' @export
deg.mast <- function(sce, sample.col, pathology, case, control, covariates=NULL,
                     method, prefix="MAST") {
    require(MAST)
    M = SummarizedExperiment::assays(sce)$counts
    M = M %*% Matrix::Diagonal(x=1e4 / Matrix::colSums(M))
    M = log1p(M) / log(2)
    cd = as.data.frame(SummarizedExperiment::colData(sce))
    cd[[pathology]] = as.factor(as.character(cd[[pathology]]))
    cd[[pathology]] = relevel(cd[[pathology]], control)
    sca = MAST::FromMatrix(as.matrix(M), cd, SummarizedExperiment::rowData(sce), check_sanity=TRUE)
    remove("M")
    covariates = covariates[covariates %in% colnames(deg.filter.design(cd[c(covariates)]))]
    ### Use design to find what WOULD have been filtered if ZLM allowed a design matrix
    design = model.matrix(as.formula(paste0(c("~0", pathology, covariates), collapse="+")), data=cd)
    design = deg.filter.design(design, rename=FALSE)
    covariates = deg.extract.covariates.from.design(design, covariates)
    if (method=="RE") {
        cat("Running MAST+Random Effect\n")
        re.cov = paste0("(1|", sample.col, ")")
        zlmCond = MAST::zlm(as.formula(paste0("~", paste0(c(pathology, covariates, re.cov), collapse="+"))), sca, method="glmer", ebayes=FALSE, strictConvergence=FALSE)
    } else {
        cat("Running MAST Hurdle model\n")
        zlmCond = MAST::zlm(as.formula(paste0("~", paste0(c(pathology, covariates), collapse="+"))), sca)
    }
    ## saveRDS(zlmCond, "zlmCond.rds")
    summaryCond = MAST::summary(zlmCond, doLRT=paste0(pathology, case))
    ## saveRDS(summaryCond, "summaryCond.rds")
    summaryDt = summaryCond$datatable
    fcHurdle = merge(summaryDt[summaryDt$contrast==paste0(pathology, case) & summaryDt$component=='H',c("primerid", "Pr(>Chisq)")], #hurdle P values
                     summaryDt[summaryDt$contrast==paste0(pathology, case) & summaryDt$component=='logFC', c("primerid", "coef", "ci.hi", "ci.lo")], by='primerid') #logFC coefficients
    fcHurdle$FDR = p.adjust(fcHurdle$`Pr(>Chisq)`, "fdr")
    colnames(fcHurdle) = gsub("^coef$", "log2FC", colnames(fcHurdle))
    colnames(fcHurdle) = gsub("^primerid$", "gene", colnames(fcHurdle))
    rownames(fcHurdle) = fcHurdle$gene
    rd = as.data.frame(SummarizedExperiment::rowData(sce))
    for (suffix in intersect(c("log2FC", "FDR"), colnames(fcHurdle))) {
        rd[[paste0(prefix, "_", method, "_", suffix)]] = NA
        rd[rownames(fcHurdle), paste0(prefix, "_", method, "_", suffix)] = fcHurdle[[suffix]]
    }
    SummarizedExperiment::rowData(sce) = rd
    return(sce)
}
