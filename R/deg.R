
#' Run edgeR on a SummarizedExperiment
#'
#' @param se SummarizedExperiment to run
#' @param pathology Column in colData(se) that contains case and control
#' @param case Case string within pathology column. Example "AD"
#' @param control Control string within pathology column. Example "control"
#' @param covariates Vector of covariates to add to design matrix
#' @param method Method used for edgeR. Either LRT or QL
#' @param assay Assay from se to utilize
#' @param cpm.cutoff Cutoff of CPM to utilize
#' @param cpm.count Number of observations passing CPM cutoff for filter
#' @export
deg.edger <- function(se, pathology, case, control, covariates,
                      method="LRT", assay=NULL, cpm.cutoff=1, cpm.count=3) {
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(se)[[1]]
    } else {
        X = SummarizedExperiment::assays(se)[[assay]]
    }
    cd = SummarizedExperiment::colData(se)
    if (!is.character(cd[[pathology]]) | !is.factor(cd[[pathology]])) {
        warning(paste0("Converting ", pathology, " to factor"))
    }
    cd[[pathology]] = as.factor(as.character(cd[[pathology]]))
    dgel = edgeR::DGEList(X, group=cd[[pathology]], remove.zeros=TRUE)
    to_keep = rowSums(edgeR::cpm(dgel) > cpm.cutoff) >= cpm.count
    dgel = dgel[to_keep,,keep.lib.sizes=FALSE]
    dgel = edgeR::calcNormFactors(dgel, method="TMM")
    design = model.matrix(as.formula(paste0("~0 + ", c(pathology, covariates), collapse=" + ")),
                          data=cd)
    colnames(design) = make.names(colnames(design)) ### need spaces to be removed
    contrasts = limma::makeContrasts(contrasts=paste0(make.names(paste0(pathology, case)),
                                                      "-",
                                                      make.names(paste0(pathology, control))),
                                     levels=design)
    if (!(method %in% c("QL", "LRT"))) {
        method = "QL"
        warning(paste0("Setting method to ", method))
    }
    if (method == "QL") {
        dgel = edgeR::estimateDisp(dgel, design, robust=TRUE)
        fit <- edgeR::glmQLFit(dgel, design, robust=TRUE)
        res <- edgeR::glmQLFTest(fit, contrast=contrasts)
    } else if (method == "LRT") {
        dgel = edgeR::estimateGLMCommonDisp(dgel, design)
        dgel = edgeR::estimateGLMTagwiseDisp(dgel, design)
        fit = edgeR::glmFit(dgel, design, robust=TRUE)
        res = edgeR::glmLRT(fit, contrast=contrasts)
    }
    tbl = edgeR::topTags(res, n=Inf, sort.by="none")$table
    colnames(tbl) = gsub("^logFC$", "log2FC", colnames(tbl))
    tbl$gene = rownames(tbl)
    tbl$method = method
    return(tbl[order(tbl$FDR),])
}

#' Add RUVseq columns to SummarizedExperiment
#'
#' @param sce SummarizedExperiment object
#' @param sample Sample on which to pseudobulk
#' @param pathology Main variable of interest
#' @param NRUV Number of RUV components to generate. Non-variable RUV components are removed.
#' @param assay Assay from SummarizedExperiment to use
#' @param cpm.cutoff Cutoff of CPM to utilize
#' @param cpm.count Number of observations passing CPM cutoff for filter
#' @param norm edgeR norm method for calcNormFactors
#' @export
ruvseq <- function(sce, sample, pathology, NRUV=10, assay=NULL, cpm.cutoff=1, cpm.count=3, norm="TMM") {
    pb = se_make_pseudobulk(sce, sample)
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(pb)[[1]]
    } else {
        X = SummarizedExperiment::assays(pb)[[assay]]
    }
    cd = SummarizedExperiment::colData(pb)
    dgel = edgeR::DGEList(X, remove.zeros=TRUE, group=cd[[pathology]])
    to_keep = rowSums(edgeR::cpm(dgel) > cpm.cutoff) >= cpm.count
    dgel = dgel[to_keep,,keep.lib.sizes=FALSE]
    design = model.matrix(~group, dgel$samples)
### Use LRT workflow
    dgel = edgeR::calcNormFactors(dgel, method=norm)
    dgel = edgeR::estimateGLMCommonDisp(dgel, design)
    dgel = edgeR::estimateGLMTagwiseDisp(dgel, design)
    fit1 = edgeR::glmFit(dgel, design)
    res1 = residuals(fit1, type="deviance")
    ruv = RUVSeq::RUVr(dgel$counts,
                        cIdx=as.character(rownames(dgel$counts)),
                        k=NRUV, residuals=res1)
    ruvw = ruv$W
    ruvw = ruvw[,apply(ruvw, 2, sd) > 0]
    colnames(ruvw) = paste0("RUV", colnames(ruvw))
    rownames(ruvw) = colnames(ruv$normalizedCounts)
    S = as.character(SummarizedExperiment::colData(sce)[[sample]])
    for (cn in colnames(ruvw)) {
        SummarizedExperiment::colData(sce)[[cn]] = ruvw[S, cn]
    }
    return(sce)
}
