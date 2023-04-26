
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
#' @param cpm.frac Fraction of samples passing CPM cutoff for filter
#' @export
deg.edger <- function(se, pathology, case, control, covariates=c(),
                      method="LRT", assay=NULL, cpm.cutoff=10, cpm.frac=0.25,
                      filter_only_case_control=FALSE) {
    if (filter_only_case_control) {
        se = se[,SummarizedExperiment::colData(se)[[pathology]] %in% c(case, control)]
    }

    if (is.null(assay)) {
        X = SummarizedExperiment::assays(se)$counts
    } else {
        X = SummarizedExperiment::assays(se)[[assay]]
    }
    if ("matrix" %in% class(X)) {
        round_diff = all(0 == zapsmall(abs(round(X) - X)))
    } else {
        round_diff = all(0 == zapsmall(abs(round(X@x) - X@x)))
    }
    if (!round_diff) {
        warning(paste0("Counts may not be integer: The difference between assay and rounded assay is ", round_diff))
    }
    if (S4Vectors::ncol(X) > 1000) {
        warning("More than 1000 samples are being called to a function under the assumption data passed in is pseudobulk")
    }

    cd = SummarizedExperiment::colData(se)
    if (!is.character(cd[[pathology]]) | !is.factor(cd[[pathology]])) {
        warning(paste0("Converting ", pathology, " to factor"))
    }
    cd[[pathology]] = as.factor(as.character(cd[[pathology]]))
    dgel = edgeR::DGEList(X, group=cd[[pathology]], remove.zeros=TRUE)
    to_keep = Matrix::rowMeans(edgeR::cpm(dgel) > cpm.cutoff) >= cpm.frac
    dgel = dgel[to_keep,,keep.lib.sizes=FALSE]
    dgel = edgeR::calcNormFactors(dgel, method="TMM")
    design = model.matrix(as.formula(paste0("~0 + ", paste0(c(pathology, covariates), collapse=" + "))),
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
    tbl$case = case
    tbl$control = control
    tbl$logCPM = edgeR::cpm(Matrix::rowSums(X), log=TRUE)[tbl$gene,]
    return(tbl[order(tbl$FDR),])
}

#' Add RUVseq columns to SummarizedExperiment
#' Based on Carles Boix's RUVSeq workflow
#'
#' @param sce SummarizedExperiment object
#' @param sample Sample on which to pseudobulk
#' @param pathology Main variable of interest
#' @param NRUV Number of RUV components to generate. Non-variable RUV components are removed.
#' @param assay Assay from SummarizedExperiment to use
#' @param cpm.cutoff Cutoff of CPM to utilize
#' @param cpm.frac Number of observations passing CPM cutoff for filter
#' @param norm edgeR norm method for calcNormFactors
#' @export
deg.ruvseq <- function(sce, sample, pathology, covariates=NULL, NRUV=3, assay=NULL, cpm.cutoff=10, cpm.frac=0.25, norm="TMM", filter_only_case_control=TRUE) {
    require(dplyr)
    pb = se_make_pseudobulk(sce, sample)
    if (filter_only_case_control) {
        pb = pb[,SummarizedExperiment::colData(pb)[[pathology]] %in% c(case, control)]
    }
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(pb)$counts
    } else {
        X = SummarizedExperiment::assays(pb)[[assay]]
    }
        if ("matrix" %in% class(X)) {
        round_diff = all(0 == zapsmall(abs(round(X) - X)))
    } else {
        round_diff = all(0 == zapsmall(abs(round(X@x) - X@x)))
    }
    if (!round_diff) {
        warning(paste0("Counts may not be integer: The difference between assay and rounded assay is ", round_diff))
    }
    cd = SummarizedExperiment::colData(pb)
    if (!is.character(cd[[pathology]]) | !is.factor(cd[[pathology]])) {
        warning(paste0("Converting ", pathology, " to factor"))
    }
    cd %>% mutate(across(where(is.factor), as.character)) -> cd
    cd[[pathology]] = as.factor(as.character(cd[[pathology]]))
    dgel = edgeR::DGEList(X, group=cd[[pathology]], remove.zeros=TRUE)
    to_keep = Matrix::rowMeans(edgeR::cpm(dgel) > cpm.cutoff) >= cpm.frac
    dgel = dgel[to_keep,,keep.lib.sizes=FALSE]
    covariates = covariates[covariates %in% colnames(cd)]
    design = model.matrix(as.formula(paste0("~", paste0(c(pathology, covariates), collapse="+"))), data=cd)
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
    colnames(ruvw) = paste0("RUV", colnames(ruvw))
    rownames(ruvw) = colnames(ruv$normalizedCounts)
    S = as.character(SummarizedExperiment::colData(sce)[[sample]])
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
deg.nebula <- function(sce, sample, pathology, case, control, covariates=c(),
                       offset="total_counts", assay=NULL, model="NBGMM",
                       filter_only_case_control=FALSE, factorize_pathology=TRUE,
                       ruv.remove.subjectlevel=TRUE,
                       cpc=0.01, NRUV=10, reml=1) {
    if (filter_only_case_control) {
        sce = sce[,SummarizedExperiment::colData(sce)[[pathology]] %in% c(case, control)]
    }
    if (NRUV > 0) {
        sce = deg.ruvseq(sce,
                         sample=sample,
                         pathology=pathology,
                         NRUV=NRUV,
                         covariates=covariates,
                         assay=assay)
        cd = SummarizedExperiment::colData(sce)
        cd.pb = SummarizedExperiment::colData(se_make_pseudobulk(sce, sample))
        if (ruv.remove.subjectlevel) {
            covariates = setdiff(covariates, colnames(cd.pb))
        }
        covariates = c(covariates, colnames(cd)[grep("^RUV", colnames(cd))])
    }
    sce[[sample]] = as.factor(as.character(sce[[sample]]))
    sce = sce[,order(sce[[sample]])]
    cd = SummarizedExperiment::colData(sce)
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(sce)$counts
    } else {
        X = SummarizedExperiment::assays(sce)[[assay]]
    }
    if ("matrix" %in% class(X)) {
        round_diff = all(0 == zapsmall(abs(round(X) - X)))
    } else {
        round_diff = all(0 == zapsmall(abs(round(X@x) - X@x)))
    }
    if (!round_diff) {
        warning(paste0("Counts may not be integer: The difference between assay and rounded assay is ", round_diff))
    }
    if (factorize_pathology) {
        cd[[pathology]] = as.factor(as.character(cd[[pathology]]))
    }
    design = model.matrix(as.formula(paste0("~",
                                            paste0(c(pathology, covariates), collapse="+"))),
                          data=cd)
    A = apply(design, 2, sd)
    design = design[,names(A)[(A>0) | (names(A) == '(Intercept)')]]
    if (is.null(offset) | !(offset %in% colnames(cd))) {
        warning("Offset not present. Using number of counts")
        offset = Matrix::colSums(X)
    } else {
        print(paste0("Using offset ", offset))
        offset = cd[[offset]]
    }
    neb = nebula::nebula(X,
                         cd[[sample]],
                         design,
                         offset=offset,
                         model=model, reml=reml, cpc=cpc)
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
    neb_df = cbind(neb_df, ovr)
    neb_df$convergence = neb$convergence
    neb_df$algorithm = neb$algorithm
    neb_df$case = case
    neb_df$control = control
    neb_df$logCPM = edgeR::cpm(Matrix::rowSums(X), log=TRUE)[neb_df$gene,]
    return(neb_df[order(neb_df$FDR),])
}

#' Run DESeq2 on pseudobulked data
#'
#' @param se SummarizedExperiment data, expected pseudobulk
#' @param pathology Column in colData(se) from which to study differential changes
#' @param case Case value within pathology column
#' @param control Control value within pathology column
#' @param covariates Covariates to pass to DESeq2
#' @param assay Assay to perform differential changes upon. Default is "counts"
#' @param filter_only_case_control logical telling whether to filter sce to only case and control values within pathology
#' @export
deg.deseq2 <- function(se,
                       pathology,
                       case,
                       control,
                       covariates=c(),
                       assay=NULL,
                       filter_only_case_control=FALSE) {
    if (filter_only_case_control) {
        se = se[,SummarizedExperiment::colData(se)[[pathology]] %in% c(case, control)]
    }
    if (is.null(assay)) {
        X = SummarizedExperiment::assays(se)$counts
    } else {
        X = SummarizedExperiment::assays(se)[[assay]]
    }
    if ("matrix" %in% class(X)) {
        round_diff = all(0 == zapsmall(abs(round(X) - X)))
    } else {
        round_diff = all(0 == zapsmall(abs(round(X@x) - X@x)))
    }
    if (!round_diff) {
        warning(paste0("Counts may not be integer: The difference between assay and rounded assay is ", round_diff))
    }
    if (S4Vectors::ncol(X) > 1000) {
        warning("More than 1000 samples are being called to a function under the assumption data passed in is pseudobulk")
    }
    formula = paste0("~", c(pathology, covariates), collapse=" + ")
    dds = DESeq2::DESeqDataSetFromMatrix(
                      X,
                      colData=SummarizedExperiment::colData(se),
                      design=as.formula(formula))
    out = DESeq2::DESeq(dds)
    df = DESeq2::results(out, contrast=c(pathology, case, control))
    df = df[!is.na(df$padj),]
    df = df[order(df$padj),]
    df$gene = rownames(df)
    colnames(df) = msub(c("^log2FoldChange$", "^padj$"),
                        c("log2FC", "FDR"), colnames(df))
    df$case = case
    df$control = control
    df$logCPM = edgeR::cpm(Matrix::rowSums(X), log=TRUE)[df$gene,]
    return(as.data.frame(df[order(df$FDR),]))
}
