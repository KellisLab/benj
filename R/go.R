

.proc.go.result <- function(obj) {
    if (is.null(obj)) {
        return(data.frame())
    }
    if ("result" %in% slotNames(obj)) {
        df = slot(obj, "result")
        if (nrow(df) == 0) {
            return(data.frame())
        } else {
            return(df)
        }
    } else {
        return(data.frame())
    }
}
#' Run clusterProfiler::enrichGO with better options and defaults
#'
#' @param gene List of genes
#' @param OrgDb OrgDb e.g. org.Hs.eg.db::org.Hs.eg.db or org.Mm.eg.db::org.Mm.eg.db
#' @param keyType Set default keyType to SYMBOL instead of ENTREZID for clusterProfiler::enrichGO
#' @param ont Set default ontology to ALL for clusterProfiler::enrichGO
#' @param minGSSize Minimum gene set size for a GO term
#' @param maxGSSize Maximum gene set size for a GO term
#' @export
enrichGO <- function(gene, OrgDb, keyType="SYMBOL", ont="ALL", minGSSize=10, maxGSSize=100, ...) {
    ego = clusterProfiler::enrichGO(gene=gene, OrgDb=OrgDb, keyType=keyType, ont=ont, minGSSize=minGSSize, maxGSSize=maxGSSize, ...)
    return(.proc.go.result(ego))
}


#' Run clusterProfiler::gseGO with better options and defaults
#'
#' @param geneList List of genes
#' @param ranks Ranked vector for GSEA
#' @param OrgDb OrgDb e.g. org.Hs.eg.db::org.Hs.eg.db or org.Mm.eg.db::org.Mm.eg.db
#' @param filter Optional logical filter to subset genes (for ease of use)
#' @param keyType Set default keyType to SYMBOL instead of ENTREZID for clusterProfiler::enrichGO
#' @param ont Set default ontology to ALL for clusterProfiler::enrichGO
#' @param minGSSize Minimum gene set size for a GO term
#' @param maxGSSize Maximum gene set size for a GO term
#' @export
gseGO <- function(geneList, ranks, OrgDb, filter=NULL, keyType="SYMBOL", ont="ALL", minGSSize=10, maxGSSize=100, ...) {
    if (!is.null(filter)) {
        geneList = setNames(ranks[filter], geneList[filter])
    } else {
        geneList = setNames(ranks, geneList)
    }
    if (length(geneList) < 2) {
        ego = NULL
    } else {
        ego = clusterProfiler::gseGO(geneList=geneList,
                                 OrgDb=OrgDb,
                                 keyType=keyType,
                                 ont=ont,
                                 minGSSize=minGSSize,
                                 maxGSSize=maxGSSize,
                                 ...)
    }
    return(.proc.go.result(ego))
}

goMatrix <- function() {
    L = as.list(org.Hs.eg.db::org.Hs.egGO2ALLEG)
    X = lapply(L, function(gl) {
        AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys=gl, column="SYMBOL", keytype="ENTREZID", multiVals="first")
    })
    allg = Reduce(union, X)
}

keggMatrix_fromGenes <- function(geneList, keytype="SYMBOL") {
    keggTerms = as.list(KEGG.db::KEGGPATHID2EXTID)
}
goMatrix_fromGenes <- function(geneList, keytype="SYMBOL") {

}
