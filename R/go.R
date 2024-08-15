

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

#' @export
enrichALL <- function(gene, OrgDb, keyType="SYMBOL", minGSSize=10, maxGSSize=100, ...) {
  if (length(gene)==0) {
    return(data.frame())
  }
  go = enrichGO(gene, OrgDb, keyType=keyType, minGSSize=minGSSize, maxGSSize=maxGSSize, ...)
  kegg = enrichKEGG(gene, OrgDb, keyType=keyType, minGSSize=minGSSize, maxGSSize=maxGSSize, ...)
  wp = enrichWP(gene, OrgDb, keyType=keyType, minGSSize=minGSSize, maxGSSize=maxGSSize, ...)
  reac = enrichReactome(gene, OrgDb, keyType=keyType, minGSSize=minGSSize, maxGSSize=maxGSSize, ...)
  do = enrichDO(gene, OrgDb, keyType=keyType, minGSSize=minGSSize, maxGSSize=maxGSSize, ...)
  dplyr::bind_rows(list(go, kegg, wp, reac))
}

#' @export
enrichALL_Excel <- function(xlsx, new_sheet, OrgDb, method,
                            min_log2FC=0, max_FDR=0.05, minGSSize=10, maxGSSize=100, ...) {
    require(dplyr)
    df <- readxl::read_xlsx(xlsx)
    sheet_names <- readxl::excel_sheets(xlsx)
    df$log2FC <- df[[paste0(method, "_log2FC")]]
    df$FDR <- df[[paste0(method, "_FDR")]]
    universe <- df$gene
    df <- df %>% filter(abs(log2FC) >= min_log2FC & FDR < max_FDR)
    wb <- openxlsx::loadWorkbook(xlsx)
    if (!(paste0(new_sheet, " genes") %in% sheet_names)) {
      openxlsx::addWorksheet(wb, paste0(new_sheet, " genes"))
    }
    openxlsx::writeData(wb, paste0(new_sheet, " genes"), df)
    go_pos <- benj::enrichALL(df$gene[df$log2FC > 0], OrgDb=OrgDb, universe=universe, minGSSize=minGSSize, maxGSSize=maxGSSize)
    go_neg <- benj::enrichALL(df$gene[df$log2FC < 0], OrgDb=OrgDb, universe=universe, minGSSize=minGSSize, maxGSSize=maxGSSize)
    df <- dplyr::bind_rows(list(Up=go_pos, Down=go_neg))
    df$log2FC=min_log2FC
    df$FDR=max_FDR
    meth_path <- paste0(method, "_pathways")
    if ((meth_path %in% sheet_names)&&(nrow(df)>0)) {
      mpdf <- openxlsx::readWorkbook(wb, sheet=meth_path)
      I <- intersect(colnames(mpdf), colnames(df))
      if (length(I) == ncol(df)) {
        df <- rbind(mpdf[,I], df[,I])
      } else {
        warning(paste0("Columns do not match:", paste0(colnames(mpdf),collapse=",") ," vs", paste0(colnames(df), collapse=",")))
      }
    } else {
      openxlsx::addWorksheet(wb, meth_path)
    }
    if (nrow(df) > 0) {
      openxlsx::writeData(wb, meth_path, df)
    }
    openxlsx::saveWorkbook(wb, xlsx, overwrite=TRUE)
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

#' @export
enrichKEGG <- function(gene, OrgDb, keyType="SYMBOL", minGSSize=10, maxGSSize=100, universe=NULL, ...) {
  if (keyType != "kegg") {
    gene = AnnotationDbi::mapIds(OrgDb, keys=gene, keytype=keyType, column="ENTREZID")
    if (!is.null(universe)) {
      universe = Reduce(c, AnnotationDbi::mapIds(OrgDb, keys=universe, keytype=keyType, column="ENTREZID", multiVals=list))
    }
    keyType = "ncbi-geneid"
  }
  if (is.null(universe)) {
    res = clusterProfiler::enrichKEGG(gene, keyType=keyType, minGSSize=minGSSize, maxGSSize=maxGSSize, ...)
  } else {
    res = clusterProfiler::enrichKEGG(gene, keyType=keyType, minGSSize=minGSSize, maxGSSize=maxGSSize, universe=universe, ...)
  }
  xf = .proc.go.result(res)
  if (nrow(xf) > 0) {
      xf$ONTOLOGY = "KEGG"
      if (keyType == "ncbi-geneid") {
          mapping = setNames(names(gene), gene)
          xf$geneID = sapply(strsplit(xf$geneID, "/"), function(gid) {
              paste0(mapping[gid], collapse="/")
          })
      }
  }
  return(xf)
}

#' @export
enrichMKEGG <- function(gene, OrgDb, keyType="SYMBOL", minGSSize=10, maxGSSize=100, universe=NULL, ...) {
  if (keyType != "kegg") {
    gene = AnnotationDbi::mapIds(OrgDb, keys=gene, keytype=keyType, column="ENTREZID")
    if (!is.null(universe)) {
      universe = Reduce(c, AnnotationDbi::mapIds(OrgDb, keys=universe, keytype=keyType, column="ENTREZID", multiVals=list))
    }
    keyType = "ncbi-geneid"
  }
  print(str(gene))
  if (is.null(universe)) {
    res = clusterProfiler::enrichMKEGG(as.integer(gene), keyType=keyType, minGSSize=minGSSize, maxGSSize=maxGSSize, ...)
  } else {
    res = clusterProfiler::enrichMKEGG(gene, keyType=keyType, minGSSize=minGSSize, maxGSSize=maxGSSize, universe=universe, ...)
  }
  xf = .proc.go.result(res)
  if (nrow(xf) > 0) {
      xf$ONTOLOGY = "MKEGG"
      if (keyType == "ncbi-geneid") {
          mapping = setNames(names(gene), gene)
          xf$geneID = sapply(strsplit(xf$geneID, "/"), function(gid) {
              paste0(mapping[gid], collapse="/")
          })
      }
  }
  return(xf)
}

#' @export
enrichWP <- function(gene, OrgDb, keyType="SYMBOL", minGSSize=10, maxGSSize=100, universe=NULL, ...) {
  mf = S4Vectors::metadata(OrgDb)
  mf = setNames(mf$value, mf$name)
  organism = mf[["ORGANISM"]]
  if (keyType != "ENTREZID") {
      gene = AnnotationDbi::mapIds(OrgDb, keys=gene, keytype=keyType, column="ENTREZID")
      if (!is.null(universe)) {
        universe = Reduce(c, AnnotationDbi::mapIds(OrgDb, keys=universe, keytype=keyType, column="ENTREZID", multiVals=list))
      }
  }
  if (is.null(universe)) {
    res = clusterProfiler::enrichWP(gene, organism, minGSSize=minGSSize, maxGSSize=maxGSSize, ...)
  } else {
    res = clusterProfiler::enrichWP(gene, organism, minGSSize=minGSSize, maxGSSize=maxGSSize, universe=universe, ...)
  }
  xf = .proc.go.result(res)
  if (nrow(xf) > 0) {
      xf$ONTOLOGY = "WP"
      if (keyType != "ENTREZID") {
          mapping = setNames(names(gene), gene)
          xf$geneID = sapply(strsplit(xf$geneID, "/"), function(gid) {
              paste0(mapping[gid], collapse="/")
          })
      }
  }
  return(xf)
}

#' @export
enrichReactome <- function(gene, OrgDb, keyType="SYMBOL", minGSSize=10, maxGSSize=100, universe=NULL, ...) {
  mf = S4Vectors::metadata(OrgDb)
  mf = setNames(mf$value, mf$name)
  species = tolower(mf[["SPECIES"]])
  if (keyType != "ENTREZID") {
      gene = AnnotationDbi::mapIds(OrgDb, keys=gene, keytype=keyType, column="ENTREZID")
      if (!is.null(universe)) {
        universe = Reduce(c, AnnotationDbi::mapIds(OrgDb, keys=universe, keytype=keyType, column="ENTREZID", multiVals=list))
      }
  }
  if (is.null(universe)) {
    res = ReactomePA::enrichPathway(gene, species, minGSSize=minGSSize, maxGSSize=maxGSSize, ...)
  } else {
    res = ReactomePA::enrichPathway(gene, species, minGSSize=minGSSize, maxGSSize=maxGSSize, universe=universe, ...)
  }
  xf = .proc.go.result(res)
  if (nrow(xf) > 0) {
      xf$ONTOLOGY = "REACTOME"
      if (keyType != "ENTREZID") {
          mapping = setNames(names(gene), gene)
          xf$geneID = sapply(strsplit(xf$geneID, "/"), function(gid) {
              paste0(mapping[gid], collapse="/")
          })
      }
  }
  return(xf)
}

#' @export
enrichDO <- function(gene, OrgDb, keyType="SYMBOL", minGSSize=10, maxGSSize=100, universe=NULL, ...) {
  if (keyType != "ENTREZID") {
    gene = AnnotationDbi::mapIds(OrgDb, keys=gene, keytype=keyType, column="ENTREZID")
    if (!is.null(universe)) {
      universe = Reduce(c, AnnotationDbi::mapIds(OrgDb, keys=universe, keytype=keyType, column="ENTREZID", multiVals=list))
    }
  }
  if (is.null(universe)) {
    res = DOSE::enrichDO(gene, minGSSize=minGSSize, maxGSSize=maxGSSize, ...)
  } else {
    res = DOSE::enrichDO(gene, minGSSize=minGSSize, maxGSSize=maxGSSize, universe=universe, ...)
  }
    xf = .proc.go.result(res)
  if (nrow(xf) > 0) {
      xf$ONTOLOGY = "DOSE"
      if (keyType != "ENTREZID") {
          mapping = setNames(names(gene), gene)
          xf$geneID = sapply(strsplit(xf$geneID, "/"), function(gid) {
              paste0(mapping[gid], collapse="/")
          })
      }
  }
  return(xf)
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

gseKEGG <- function(geneList, ranks, OrgDb, filter=NULL, keyType="SYMBOL", minGSSize=10, maxGSSize=100, universe=NULL, ...) {
    if (keyType != "kegg") {
    gene = AnnotationDbi::mapIds(OrgDb, keys=gene, keytype=keyType, column="ENTREZID")
    if (!is.null(universe)) {
      universe = Reduce(c, AnnotationDbi::mapIds(OrgDb, keys=universe, keytype=keyType, column="ENTREZID", multiVals=list))
    }
    keyType = "ncbi-geneid"
  }

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
