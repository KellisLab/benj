
#' Get SRA ID and GEO accession from GEO id
#'
#' @export
sra_from_geo <- function(geo) {
    L = unlist(lapply(geo_list, GEOquery::getGEO))
    pf = dplyr::bind_rows(lapply(L, Biobase::pData))
    pf = pf[c("geo_accession", "relation", "relation.1")]
    pf$sra = ifelse(1:nrow(pf) %in% grep("^SRA: https://www.ncbi.nlm.nih.gov/sra[?]term=", pf$relation), pf$relation, pf$relation.1)
    pf$sra = gsub("^SRA: https://www.ncbi.nlm.nih.gov/sra[?]term=", "", pf$sra)
    return(data.frame(sra=pf$sra, geo_accession=pf$geo_accession))
}
