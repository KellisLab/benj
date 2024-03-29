

# Color schemes for cell types
#' @export
colors = list(celltypes=c("Ast"="#D6242A",
                          "Exc"="#45A23A",
                          "Inh"="#9ED97A",
                          "Mic"="#BE9FCE",
                          "Oli"="#F4A85E",
                          "Opc"="#A14E2D",
                          "Msn"="#F59695",
                          "Vas"="#FED9A6"),
              regions=c("AG"="#8DD3C7",
                        "MT"="#80B1D3",
                        "NAc"="#32768C",
                        "MB"="#FAC0DF",
                        "PFC"="#FDB462",
                        "EC"="#FFED6F",
                        "HC"="#BEBADA",
                        "DS"="",
                        "VS"="#F59695",
                        "TH"="#B3DE69"))


#' Get named color categories, subsetted by 'select'.
#'
#' @param name Name of color palette
#' @param select Subset of palette to select. Use names to rename
#' @param ordered Use specific order
get.colors <- function(name, select, ordered=FALSE) {

}
