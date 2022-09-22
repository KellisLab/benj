
#' Set default heatmap options
#'
#' @param htsc Padding ratio
#' @export
set_ht_opt <- function(htsc=2.5) {
    ComplexHeatmap::ht_opt(
        heatmap_row_names_gp = grid::gpar(fontsize = 5),
        heatmap_column_names_gp = grid::gpar(fontsize = 5),
        heatmap_row_title_gp = grid::gpar(fontsize = 5.5, fontface="bold"),
        heatmap_column_title_gp = grid::gpar(fontsize = 5.5, fontface="bold"),
        legend_title_gp = grid::gpar(fontsize = 5.5, font=2),
        legend_labels_gp = grid::gpar(fontsize = 5),
        legend_grid_height = grid::unit(2.5, 'mm'),
        legend_grid_width = grid::unit(2.5, 'mm'),
        DENDROGRAM_PADDING = grid::unit(.5 / htsc, 'mm'),
        DIMNAME_PADDING = grid::unit(1 / htsc, 'mm'),
        COLUMN_ANNO_PADDING = grid::unit(1 / htsc, 'mm'),
        ROW_ANNO_PADDING = grid::unit(1 / htsc, 'mm'),
        HEATMAP_LEGEND_PADDING = grid::unit(2 / htsc, 'mm'),
        ANNOTATION_LEGEND_PADDING = grid::unit(2 / htsc, 'mm'))
}

#' Save ComplexHeatmap
#'
#' @param ht ComplexHeatmap heatmap
#' @param pltprefix Prefix for images
#' @param w Width in inches
#' @param h Height in inches
#' @param dpi DPI of PNG
#' @export
saveHeatmap <- function(ht, pltprefix, w, h, dpi=600, ...) {
    pdf(paste0(pltprefix, ".pdf"), width=w, height=h)
    ComplexHeatmap::draw(ht, ht_gap=grid::unit(0.5, "mm"), ...)
    dev.off()
    png(paste0(pltprefix, ".png"), res=dpi, width=w, height=h, units="in")
    ComplexHeatmap::draw(ht, ht_gap=grid::unit(0.5, "mm"), ...)
    dev.off()
    print(pltprefix)
}

#' Save ggplot
#'
#' @param gp GGplot object
#' @param pltprefix Prefix for images
#' @param w Width in inches
#' @param h Height in inches
#' @param dpi DPI of PNG
#' @export
saveGGplot <- function(gp, pltprefix, w, h, dpi=600) {
    ggplot2::ggsave(paste0(pltprefix, ".pdf"), gp, units="in", dpi=dpi, width=w, height=h)
    ggplot2::ggsave(paste0(pltprefix, ".png"), gp, units="in", dpi=dpi, width=w, height=h)
    print(pltprefix)
}
