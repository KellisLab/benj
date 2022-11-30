
#' Set default heatmap options, modified from Carles Boix
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


htSortMatrix <- function(M, method="euclidean", ratio=0.5, cutoff=cutoff) {
    M = order.tsp(M, rows=TRUE, method=method)
    M = order.tsp(M, rows=FALSE, method=method)
    M = diag.mat3(M, rows=TRUE, ratio=ratio, cutoff=cutoff)
    M = diag.mat3(M, rows=FALSE, ratio=ratio, cutoff=cutoff)
    return(M)
}
#' Automatic heatmap
#'
#' @param M matrix to pass
#' @param ux millimeters per grid
#' @export
autoHeatmap <- function(M, ux=1.5, sort=TRUE, method="euclidean", dimname_fontsize=3.5, ratio=0.5, cutoff=0.25, ...) {
    if (sort) {
        M = htSortMatrix(M, method=method, ratio=ratio, cutoff=cutoff)
    }
    return(ComplexHeatmap::Heatmap(
        M, cluster_rows=FALSE, cluster_columns=FALSE,
        width = ncol(M)*grid::unit(ux, "mm"),
        height = nrow(M)*grid::unit(ux, "mm"),
        row_names_gp=grid::gpar(fontsize=dimname_fontsize),
        column_names_gp=grid::gpar(fontsize=dimname_fontsize),
        ...
    ))
}
#' Save ComplexHeatmap, modified from Carles Boix
#'
#' @param ht ComplexHeatmap heatmap
#' @param pltprefix Prefix for images
#' @param w Width in inches
#' @param h Height in inches
#' @param dpi DPI of PNG
#' @export
saveHeatmap <- function(ht, pltprefix, w, h, dpi=600, extra=NULL, ...) {
    pdf(paste0(pltprefix, ".pdf"), width=w, height=h)
    ComplexHeatmap::draw(ht, ht_gap=grid::unit(0.5, "mm"), ...)
    if (is.function(extra)) {
        extra()
    }
    dev.off()
    png(paste0(pltprefix, ".png"), res=dpi, width=w, height=h, units="in")
    ComplexHeatmap::draw(ht, ht_gap=grid::unit(0.5, "mm"), ...)
    if (is.function(extra)) {
        extra()
    }
    dev.off()
    print(pltprefix)
}

#' Save ggplot, modified from Carles Boix
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

#' Draw triangles in heatmap
#'
#' @export
ht_triangle_split <- function(mat.ul, mat.lr, col.ul, col.lr, lwd=0, ...) {
    return(function(j, i, x, y, width, height, fill) {
        ### start with UL, UR, LR, LL
        corners = data.frame(x=as.numeric(width) * c(-1/2, 1/2, 1/2, -1/2),
                             y=as.numeric(height) * c(-1/2, -1/2, 1/2, 1/2),
                             row.names=c("LL", "UL", "UR", "LR"))
        ### upper left
        grid::grid.polygon(x=as.numeric(x) + corners[c("UL", "UR", "LL", "UL"),"x"],
                           y=as.numeric(y) + corners[c("UL", "UR", "LL", "UL"),"y"],
                           gp=grid::gpar(fill=col.ul(mat.ul[i, j]), lwd=lwd, ...))
        grid::grid.polygon(x=as.numeric(x) + corners[c("UR", "LR", "LL", "UR"),"x"],
                           y=as.numeric(y) + corners[c("UR", "LR", "LL", "UR"),"y"],
                           gp=grid::gpar(fill=col.lr(mat.lr[i, j]), lwd=lwd, ...))
    })
}

#ht_dotplot <- function(pct,
volcano <- function(df, min.FDR=0.99999) {

}
## manhattan <- function(gr, chrom.sizes) {

## }

#' Get ggrepel locations
#' @export
repel_text <- function(x, y, text, cex=1, ord=NULL, max.overlaps=NULL, font=NULL, vfont=NULL, q=0.5, xlim=NULL, ylim=NULL, hjust=NULL, vjust=NULL) {
    df = data.frame(point.x=x, point.y=y, text=text)
    if (!is.null(ord)) {
        df = df[order(ord),]
    }
    bb_w = graphics::strwidth(text, cex=cex, font=font, vfont=vfont)
    bb_h = graphics::strheight(text, cex=cex, font=font, vfont=vfont)
    bb_px = .1 * quantile(bb_h, q)
    bb_py = .1 * quantile(bb_w, q)
    rb = ggrepel:::repel_boxes2(data_points=data_points,
                                point_size=point_size,
                                point_padding_x=point_padding_x,
                                point_padding_y=point_padding_y,
                                boxes=do.call(rbind, boxes),
                                xlim=range(xlim),
                                ylim=range(ylim),
                                hjust=hjust,
                                vjust=vjust,
                                max_overlaps=max_overlaps)


}

#' @export
colorRamp2 <- function(breaks, colors, transparency=0, space="LAB", hcl_palette=NULL, reverse=FALSE, transform=NULL, transform.length=100) {
    cr = circlize::colorRamp2(breaks, colors, transparency=transparency, space=space, hcl_palette=hcl_palette, reverse=reverse)
    if (is.function(transform)) {
        S = seq(min(breaks), max(breaks), length.out=transform.length)
        colors = cr(S)
        return(circlize::colorRamp2(transform(S), cr(S),
                                    transparency=transparency,
                                    space=space, hcl_palette=hcl_palette,
                                    reverse=reverse))
    } else {
        return(cr)
    }
}
