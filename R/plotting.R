
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

savePlot <- function(p, pltprefix, w=7, h=7, dpi=600, ...){
  p <- substitute(p)
  pdf(paste0(pltprefix, ".pdf"), res=dpi, width=w, height=h, units="in")
  eval(p)
  dev.off()

  png(paste0(pltprefix, ".png"), res=dpi, width=w, height=h, units="in")
  eval(p)
  dev.off()
  print(pltprefix)
}
#' @export
htSortMatrix <- function(M, method="euclidean", ratio=0.5, cutoff=0.25, sort=c(1,2)) {
    if (is.logical(sort)) {
        if (sort) {
            sort = c(1,2)
        } else {
            sort = c()
        }
    }
    rows = 1 %in% sort
    columns = 2 %in% sort
    if (rows) {
        M = order.tsp(M, rows=TRUE, method=method)
    }
    if (columns) {
        M = order.tsp(M, rows=FALSE, method=method)
    }
    if (rows) {
        M = diag.mat3(M, rows=TRUE, ratio=ratio, cutoff=cutoff)
    }
    if (columns) {
        M = diag.mat3(M, rows=FALSE, ratio=ratio, cutoff=cutoff)
    }
    return(M)
}

#' Automatic heatmap
#'
#' @param M matrix to pass
#' @param ux millimeters per grid
#' @param sort Axes to sort upon. Be careful as the matrix is re-sorted before plotting with htSortMatrix
#' @export
autoHeatmap <- function(M, ux=1.5, sort=c(1, 2), method="euclidean",
                        dimname_fontsize=3.5, ratio=0.5, cutoff=0.25,
                        aspect_ratio=1, cell_text=NULL, cell_text_fontsize=10,
                        cluster_rows=FALSE, cluster_columns=FALSE, ...) {
    if (is.logical(cluster_rows)) {
      cluster_rows = cluster_rows || -1 %in% sort
    }
    if (is.logical(cluster_columns)) {
      cluster_columns = cluster_columns || -2 %in% sort
    }
    M = htSortMatrix(M, method=method, ratio=ratio, cutoff=cutoff, sort=sort)
    cell_fun=NULL
    if (!is.null(cell_text)) {
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf(cell_text, M[i, j]), x, y, gp = gpar(fontsize = cell_text_fontsize)) }
      return(ComplexHeatmap::Heatmap(
                               M, cluster_rows=cluster_rows, cluster_columns=cluster_columns,
                               width = ncol(M)*grid::unit(aspect_ratio * ux, "mm"),
                               height = nrow(M)*grid::unit(ux, "mm"),
                               row_names_gp=grid::gpar(fontsize=dimname_fontsize),
                               column_names_gp=grid::gpar(fontsize=dimname_fontsize),
                               cell_fun=cell_fun,
                               ...
                             ))
    } else {
      return(ComplexHeatmap::Heatmap(
                               M, cluster_rows=cluster_rows, cluster_columns=cluster_columns,
                               width = ncol(M)*grid::unit(aspect_ratio * ux, "mm"),
                               height = nrow(M)*grid::unit(ux, "mm"),
                               row_names_gp=grid::gpar(fontsize=dimname_fontsize),
                               column_names_gp=grid::gpar(fontsize=dimname_fontsize),
                               ...
                             ))
    }
}
#' Save ComplexHeatmap, modified from Carles Boix
#'
#' @param ht ComplexHeatmap heatmap
#' @param pltprefix Prefix for images
#' @param w Width in inches
#' @param h Height in inches
#' @param extra Function to call after draw, before figure is saved
#' @param dpi DPI of PNG
#' @export
saveHeatmap <- function(ht, pltprefix, w=7, h=7, dpi=600, extra=NULL, ...) {
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
saveGGplot <- function(gp, pltprefix, w=7, h=7, dpi=600) {
    ggplot2::ggsave(paste0(pltprefix, ".pdf"), gp, units="in", dpi=dpi, width=w, height=h)
    ggplot2::ggsave(paste0(pltprefix, ".png"), gp, units="in", dpi=dpi, width=w, height=h)
    print(pltprefix)
}

#' Draw triangles in heatmap
#' Pass cell_fun=ht_triangle_split(...) to Heatmap()
#' @param mat.ul Upper left matrix
#' @param mat.lr Lower right matrix
#' @param col.ul Upper left color
#' @param col.lr Lower right color
#' @param lwd lwd for grid::gpar
#' @param ... Arguments to grid::gpar
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

#' Draw asterisks
#' Pass cell_fun=ht_asterisks(...) to Heatmap()
#' @param p.matrix P value matrix, unsorted
#' @param gp graphic parameters from grid::gpar(). Recommended to set fontsize
#' @param ... Arguments to grid::grid.text
#' @export
ht_asterisks <- function(p.matrix, gp, ...) {
    return(function(j, i, x, y, w, h, fill) {
        if (is.na(p.matrix[i, j]) | is.nan(p.matrix[i, j])) {
        } else if (p.matrix[i, j] < 0.001) {
            grid::grid.text("***", x, y, gp=gp, ...)
        } else if (p.matrix[i, j] < 0.01) {
            grid::grid.text("**", x, y, gp=gp, ...)
        } else if (p.matrix[i, j] < 0.05) {
            grid::grid.text("*", x, y, gp=gp, ...)
        }
    })
}

#' Draw volcano plot
#' @param df dataframe to pass
#' @param label Text label for points
#' @param title Tile for plot
#' @param threshold.FDR FDR threshold for plotting text
#' @param threshold.log2FC Log2FC threshold for plotting text
#' @param force Force parameter for ggrepel:::geom_text_repel()
#' @param max.overlaps Max overlaps for ggrepel::geom_text_repel()
#' @param quantile.log2FC If less than 1, quantile clip log2FC to fit inside smaller plot area
#' @export
volcano <- function(df, label="gene", title="Volcano plot of",
                    threshold.FDR=0.99999, threshold.log2FC=0.5,
                    force=1, max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
                    repel=TRUE,
                    prefix="",
                    quantile.log2FC=1) {
    require(ggplot2)
    df$log2FC = qclip(df[[paste0(prefix, "log2FC")]], quantile.log2FC)
    df$FDR = vclip(df[[paste0(prefix, "FDR")]], min=1e-300)
    df$score = -log10(df$FDR) * abs(df$log2FC)
    if (all(label %in% colnames(df))) {
        df$label = apply(df[label], 1, function(x) { paste0(x, collapse=" ") })
    } else {
        df$label = ""
    }
    df$label = ifelse(df$FDR <= threshold.FDR, df$label, rep("", nrow(df)))
    df$label = ifelse(abs(df$log2FC) >= threshold.log2FC, df$label, rep("", nrow(df)))
    g = ggplot(df, aes(x=log2FC, y=-log10(FDR), label=label))
    if (repel) {
        g = g + geom_point() + ggrepel::geom_text_repel(force=force, max.overlaps=max.overlaps)
    } else {
        g = g + geom_text()
    }
    g = g + ggpubr::theme_pubr() + ggtitle(title)
    return(g)
}

#' @export
pbviolin <- function(sce, gene, groupby, sample.col="Sample", target_sum=10000, plot=TRUE) {
    require(ggplot2)
    X = SummarizedExperiment::assays(sce)$counts
    cd = as.data.frame(SummarizedExperiment::colData(sce))
    if ("total_counts" %in% colnames(cd)) {
        D = Matrix::Diagonal(x=target_sum / cd$total_counts)
    } else {
        D = Matrix::Diagonal(x=target_sum / Matrix::colSums(X))
    }
    X = as.matrix(X[gene,,drop=FALSE] %*% D)
    X = as.data.frame(log1p(t(X)))
    cd = cbind(cd, X)
    ## pseudobulk
    pb = se_make_pseudobulk(sce, sample.col)
    pd = as.data.frame(SummarizedExperiment::colData(pb))
    PX = SummarizedExperiment::assays(pb)$counts
    D = Matrix::Diagonal(x=target_sum / Matrix::colSums(PX))
    PX = as.matrix(PX[gene,,drop=FALSE] %*% D)
    PX = as.data.frame(log1p(t(PX)))
    PX$modality="pseudobulk"
    pd = cbind(pd, PX)
    if (!plot) {
        return(list(sc=cd, pb=pd))
    } else {
        L = lapply(setNames(gene, gene), function(gx) {
            ggplot(cd, aes_string(x=groupby, y=gx, fill=groupby)) + geom_violin() + geom_point(data=pd, position = position_jitter(width = 0.2)) + ggpubr::theme_pubr() + ylab(gx) + xlab(groupby) + ggtitle(gx)
        })
        if (length(L) == 1) {
            return(L[[1]])
        } else {
            return(L)
        }
    }
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
    point_padding_x <- graphics::strheight('a', cex=pt.cex)
    point_padding_y <- graphics::strheight('a', cex=pt.cex)

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

#' @export
convert_pvalue_to_star <- function(pvalue) {
    stars = rep("", length(pvalue))
    stars[pvalue < 0.05] = "*"
    stars[pvalue < 0.01] = "**"
    stars[pvalue < 0.001] = "***"
    return(setNames(stars, names(pvalue)))
}
