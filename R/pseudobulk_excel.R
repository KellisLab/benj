#' @export
dump.excel <- function(se, output) {
    stopifnot(dir.exists(dirname(output)))
    md = S4Vectors::metadata(se)
    cd = SummarizedExperiment::colData(se)
    rd = SummarizedExperiment::rowData(se)
    A = list(gene_info=cbind(data.frame(gene=rownames(se)), rd),
             sample_metadata=cbind(data.frame(name=colnames(se)), cd),
             deg_metadata=as.data.frame(md$deg))
    if ("subset" %in% names(md)) {
        sf = as.data.frame(md$subset)
        if (nrow(sf) == 0) {
            sf = data.frame(h5ad=md$h5ad)
        } else {
            sf$h5ad = md$h5ad
        }
        A$subset = sf
    }
    writexl::write_xlsx(A, output)
}


#' @export
load.persample.excel <- function(xlsx) {
    sheets = readxl::excel_sheets(xlsx)
    obs = readxl::read_xlsx(xlsx, "sample_metadata")
    if ("genes" %in% sheets) {
        var = readxl::read_xlsx(xlsx, "gene_info")
    } else {
        var = NULL
    }
    qs = c("counts", sheets[grep("^Quantile ", sheets)])
    qs = qs[qs %in% sheets]
    return(SummarizedExperiment::SummarizedExperiment(lapply(setNames(qs, qs), function(qsheet) {
        df = as.data.frame(readxl::read_xlsx(xlsx, qsheet))
        rownames(df) = df[[1]]
        return(as.matrix(df[,-1]))
    }), colData=obs, rowData=var))
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
ht_boxplot <- function(se, box_width=0.6, direction=c("vertical", "horizontal"), ...) {
    direction = match.arg(direction)[1]
    boxplot_stats = as.list(SummarizedExperiment::assays(se))
    boxplot_stats = boxplot_stats[sprintf("Quantile %.2f", seq(0, 1, length.out=5))]
    boxplot_stats = array(unlist(boxplot_stats), dim=c(nrow(boxplot_stats[[1]]),
                                                       ncol(boxplot_stats[[1]]),
                                                       length(boxplot_stats)))
    extremes = list(min=apply(boxplot_stats, 1, min),
                    max=apply(boxplot_stats, 1, max))
    extra_args = list(...)
    print(str(extra_args))
    boxplot_stats = aperm(apply(boxplot_stats, c(2,3), function(vec) {
        ## scale each gene to go from -0.5 to 0.5
        denom = extremes$max - extremes$min
        denom[denom < 1e-100] = mean(denom)
        (vec - extremes$min) / denom - 0.5
    }), c(3, 1, 2))
### https://github.com/jokergoo/ComplexHeatmap/blob/ae0ec42cd2e4e0446c114d23dcf43cf2c2f585c8/R/utils.R#L945
    if(direction == "vertical") {
        return(function(j, i, x, y, width, height, fill) {
### TODO where is X and Y? Should set baselines
            gpar_args = as.list(extra_args)
            gpar_args$fill = fill
            gp = do.call(grid::gpar, gpar_args)
            ### box
            grid::grid.rect(x = x,
                            y = y + height * boxplot_stats[2, i, j],
                            height = height * (boxplot_stats[4, i, j] - boxplot_stats[2, i, j]),
                            width = width*box_width,
                            just = "bottom",
                            default.units = "native", gp = gp)
            ### top whisker
            grid::grid.segments(x - 0.5*box_width * width,
                                y + height * boxplot_stats[5, i, j],
                                x + 0.5*box_width * width,
                                y + height * boxplot_stats[5, i, j],
                                default.units = "native", gp = gp)
            ### top whisker connector
            grid::grid.segments(x, y + height * boxplot_stats[5, i, j],
                                x, y + height * boxplot_stats[4, i, j],
                                default.units = "native", gp = gp)
            ### bottom whisker connector
            grid::grid.segments(x, y + height * boxplot_stats[1, i, j],
                                x, y + height * boxplot_stats[2, i, j],
                                default.units = "native", gp = gp)
            ### bottom whisker
            grid::grid.segments(x - 0.5*box_width*width,
                                y + height * boxplot_stats[1, i, j],
                                x + 0.5*box_width*width,
                                y + height * boxplot_stats[1, i, j],
                                default.units = "native", gp = gp)
            ### median whisker
            grid::grid.segments(x - 0.5*box_width*width,
                                y + height * boxplot_stats[3, i, j],
                                x + 0.5*box_width*width,
                                y + height * boxplot_stats[3, i, j],
                                default.units = "native", gp = gp)
        })
    }
}
