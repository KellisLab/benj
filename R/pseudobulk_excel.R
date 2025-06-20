#' @export
dump.excel <- function(se, output, use.openxlsx=TRUE) {
    stopifnot(dir.exists(dirname(output)))
    md = S4Vectors::metadata(se)
    cd = SummarizedExperiment::colData(se)
    rd = SummarizedExperiment::rowData(se)
    A = list(gene_info=cbind(data.frame(gene=rownames(se)), rd),
             sample_metadata=cbind(data.frame(name=colnames(se)), cd),
             deg_metadata=as.data.frame(md$deg))
    if ("subset" %in% names(md)) {
        A$subset = as.data.frame(md$subset)
    }
    if ("h5ad" %in% names(md)) {
        A$files = data.frame(h5ad=normalizePath(md$h5ad))
        if (!is.null(md$annotation)) {
            A$files$annotation = normalizePath(md$annotation)
        }
    }
    gic = colnames(A$gene_info)
    if (use.openxlsx) {
        wb = openxlsx::createWorkbook()
        # Loop over the list and add each dataframe as a new sheet
        for (sheet_name in names(A)) {
            openxlsx::addWorksheet(wb, sheet_name)
            openxlsx::writeData(wb, sheet_name, A[[sheet_name]], withFilter=TRUE)
        }

        openxlsx::conditionalFormatting(wb=wb, sheet="gene_info",
                                        type="colourScale",
                                        style=c("#FF0000", "#FFFFFF", "#0000FF"),
                                        rows=1:nrow(rd), cols=grep("^logCPM", gic))
        for (col in grep("log2FC$", gic)) {
            openxlsx::conditionalFormatting(wb=wb, sheet="gene_info",
                                            type="colourScale",
                                            style=c("#FF0000", "#FFFF00", "#00FF00"),
                                            rows=1:nrow(rd), cols=col)
        }
        for (col in grep("FDR$", gic)) {
            openxlsx::conditionalFormatting(wb=wb, sheet="gene_info",
                                            type="colourScale",
                                            style=c("#FEFEFF", "#AAAAFF", "#4444FF"),
                                            rows=1:nrow(rd), cols=col)
        }
        openxlsx::freezePane(wb, sheet="gene_info", firstRow=TRUE, firstCol=TRUE)
        openxlsx::saveWorkbook(wb, output, overwrite=TRUE)
    } else {
        writexl::write_xlsx(A, output)
    }
}


#' @export
load.persample.excel <- function(xlsx) {
    sheets = readxl::excel_sheets(xlsx)
    md = list()
    obs = as.data.frame(readxl::read_xlsx(xlsx, "sample_metadata"))
    rownames(obs) = obs[[1]]
    if ("gene_info" %in% sheets) {
        var = as.data.frame(readxl::read_xlsx(xlsx, "gene_info"))
        rownames(var) = var[[1]]
    } else {
        var = NULL
    }
    if ("deg_metadata" %in% sheets) {
        md$deg = as.list(readxl::read_xlsx(xlsx, "deg_metadata"))
    }
    if ("subset" %in% sheets) {
        md$subset = as.list(readxl::read_xlsx(xlsx, "subset"))
    }
    qs = c("counts", sheets[grep("^Quantile ", sheets)])
    qs = qs[qs %in% sheets]
    return(SummarizedExperiment::SummarizedExperiment(lapply(setNames(qs, qs), function(qsheet) {
        df = as.data.frame(readxl::read_xlsx(xlsx, qsheet))
        rownames(df) = df[[1]]
        return(as.matrix(df[,-1]))
    }), colData=obs, rowData=var, metadata=md))
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
    }), c(3, 2, 1))
    print(str(boxplot_stats))
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
