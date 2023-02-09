
make_aggregation_table <- function(dir.list) {
    dirlist = do.call(rbind, lapply(dir.list, function(dir) {
        file = list.files(dir, full.names=TRUE)
        data.frame(file=file, dir=dir)
    }))
    is_atac = file.exists(paste0(dirlist$file, "/outs/atac_fragments.tsv.gz"))
    dirlist = dirlist[is_atac,]
    df = data.frame(library_id=make.unique(basename(dirlist$file), sep="_"),
                    atac_fragments=paste0(dirlist$file, "/outs/atac_fragments.tsv.gz"),
                    per_barcode_metrics=paste0(dirlist$file, "/outs/per_barcode_metrics.csv"),
                    gex_molecule_info=paste0(dirlist$file, "/outs/gex_molecule_info.h5"))
    if (!is.null(names(dir.list))) {
        df$batch = names(dir.list)[match(dirlist$dir, dir.list)]
    }
    return(df)
}
