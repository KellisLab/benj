#' Install extra packages through devtools
#' @export
installExtraPackages <- function(upgrade="never", seurat=TRUE, ...) {
### Reinstall seurat
    if (seurat) {
        remotes::install_github("satijalab/seurat", "seurat5", upgrade=upgrade, ...)
        remotes::install_github("satijalab/seurat-data", "seurat5", upgrade=upgrade, ...)
        remotes::install_github("satijalab/seurat-wrappers", "seurat5", upgrade=upgrade, ...)
        remotes::install_github("stuart-lab/signac", upgrade=upgrade, ...)
        remotes::install_github("satijalab/azimuth", "seurat5", upgrade=upgrade, ...)
    }
### Tools
    devtools::install_github('lhe17/nebula', upgrade=upgrade, ...)
    devtools::install_github('cole-trapnell-lab/cicero-release', ref='monocle3', upgrade=upgrade, ...)
    devtools::install_github('scverse/anndataR', upgrade=upgrade, ...)
    devtools::install_github('bnprks/BPCells', upgrade=upgrade, ...)
    devtools::install_github('chrchang/plink-ng', subdir='2.0/pgenlibr', upgrade=upgrade, ...)
    devtools::install_github('GreenleafLab/ArchR', ref='master', repos=BiocManager::repositories(), upgrade=upgrade, ...)
    devtools::install_github("caleblareau/gchromVAR", upgrade=upgrade, ...)
    devtools::install_github("sankaranlab/SCAVENGE", upgrade=upgrade, ...)
    BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
    ArchR::installExtraPackages()
}
