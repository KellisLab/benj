#' Install extra packages through devtools
#' @export
installExtraPackages <- function() {
    devtools::install_github('lhe17/nebula', upgrade='never')
    devtools::install_github('cole-trapnell-lab/cicero-release', ref='monocle3', upgrade='never')
    devtools::install_github('scverse/anndataR', upgrade='never')
    devtools::install_github('bnprks/BPCells', upgrade='never')
    devtools::install_github('chrchang/plink-ng', subdir='2.0/pgenlibr', upgrade='never')
    devtools::install_github('GreenleafLab/ArchR', ref='master', repos=BiocManager::repositories(), upgrade='never')
    BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
    ArchR::installExtraPackages()
}
