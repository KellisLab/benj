
### if you want to test your own H5AD just assign "h5ad" and delete the download.file(), unlink() lines
h5ad=tempfile()
download.file('http://falexwolf.de/data/pbmc3k_raw.h5ad', h5ad, quiet=TRUE)

C = read_h5ad(h5ad)
obs = SingleCellExperiment::colData(C)
var = SingleCellExperiment::rowData(C)
C = SingleCellExperiment::counts(C)

testthat::test_that("AnnData complete", {
    A = C - Matrix::t(anndata::read_h5ad(h5ad)$X)
    testthat::expect_true(zapsmall(max(A)) == 0)
    testthat::expect_true(zapsmall(min(A)) == 0)
})

set.seed(0)
for (i in seq(1, 5)) {
    OI = sample(rownames(obs),
                sample(seq(2, nrow(obs)), 1),
                replace=FALSE)
    VI = sample(rownames(var),
                sample(seq(2, nrow(var)), 1),
                replace=FALSE)
    testthat::test_that(".obs selection", {
        M = SingleCellExperiment::counts(read_h5ad(h5ad, obs=OI)) - C[,OI]
        testthat::expect_true(zapsmall(min(M)) == 0)
        testthat::expect_true(zapsmall(max(M)) == 0)
    })
    testthat::test_that(".var selection", {
        M = SingleCellExperiment::counts(read_h5ad(h5ad, var=VI)) - C[VI,]
        testthat::expect_true(zapsmall(min(M)) == 0)
        testthat::expect_true(zapsmall(max(M)) == 0)
    })
    testthat::test_that(".obs .var selection", {
        M = SingleCellExperiment::counts(read_h5ad(h5ad, obs=OI, var=VI)) - C[VI,OI]
        testthat::expect_true(zapsmall(min(M)) == 0)
        testthat::expect_true(zapsmall(max(M)) == 0)
    })
}

unlink(h5ad)
