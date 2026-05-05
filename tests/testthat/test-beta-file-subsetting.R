suppressPackageStartupMessages({
    library(testthat)
    library(CMEnt)
})
options("CMEnt.verbose" = 0)

test_that(".subsetBetaFile preserves every sample column", {
    beta_file <- tempfile(fileext = ".tsv")
    withr::defer(unlink(beta_file), teardown_env())

    writeLines(
        c(
            "\tS1\tS2\tS3",
            "cg1\t0.1\t0.2\t0.3",
            "cg2\t0.4\t0.5\t0.6"
        ),
        con = beta_file
    )

    beta_subset <- CMEnt:::.subsetBetaFile(
        beta_file,
        sites = c("cg1", "cg2")
    )

    expect_equal(dim(beta_subset), c(2, 3))
    expect_equal(colnames(beta_subset), c("S1", "S2", "S3"))
    expect_equal(
        unname(as.matrix(beta_subset)),
        matrix(
            c(0.1, 0.4, 0.2, 0.5, 0.3, 0.6),
            nrow = 2
        )
    )
})
