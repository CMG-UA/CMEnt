suppressPackageStartupMessages({
    library(testthat)
    library(CMEnt)
})
options("CMEnt.verbose" = 0)

test_that("getRegistry creates an empty HDF5-backed registry for header-only input", {
    input_file <- tempfile(fileext = ".tsv.gz")
    con <- gzfile(input_file, open = "wt")
    writeLines("#chrom\tstart\tend\tid", con)
    close(con)
    withr::defer(unlink(input_file), teardown_env())

    output_h5file <- tempfile(fileext = ".h5")
    withr::defer(unlink(output_h5file), teardown_env())

    registry <- CMEnt:::getRegistry(
        input_file,
        select = c("#chrom", "start"),
        output_h5file = output_h5file
    )

    expect_s4_class(registry, "DelayedDataFrame")
    expect_equal(nrow(registry), 0)
    expect_equal(ncol(registry), 2)
    expect_equal(colnames(registry), c("#chrom", "start"))

    h5_contents <- rhdf5::h5ls(output_h5file)
    expect_true(all(c("data", "data_colnames", "data_nrows") %in% h5_contents$name))
    expect_equal(as.integer(rhdf5::h5read(output_h5file, "data_nrows")), 0)
})
