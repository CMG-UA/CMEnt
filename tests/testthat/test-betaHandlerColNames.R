suppressPackageStartupMessages({
    library(testthat)
    library(CMEnt)
})
options("CMEnt.verbose" = 0)

test_that("BetaHandler preserves sample names for header-only beta files loaded in memory", {
    beta_file <- tempfile(fileext = ".tsv")
    withr::defer(unlink(beta_file))

    writeLines(
        c(
            "SampleA\tSampleB\tSampleC",
            "cg1\t0.1\t0.2\t0.3",
            "cg2\t0.4\t0.5\t0.6"
        ),
        beta_file
    )

    sorted_locs <- data.frame(
        chr = c("chr1", "chr1"),
        start = c(100L, 200L),
        end = c(101L, 201L),
        row.names = c("cg1", "cg2"),
        stringsAsFactors = FALSE
    )

    beta_handler <- getBetaHandler(beta = beta_file, sorted_locs = sorted_locs)

    expect_equal(beta_handler$getBetaColNames(), c("SampleA", "SampleB", "SampleC"))
    expect_equal(colnames(beta_handler$getBeta()), c("SampleA", "SampleB", "SampleC"))
})

test_that("sortBetaFileByCoordinates preserves sample names for header-only beta files", {
    beta_file <- tempfile(fileext = ".tsv")
    sorted_file <- tempfile(fileext = ".tsv.gz")
    withr::defer(unlink(c(beta_file, sorted_file)))

    writeLines(
        c(
            "SampleA\tSampleB\tSampleC",
            "cg2\t0.4\t0.5\t0.6",
            "cg1\t0.1\t0.2\t0.3"
        ),
        beta_file
    )

    sorted_locs <- data.frame(
        chr = c("chr1", "chr1"),
        start = c(100L, 200L),
        end = c(101L, 201L),
        row.names = c("cg1", "cg2"),
        stringsAsFactors = FALSE
    )

    sorted_beta_file <- sortBetaFileByCoordinates(
        beta_file = beta_file,
        genomic_locs = sorted_locs,
        output_file = sorted_file,
        overwrite = TRUE
    )
    sorted_data <- read.delim(
        gzfile(sorted_beta_file),
        header = TRUE,
        check.names = FALSE
    )

    expect_equal(colnames(sorted_data), c("ID", "SampleA", "SampleB", "SampleC"))
    expect_equal(sorted_data$ID, c("cg1", "cg2"))
})

test_that("BetaHandler extracts sample names from minimal tabix headers", {
    tabix_file <- tempfile(fileext = ".bed.gz")
    withr::defer(unlink(c(tabix_file, paste0(tabix_file, ".tbi"))))

    con <- gzfile(tabix_file, "w")
    writeLines(
        c(
            "#chrom\tstart\tSample1\tSample2",
            "chr1\t100\t0.1\t0.2",
            "chr1\t200\t0.3\t0.4"
        ),
        con
    )
    close(con)
    file.create(paste0(tabix_file, ".tbi"))

    beta_handler <- getBetaHandler(beta = tabix_file)

    expect_equal(beta_handler$getBetaColNames(), c("Sample1", "Sample2"))
})
