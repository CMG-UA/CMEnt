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

test_that("BetaHandler tabix backend supports repeated row requests", {
    skip_if_not(nzchar(Sys.which("tabix")))
    skip_if_not(nzchar(Sys.which("bgzip")))
    withr::local_options(list(CMEnt.beta_in_mem_threshold_mb = 0))

    beta_file <- tempfile(fileext = ".tsv")
    withr::defer(unlink(beta_file))
    writeLines(
        c(
            "ID\tSample1\tSample2",
            "cg1\t0.1\t0.2",
            "cg2\t0.3\t0.4",
            "cg3\t0.5\t0.6"
        ),
        beta_file
    )
    sorted_locs <- data.frame(
        chr = rep("chr1", 3),
        start = c(100L, 200L, 300L),
        end = c(101L, 201L, 301L),
        row.names = c("cg1", "cg2", "cg3")
    )

    beta_handler <- getBetaHandler(beta = beta_file, sorted_locs = sorted_locs)
    beta_subset <- beta_handler$getBeta(
        row_names = c("cg1", "cg2", "cg2", "cg3"),
        col_names = c("Sample1", "Sample2")
    )

    expect_equal(nrow(beta_subset), 4)
    expect_equal(unname(unlist(beta_subset[2, ])), unname(unlist(beta_subset[3, ])))
    expect_equal(rownames(beta_subset), c("cg1", "cg2", "cg2.1", "cg3"))
})

test_that("BetaHandler file fallback refreshes row order after sorting", {
    withr::local_options(list(CMEnt.beta_in_mem_threshold_mb = 0))
    withr::local_envvar(list(PATH = tempfile("no-tabix-")))

    beta_file <- tempfile(fileext = ".tsv")
    withr::defer(unlink(beta_file))
    writeLines(
        c(
            "ID\tSample1\tSample2",
            "cg2\t0.3\t0.4",
            "cg1\t0.1\t0.2"
        ),
        beta_file
    )
    sorted_locs <- data.frame(
        chr = rep("chr1", 2),
        start = c(100L, 200L),
        end = c(101L, 201L),
        row.names = c("cg1", "cg2")
    )

    beta_handler <- suppressWarnings(
        getBetaHandler(beta = beta_file, sorted_locs = sorted_locs)
    )

    expect_equal(rownames(beta_handler$getBetaLocs()), c("cg1", "cg2"))
    expect_equal(
        beta_handler$getBeta(row_names = c("cg2", "cg1")),
        matrix(
            c(0.3, 0.1, 0.4, 0.2),
            nrow = 2,
            dimnames = list(c("cg2", "cg1"), c("Sample1", "Sample2"))
        )
    )
    expect_error(beta_handler$getBeta(row_names = "missing"), "not found")
})

test_that("BetaHandler subset uses exact row-name matching for genomic locations", {
    beta <- matrix(
        c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
        nrow = 3,
        dimnames = list(c("cg1", "cg10", "cg2"), c("Sample1", "Sample2"))
    )
    sorted_locs <- data.frame(
        chr = rep("chr1", 3),
        start = c(100L, 200L, 300L),
        end = c(101L, 201L, 301L),
        row.names = c("cg1", "cg10", "cg2")
    )

    beta_handler <- getBetaHandler(beta = beta, sorted_locs = sorted_locs)
    subset_handler <- beta_handler$subset(row_names = c("cg1", "cg10"))

    expect_equal(rownames(subset_handler$getBetaLocs()), c("cg1", "cg10"))
})
