test_that("plotDMRsManhattan returns a ggplot object", {
    skip_if_not_installed("ggplot2")

    dmrs_path <- system.file("extdata/example_output.rds", package = "DMRsegal")
    if (!nzchar(dmrs_path) || !file.exists(dmrs_path)) {
        skip("Benchmark DMRs not available")
    }
    dmrs <- readRDS(dmrs_path)
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("Benchmark DMRs not available")
    }

    p <- plotDMRsManhattan(dmrs, score_col = "score")
    expect_s3_class(p, "ggplot")
})

test_that("plotDMRsManhattan validates score column", {
    dmrs_path <- system.file("extdata/example_output.rds", package = "DMRsegal")
    if (!nzchar(dmrs_path) || !file.exists(dmrs_path)) {
        skip("Benchmark DMRs not available")
    }
    dmrs <- readRDS(dmrs_path)
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("Benchmark DMRs not available")
    }

    expect_error(
        plotDMRsManhattan(dmrs, score_col = "missing_score_col"),
        "not found"
    )
})

test_that(".normalizeCircosRegion parses region strings", {
    region_df <- DMRsegal:::.normalizeCircosRegion("chr7:1000-2500")
    expect_equal(nrow(region_df), 1)
    expect_equal(region_df$chr[1], "chr7")
    expect_equal(region_df$start[1], 1000)
    expect_equal(region_df$end[1], 2500)
})

test_that(".selectCircosInteractions prioritizes JASPAR-matched interactions", {
    link_data <- data.frame(
        component_id = c(1, 1, 2, 2, 3),
        sim = c(0.91, 0.90, 0.89, 0.88, 0.87),
        has_jaspar_match = c(TRUE, FALSE, TRUE, FALSE, FALSE),
        component_best_rank = c(1, 1, 2, 2, 3),
        stringsAsFactors = FALSE
    )

    selected <- DMRsegal:::.selectCircosInteractions(link_data, max_components = 3)
    expect_equal(nrow(selected), 3)
    expect_equal(sum(selected$has_jaspar_match), 2)
    expect_true(all(selected$has_jaspar_match[seq_len(2)]))
})
