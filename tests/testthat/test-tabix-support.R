# Test tabix file support
library(testthat)
library(Rsamtools)

test_that("DMR finding works with tabix files", {
    skip_if_not_installed("Rsamtools")

    # Create test data
    test_data <- create_test_data(n_cpgs = 100, n_dmps = 10, n_samples = 10)

    # Create bgzipped and tabix-indexed file
    beta_bgz_file <- paste0(test_data$beta_file, ".bgz")
    system2("bgzip", c("-c", test_data$beta_file), stdout = beta_bgz_file)
    system2("tabix", c("-s", "1", "-b", "2", "-e", "2", beta_bgz_file))

    # Create DMPs
    dmps <- data.frame(
        dmp = test_data$cpg_ids[1:50],
        chr = rep("chr1", 50),
        pos = seq(1000, 10000, length.out = 50),
        pval = runif(50, 0, 0.01),
        pval_adj = runif(50, 0, 0.01),
        delta_beta = rep(0.4, 50)
    )

    dmps_file <- tempfile(fileext = ".txt")
    write.table(dmps, file = dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)

    # Run with tabix file
    result_tabix <- findDMRsFromSeeds(
        tabix_file = beta_bgz_file,
        beta_file = NULL,
        dmps_file = dmps_file,
        pheno = test_data$pheno,
        sample_group_col = "group",
        min_dmps = 2,
        min_cpgs = 2
    )

    # Run with regular file for comparison
    result_regular <- findDMRsFromSeeds(
        beta_file = test_data$beta_file,
        dmps_file = dmps_file,
        pheno = test_data$pheno,
        sample_group_col = "group",
        min_dmps = 2,
        min_cpgs = 2
    )

    # Check that results are consistent between tabix and regular file
    expect_equal(length(result_tabix), length(result_regular))
    expect_equal(
        sort(mcols(result_tabix)$mean_delta_beta),
        sort(mcols(result_regular)$mean_delta_beta),
        tolerance = 1e-6
    )

    # Test tabix specific features
    # Test querying specific regions
    region_result <- findDMRsFromSeeds(
        tabix_file = beta_bgz_file,
        dmps_file = dmps_file,
        pheno = test_data$pheno,
        sample_group_col = "group",
        min_dmps = 2,
        min_cpgs = 2,
        region = "chr1:1000-5000"
    )

    expect_true(all(start(region_result) >= 1000))
    expect_true(all(end(region_result) <= 5000))
})
