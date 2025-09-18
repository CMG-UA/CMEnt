# Test parallel processing functionality
library(testthat)

test_that("Parallel processing produces same results as sequential", {
    # Create test data
    test_data <- create_test_data(n_cpgs = 100, n_samples = 20)
    beta_file <- test_data$beta_file
    pheno <- test_data$pheno
    
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
    
    # Run with single core
    result_single <- findDMRsFromDMPs(
        beta_file = beta_file,
        dmps_tsv_file = dmps_file,
        pheno = pheno,
        sample_group_col = "group",
        min_dmps = 2,
        min_cpgs = 2,
        njobs = 1
    )
    
    # Run with multiple cores
    result_parallel <- findDMRsFromDMPs(
        beta_file = beta_file,
        dmps_tsv_file = dmps_file,
        pheno = pheno,
        sample_group_col = "group",
        min_dmps = 2,
        min_cpgs = 2,
        njobs = 2
    )
    
    # Check that results are identical
    expect_equal(length(result_single), length(result_parallel))
    expect_equal(
        sort(mcols(result_single)$mean_delta_beta), 
        sort(mcols(result_parallel)$mean_delta_beta),
        tolerance = 1e-6
    )
    expect_equal(
        sort(mcols(result_single)$min_pval),
        sort(mcols(result_parallel)$min_pval),
        tolerance = 1e-6
    )
})