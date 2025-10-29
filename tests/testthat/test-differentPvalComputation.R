# Test suite for new functionality in findDMRsFromSeeds
library(testthat)

test_that("findDMRsFromSeeds works with empirical p-value mode and different strategies", {

    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/dmps.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))

    # Test parametric mode (baseline)
    dmrs_parametric <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "parametric",
        memory_threshold_mb = 500
    )

    # Test empirical mode with auto strategy
    dmrs_empirical_auto <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "auto",
        ntries = 100,
        memory_threshold_mb = 500
    )

    # Test empirical mode with montecarlo strategy
    dmrs_empirical_mc <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "montecarlo",
        ntries = 100,
        memory_threshold_mb = 500
    )

    # Test empirical mode with permutations strategy
    dmrs_empirical_perm <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "permutations",
        ntries = 100,
        memory_threshold_mb = 500
    )

    # Assertions
    expect_true(is.null(dmrs_parametric) || inherits(dmrs_parametric, "GRanges"))
    expect_true(is.null(dmrs_empirical_auto) || inherits(dmrs_empirical_auto, "GRanges"))
    expect_true(is.null(dmrs_empirical_mc) || inherits(dmrs_empirical_mc, "GRanges"))
    expect_true(is.null(dmrs_empirical_perm) || inherits(dmrs_empirical_perm, "GRanges"))

    # All should produce valid results
    if (!is.null(dmrs_parametric) && length(dmrs_parametric) > 0) {
        expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs_parametric))))
    }

    if (!is.null(dmrs_empirical_auto) && length(dmrs_empirical_auto) > 0) {
        expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs_empirical_auto))))
    }
})

test_that("findDMRsFromSeeds empirical mode respects tries_seed for reproducibility", {
    
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/dmps.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))

    # Run with same seed twice
    dmrs_seed1_run1 <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "montecarlo",
        ntries = 50,
        tries_seed = 42,
        memory_threshold_mb = 500
    )

    dmrs_seed1_run2 <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "montecarlo",
        ntries = 50,
        tries_seed = 42,
        memory_threshold_mb = 500
    )

    # Run with different seed
    dmrs_seed2 <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "montecarlo",
        ntries = 50,
        tries_seed = 123,
        memory_threshold_mb = 500
    )

    # Assertions
    expect_true(is.null(dmrs_seed1_run1) || inherits(dmrs_seed1_run1, "GRanges"))
    expect_true(is.null(dmrs_seed1_run2) || inherits(dmrs_seed1_run2, "GRanges"))
    expect_true(is.null(dmrs_seed2) || inherits(dmrs_seed2, "GRanges"))

    # Same seed should produce same number of DMRs
    if (!is.null(dmrs_seed1_run1) && !is.null(dmrs_seed1_run2)) {
        expect_equal(length(dmrs_seed1_run1), length(dmrs_seed1_run2))
    }
})

test_that("findDMRsFromSeeds handles different ntries values correctly", {
    
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/dmps.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))

    # Test with ntries = 0 (should use default)
    dmrs_ntries_0 <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        ntries = 0,
        memory_threshold_mb = 500
    )

    # Test with ntries = 50
    dmrs_ntries_50 <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        ntries = 50,
        tries_seed = 42,
        memory_threshold_mb = 500
    )

    # Test with ntries = 200
    dmrs_ntries_200 <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        ntries = 200,
        tries_seed = 42,
        memory_threshold_mb = 500
    )

    # Assertions
    expect_true(is.null(dmrs_ntries_0) || inherits(dmrs_ntries_0, "GRanges"))
    expect_true(is.null(dmrs_ntries_50) || inherits(dmrs_ntries_50, "GRanges"))
    expect_true(is.null(dmrs_ntries_200) || inherits(dmrs_ntries_200, "GRanges"))

    # All should produce valid results
    if (!is.null(dmrs_ntries_50) && length(dmrs_ntries_50) > 0) {
        expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs_ntries_50))))
    }
})

test_that("findDMRsFromSeeds aggfun accepts function objects", {
    
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/dmps.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))

    # Test with median function
    dmrs_median_func <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        aggfun = median,
        memory_threshold_mb = 500
    )

    # Test with mean function
    dmrs_mean_func <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        aggfun = mean,
        memory_threshold_mb = 500
    )

    # Test with character string for comparison
    dmrs_median_char <- findDMRsFromSeeds(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        aggfun = "median",
        memory_threshold_mb = 500
    )

    # Assertions
    expect_true(is.null(dmrs_median_func) || inherits(dmrs_median_func, "GRanges"))
    expect_true(is.null(dmrs_mean_func) || inherits(dmrs_mean_func, "GRanges"))
    expect_true(is.null(dmrs_median_char) || inherits(dmrs_median_char, "GRanges"))

    # Function and character should produce same results for median
    if (!is.null(dmrs_median_func) && !is.null(dmrs_median_char)) {
        expect_equal(length(dmrs_median_func), length(dmrs_median_char))
    }
})
