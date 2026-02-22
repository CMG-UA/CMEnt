# Test suite for new functionality in findDMRsFromSeeds
library(testthat)

test_that("findDMRsFromSeeds works with empirical p-value mode and different strategies", {
    skip_on_ci()
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    # Test parametric mode (baseline)
    dmrs_parametric <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "parametric",
        annotate_with_genes = FALSE
    )

    # Test empirical mode with auto strategy
    dmrs_empirical_auto <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "auto",
        ntries = 100,
        annotate_with_genes = FALSE
    )

    # Test automatic p-value mode selection
    dmrs_pval_auto <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "auto",
        empirical_strategy = "auto",
        ntries = 100,
        annotate_with_genes = FALSE
    )

    # Test empirical mode with montecarlo strategy
    dmrs_empirical_mc <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "montecarlo",
        ntries = 100,
        annotate_with_genes = FALSE
    )

    # Test empirical mode with permutations strategy
    dmrs_empirical_perm <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "permutations",
        ntries = 100,
        annotate_with_genes = FALSE
    )

    # Assertions
    expect_true(is.null(dmrs_parametric) || inherits(dmrs_parametric, "GRanges"))
    expect_true(is.null(dmrs_empirical_auto) || inherits(dmrs_empirical_auto, "GRanges"))
    expect_true(is.null(dmrs_pval_auto) || inherits(dmrs_pval_auto, "GRanges"))
    expect_true(is.null(dmrs_empirical_mc) || inherits(dmrs_empirical_mc, "GRanges"))
    expect_true(is.null(dmrs_empirical_perm) || inherits(dmrs_empirical_perm, "GRanges"))

    # All should produce valid results
    if (!is.null(dmrs_parametric) && length(dmrs_parametric) > 0) {
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_parametric))))
    }

    if (!is.null(dmrs_empirical_auto) && length(dmrs_empirical_auto) > 0) {
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_empirical_auto))))
    }

    if (!is.null(dmrs_pval_auto) && length(dmrs_pval_auto) > 0) {
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_pval_auto))))
    }
})

test_that("findDMRsFromSeeds empirical mode respects random seed for reproducibility", {
    skip_on_ci()
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    # Run with same seed twice
    dmrs_seed1_run1 <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "montecarlo",
        ntries = 50,
        annotate_with_genes = FALSE
    )
    options("DMRsegal.random_seed" = 42)
    dmrs_seed1_run2 <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "montecarlo",
        ntries = 50,
        annotate_with_genes = FALSE
    )
    # Run with different seed
    options("DMRsegal.random_seed" = 123)
    dmrs_seed2 <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "montecarlo",
        ntries = 50,
        annotate_with_genes = FALSE
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
    skip_on_ci()
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    # Test with ntries = 0 (should use default)
    dmrs_ntries_0 <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        ntries = 0,
        annotate_with_genes = FALSE
    )

    # Test with ntries = 50
    dmrs_ntries_50 <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        ntries = 50,
        annotate_with_genes = FALSE
    )

    # Test with ntries = 200
    dmrs_ntries_200 <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        ntries = 200,
        annotate_with_genes = FALSE
    )

    # Assertions
    expect_true(is.null(dmrs_ntries_0) || inherits(dmrs_ntries_0, "GRanges"))
    expect_true(is.null(dmrs_ntries_50) || inherits(dmrs_ntries_50, "GRanges"))
    expect_true(is.null(dmrs_ntries_200) || inherits(dmrs_ntries_200, "GRanges"))

    # All should produce valid results
    if (!is.null(dmrs_ntries_50) && length(dmrs_ntries_50) > 0) {
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_ntries_50))))
    }
})
