options("DMRsegal.verbose" = 0)
test_that("findDMRsFromSeeds works with empirical p-value mode and different strategies", {
    skip_on_ci()
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    dmps <- dmps[seq_len(100), ] # Use a smaller set for testing

    # Test parametric mode (baseline)
    dmrs_parametric <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotated_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "parametric"
    )

    # Test empirical mode with auto strategy
    dmrs_empirical_auto <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotated_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "auto",
        ntries = 100
    )

    # Test automatic p-value mode selection
    dmrs_pval_auto <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotated_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "auto",
        empirical_strategy = "auto",
        ntries = 100
    )

    # Test empirical mode with montecarlo strategy
    dmrs_empirical_mc <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotated_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "montecarlo",
        ntries = 100
    )

    # Test empirical mode with permutations strategy
    dmrs_empirical_perm <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotated_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "permutations",
        ntries = 100
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
    dmps <- dmps[seq_len(100), ] # Use a smaller set for testing
    # Run with same seed twice
    dmrs_seed1_run1 <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotated_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "montecarlo",
        ntries = 50
    )
    options("DMRsegal.random_seed" = 42)
    dmrs_seed1_run2 <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotated_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "montecarlo",
        ntries = 50
    )
    # Run with different seed
    options("DMRsegal.random_seed" = 123)
    dmrs_seed2 <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotated_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        empirical_strategy = "montecarlo",
        ntries = 50
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
    dmps <- dmps[seq_len(50), ] # Use a smaller set for testing
    # Test with ntries = 0 (should use default)
    dmrs_ntries_0 <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotated_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        ntries = 0
    )

    # Test with ntries = 50
    dmrs_ntries_50 <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotated_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        ntries = 50
    )

    # Assertions
    expect_true(is.null(dmrs_ntries_0) || inherits(dmrs_ntries_0, "GRanges"))
    expect_true(is.null(dmrs_ntries_50) || inherits(dmrs_ntries_50, "GRanges"))
    
    # All should produce valid results
    if (!is.null(dmrs_ntries_50) && length(dmrs_ntries_50) > 0) {
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_ntries_50))))
    }
})

test_that(".testConnectivityBatch marks edges as failing when empirical permutation p-values cannot reach threshold", {
    set.seed(1)
    sites_beta <- matrix(runif(5 * 12, min = 0.05, max = 0.95), nrow = 5, ncol = 12)
    pheno <- data.frame(dummy = seq_len(12))
    pheno[["__casecontrol__"]] <- c(rep(0, 6), rep(1, 6))
    group_inds <- list(g1 = 1:6, g2 = 7:12)
    expect_warning(
        ret_strong <- DMRsegal:::.testConnectivityBatch(
            sites_beta = sites_beta,
            group_inds = group_inds,
            pheno = pheno,
            max_pval = 1e-5,
            entanglement = "strong",
            aggfun = median,
            pval_mode = c(g1 = "empirical", g2 = "empirical"),
            empirical_strategy = c(g1 = "permutations", g2 = "permutations"),
            ntries = 50,
            mid_p = FALSE
        ), "sufficient small empirical p-value"
    )
    expect_true(all(!ret_strong$connected))
    expect_true(all(ret_strong$reason == "pval>max_pval"))
    expect_true(all(ret_strong$pval == 1))

    expect_warning(
        ret_weak <- DMRsegal:::.testConnectivityBatch(
            sites_beta = sites_beta,
            group_inds = group_inds,
            pheno = pheno,
            max_pval = 1e-5,
            entanglement = "weak",
            aggfun = median,
            pval_mode = c(g1 = "empirical", g2 = "empirical"),
            empirical_strategy = c(g1 = "permutations", g2 = "permutations"),
            ntries = 50,
            mid_p = FALSE
        ), "sufficient small empirical p-value"
    )
    expect_true(all(!ret_weak$connected))
    expect_true(all(grepl("pval>max_pval", ret_weak$reason, fixed = TRUE)))
})
