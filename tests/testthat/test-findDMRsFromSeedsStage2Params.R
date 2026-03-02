
library(testthat)

test_that("findDMRsFromSeeds with expansion_window and max_bridge_seeds_gaps parameters", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    # Test with expansion_window and max_bridge_seeds_gaps
    expect_message(
        dmrs_expanded <- findDMRsFromSeeds(
            rank_dmrs = FALSE,
            beta = beta,
            seeds = dmps,
            pheno = pheno,
            sample_group_col = "Sample_Group",
            min_seeds = 2,
            min_cpgs = 3,
            max_lookup_dist = 1000,
            expansion_window = 1, # Expand DMRs by 1bp
            max_bridge_seeds_gaps = 2, # Allow bridging up to 2 seeds apart
            annotate_with_genes = FALSE,
            verbose = 2
        ), "Stage 2 connectivity restricted to 383 seed-derived windows" # Expect the windows to be the same number as the DMRs at that point
    )

    # Assertions
    expect_s4_class(dmrs_expanded, "GRanges")
    if (!is.null(dmrs_expanded)) {
        expect_true(length(dmrs_expanded) >= 0)
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_expanded))))
    }
})

test_that("findDMRsFromSeeds handles min_cpg_delta_beta filtering", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    # Test with no delta beta filtering
    dmrs_no_filter <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        min_cpg_delta_beta = 0,
        max_lookup_dist = 1000,
        annotate_with_genes = FALSE
    )

    # Test with delta beta filtering
    dmrs_with_filter <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        min_cpg_delta_beta = 0.1, # Filter out small changes
        max_lookup_dist = 1000,
        annotate_with_genes = FALSE
    )


    # Assertions
    expect_true(is.null(dmrs_no_filter) || inherits(dmrs_no_filter, "GRanges"))
    expect_true(is.null(dmrs_with_filter) || inherits(dmrs_with_filter, "GRanges"))

    # Filtered results should have fewer or equal DMRs
    if (!is.null(dmrs_no_filter) && !is.null(dmrs_with_filter)) {
        expect_true(length(dmrs_with_filter) <= length(dmrs_no_filter))
    }
})