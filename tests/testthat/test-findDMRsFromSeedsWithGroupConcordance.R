# Test suite for group_concordance_strategy parameter
library(testthat)
library(DMRsegal)

test_that("findDMRsFromSeeds works with relaxed group concordance strategy", {
    beta <- loadExampleInputData("beta")
    dmps <- loadExampleInputData("dmps")
    pheno <- loadExampleInputData("pheno")

    dmrs_relaxed <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        group_concordance_strategy = "relaxed",
        pval_mode = "parametric",
        annotate_with_genes = FALSE,
        verbose = 2
    )

    expect_true(is.null(dmrs_relaxed) || inherits(dmrs_relaxed, "GRanges"))
    if (!is.null(dmrs_relaxed) && length(dmrs_relaxed) > 0) {
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_relaxed))))
    }
})

test_that("relaxed strategy produces more or equal DMRs than strict strategy", {
    beta <- loadExampleInputData("beta")
    dmps <- loadExampleInputData("dmps")
    pheno <- loadExampleInputData("pheno")

    dmrs_strict <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        group_concordance_strategy = "strict",
        pval_mode = "parametric",
        annotate_with_genes = FALSE,
        verbose = 2
    )

    dmrs_relaxed <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        group_concordance_strategy = "relaxed",
        pval_mode = "parametric",
        annotate_with_genes = FALSE,
        verbose = 2
    )

    strict_count <- if (is.null(dmrs_strict)) 0 else length(dmrs_strict)
    relaxed_count <- if (is.null(dmrs_relaxed)) 0 else length(dmrs_relaxed)

    expect_true(relaxed_count >= strict_count)
})

test_that("group_concordance_strategy parameter validates correctly", {
    beta <- loadExampleInputData("beta")
    dmps <- loadExampleInputData("dmps")
    pheno <- loadExampleInputData("pheno")

    expect_error(
        findDMRsFromSeeds(
            rank_dmrs = FALSE,
            beta = beta,
            seeds = dmps,
            pheno = pheno,
            sample_group_col = "Sample_Group",
            group_concordance_strategy = "invalid"
        ),
        "is not a prefix"
    )
})
