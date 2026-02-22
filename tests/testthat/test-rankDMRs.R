# Test suite for rankDMRs function
library(testthat)

test_that("rankDMRs adds score column to DMRs", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")

    example_output_path <- system.file("extdata", "example_outputChr5And11.rds", package = "DMRsegal")
    dmrs <- readRDS(example_output_path)
    mcols(dmrs)$score <- NULL

    expect_false("score" %in% names(mcols(dmrs)))

    ranked_dmrs <- rankDMRs(
        dmrs = dmrs,
        beta = beta,
        pheno = pheno,
        sample_group_col = "Sample_Group"
    )

    expect_true("score" %in% names(mcols(ranked_dmrs)))
    expect_true(all(mcols(ranked_dmrs)$score >= 0))
    expect_true(all(mcols(ranked_dmrs)$score <= 1))
    expect_equal(length(dmrs), length(ranked_dmrs))
})

test_that("rankDMRs works when called from findDMRsFromSeeds with rank_dmrs=TRUE", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    dmrs <- findDMRsFromSeeds(
        rank_dmrs = TRUE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        annotate_with_genes = FALSE,
        verbose = 0
    )

    expect_true("score" %in% names(mcols(dmrs)))
    expect_true(all(mcols(dmrs)$score >= 0))
    expect_true(all(mcols(dmrs)$score <= 1))
})

test_that("ignored_sample_groups affects detection only, not downstream ranking", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    dmrs <- expect_no_error(findDMRsFromSeeds(
        rank_dmrs = TRUE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        ignored_sample_groups = "cancer",
        min_seeds = 2,
        min_cpgs = 2,
        max_lookup_dist = 2000,
        annotate_with_genes = FALSE,
        verbose = 3,
        njobs = 1
    ))

    expect_s4_class(dmrs, "GRanges")
    expect_true("score" %in% names(mcols(dmrs)))
    expect_true("cv_accuracy" %in% names(mcols(dmrs)))
})

test_that("rankDMRs score values are meaningful", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")

    example_output_path <- system.file("extdata", "example_outputChr5And11.rds", package = "DMRsegal")
    dmrs <- readRDS(example_output_path)
    mcols(dmrs)$score <- NULL

    ranked_dmrs <- rankDMRs(
        dmrs = dmrs,
        beta = beta,
        pheno = pheno,
        sample_group_col = "Sample_Group"
    )

    expect_true(length(unique(mcols(ranked_dmrs)$score)) > 1)
    expect_true(any(mcols(ranked_dmrs)$score > 0.5))
})

test_that("rankDMRs works with different nfold values", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")

    example_output_path <- system.file("extdata", "example_outputChr5And11.rds", package = "DMRsegal")
    dmrs <- readRDS(example_output_path)
    mcols(dmrs)$score <- NULL

    options(DMRsegal.ranking_nfold = 3)
    ranked_dmrs_3fold <- rankDMRs(
        dmrs = dmrs,
        beta = beta,
        pheno = pheno,
        sample_group_col = "Sample_Group"
    )

    options(DMRsegal.ranking_nfold = 5)
    ranked_dmrs_5fold <- rankDMRs(
        dmrs = dmrs,
        beta = beta,
        pheno = pheno,
        sample_group_col = "Sample_Group"
    )

    expect_true("score" %in% names(mcols(ranked_dmrs_3fold)))
    expect_true("score" %in% names(mcols(ranked_dmrs_5fold)))
    expect_equal(length(ranked_dmrs_3fold), length(ranked_dmrs_5fold))
})

test_that("rankDMRs returns DMRs ordered by p-value", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")

    example_output_path <- system.file("extdata", "example_outputChr5And11.rds", package = "DMRsegal")
    dmrs <- readRDS(example_output_path)
    mcols(dmrs)$score <- NULL

    ranked_dmrs <- rankDMRs(
        dmrs = dmrs,
        beta = beta,
        pheno = pheno,
        sample_group_col = "Sample_Group"
    )

    pvals <- mcols(ranked_dmrs)$pval
    expect_true(all(diff(pvals) >= 0))
})
