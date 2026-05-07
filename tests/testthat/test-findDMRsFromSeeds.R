options("CMEnt.verbose" = 0)

test_that("findDMRsFromSeeds accepts covariates without changing synthetic DMR calls", {
    fixture <- makeSyntheticFindDMRsFixture()

    baseline <- runSyntheticFindDMRs(fixture)
    adjusted <- runSyntheticFindDMRs(
        fixture,
        covariates = c("Age", "Gender")
    )

    expect_s4_class(adjusted, "GRanges")
    expect_identical(S4Vectors::mcols(adjusted)$id, S4Vectors::mcols(baseline)$id)
    expect_identical(as.integer(S4Vectors::mcols(adjusted)$seeds_num), 3L)
})

test_that("findDMRsFromSeeds normalizes scalar list columns in pheno", {
    fixture <- makeSyntheticFindDMRsFixture()
    fixture$pheno$Sample_Group <- as.list(fixture$pheno$Sample_Group)
    fixture$pheno$casecontrol <- as.list(fixture$pheno$casecontrol)

    dmrs <- runSyntheticFindDMRs(fixture)

    expect_s4_class(dmrs, "GRanges")
    expect_identical(S4Vectors::mcols(dmrs)$id, "chr1:cgA-cgC")
})

test_that("findDMRsFromSeeds filters synthetic DMRs with min_seeds and max_pval", {
    fixture <- makeSyntheticFindDMRsFixture()

    baseline <- runSyntheticFindDMRs(fixture)
    strict_seed_filter <- runSyntheticFindDMRs(fixture, min_seeds = 3)
    strict_pval_filter <- runSyntheticFindDMRs(fixture, max_pval = 0.01)

    expect_identical(
        S4Vectors::mcols(baseline)$id,
        "chr1:cgA-cgC"
    )
    expect_identical(S4Vectors::mcols(strict_seed_filter)$id, "chr1:cgA-cgC")
    expect_identical(S4Vectors::mcols(strict_pval_filter)$id, "chr1:cgA-cgC")
})

test_that("DMR beta aggregation changes with aggfun while preserving effect direction", {
    beta_stats <- makeSyntheticDMRBetaStats()

    agg_mean <- CMEnt:::.aggregateDMRBetaStats(beta_stats, aggfun = mean)
    agg_median <- CMEnt:::.aggregateDMRBetaStats(beta_stats, aggfun = stats::median)

    expect_lt(agg_mean$cases_beta[agg_mean$dmr_id == 1L], agg_median$cases_beta[agg_median$dmr_id == 1L])
    expect_gt(agg_mean$controls_beta[agg_mean$dmr_id == 1L], agg_median$controls_beta[agg_median$dmr_id == 1L])
    expect_true(all(sign(agg_mean$cases_beta) == sign(agg_median$cases_beta)))
    expect_true(all(sign(agg_mean$controls_beta) == sign(agg_median$controls_beta)))
})

test_that("findDMRsFromSeeds preserves non-tabular columns in TSV outputs", {
    fixture <- makeSyntheticFindDMRsFixture()
    output_prefix <- tempfile("cment-synth-")
    withr::defer(unlink(Sys.glob(paste0(output_prefix, "*"))))

    dmrs <- expect_no_error(runSyntheticFindDMRs(
        fixture,
        covariates = c("Age", "Gender"),
        entanglement = "weak",
        max_bridge_extension_gaps = 1,
        max_bridge_seeds_gaps = 1,
        output_prefix = output_prefix
    ))

    dmrs_file <- paste0(output_prefix, ".dmrs.tsv.gz")
    beta_file <- paste0(output_prefix, ".seeds_beta.tsv.gz")
    expect_true(file.exists(dmrs_file))
    expect_true(file.exists(beta_file))

    dmrs_df <- read.delim(gzfile(dmrs_file), check.names = FALSE)
    saved_beta <- read.delim(gzfile(beta_file), row.names = 1, check.names = FALSE)
    supporting_sites <- unique(unlist(lapply(as.character(S4Vectors::mcols(dmrs)$sites), CMEnt:::.splitCsvValues), use.names = FALSE))

    expect_true(all(c("sites", "seeds", "stop_connection_reason") %in% colnames(dmrs_df)))
    expect_setequal(rownames(saved_beta), supporting_sites)
    expect_true(all(nzchar(dmrs_df$sites)))
    expect_true(all(nzchar(dmrs_df$seeds)))
})
