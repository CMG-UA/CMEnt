options("CMEnt.verbose" = 0)

test_that("buildDMRs accepts covariates without changing synthetic DMR calls", {
    fixture <- makeSyntheticBuildDMRsFixture()

    baseline <- runSyntheticBuildDMRs(fixture)
    adjusted <- runSyntheticBuildDMRs(
        fixture,
        covariates = c("Age", "Gender")
    )

    expect_s4_class(adjusted, "GRanges")
    expect_identical(S4Vectors::mcols(adjusted)$id, S4Vectors::mcols(baseline)$id)
    expect_identical(as.integer(S4Vectors::mcols(adjusted)$seeds_num), 3L)
})

test_that("buildDMRs reuses prepared covariate models across connectivity passes", {
    fixture <- makeSyntheticBuildDMRsFixture()
    original_prepare <- CMEnt:::.prepareGroupCovariateModels
    prepare_calls <- 0L
    local_mocked_bindings(
        .prepareGroupCovariateModels = function(...) {
            prepare_calls <<- prepare_calls + 1L
            original_prepare(...)
        },
        .package = "CMEnt"
    )

    dmrs <- runSyntheticBuildDMRs(
        fixture,
        covariates = c("Age", "Gender"),
        entanglement = "weak",
        max_bridge_extension_gaps = 1L,
        max_bridge_seeds_gaps = 1L
    )

    expect_s4_class(dmrs, "GRanges")
    expect_identical(prepare_calls, 1L)
})

test_that("buildDMRs normalizes scalar list columns in pheno", {
    fixture <- makeSyntheticBuildDMRsFixture()
    fixture$pheno$Sample_Group <- as.list(fixture$pheno$Sample_Group)
    fixture$pheno$casecontrol <- as.list(fixture$pheno$casecontrol)

    dmrs <- runSyntheticBuildDMRs(fixture)

    expect_s4_class(dmrs, "GRanges")
    expect_identical(S4Vectors::mcols(dmrs)$id, "chr1:cgA-cgC")
})

test_that("buildDMRs reads phenotype TSV sample IDs from the first column", {
    fixture <- makeSyntheticBuildDMRsFixture()
    pheno_file <- tempfile(fileext = ".tsv")
    withr::defer(unlink(pheno_file))

    pheno_out <- cbind(sample_id = rownames(fixture$pheno), fixture$pheno)
    write.table(pheno_out, pheno_file, sep = "\t", row.names = FALSE, quote = FALSE)

    dmrs <- runSyntheticBuildDMRs(fixture, pheno = pheno_file)

    expect_s4_class(dmrs, "GRanges")
    expect_identical(S4Vectors::mcols(dmrs)$id, "chr1:cgA-cgC")
})

test_that("buildDMRs derives case/control after aligning pheno to beta samples", {
    fixture <- makeSyntheticBuildDMRsFixture()
    fixture$pheno$Sample_Group <- rep(c("control", "case"), each = 3)
    fixture$pheno$casecontrol <- NULL

    unused_sample <- fixture$pheno[1, , drop = FALSE]
    rownames(unused_sample) <- "unused_sample"
    unused_sample$Sample_Group <- "aaa"
    fixture$pheno <- rbind(unused_sample, fixture$pheno)

    dmrs <- runSyntheticBuildDMRs(fixture, casecontrol_col = NULL)
    dmr_stats <- S4Vectors::mcols(dmrs)

    expect_s4_class(dmrs, "GRanges")
    expect_true(all(is.finite(dmr_stats$cases_beta)))
    expect_true(all(is.finite(dmr_stats$controls_beta)))
    expect_true(all(is.finite(dmr_stats$delta_beta)))
})

test_that("buildDMRs filters synthetic DMRs with min_seeds and max_pval", {
    fixture <- makeSyntheticBuildDMRsFixture()

    baseline <- runSyntheticBuildDMRs(fixture)
    strict_seed_filter <- runSyntheticBuildDMRs(fixture, min_seeds = 3)
    strict_pval_filter <- runSyntheticBuildDMRs(fixture, max_pval = 0.01)

    expect_identical(
        S4Vectors::mcols(baseline)$id,
        "chr1:cgA-cgC"
    )
    expect_identical(S4Vectors::mcols(strict_seed_filter)$id, "chr1:cgA-cgC")
    expect_identical(S4Vectors::mcols(strict_pval_filter)$id, "chr1:cgA-cgC")
})

test_that("buildDMRs adds Stouffer meta p-values and FDR q-values", {
    fixture <- makeSyntheticBuildDMRsFixture()
    dmrs <- runSyntheticBuildDMRs(fixture)
    dmr_stats <- S4Vectors::mcols(dmrs)
    seed_pvals <- fixture$seeds$pval[match(
        CMEnt:::.splitCsvValues(dmr_stats$seeds[[1L]]),
        rownames(fixture$seeds)
    )]
    expected_pval <- stats::pnorm(
        sum(stats::qnorm(seed_pvals, lower.tail = FALSE)) / sqrt(length(seed_pvals)),
        lower.tail = FALSE
    )

    expect_true(all(c("pval", "qval") %in% colnames(dmr_stats)))
    expect_equal(dmr_stats$pval[[1L]], expected_pval)
    expect_equal(dmr_stats$qval[[1L]], expected_pval)
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

test_that("DMR beta aggregation ignores missing per-site beta statistics", {
    beta_stats <- data.frame(
        dmr_id = c(1L, 1L, 1L, 2L, 2L),
        cases_beta = c(0.8, NA_real_, 0.6, NA_real_, NA_real_),
        controls_beta = c(0.2, 0.3, NA_real_, 0.4, NA_real_),
        cases_beta_sd = c(0.02, NA_real_, 0.03, NA_real_, NA_real_),
        controls_beta_sd = c(0.01, 0.04, NA_real_, 0.05, NA_real_)
    )

    agg <- CMEnt:::.aggregateDMRBetaStats(beta_stats, aggfun = mean)

    expect_equal(agg$cases_beta[agg$dmr_id == 1L], 0.7)
    expect_equal(agg$controls_beta[agg$dmr_id == 1L], 0.25)
    expect_equal(agg$cases_beta_min[agg$dmr_id == 1L], 0.6)
    expect_equal(agg$cases_beta_max[agg$dmr_id == 1L], 0.8)
    expect_true(is.na(agg$cases_beta[agg$dmr_id == 2L]))
    expect_equal(agg$controls_beta[agg$dmr_id == 2L], 0.4)
    expect_true(is.na(agg$cases_beta_min[agg$dmr_id == 2L]))
    expect_true(is.na(agg$cases_beta_max[agg$dmr_id == 2L]))
})

test_that("buildDMRs preserves non-tabular columns in TSV outputs", {
    fixture <- makeSyntheticBuildDMRsFixture()
    output_prefix <- tempfile("cment-synth-")
    withr::defer(unlink(Sys.glob(paste0(output_prefix, "*"))))

    dmrs <- expect_no_error(runSyntheticBuildDMRs(
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
