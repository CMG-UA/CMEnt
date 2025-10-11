# Test array platform support
library(testthat)

test_that("DMR finding works with both 450k and EPIC array platforms", {
    # 450K test data (already sorted and with internally consistent DMPs)
    test_450k <- create_test_data(n_cpgs = 50, n_dmps = 10, n_samples = 10, platform = "450K")
    # Use a subset of existing helper-generated dmps (ensures matching chr/pos)
    dmps_450k <- test_450k$dmps[1:10, ]
    # Ensure qval column exists (core function expects it when computing summaries)
    if (!"qval" %in% colnames(dmps_450k)) dmps_450k$qval <- dmps_450k$pval_adj
    test_450k$dmps_file <- tempfile(fileext = ".txt")
    write.table(dmps_450k, file = test_450k$dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)

    # EPIC test data
    test_epic <- create_test_data(n_cpgs = 50, n_dmps = 10, n_samples = 10, platform = "EPIC")
    dmps_epic <- test_epic$dmps[1:10, ]
    if (!"qval" %in% colnames(dmps_epic)) dmps_epic$qval <- dmps_epic$pval_adj
    test_epic$dmps_file <- tempfile(fileext = ".txt")
    write.table(dmps_epic, file = test_epic$dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)

    # Run DMR finding for 450k
    result_450k <- findDMRsFromSeeds(
        beta_file = test_450k$beta_file,
        dmps_file = test_450k$dmps_file,
        pheno = test_450k$pheno,
        sample_group_col = "Sample_Group",
        array = "450K",
        output.dir = tempdir(),
        output.id = paste0("aps_450k_", as.integer(Sys.time()))
    )

    # Run DMR finding for EPIC
    result_epic <- findDMRsFromSeeds(
        beta_file = test_epic$beta_file,
        dmps_file = test_epic$dmps_file,
        pheno = test_epic$pheno,
        sample_group_col = "Sample_Group",
        array = "EPIC",
        output.dir = tempdir(),
        output.id = paste0("aps_epic_", as.integer(Sys.time()))
    )

    # Check that results have expected structure
    expect_s4_class(result_450k, "GRanges")
    expect_s4_class(result_epic, "GRanges")

    # Check additional metadata columns
    epic_metadata_cols <- c(
        "n_cpgs", "n_dmps", "mean_delta_beta", "max_delta_beta",
        "min_pval", "cases_beta", "controls_beta", "corr_pval"
    )

    expect_true(all(epic_metadata_cols %in% colnames(mcols(result_epic))))
    expect_true(all(epic_metadata_cols %in% colnames(mcols(result_450k))))
})
