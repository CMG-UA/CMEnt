# Test array platform handling
library(testthat)

# Replace local create_test_data in this file by importing from helper if available
if (exists("create_test_data", where = asNamespace("testthat"))) {
    # no-op
}

# Adjust tests to use helper-generated DMPs directly to guarantee presence
test_that("Package handles 450k probes", {
    # Create 450k-style test data
    test_data <- create_test_data(platform = "450k")

    # Should work with 450k annotations
    result <- findDMRsFromSeeds(
        beta_file = test_data$beta_file,
        dmps_file = test_data$dmps_file,
        pheno = test_data$pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 1,
        min_cpgs = 0
    )

    expect_s4_class(result, "GRanges")
    expect_true(length(result) > 0)
})

test_that("Package handles EPIC probes", {
    # Create EPIC-style test data
    test_data <- create_test_data(platform = "EPIC")

    # Should work with EPIC annotations
    result <- findDMRsFromSeeds(
        beta_file = test_data$beta_file,
        dmps_file = test_data$dmps_file,
        pheno = test_data$pheno,
        sample_group_col = "Sample_Group",
        array = "EPIC",
        min_dmps = 1,
        min_cpgs = 0
    )

    expect_s4_class(result, "GRanges")
    expect_true(length(result) > 0)
})
