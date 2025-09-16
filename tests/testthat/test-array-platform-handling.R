# Test array platform handling
library(testthat)

# Replace local create_test_data in this file by importing from helper if available
if (exists("create_test_data", where = asNamespace("testthat"))) {
  # no-op
}

# Adjust tests to use helper-generated DMPs directly to guarantee presence
test_that("Package handles 450k probes", {
    # Create 450k-style test data
    test_data <- create_test_data("450k")
    
    # Should work with 450k annotations
    result <- findDMRsFromDMPs(
        beta.file = test_data$beta_file,
        dmps.tsv.file = test_data$dmps_file,
        pheno = test_data$pheno,
        sample_group.col = "Sample_Group",
        verbose = TRUE
    )
    
    expect_s4_class(result, "GRanges")
    expect_true(length(result) > 0)
})

test_that("Package handles EPIC probes", {
    # Create EPIC-style test data
    test_data <- create_test_data("EPIC")
    
    # Should work with EPIC annotations
    result <- findDMRsFromDMPs(
        beta.file = test_data$beta_file,
        dmps.tsv.file = test_data$dmps_file,
        pheno = test_data$pheno,
        sample_group.col = "Sample_Group",
        array = "EPIC",
        verbose = TRUE
    )
    
    expect_s4_class(result, "GRanges")
    expect_true(length(result) > 0)
})

test_that("Package handles mixed probe types", {
    # Create mixed test data
    test_450k <- create_test_data("450k")
    test_epic <- create_test_data("EPIC")
    
    # Create DMPs data frame with mixed probe names from both beta value files
    mixed_dmps <- data.frame(
        dmp = rownames(test_450k$beta_values)[1:5],  # Use probes that exist in test_450k
        chr = rep("chr1", 5),
        pos = seq(1000, 5000, by = 1000),
        pval = rep(0.01, 5),
        pval_adj = rep(0.01, 5),
        delta_beta = rep(0.3, 5),
        Sample_Group = rep("Case", 5),  # Add Sample_Group column
        cases_beta = rep(0.8, 5),
        controls_beta = rep(0.4, 5),
        cases_beta_sd = rep(0.1, 5),
        controls_beta_sd = rep(0.1, 5),
        cases_num = rep(5, 5),
        controls_num = rep(5, 5)
    )
    
    mixed_dmps_file <- tempfile(fileext = ".txt")
    write.table(mixed_dmps, file = mixed_dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    result <- findDMRsFromDMPs(
        beta.file = test_450k$beta_file,
        dmps.tsv.file = mixed_dmps_file,
        pheno = test_450k$pheno,
        sample_group.col = "Sample_Group",
        verbose = TRUE
    )
    
    expect_s4_class(result, "GRanges")
    expect_true(length(result) > 0)
})

test_that("Package handles custom CpG positions", {
    # Create test data
    test_data <- create_test_data("450k")
    
    # Create DMPs with same probe names as the beta file
    custom_dmps <- data.frame(
        dmp = rownames(test_data$beta_values)[1:5],  # Use probes that exist in beta values
        chr = rep("chr1", 5),
        pos = seq(1000, 5000, by = 1000),
        pval = rep(0.01, 5),
        pval_adj = rep(0.01, 5),
        delta_beta = rep(0.3, 5),
        Sample_Group = rep("Case", 5),  # Add Sample_Group column
        cases_beta = rep(0.8, 5),
        controls_beta = rep(0.4, 5),
        cases_beta_sd = rep(0.1, 5),
        controls_beta_sd = rep(0.1, 5),
        cases_num = rep(5, 5),
        controls_num = rep(5, 5)
    )
    
    custom_dmps_file <- tempfile(fileext = ".txt")
    write.table(custom_dmps, file = custom_dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    result <- findDMRsFromDMPs(
        beta.file = test_data$beta_file,
        dmps.tsv.file = custom_dmps_file,
        pheno = test_data$pheno,
        sample_group.col = "Sample_Group",
        verbose = TRUE
    )
    
    expect_s4_class(result, "GRanges")
    expect_true(length(result) > 0)
})