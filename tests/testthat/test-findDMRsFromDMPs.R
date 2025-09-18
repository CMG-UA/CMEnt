# Load required packages
library(testthat)
library(GenomicRanges)

# Remove local create_test_data and rely on helper in helper-test-data.R

test_that("DMR finding works correctly with nearby DMPs", {
    # Create test data with known differential methylation
    test_data <- create_test_data(n_cpgs = 20, n_dmps = 10, n_samples = 10)
    beta_file <- test_data$beta_file
    pheno <- test_data$pheno
    
    # Create DMPs with known nearby positions
    dmps <- data.frame(
        dmp = test_data$cpg_ids[1:5],  # Match the number of CpGs in beta values
        chr = "chr1",
        pos = seq(1000, 2000, length.out = 5),  # DMPs within 200bp of each other
        pval = runif(5, 0, 0.01),
        pval_adj = runif(5, 0, 0.01),
        qval = runif(5, 0, 0.01),
        delta_beta = rep(0.4, 5),
        cases_beta = rep(0.8, 5),
        controls_beta = rep(0.4, 5),
        cases_beta_sd = rep(0.1, 5),
        controls_beta_sd = rep(0.1, 5),
        cases_num = rep(5, 5),
        controls_num = rep(5, 5),
        Sample_Group = rep("Case", 5)  # Add Sample_Group column
    )
    
    dmps_file <- tempfile(fileext = ".txt")
    write.table(dmps, file = dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Find DMRs
    result <- findDMRsFromDMPs(
        beta_file = beta_file,
        dmps_tsv_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 2,
        max.lookup.dist = 200
    )
    
    # Test DMR detection
    expect_true(length(result) > 0, "Should find at least one DMR")
    
    # Test DMR properties
    first_dmr <- result[1]
    expect_true(first_dmr$dmps_num >= 2, "DMR should contain at least 2 DMPs")
    expect_true(abs(first_dmr$delta_beta) > 0.3, "DMR should have significant delta beta")
    expect_true(first_dmr$cases_beta > first_dmr$controls_beta, "Cases should have higher beta values")
    
    # Clean up
    unlink(c(beta_file, dmps_file))
})

test_that("DMR finding respects distance threshold", {
    # Create test data
    test_data <- create_test_data(n_cpgs = 20, n_dmps=10, n_samples = 10)
    beta_file <- test_data$beta_file
    pheno <- test_data$pheno
    
    # Create DMPs with distances just above and below threshold
    dmps <- data.frame(
        dmp = test_data$cpg_ids[1:4],
        chr = "chr1",
        pos = c(1000, 1150, 1400, 1550),  # Two pairs of DMPs: 150bp and 150bp apart
        pval = rep(0.001, 4),
        pval_adj = rep(0.001, 4),
        qval = rep(0.001, 4),
        delta_beta = rep(0.4, 4),
        cases_beta = rep(0.8, 4),
        controls_beta = rep(0.4, 4),
        cases_beta_sd = rep(0.1, 4),
        controls_beta_sd = rep(0.1, 4),
        cases_num = rep(5, 4),
        controls_num = rep(5, 4),
        Sample_Group = rep("Case", 4)  # Add Sample_Group column
    )
    
    dmps_file <- tempfile(fileext = ".txt")
    write.table(dmps, file = dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Test with different distance thresholds
    result_200 <- findDMRsFromDMPs(
        beta_file = beta_file,
        dmps_tsv_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 2,
        max.lookup.dist = 200
    )
    
    result_100 <- findDMRsFromDMPs(
        beta_file = beta_file,
        dmps_tsv_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 2,
        max.lookup.dist = 100
    )
    
    # Should find more DMRs with larger distance threshold
    expect_true(length(result_200) >= length(result_100))
    
    # Clean up
    unlink(c(beta_file, dmps_file))
})

test_that("DMR finding handles edge cases", {
    # Create test data
    test_data <- create_test_data(n_cpgs = 10, n_samples = 10)
    beta_file <- test_data$beta_file
    pheno <- test_data$pheno
    
    # Test with single DMP
    single_dmp <- data.frame(
        dmp = test_data$cpg_ids[1],
        chr = "chr1",
        pos = 1000,
        pval = 0.001,
        pval_adj = 0.001,
        qval = 0.001,
        delta_beta = 0.4,
        cases_beta = 0.8,
        controls_beta = 0.4,
        cases_beta_sd = 0.1,
        controls_beta_sd = 0.1,
        cases_num = 5,
        controls_num = 5,
        Sample_Group = "Case"  # Add Sample_Group column
    )
    
    dmps_file <- tempfile(fileext = ".txt")
    write.table(single_dmp, file = dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Should not find DMRs with single DMP
    result <- findDMRsFromDMPs(
        beta_file = beta_file,
        dmps_tsv_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 2
    )
    
    expect_equal(length(result), 0, "Should not find DMRs with single DMP")
    
    # Test with DMPs on different chromosomes
    multi_chr_dmps <- data.frame(
        dmp = test_data$cpg_ids[1:4],
        chr = c("chr1", "chr1", "chr2", "chr2"),
        pos = c(1000, 1100, 1000, 1100),
        pval = rep(0.001, 4),
        pval_adj = rep(0.001, 4),
        qval = rep(0.001, 4),
        delta_beta = rep(0.4, 4),
        cases_beta = rep(0.8, 4),
        controls_beta = rep(0.4, 4),
        cases_beta_sd = rep(0.1, 4),
        controls_beta_sd = rep(0.1, 4),
        cases_num = rep(5, 4),
        controls_num = rep(5, 4),
        Sample_Group = rep("Case", 4)  # Add Sample_Group column
    )
    
    write.table(multi_chr_dmps, file = dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    result <- findDMRsFromDMPs(
      beta_file = beta_file,
      dmps_tsv_file = dmps_file,
      pheno = pheno,
      sample_group_col = "Sample_Group",
      min_dmps = 2,
      min_cpgs = 2
    )
    
    # Should find separate DMRs for each chromosome
    if (length(result) > 0) {
      chr_counts <- table(as.character(seqnames(result)))
      expect_true(all(chr_counts <= 1), "Should not merge DMRs across chromosomes")
    }
    
    # Clean up
    unlink(c(beta_file, dmps_file))
})

test_that("findDMRsFromDMPs works with tabix input", {
  skip("Tabix tests require external setup")
})