context("CpG ID handling")

# Create test data with mixed 450k and EPIC IDs
create_mixed_test_data <- function() {
    # Create beta matrix with mixed IDs
    n_samples <- 10
    n_cpgs <- 10
    mixed_cpg_ids <- c(
        paste0("cg", sprintf("%07d", 1:5)),     # 450k style
        paste0("cg", sprintf("%08d", 6:10))     # EPIC style
    )
    
    beta_values <- matrix(
        runif(n_cpgs * n_samples, 0, 1),
        nrow = n_cpgs,
        dimnames = list(
            mixed_cpg_ids,
            paste0("sample", 1:n_samples)
        )
    )
    
    # Create phenotype data
    pheno <- data.frame(
        Sample_Group = factor(rep(c("Case", "Control"), each = 5)),
        casecontrol = rep(c(1, 0), each = 5),
        row.names = colnames(beta_values)
    )
    
    # Create beta file
    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_values), beta_values),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    
    # Create DMPs with mixed IDs
    dmps <- data.frame(
        dmp = mixed_cpg_ids[1:6],  # Include both 450k and EPIC IDs
        chr = rep("chr1", 6),
        pos = seq(1000, 6000, by = 1000),
        pval = rep(0.01, 6),
        pval_adj = rep(0.01, 6),
        delta_beta = rep(0.3, 6),
        Sample_Group = rep("Case", 6),
        cases_beta = rep(0.8, 6),
        controls_beta = rep(0.4, 6),
        cases_beta_sd = rep(0.1, 6),
        controls_beta_sd = rep(0.1, 6),
        cases_num = rep(5, 6),
        controls_num = rep(5, 6),
        stringsAsFactors = FALSE
    )
    
    # Create DMPs file
    dmps_file <- tempfile(fileext = ".txt")
    write.table(dmps, file = dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    list(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        beta_values = beta_values,
        dmps = dmps,
        mixed_cpg_ids = mixed_cpg_ids
    )
}

test_that("Initial regions preserve CpG IDs", {
    # Create test data with mixed CpG IDs
    test_data <- create_mixed_test_data()
    
    # Get initial DMRs
    dmrs <- findDMRsFromDMPs(
        beta.file = test_data$beta_file,
        dmps.tsv.file = test_data$dmps_file,
        pheno = test_data$pheno,
        sample_group.col = "Sample_Group",
        verbose = TRUE
    )
    
    # Verify CpG IDs are preserved
    expect_s4_class(dmrs, "GRanges")
    expect_true(all(dmrs$start_dmp %in% test_data$mixed_cpg_ids))
    expect_true(all(dmrs$end_dmp %in% test_data$mixed_cpg_ids))
})

test_that("DMR expansion maintains CpG IDs", {
    # Create test data
    test_data <- create_mixed_test_data()
    
    # Get expanded DMRs
    dmrs <- findDMRsFromDMPs(
        beta.file = test_data$beta_file,
        dmps.tsv.file = test_data$dmps_file,
        pheno = test_data$pheno,
        sample_group.col = "Sample_Group",
        expansion.step = 50,  # Small step for testing
        max.lookup.dist = 2000,  # Large enough to find neighbors
        min.cpgs = 2,
        verbose = TRUE
    )
    
    # Verify CpG IDs after expansion
    expect_s4_class(dmrs, "GRanges")
    expect_true(length(dmrs) > 0)
    expect_true(all(dmrs$start_cpg %in% test_data$mixed_cpg_ids))
    expect_true(all(dmrs$end_cpg %in% test_data$mixed_cpg_ids))
})

test_that("findDMRsFromDMPs maintains CpG ID integrity with custom IDs", {
    # Create test data with custom CpG IDs
    test_data <- create_mixed_test_data()
    custom_dmps <- test_data$dmps
    custom_dmps$dmp <- paste0("custom_", seq_len(nrow(custom_dmps)))
    
    custom_dmps_file <- tempfile(fileext = ".txt")
    write.table(custom_dmps, file = custom_dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Find DMRs with custom IDs
    dmrs <- findDMRsFromDMPs(
        beta.file = test_data$beta_file,
        dmps.tsv.file = custom_dmps_file,
        pheno = test_data$pheno,
        sample_group.col = "Sample_Group",
        verbose = TRUE
    )
    
    # Verify custom CpG IDs are preserved
    expect_s4_class(dmrs, "GRanges")
    expect_true(length(dmrs) > 0)
    expect_true(all(startsWith(dmrs$start_dmp, "custom_")))
    expect_true(all(startsWith(dmrs$end_dmp, "custom_")))
    
    # Clean up
    unlink(c(custom_dmps_file))
})