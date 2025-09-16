# Load required packages
library(testthat)
library(GenomicRanges)

# Load required packages
library(testthat)
library(GenomicRanges)

#' Create test data for DMR finding tests
#' @param n_cpgs Number of CpGs to generate
#' @param n_samples Number of samples to generate
#' @param seed Random seed for reproducibility
#' @return List containing beta_file, beta_values, cpg_ids, and pheno
create_test_data <- function(n_cpgs = 10, n_samples = 10, seed = 42, platform = "450k") {
    set.seed(seed)
    # Create mix of 450k and EPIC style IDs
    if (platform == "mixed") {
        cpg_ids <- c(
            paste0("cg", sprintf("%07d", 1:(n_cpgs/2))),  # 450k style
            paste0("cg", sprintf("%08d", (n_cpgs/2 + 1):n_cpgs))  # EPIC style
        )
    } else {
        cpg_ids <- paste0("cg", sprintf("%07d", 1:n_cpgs))
    }
    
    # Create beta values with clear differential methylation pattern
    n_half <- n_samples %/% 2
    beta_values <- matrix(
        c(
            # First half of CpGs highly methylated in cases
            matrix(runif(floor(n_cpgs/2) * n_half, 0.7, 0.9), nrow = floor(n_cpgs/2)), # Cases
            matrix(runif(floor(n_cpgs/2) * n_half, 0.2, 0.4), nrow = floor(n_cpgs/2)), # Controls
            # Last half of CpGs no differential methylation
            matrix(runif(ceiling(n_cpgs/2) * n_samples, 0.4, 0.6), nrow = ceiling(n_cpgs/2))
        ),
        nrow = n_cpgs,
        dimnames = list(cpg_ids, paste0("sample", 1:n_samples))
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
    
    # Create phenotype data
    pheno <- data.frame(
        row.names = colnames(beta_values),
        group = factor(rep(c("Case", "Control"), c(n_half, n_samples - n_half))),
        casecontrol = rep(c(TRUE, FALSE), c(n_half, n_samples - n_half))
    )
    
    # Return all created test data
    list(
        beta_file = beta_file,
        beta_values = beta_values,
        cpg_ids = cpg_ids,
        pheno = pheno
    )
}

test_that("DMR finding works correctly with nearby DMPs", {
    # Create test data with known differential methylation
    test_data <- create_test_data(n_cpgs = 20, n_samples = 10)
    beta_file <- test_data$beta_file
    pheno <- test_data$pheno
    
    # Create DMPs with known nearby positions
    dmps <- data.frame(
        dmp = test_data$cpg_ids[1:10],
        chr = "chr1",
        pos = seq(1000, 2800, length.out = 10),  # DMPs within 200bp of each other
        pval = runif(10, 0, 0.01),
        pval_adj = runif(10, 0, 0.01),
        qval = runif(10, 0, 0.01),
        delta_beta = rep(0.4, 10),
        cases_beta = rep(0.8, 10),
        controls_beta = rep(0.4, 10),
        cases_beta_sd = rep(0.1, 10),
        controls_beta_sd = rep(0.1, 10),
        cases_num = rep(5, 10),
        controls_num = rep(5, 10)
    )
    
    dmps_file <- tempfile(fileext = ".txt")
    write.table(dmps, file = dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Find DMRs
    result <- findDMRsFromDMPs(
        beta.file = beta_file,
        dmps.tsv.file = dmps_file,
        pheno = pheno,
        sample_group.col = "group",
        min.dmps = 2,
        min.cpgs = 2,
        max.dist = 200
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
    test_data <- create_test_data(n_cpgs = 20, n_samples = 10)
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
        controls_num = rep(5, 4)
    )
    
    dmps_file <- tempfile(fileext = ".txt")
    write.table(dmps, file = dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Test with different distance thresholds
    result_200 <- findDMRsFromDMPs(
        beta.file = beta_file,
        dmps.tsv.file = dmps_file,
        pheno = pheno,
        sample_group.col = "group",
        min.dmps = 2,
        min.cpgs = 2,
        max.dist = 200
    )
    
    result_100 <- findDMRsFromDMPs(
        beta.file = beta_file,
        dmps.tsv.file = dmps_file,
        pheno = pheno,
        sample_group.col = "group",
        min.dmps = 2,
        min.cpgs = 2,
        max.dist = 100
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
        controls_num = 5
    )
    
    dmps_file <- tempfile(fileext = ".txt")
    write.table(single_dmp, file = dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Should not find DMRs with single DMP
    result <- findDMRsFromDMPs(
        beta.file = beta_file,
        dmps.tsv.file = dmps_file,
        pheno = pheno,
        sample_group.col = "group",
        min.dmps = 2,
        min.cpgs = 2
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
        controls_num = rep(5, 4)
    )
    
    write.table(multi_chr_dmps, file = dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    result <- findDMRsFromDMPs(
        beta.file = beta_file,
        dmps.tsv.file = dmps_file,
        pheno = pheno,
        sample_group.col = "group",
        min.dmps = 2,
        min.cpgs = 2
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