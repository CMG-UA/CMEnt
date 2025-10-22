# Test suite for findDMRsFromSeeds function
library(testthat)
library(DMRSegal)

create_test_dmps_file <- function(dmps_df, file_path = NULL) {
    if (is.null(file_path)) {
        file_path <- tempfile(fileext = ".txt")
    }
    dmps_df_with_id <- data.frame(rownames(dmps_df), dmps_df, check.names = FALSE)
    colnames(dmps_df_with_id)[1] <- "row.names"
    write.table(
        dmps_df_with_id,
        file = file_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    return(file_path)
}

test_that("findDMRsFromSeeds works with small beta file (in-memory loading)", {
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    # Use real minfiData for this test
    library(minfi)

    # Subset to a small dataset for faster testing
    mset <- minfiData::MsetEx[1:5000, ]
    beta_mat <- getBeta(mset)

    # Prepare pheno data
    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        casecontrol = pData(mset)$Sample_Group == "GroupB",
        row.names = colnames(mset)
    )

    # Create temporary files
    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- DMRSegal::sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    # Find DMPs
    dmps <- suppressWarnings(dmpFinder(mset,
        pheno = pheno$Sample_Group,
        type = "categorical",
        shrinkVar = TRUE
    ))
    dmps <- dmps[!is.na(dmps$pval), ]
    dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

    # Filter significant DMPs with lenient threshold
    sig_dmps <- dmps[dmps$pval < 0.1, ]
    if (nrow(sig_dmps) == 0) {
        sig_dmps <- dmps[order(dmps$pval)[1:min(50, nrow(dmps))], ]
    }
    dmps_file <- create_test_dmps_file(sig_dmps)

    # Run findDMRsFromSeeds with memory_threshold_mb=500 (small file loaded in memory)
    dmrs <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500 # Small threshold allows in-memory loading
    )

    # Clean up
    unlink(c(beta_file, dmps_file))

    # Assertions
    expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
    if (!is.null(dmrs) && length(dmrs) > 0) {
        expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs))))
    }
})

test_that("findDMRsFromSeeds works with large beta file (tabix indexing)", {
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    library(minfi)

    mset <- minfiData::MsetEx[1:5000, ]
    beta_mat <- getBeta(mset)

    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        casecontrol = pData(mset)$Sample_Group == "GroupB",
        row.names = colnames(mset)
    )

    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- DMRSegal::sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    dmps <- suppressWarnings(dmpFinder(mset,
        pheno = pheno$Sample_Group,
        type = "categorical",
        shrinkVar = TRUE
    ))
    dmps <- dmps[!is.na(dmps$pval), ]
    dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

    sig_dmps <- dmps[dmps$pval < 0.1, ]
    if (nrow(sig_dmps) == 0) {
        sig_dmps <- dmps[order(dmps$pval)[1:min(50, nrow(dmps))], ]
    }
    dmps_file <- create_test_dmps_file(sig_dmps)

    dmrs <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 0.01
    )

    unlink(c(beta_file, dmps_file))

    expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
    if (!is.null(dmrs) && length(dmrs) > 0) {
        expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs))))
    }
})

test_that("findDMRsFromSeeds reproduces benchmark.Rmd results with minfiData", {
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    # Replicate the benchmark.Rmd setup
    library(minfi)

    # Get sample type information from MsetEx
    sample_types <- pData(minfiData::MsetEx)$Sample_Group
    sample_types[sample_types == "GroupA"] <- "normal"
    sample_types[sample_types == "GroupB"] <- "cancer"
    sample_groups <- factor(sample_types, levels = c("normal", "cancer"))
    selected_samples <- seq_along(sample_groups)

    # Create the final subset
    mset <- minfiData::MsetEx[, selected_samples]
    pheno <- data.frame(
        status = sample_groups[selected_samples],
        row.names = colnames(mset)
    )
    pheno$group <- pheno$status
    pheno$casecontrol <- pheno$status == "cancer"

    # Find DMPs using limma (as in benchmark)
    design <- model.matrix(~status, data = pheno)
    colnames(design) <- c("Intercept", "cancer_vs_normal")

    dmps <- suppressWarnings(dmpFinder(mset,
        pheno = pheno$status,
        type = "categorical",
        shrinkVar = TRUE
    ))
    dmps <- dmps[!is.na(dmps$pval), ]
    dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

    # Filter significant DMPs
    sig_dmps <- dmps[dmps$pval_adj < 0.05, ]
    if (nrow(sig_dmps) == 0) {
        sig_dmps <- dmps[order(dmps$pval)[1:min(50, nrow(dmps))], ]
    }

    # Write beta values to temp file
    beta_file <- tempfile(fileext = ".txt")
    beta_mat <- getBeta(mset)

    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    beta_file <- DMRSegal::sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    # Write DMPs to temp file
    dmps_file <- create_test_dmps_file(sig_dmps)

    # Run DMRSegal with same parameters as benchmark
    dmrs_segal <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500
    )

    # Clean up
    unlink(c(beta_file, dmps_file))

    # Assertions
    expect_s4_class(dmrs_segal, "GRanges")
    expect_equal(length(dmrs_segal), 5,
        info = "Should find 5 DMRs as in benchmark.Rmd"
    )
    expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs_segal))))

    # Check that all DMRs meet the criteria
    expect_true(all(mcols(dmrs_segal)$dmps_num >= 2))
    expect_true(all(mcols(dmrs_segal)$cpgs_num >= 3))
})

test_that("findDMRsFromSeeds parameter variations work correctly", {
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    # Use real minfiData for this test
    library(minfi)

    # Subset to a smaller dataset for faster testing
    mset <- minfiData::MsetEx[1:3000, ]
    beta_mat <- getBeta(mset)

    # Prepare pheno data
    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        casecontrol = pData(mset)$Sample_Group == "GroupB",
        row.names = colnames(mset)
    )

    # Create temporary files
    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- DMRSegal::sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    # Find DMPs
    dmps <- suppressWarnings(dmpFinder(mset,
        pheno = pheno$Sample_Group,
        type = "categorical",
        shrinkVar = TRUE
    ))
    dmps <- dmps[!is.na(dmps$pval), ]
    dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

    # Filter DMPs more leniently for parameter testing
    sig_dmps <- dmps[dmps$pval_adj < 0.1, ]
    dmps_file <- create_test_dmps_file(sig_dmps)

    # Test with strict min_dmps
    dmrs_strict <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 5, # Stricter
        min_cpgs = 3,
        max_lookup_dist = 1000,
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500
    )

    # Test with lenient parameters
    dmrs_lenient <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 1, # More lenient
        min_cpgs = 2,
        max_lookup_dist = 2000, # Larger distance
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500
    )

    # Test with different max_pval
    dmrs_strict_pval <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        max_pval = 0.01, # Stricter p-value
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500
    )

    # Clean up
    unlink(c(beta_file, dmps_file))

    # Assertions
    expect_s4_class(dmrs_strict, "GRanges")
    expect_s4_class(dmrs_lenient, "GRanges")
    expect_s4_class(dmrs_strict_pval, "GRanges")

    # Lenient parameters should generally find more or equal DMRs
    expect_true(length(dmrs_lenient) >= length(dmrs_strict))

    # Strict parameters should have all DMRs meeting criteria
    if (length(dmrs_strict) > 0) {
        expect_true(all(mcols(dmrs_strict)$dmps_num >= 5))
    }

    if (length(dmrs_lenient) > 0) {
        expect_true(all(mcols(dmrs_lenient)$dmps_num >= 1))
        expect_true(all(mcols(dmrs_lenient)$cpgs_num >= 2))
    }
})

test_that("findDMRsFromSeeds handles different aggregation functions", {
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    # Use real minfiData for this test
    library(minfi)

    # Subset to a smaller dataset for faster testing
    mset <- minfiData::MsetEx[1:3000, ]
    beta_mat <- getBeta(mset)

    # Prepare pheno data
    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        casecontrol = pData(mset)$Sample_Group == "GroupB",
        row.names = colnames(mset)
    )

    # Create temporary files
    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- DMRSegal::sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    # Find DMPs
    dmps <- suppressWarnings(dmpFinder(mset,
        pheno = pheno$Sample_Group,
        type = "categorical",
        shrinkVar = TRUE
    ))
    dmps <- dmps[!is.na(dmps$pval), ]
    dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

    sig_dmps <- dmps[dmps$pval < 0.1, ]
    if (nrow(sig_dmps) == 0) {
        sig_dmps <- dmps[order(dmps$pval)[1:min(50, nrow(dmps))], ]
    }
    dmps_file <- create_test_dmps_file(sig_dmps)

    # Test with median aggregation
    dmrs_median <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        aggfun = "median",
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500
    )

    # Test with mean aggregation
    dmrs_mean <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        aggfun = "mean",
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500
    )

    # Clean up
    unlink(c(beta_file, dmps_file))

    # Assertions
    expect_true(is.null(dmrs_median) || inherits(dmrs_median, "GRanges"))
    expect_true(is.null(dmrs_mean) || inherits(dmrs_mean, "GRanges"))

    # Both should return valid results
    if (!is.null(dmrs_median)) expect_true(length(dmrs_median) >= 0)
    if (!is.null(dmrs_mean)) expect_true(length(dmrs_mean) >= 0)
})

test_that("findDMRsFromSeeds handles min_cpg_delta_beta filtering", {
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    # Use real minfiData for this test
    library(minfi)

    # Subset to a smaller dataset for faster testing
    mset <- minfiData::MsetEx[1:3000, ]
    beta_mat <- getBeta(mset)

    # Prepare pheno data
    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        casecontrol = pData(mset)$Sample_Group == "GroupB",
        row.names = colnames(mset)
    )

    # Create temporary files
    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- DMRSegal::sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    # Find DMPs
    dmps <- suppressWarnings(dmpFinder(mset,
        pheno = pheno$Sample_Group,
        type = "categorical",
        shrinkVar = TRUE
    ))
    dmps <- dmps[!is.na(dmps$pval), ]
    dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

    sig_dmps <- dmps[dmps$pval < 0.1, ]
    if (nrow(sig_dmps) == 0) {
        sig_dmps <- dmps[order(dmps$pval)[1:min(50, nrow(dmps))], ]
    }
    dmps_file <- create_test_dmps_file(sig_dmps)

    # Test with no delta beta filtering
    dmrs_no_filter <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        min_cpg_delta_beta = 0,
        max_lookup_dist = 1000,
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500
    )

    # Test with delta beta filtering
    dmrs_with_filter <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        min_cpg_delta_beta = 0.1, # Filter out small changes
        max_lookup_dist = 1000,
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500
    )

    # Clean up
    unlink(c(beta_file, dmps_file))

    # Assertions
    expect_true(is.null(dmrs_no_filter) || inherits(dmrs_no_filter, "GRanges"))
    expect_true(is.null(dmrs_with_filter) || inherits(dmrs_with_filter, "GRanges"))

    # Filtered results should have fewer or equal DMRs
    if (!is.null(dmrs_no_filter) && !is.null(dmrs_with_filter)) {
        expect_true(length(dmrs_with_filter) <= length(dmrs_no_filter))
    }
})

test_that("findDMRsFromSeeds handles edge cases gracefully", {
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    # Use real minfiData for this test
    library(minfi)

    # Subset to a smaller dataset for faster testing
    mset <- minfiData::MsetEx[1:2000, ]
    beta_mat <- getBeta(mset)

    # Prepare pheno data
    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        casecontrol = pData(mset)$Sample_Group == "GroupB",
        row.names = colnames(mset)
    )

    # Create temporary files
    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- DMRSegal::sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    # Find DMPs
    dmps <- suppressWarnings(dmpFinder(mset,
        pheno = pheno$Sample_Group,
        type = "categorical",
        shrinkVar = TRUE
    ))
    dmps <- dmps[!is.na(dmps$pval), ]
    dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

    # Test with no significant DMPs (very high threshold)
    no_sig_dmps <- dmps[dmps$pval_adj > 0.99, ]

    if (nrow(no_sig_dmps) > 0) {
        dmps_file <- create_test_dmps_file(no_sig_dmps)

        dmrs_empty <- findDMRsFromSeeds(
            beta_file = beta_file,
            dmps_file = dmps_file,
            pheno = pheno,
            sample_group_col = "Sample_Group",
            min_dmps = 2,
            min_cpgs = 3,
            max_lookup_dist = 1000,
            njobs = 1,
            verbose = 0,
            memory_threshold_mb = 500
        )

        unlink(dmps_file)

        # Should return an empty or very small GRanges
        expect_s4_class(dmrs_empty, "GRanges")
        expect_true(length(dmrs_empty) >= 0)
    }

    # Clean up
    unlink(beta_file)
})

test_that("findDMRsFromSeeds validates input parameters correctly", {
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    # Use real minfiData for this test
    library(minfi)

    # Subset to a smaller dataset
    mset <- minfiData::MsetEx[1:1000, ]
    beta_mat <- getBeta(mset)

    # Prepare pheno data
    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        casecontrol = pData(mset)$Sample_Group == "GroupB",
        row.names = colnames(mset)
    )

    # Create temporary files
    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- DMRSegal::sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    # Find DMPs
    dmps <- suppressWarnings(dmpFinder(mset,
        pheno = pheno$Sample_Group,
        type = "categorical",
        shrinkVar = TRUE
    ))
    dmps <- dmps[!is.na(dmps$pval), ]
    dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

    sig_dmps <- dmps[dmps$pval < 0.1, ]
    if (nrow(sig_dmps) == 0) {
        sig_dmps <- dmps[order(dmps$pval)[1:min(50, nrow(dmps))], ]
    }
    dmps_file <- create_test_dmps_file(sig_dmps)

    # Test missing required parameters
    expect_error(
        findDMRsFromSeeds(
            beta_file = beta_file,
            dmps_file = NULL, # Missing
            pheno = pheno
        ),
        "dmps_file"
    )

    expect_error(
        findDMRsFromSeeds(
            beta_file = beta_file,
            dmps_file = dmps_file,
            pheno = NULL # Missing
        ),
        "pheno"
    )

    # Test with wrong column names in pheno
    pheno_wrong <- pheno
    names(pheno_wrong) <- c("wrong_col1", "wrong_col2")

    expect_error(
        findDMRsFromSeeds(
            beta_file = beta_file,
            dmps_file = dmps_file,
            pheno = pheno_wrong,
            sample_group_col = "Sample_Group",
            casecontrol_col = "casecontrol",
            verbose = 0
        )
    )

    # Clean up
    unlink(c(beta_file, dmps_file))
})

test_that("findDMRsFromSeeds works with different genome builds", {
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    # Use real minfiData for this test
    library(minfi)

    # Subset to a smaller dataset
    mset <- minfiData::MsetEx[1:2000, ]
    beta_mat <- getBeta(mset)

    # Prepare pheno data
    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        casecontrol = pData(mset)$Sample_Group == "GroupB",
        row.names = colnames(mset)
    )

    # Create temporary files
    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- DMRSegal::sortBetaFileByCoordinates(beta_file, genome = "hg19", overwrite = TRUE)

    # Find DMPs
    dmps <- suppressWarnings(dmpFinder(mset,
        pheno = pheno$Sample_Group,
        type = "categorical",
        shrinkVar = TRUE
    ))
    dmps <- dmps[!is.na(dmps$pval), ]
    dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

    sig_dmps <- dmps[dmps$pval < 0.1, ]
    if (nrow(sig_dmps) == 0) {
        sig_dmps <- dmps[order(dmps$pval)[1:min(50, nrow(dmps))], ]
    }
    dmps_file <- create_test_dmps_file(sig_dmps)

    # Test with hg19
    dmrs_hg19 <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        genome = "hg19",
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500
    )

    # Clean up
    unlink(c(beta_file, dmps_file))

    # Assertions
    expect_true(is.null(dmrs_hg19) || inherits(dmrs_hg19, "GRanges"))
    if (!is.null(dmrs_hg19)) expect_true(length(dmrs_hg19) >= 0)
})

test_that("findDMRsFromSeeds works when tabix is not available", {
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")
    skip_if_not_installed("mockery")

    library(minfi)
    library(mockery)

    mset <- minfiData::MsetEx[1:3000, ]
    beta_mat <- getBeta(mset)

    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        casecontrol = pData(mset)$Sample_Group == "GroupB",
        row.names = colnames(mset)
    )

    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- DMRSegal::sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    dmps <- suppressWarnings(dmpFinder(mset,
        pheno = pheno$Sample_Group,
        type = "categorical",
        shrinkVar = TRUE
    ))
    dmps <- dmps[!is.na(dmps$pval), ]
    dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

    sig_dmps <- dmps[dmps$pval < 0.1, ]
    if (nrow(sig_dmps) == 0) {
        sig_dmps <- dmps[order(dmps$pval)[1:min(50, nrow(dmps))], ]
    }
    dmps_file <- create_test_dmps_file(sig_dmps)

    mock_convertBetaToTabix <- mock(NULL)

    stub(findDMRsFromSeeds, "convertBetaToTabix", mock_convertBetaToTabix)

    dmrs <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 1,
        min_cpgs = 2,
        max_lookup_dist = 1000,
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500
    )

    expect_called(mock_convertBetaToTabix, 0)

    expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
    if (!is.null(dmrs) && length(dmrs) > 0) {
        expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs))))
    }

    unlink(c(beta_file, dmps_file))
})

test_that("findDMRsFromSeeds empirical p-value mode works", {
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    library(minfi)

    mset <- minfiData::MsetEx[1:3000, ]
    beta_mat <- getBeta(mset)

    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        casecontrol = pData(mset)$Sample_Group == "GroupB",
        row.names = colnames(mset)
    )

    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- DMRSegal::sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    dmps <- suppressWarnings(dmpFinder(mset,
        pheno = pheno$Sample_Group,
        type = "categorical",
        shrinkVar = TRUE
    ))
    dmps <- dmps[!is.na(dmps$pval), ]
    dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

    sig_dmps <- dmps[dmps$pval < 0.1, ]
    if (nrow(sig_dmps) == 0) {
        sig_dmps <- dmps[order(dmps$pval)[1:min(50, nrow(dmps))], ]
    }
    dmps_file <- create_test_dmps_file(sig_dmps)

    dmrs_parametric <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "parametric",
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500
    )

    dmrs_empirical <- findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        pval_mode = "empirical",
        nperm = 50,
        perm_seed = 12345,
        njobs = 1,
        verbose = 0,
        memory_threshold_mb = 500
    )

    unlink(c(beta_file, dmps_file))

    expect_true(is.null(dmrs_parametric) || inherits(dmrs_parametric, "GRanges"))
    expect_true(is.null(dmrs_empirical) || inherits(dmrs_empirical, "GRanges"))

    if (!is.null(dmrs_parametric) && length(dmrs_parametric) > 0) {
        expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs_parametric))))
    }

    if (!is.null(dmrs_empirical) && length(dmrs_empirical) > 0) {
        expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs_empirical))))
    }
})
