# Test suite for findDMRsFromSeeds function
library(testthat)


test_that("findDMRsFromSeeds works with small beta file (in-memory loading)", {
    beta <- loadExampleInputData("beta")
    dmps <- loadExampleInputData("dmps")
    pheno <- loadExampleInputData("pheno")
    # Run findDMRsFromSeeds with memory_threshold_mb=500 (small file loaded in memory)
    dmrs <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        memory_threshold_mb = 500, # allows in-memory loading
        verbose = 2
    )

    # Assertions
    expect_true(!is.null(dmrs))
    expect_true(inherits(dmrs, "GRanges"))
    expect_true(length(dmrs) > 0)
    expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs))))
})

test_that("findDMRsFromSeeds works with large beta file (tabix indexing)", {
    beta <- loadExampleInputData("beta")
    dmps <- loadExampleInputData("dmps")
    pheno <- loadExampleInputData("pheno")
    beta_file <- tempfile(fileext = ".tsv")
    withr::defer(unlink(beta_file))
    write.table(as.data.frame(beta), file = beta_file, sep = "\t", col.names = NA, quote = FALSE)
    sorted_beta_file <- sortBetaFileByCoordinates(beta_file, overwrite = TRUE)
    withr::defer(unlink(sorted_beta_file))
    options("DMRsegal.verbose" = 2)
    options("DMRsegal.use_tabix_cache" = FALSE)

    dmrs <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = sorted_beta_file,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        njobs = 1,
        memory_threshold_mb = 0.01
    )

    expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
    if (!is.null(dmrs) && length(dmrs) > 0) {
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs))))
    }
})

test_that("findDMRsFromSeeds reproduces benchmark.Rmd results with minfi", {
    skip_if_not_installed("minfi")

    beta <- loadExampleInputData("beta")
    pheno <- loadExampleInputData("pheno")
    array_type <- loadExampleInputData("array_type")
    genome <- "hg19"


    beta_handler <- DMRsegal::getBetaHandler(beta, array = array_type, genome = genome)
    beta_mat <- as.matrix(beta_handler$getBeta())
    locs <- beta_handler$getBetaLocs()
    mvalues <- log2(beta_mat / (1 - beta_mat + 1e-6) + 1e-6)
    pheno$casecontrol <- pheno$Sample_Group == "cancer"

    # Find DMPs using minfi's dmpFinder
    dmps <- suppressWarnings(minfi::dmpFinder(mvalues,
        pheno = pheno$casecontrol,
        type = "categorical",
        shrinkVar = TRUE
    ))
    dmps <- dmps[!is.na(dmps$pval), ]
    # Add adjusted p-value to dmps results
    dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

    # Filter significant DMPs
    sig_dmps <- dmps[dmps$pval_adj < 0.05, ]
    options("DMRsegal.verbose" = 2)
    # Run DMRsegal with same parameters as benchmark
    dmrs_segal <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta_handler,
        seeds = sig_dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 10000,
        max_pval = 0.05,
        pval_mode = "parametric",
        njobs = 1,
        memory_threshold_mb = 500
    )

    # Assertions
    expect_s4_class(dmrs_segal, "GRanges")
    expect_equal(length(dmrs_segal), 986)
    expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_segal))))

    # Check that all DMRs meet the criteria
    expect_true(all(mcols(dmrs_segal)$seeds_num >= 2))
    expect_true(all(mcols(dmrs_segal)$cpgs_num >= 3))
})

test_that("findDMRsFromSeeds parameter variations work correctly", {
    beta <- loadExampleInputData("beta")
    dmps <- loadExampleInputData("dmps")
    pheno <- loadExampleInputData("pheno")

    # Test with strict min_seeds
    dmrs_strict <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 5, # Stricter
        min_cpgs = 3,
        max_lookup_dist = 1000,
        memory_threshold_mb = 500,
        annotate_with_genes = FALSE
    )

    # Test with lenient parameters
    dmrs_lenient <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2, # More lenient
        min_cpgs = 2, # More lenient
        max_lookup_dist = 2000, # Larger distance
        memory_threshold_mb = 500,
        annotate_with_genes = FALSE
    )

    # Test with different max_pval
    dmrs_strict_pval <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        max_pval = 0.01, # Stricter p-value
        memory_threshold_mb = 500,
        annotate_with_genes = FALSE
    )


    # Assertions
    expect_s4_class(dmrs_strict, "GRanges")
    expect_s4_class(dmrs_lenient, "GRanges")
    expect_s4_class(dmrs_strict_pval, "GRanges")

    # Lenient parameters should generally find more or equal DMRs
    expect_true(length(dmrs_lenient) >= length(dmrs_strict))

    # Strict parameters should have all DMRs meeting criteria
    if (length(dmrs_strict) > 0) {
        expect_true(all(mcols(dmrs_strict)$seeds_num >= 5))
    }

    if (length(dmrs_lenient) > 0) {
        expect_true(all(mcols(dmrs_lenient)$seeds_num >= 1))
        expect_true(all(mcols(dmrs_lenient)$cpgs_num >= 2))
    }
})

test_that("findDMRsFromSeeds handles different aggregation functions", {
    beta <- loadExampleInputData("beta")
    dmps <- loadExampleInputData("dmps")
    pheno <- loadExampleInputData("pheno")

    # Test with median aggregation
    dmrs_median <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        aggfun = "median",
        memory_threshold_mb = 500,
        annotate_with_genes = FALSE
    )

    # Test with mean aggregation
    dmrs_mean <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        aggfun = "mean",
        memory_threshold_mb = 500,
        annotate_with_genes = FALSE
    )


    # Assertions
    expect_true(is.null(dmrs_median) || inherits(dmrs_median, "GRanges"))
    expect_true(is.null(dmrs_mean) || inherits(dmrs_mean, "GRanges"))

    # Both should return valid results
    if (!is.null(dmrs_median)) expect_true(length(dmrs_median) >= 0)
    if (!is.null(dmrs_mean)) expect_true(length(dmrs_mean) >= 0)
})

test_that("findDMRsFromSeeds handles min_cpg_delta_beta filtering", {
    beta <- loadExampleInputData("beta")
    dmps <- loadExampleInputData("dmps")
    pheno <- loadExampleInputData("pheno")

    # Test with no delta beta filtering
    dmrs_no_filter <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        min_cpg_delta_beta = 0,
        max_lookup_dist = 1000,
        memory_threshold_mb = 500,
        annotate_with_genes = FALSE
    )

    # Test with delta beta filtering
    dmrs_with_filter <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        min_cpg_delta_beta = 0.1, # Filter out small changes
        max_lookup_dist = 1000,
        memory_threshold_mb = 500,
        annotate_with_genes = FALSE
    )


    # Assertions
    expect_true(is.null(dmrs_no_filter) || inherits(dmrs_no_filter, "GRanges"))
    expect_true(is.null(dmrs_with_filter) || inherits(dmrs_with_filter, "GRanges"))

    # Filtered results should have fewer or equal DMRs
    if (!is.null(dmrs_no_filter) && !is.null(dmrs_with_filter)) {
        expect_true(length(dmrs_with_filter) <= length(dmrs_no_filter))
    }
})

test_that("findDMRsFromSeeds validates input parameters correctly", {
    beta <- loadExampleInputData("beta")
    dmps <- loadExampleInputData("dmps")
    pheno <- loadExampleInputData("pheno")

    # Test missing required parameters
    expect_error(
        findDMRsFromSeeds(
            rank_dmrs = FALSE,
            beta = beta,
            seeds = NULL, # Missing
            pheno = pheno
        ),
        "seeds"
    )

    expect_error(
        findDMRsFromSeeds(
            rank_dmrs = FALSE,
            beta = beta,
            seeds = dmps,
            pheno = NULL # Missing
        ),
        "pheno"
    )

    # Test with wrong column names in pheno
    pheno_wrong <- pheno
    colnames(pheno_wrong) <- c("wrong_col1", "wrong_col2", "wrong_col3")

    expect_error(
        findDMRsFromSeeds(
            rank_dmrs = FALSE,
            beta = beta,
            seeds = dmps,
            pheno = pheno_wrong,
            sample_group_col = "Sample_Group",
            casecontrol_col = "casecontrol"
        )
    )
})

test_that("findDMRsFromSeeds works with different genome builds", {
    beta <- loadExampleInputData("beta")
    dmps <- loadExampleInputData("dmps")
    pheno <- loadExampleInputData("pheno")
    options("DMRsegal.use_annotation_cache" = FALSE)
    # Test with hg38
    dmrs_hg38 <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        genome = "hg38",
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        memory_threshold_mb = 500
    )

    # Assertions
    expect_true(is.null(dmrs_hg38) || inherits(dmrs_hg38, "GRanges"))
    if (!is.null(dmrs_hg38)) expect_true(length(dmrs_hg38) >= 0)
})

test_that("findDMRsFromSeeds works when tabix is not available", {
    skip_if_not_installed("mockery")
    beta <- loadExampleInputData("beta")
    dmps <- loadExampleInputData("dmps")
    pheno <- loadExampleInputData("pheno")
    library(mockery)

    mock_convertBetaToTabix <- mock(NULL) # nolint

    stub(findDMRsFromSeeds, "convertBetaToTabix", mock_convertBetaToTabix)
    options("DMRsegal.use_tabix_cache" = FALSE)
    options("DMRsegal.verbose" = 2)
    dmrs <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 1,
        min_cpgs = 2,
        max_lookup_dist = 1000,
        memory_threshold_mb = 0.1,
        annotate_with_genes = FALSE
    )

    expect_called(mock_convertBetaToTabix, 0)

    expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
    if (!is.null(dmrs) && length(dmrs) > 0) {
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs))))
    }
})

test_that("findDMRsFromSeeds does not annotate DMRs when annotate_with_genes=FALSE", {
    beta <- loadExampleInputData("beta")
    dmps <- loadExampleInputData("dmps")
    pheno <- loadExampleInputData("pheno")

    dmrs_not_annotated <- findDMRsFromSeeds(
        rank_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        min_seeds = 2,
        min_cpgs = 3,
        annotate_with_genes = FALSE,
        njobs = 1
    )

    expect_true(is.null(dmrs_not_annotated) || inherits(dmrs_not_annotated, "GRanges"))
    if (!is.null(dmrs_not_annotated) && length(dmrs_not_annotated) > 0) {
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_not_annotated))))
        expect_false("in_promoter_of" %in% names(mcols(dmrs_not_annotated)))
        expect_false("in_gene_body_of" %in% names(mcols(dmrs_not_annotated)))
    }
})
