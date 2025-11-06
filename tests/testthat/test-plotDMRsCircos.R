test_that("plotDMRsCircos creates a circos plot", {
    beta <- loadExampleInputData("beta")
    pheno <- loadExampleInputData("pheno")
    array_type <- loadExampleInputData("array_type")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[seq_len(min(5, length(dmrs)))]
    options("DMRsegal.verbose" = 2)
    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group"
        )
    )
})

test_that("plotDMRsCircos works with interactions", {
    beta <- loadExampleInputData("beta")
    pheno <- loadExampleInputData("pheno")
    array_type <- loadExampleInputData("array_type")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs

    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            min_sim = 0.8
        )
    )
})

test_that("plotDMRsCircos handles BetaHandler input", {
    beta <- loadExampleInputData("beta")
    pheno <- loadExampleInputData("pheno")
    array_type <- loadExampleInputData("array_type")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }
    options("DMRsegal.verbose" = 3)
    dmrs_subset <- dmrs[seq_len(min(3, length(dmrs)))]

    beta_handler <- getBetaHandler(beta, array = array_type, genome = "hg19")

    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta_handler,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group"
        )
    )
})

test_that("plotDMRsCircos handles data frame DMRs input", {
    beta <- loadExampleInputData("beta")
    pheno <- loadExampleInputData("pheno")
    array_type <- loadExampleInputData("array_type")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[seq_len(min(3, length(dmrs)))]
    dmrs_df <- as.data.frame(dmrs_subset)
    colnames(dmrs_df)[1] <- "chr"

    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs_df,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group"
        )
    )
})

test_that("plotDMRsCircos validates inputs", {
    beta <- loadExampleInputData("beta")
    pheno <- loadExampleInputData("pheno")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[seq_len(min(3, length(dmrs)))]

    expect_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = "nonexistent_file.txt",
            pheno = pheno,
            genome = "hg19"
        ),
        "does not exist"
    )

    bad_pheno <- pheno
    colnames(bad_pheno) <- c("col1", "col2", "col3")

    expect_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = bad_pheno,
            genome = "hg19",
            sample_group_col = "Sample_Group"
        ),
        "not found in pheno"
    )
})
