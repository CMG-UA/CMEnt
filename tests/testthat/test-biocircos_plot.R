test_that("plotDMRsCircos creates a BioCircos plot", {
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))
    load(system.file("data/array_type.rda", package = "DMRsegal"))

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[seq_len(min(5, length(dmrs)))]

    plot_obj <- plotDMRsCircos(
        dmrs = dmrs_subset,
        beta = beta,
        pheno = pheno,
        array = array_type,
        genome = "hg19",
        sample_group_col = "Sample_Group",
        plot_interactions = FALSE
    )

    expect_true(!is.null(plot_obj))
})

test_that("plotDMRsCircos works with interactions", {
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))
    load(system.file("data/array_type.rda", package = "DMRsegal"))

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[1:min(5, length(dmrs))]

    plot_obj <- plotDMRsCircos(
        dmrs = dmrs_subset,
        beta = beta,
        pheno = pheno,
        array = array_type,
        genome = "hg19",
        sample_group_col = "Sample_Group",
        plot_interactions = TRUE,
        min_sim = 0.7
    )

    expect_true(!is.null(plot_obj))
})

test_that("plotDMRsCircos handles BetaHandler input", {
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))
    load(system.file("data/array_type.rda", package = "DMRsegal"))

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[1:min(3, length(dmrs))]

    beta_handler <- getBetaHandler(beta, array = array_type, genome = "hg19")

    plot_obj <- plotDMRsCircos(
        dmrs = dmrs_subset,
        beta = beta_handler,
        pheno = pheno,
        array = array_type,
        genome = "hg19",
        sample_group_col = "Sample_Group",
        plot_interactions = FALSE
    )

    expect_true(!is.null(plot_obj))
})

test_that("plotDMRsCircos handles data frame DMRs input", {
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))
    load(system.file("data/array_type.rda", package = "DMRsegal"))

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[1:min(3, length(dmrs))]
    dmrs_df <- as.data.frame(dmrs_subset)
    colnames(dmrs_df)[1] <- "chr"

    plot_obj <- plotDMRsCircos(
        dmrs = dmrs_df,
        beta = beta,
        pheno = pheno,
        array = array_type,
        genome = "hg19",
        sample_group_col = "Sample_Group",
        plot_interactions = FALSE
    )

    expect_true(!is.null(plot_obj))
})

test_that("plotDMRsCircos validates inputs", {
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[1:min(3, length(dmrs))]

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
