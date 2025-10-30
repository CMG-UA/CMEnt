test_that("plotDMRsBioCircos creates a BioCircos plot", {
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))
    load(system.file("data/array_type.rda", package = "DMRsegal"))
    
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))
    
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }
    
    dmrs_subset <- dmrs[1:min(5, length(dmrs))]
    
    plot_obj <- plotDMRsBioCircos(
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

test_that("plotDMRsBioCircos works with interactions", {
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))
    load(system.file("data/array_type.rda", package = "DMRsegal"))
    
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))
    
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }
    
    dmrs_subset <- dmrs[1:min(5, length(dmrs))]
    
    plot_obj <- plotDMRsBioCircos(
        dmrs = dmrs_subset,
        beta = beta,
        pheno = pheno,
        array = array_type,
        genome = "hg19",
        sample_group_col = "Sample_Group",
        plot_interactions = TRUE,
        min_fdr = 0.05
    )
    
    expect_true(!is.null(plot_obj))
})

test_that("plotDMRsBioCircos handles BetaHandler input", {
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))
    load(system.file("data/array_type.rda", package = "DMRsegal"))
    
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))
    
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }
    
    dmrs_subset <- dmrs[1:min(3, length(dmrs))]
    
    beta_handler <- getBetaHandler(beta, array = array_type, genome = "hg19")
    
    plot_obj <- plotDMRsBioCircos(
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

test_that("plotDMRsBioCircos handles data frame DMRs input", {
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
    
    plot_obj <- plotDMRsBioCircos(
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

test_that("plotDMRsBioCircos validates inputs", {
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))
    
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))
    
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }
    
    dmrs_subset <- dmrs[1:min(3, length(dmrs))]
    
    expect_error(
        plotDMRsBioCircos(
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
        plotDMRsBioCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = bad_pheno,
            genome = "hg19",
            sample_group_col = "Sample_Group"
        ),
        "not found in pheno"
    )
})
