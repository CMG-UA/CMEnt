test_that("plotDMRsCircos creates a circos plot", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[seq_len(min(5, length(dmrs))), drop = FALSE]
    options("DMRsegal.verbose" = 2)
    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            min_similarity = 0.8
        )
    )
})

test_that("plotDMRsCircos works with interactions", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[as.character(dmrs@seqnames) %in% c("chr5", "chr11"), , drop = FALSE]

    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            min_similarity = 0.8
        )
    )
})

test_that("plotDMRsCircos handles BetaHandler input", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))

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
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[seq_len(min(3, length(dmrs)))]
    dmrs_df <- convertToDataFrame(dmrs_subset)

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
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))

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

test_that("plotDMRsCircos supports chromosome and region filters", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[seq_len(min(12, length(dmrs)))]
    target_chr <- as.character(GenomicRanges::seqnames(dmrs_subset)[1])
    target_start <- GenomicRanges::start(dmrs_subset)[1]
    target_end <- GenomicRanges::end(dmrs_subset)[1]
    target_region <- sprintf("%s:%d-%d", target_chr, target_start, target_end)

    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            chromosomes = target_chr,
            region = target_region,
            max_components = 5
        )
    )

    expect_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            region = "chr5-100-200"
        ),
        "chr:start-end"
    )
})
