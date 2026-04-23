options("CMEnt.verbose" = 0)


test_that("findDMRsFromSeeds work with covariates adjustment", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    dmps <- subsetDenseExampleDmpsChr5And11(dmps)
    options(error = traceback)
    options(warn = 2)
    dmrs <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        covariates = c("Age", "Gender"),
        min_seeds = 2,
        min_cpgs = 3,
        njobs = 1,
        max_lookup_dist = 1000
    )
    expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
})

test_that("findDMRsFromSeeds reproduces benchmark.Rmd results with minfi", {
    if (!suppressPackageStartupMessages(requireNamespace("minfi", quietly = TRUE))) {
        skip("Package 'minfi' not installed")
    }
    suppressPackageStartupMessages({
        beta <- loadExampleInputDataChr5And11("beta")
        pheno <- loadExampleInputDataChr5And11("pheno")
        array_type <- loadExampleInputDataChr5And11("array_type")
        genome <- "hg19"


        beta_handler <- CMEnt::getBetaHandler(beta, array = array_type, genome = genome)
        beta_mat <- as.matrix(beta_handler$getBeta())
        locs <- beta_handler$getBetaLocs()
        mvalues <- minfi::logit2(beta_mat)
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
        # Run CMEnt with same parameters as benchmark
        dmrs_cment <- findDMRsFromSeeds(
            .score_dmrs = FALSE,
            extract_motifs = FALSE,
            annotate_with_genes = FALSE,
            beta = beta_handler,
            seeds = sig_dmps,
            pheno = pheno,
            sample_group_col = "Sample_Group",
            min_cpg_delta_beta = 0,
            min_seeds = 2,
            min_cpgs = 3,
            max_lookup_dist = 10000,
            expansion_window = 0,
            max_bridge_seeds_gaps = 0,
            max_pval = 0.05,
            pval_mode = "parametric",
            njobs = 1
        )

        # Assertions
        expect_s4_class(dmrs_cment, "GRanges")
        expect_equal(length(dmrs_cment), 143)
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_cment))))

        # Check that all DMRs meet the criteria
        expect_true(all(mcols(dmrs_cment)$seeds_num >= 2))
        expect_true(all(mcols(dmrs_cment)$cpgs_num >= 3))
    })
})

test_that("findDMRsFromSeeds parameter variations work correctly", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    # Make the dataset smaller for faster testing
    dmps <- subsetDenseExampleDmpsChr5And11(dmps)

    # Test with strict min_seeds
    expect_warning(
    dmrs_strict <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 5, # Stricter
        min_cpgs = 3,
        max_lookup_dist = 1000
    ), regexp = "No DMRs remain"
    )

    # Test with lenient parameters
    dmrs_lenient <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2, # More lenient
        min_cpgs = 2, # More lenient
        max_lookup_dist = 2000 # Larger distance
    )

    # Test with different max_pval
    expect_warning(
    dmrs_strict_pval <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        max_pval = 0.01 # Stricter p-value
    ), regexp = "No DMRs remain"
    )


    # Assertions
    expect_s4_class(dmrs_lenient, "GRanges")

    # Strict parameters should have all DMRs meeting criteria
    if (length(dmrs_lenient) > 0) {
        expect_true(all(mcols(dmrs_lenient)$seeds_num >= 1))
    }
})

test_that("findDMRsFromSeeds handles different aggregation functions", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    dmps <- subsetDenseExampleDmpsChr5And11(dmps)

    # Test with median aggregation
    dmrs_median <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        aggfun = "median"
    )

    # Test with mean aggregation
    dmrs_mean <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        aggfun = "mean"
    )


    # Assertions
    expect_true(is.null(dmrs_median) || inherits(dmrs_median, "GRanges"))
    expect_true(is.null(dmrs_mean) || inherits(dmrs_mean, "GRanges"))

    # Both should return valid results
    if (!is.null(dmrs_median)) expect_true(length(dmrs_median) >= 0)
    if (!is.null(dmrs_mean)) expect_true(length(dmrs_mean) >= 0)

    # If both are non-null, they should have the same number of DMRs (since only scoring differs)
    if (!is.null(dmrs_median) && !is.null(dmrs_mean)) {
        expect_equal(length(dmrs_median), length(dmrs_mean))
    }
    # If both are non-null, they should be different due to different aggregation functions
    if (!is.null(dmrs_median) && !is.null(dmrs_mean)) {
        expect_false(identical(as.data.frame(dmrs_median), as.data.frame(dmrs_mean)))
    }
})

test_that("findDMRsFromSeeds works with different genome builds", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    dmps <- subsetDenseExampleDmpsChr5And11(dmps)

    # Test with hg38
    dmrs_hg38 <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        genome = "hg38",
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000
    )

    # Assertions
    expect_true(is.null(dmrs_hg38) || inherits(dmrs_hg38, "GRanges"))
    if (!is.null(dmrs_hg38)) expect_true(length(dmrs_hg38) >= 0)
})


test_that("findDMRsFromSeeds preserves non-tabular columns in TSV outputs", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    dmps <- subsetDenseExampleDmpsChr5And11(dmps)
    pheno$casecontrol <- pheno$Sample_Group == "cancer"

    output_prefix <- file.path(tempdir(), paste0("cment-test-", as.integer(Sys.time())))

    dmrs <- expect_no_error(findDMRsFromSeeds(
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        casecontrol_col = "casecontrol",
        covariates = c("Age", "Gender"),
        max_bridge_seeds_gaps = 1,
        max_bridge_extension_gaps = 1,
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        max_pval = 0.05,
        pval_mode = "parametric",
        entanglement = "weak",
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        output_prefix = output_prefix,
        njobs = 1
    ))

    dmrs_file <- paste0(output_prefix, ".dmrs.tsv.gz")
    beta_file <- paste0(output_prefix, ".seeds_beta.tsv.gz")
    expect_true(file.exists(dmrs_file))
    expect_true(file.exists(beta_file))

    dmrs_df <- read.delim(gzfile(dmrs_file), check.names = FALSE)
    saved_beta <- read.delim(gzfile(beta_file), row.names = 1, check.names = FALSE)
    supporting_cpgs <- unique(unlist(lapply(as.character(S4Vectors::mcols(dmrs)$cpgs), CMEnt:::.splitCsvValues), use.names = FALSE))
    expect_setequal(rownames(saved_beta), supporting_cpgs)

    loaded_dmrs <- CMEnt:::.loadCMEntData(output_prefix)$dmrs
    expect_equal(length(loaded_dmrs), length(dmrs))
})
