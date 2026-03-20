# Test suite for findDMRsFromSeeds function
library(testthat)


test_that("findDMRsFromSeeds work with covariates adjustment", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    options(error = traceback)
    options(warn = 2)
    options("DMRsegal.verbose" = 2)
    dmrs <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        covariates = c("Age", "Gender"),
        min_seeds = 2,
        min_cpgs = 3,
        njobs = 1,
        max_lookup_dist = 1000,
    )
    expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
})

test_that("findDMRsFromSeeds reproduces benchmark.Rmd results with minfi", {
    skip_if_not_installed("minfi")

    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")
    genome <- "hg19"


    beta_handler <- DMRsegal::getBetaHandler(beta, array = array_type, genome = genome)
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
    options("DMRsegal.verbose" = 2)
    # Run DMRsegal with same parameters as benchmark
    dmrs_segal <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = beta_handler,
        seeds = sig_dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_cpg_delta_beta = 0,
        adaptive_min_cpg_delta_beta = FALSE,
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
    expect_s4_class(dmrs_segal, "GRanges")
    expect_equal(length(dmrs_segal), 165)
    expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_segal))))

    # Check that all DMRs meet the criteria
    expect_true(all(mcols(dmrs_segal)$seeds_num >= 2))
    expect_true(all(mcols(dmrs_segal)$cpgs_num >= 3))
})

test_that("findDMRsFromSeeds parameter variations work correctly", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    # Test with strict min_seeds
    dmrs_strict <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 5, # Stricter
        min_cpgs = 3,
        max_lookup_dist = 1000,
        annotate_with_genes = FALSE
    )

    # Test with lenient parameters
    dmrs_lenient <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2, # More lenient
        min_cpgs = 2, # More lenient
        max_lookup_dist = 2000, # Larger distance
        annotate_with_genes = FALSE
    )

    # Test with different max_pval
    dmrs_strict_pval <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        max_pval = 0.01, # Stricter p-value
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
    }
})

test_that("findDMRsFromSeeds handles different aggregation functions", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    # Test with median aggregation
    dmrs_median <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        aggfun = "median",
        annotate_with_genes = FALSE
    )

    # Test with mean aggregation
    dmrs_mean <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        aggfun = "mean",
        annotate_with_genes = FALSE
    )


    # Assertions
    expect_true(is.null(dmrs_median) || inherits(dmrs_median, "GRanges"))
    expect_true(is.null(dmrs_mean) || inherits(dmrs_mean, "GRanges"))

    # Both should return valid results
    if (!is.null(dmrs_median)) expect_true(length(dmrs_median) >= 0)
    if (!is.null(dmrs_mean)) expect_true(length(dmrs_mean) >= 0)
})

test_that("findDMRsFromSeeds works with different genome builds", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    options("DMRsegal.use_annotation_cache" = FALSE)
    # Test with hg38
    dmrs_hg38 <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
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


test_that("findDMRsFromSeeds does not annotate DMRs when annotate_with_genes=FALSE", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    dmrs_not_annotated <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
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

test_that("findDMRsFromSeeds preserves non-tabular columns in TSV outputs", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    pheno$casecontrol <- pheno$Sample_Group == "cancer"

    output_prefix <- file.path(tempdir(), paste0("dmrsegal-test-", as.integer(Sys.time())))

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
        annotate_with_genes = FALSE,
        output_prefix = output_prefix,
        njobs = 1
    ))

    dmrs_file <- paste0(output_prefix, ".dmrs.tsv.gz")
    beta_file <- paste0(output_prefix, ".seeds_beta.tsv.gz")
    expect_true(file.exists(dmrs_file))
    expect_true(file.exists(beta_file))

    dmrs_df <- read.delim(gzfile(dmrs_file), check.names = FALSE)
    expect_true("pwm" %in% colnames(dmrs_df))
    saved_beta <- read.delim(gzfile(beta_file), row.names = 1, check.names = FALSE)
    supporting_cpgs <- unique(unlist(lapply(as.character(S4Vectors::mcols(dmrs)$cpgs), DMRsegal:::.splitCsvValues), use.names = FALSE))
    expect_setequal(rownames(saved_beta), supporting_cpgs)

    loaded_dmrs <- DMRsegal:::.loadDMRsegalData(output_prefix)$dmrs
    expect_true("pwm" %in% names(S4Vectors::mcols(loaded_dmrs)))
    expect_true(any(vapply(S4Vectors::mcols(loaded_dmrs)$pwm, is.matrix, logical(1))))
    expect_equal(length(loaded_dmrs), length(dmrs))
})
