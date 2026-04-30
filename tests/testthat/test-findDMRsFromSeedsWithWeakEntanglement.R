options("CMEnt.verbose" = 0)

test_that("findDMRsFromSeeds works with weak entanglement", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    dmps <- subsetDenseExampleDmpsChr5And11(dmps)

     dmrs_relaxed <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotate_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_sites = 3,
        max_lookup_dist = 1000,
        entanglement = "weak",
        pval_mode = "parametric"
    )


    dmrs_relaxed <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotate_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_sites = 3,
        max_lookup_dist = 1000,
        entanglement = "weak",
        pval_mode = "parametric",
    )

    expect_true(is.null(dmrs_relaxed) || inherits(dmrs_relaxed, "GRanges"))
    if (!is.null(dmrs_relaxed) && length(dmrs_relaxed) > 0) {
        expect_true(all(c("sites_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_relaxed))))
    }
})

test_that("weak entanglement produces more or equal DMRs than strong entanglement", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    dmps <- subsetDenseExampleDmpsChr5And11(dmps)
    handler <- getBetaHandler(beta, array = "450K", genome = "hg19")

    dmrs_strict <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotate_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_sites = 3,
        max_lookup_dist = 1000,
        entanglement = "strong",
        pval_mode = "parametric"
    )

    dmrs_relaxed <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        annotate_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_sites = 3,
        max_lookup_dist = 1000,
        entanglement = "weak",
        pval_mode = "parametric"
    )

    strict_count <- if (is.null(dmrs_strict)) 0 else length(dmrs_strict)
    relaxed_count <- if (is.null(dmrs_relaxed)) 0 else length(dmrs_relaxed)

    expect_true(relaxed_count >= strict_count)
})

test_that("entanglement parameter validates correctly", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    dmps <- subsetDenseExampleDmpsChr5And11(dmps)

    expect_error(
        findDMRsFromSeeds(
            .score_dmrs = FALSE,
            annotate_with_genes = FALSE,
            extract_motifs = FALSE,
            beta = beta,
            seeds = dmps,
            pheno = pheno,
            sample_group_col = "Sample_Group",
            entanglement = "invalid"
        ),
        "is not a prefix"
    )
})
