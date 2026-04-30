options("CMEnt.verbose" = 0)

.makeOrderingInput <- function(locs) {
    site_ids <- rownames(locs)
    n_site <- length(site_ids)
    beta <- matrix(
        seq(0.1, by = 0.01, length.out = n_site * 6),
        nrow = n_site,
        byrow = TRUE
    )
    rownames(beta) <- site_ids
    colnames(beta) <- paste0("S", seq_len(6))

    pheno <- data.frame(
        Sample_Group = c("A", "A", "A", "B", "B", "B"),
        casecontrol = c(0, 0, 0, 1, 1, 1),
        row.names = colnames(beta),
        stringsAsFactors = FALSE
    )

    seeds <- data.frame(
        pval = rep(1e-6, n_site),
        row.names = site_ids,
        stringsAsFactors = FALSE
    )

    list(
        beta_handler = getBetaHandler(beta = beta, sorted_locs = locs),
        seeds = seeds,
        pheno = pheno
    )
}

test_that("findDMRsFromSeeds validates start ordering within chromosome", {
    locs <- data.frame(
        chr = c("chr1", "chr1", "chr1"),
        start = c(100, 300, 200),
        end = c(101, 301, 201),
        row.names = c("cg1", "cg2", "cg3"),
        stringsAsFactors = FALSE
    )
    input <- .makeOrderingInput(locs)

    expect_error(
        findDMRsFromSeeds(
            .score_dmrs = FALSE,
            extract_motifs = FALSE,
            annotate_with_genes = FALSE,
            beta = input$beta_handler,
            seeds = input$seeds,
            pheno = input$pheno,
            sample_group_col = "Sample_Group",
            casecontrol_col = "casecontrol",
            min_seeds = 2,
            min_sites = 2,
            ext_site_delta_beta = 0,
            max_lookup_dist = 1000,
            max_pval = 0.05,
            pval_mode = "parametric",
            njobs = 1
        ),
        regexp = "Beta locations are not sorted within chromosome chr1"
    )
})

test_that("findDMRsFromSeeds validates chromosome block contiguity", {
    locs <- data.frame(
        chr = c("chr1", "chr2", "chr1"),
        start = c(100, 100, 200),
        end = c(101, 101, 201),
        row.names = c("cg1", "cg2", "cg3"),
        stringsAsFactors = FALSE
    )
    input <- .makeOrderingInput(locs)

    expect_error(
        findDMRsFromSeeds(
            .score_dmrs = FALSE,
            extract_motifs = FALSE,
            annotate_with_genes = FALSE,
            beta = input$beta_handler,
            seeds = input$seeds,
            pheno = input$pheno,
            sample_group_col = "Sample_Group",
            casecontrol_col = "casecontrol",
            min_seeds = 2,
            min_sites = 2,
            ext_site_delta_beta = 0,
            max_lookup_dist = 1000,
            max_pval = 0.05,
            pval_mode = "parametric",
            njobs = 1
        ),
        regexp = "Beta locations are not grouped by chromosome: chr1 appears in multiple blocks"
    )
})
