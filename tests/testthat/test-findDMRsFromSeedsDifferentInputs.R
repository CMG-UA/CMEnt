suppressPackageStartupMessages({
    library(testthat)
    library(DMRsegal)
    library(bsseq)
    library(GenomicRanges)
})

create_seeds_with_chr_pos <- function(seeds, beta_mat, locs) {
    seed_row_names <- rownames(seeds)
    seed_indices <- match(seed_row_names, rownames(beta_mat))
    valid_indices <- !is.na(seed_indices)

    seeds_subset <- seeds[valid_indices, , drop = FALSE]
    seeds_subset$ID <- paste0(as.character(locs$chr[seed_indices[valid_indices]]), ":", locs$start[seed_indices[valid_indices]])
    seeds_subset <- seeds_subset[!grepl("NA", seeds_subset$ID), , drop = FALSE]

    seeds_subset
}

test_that("findDMRsFromSeeds validates input parameters correctly", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    # Test missing required parameters
    expect_error(
        findDMRsFromSeeds(
            .score_dmrs = FALSE,
            beta = beta,
            seeds = NULL, # Missing
            pheno = pheno
        ),
        "seeds"
    )

    expect_error(
        findDMRsFromSeeds(
            .score_dmrs = FALSE,
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
            .score_dmrs = FALSE,
            beta = beta,
            seeds = dmps,
            pheno = pheno_wrong,
            sample_group_col = "Sample_Group",
            casecontrol_col = "casecontrol"
        )
    )
})

skip_if_covr_expensive("Skipping expensive input-mode integration tests under covr.")

test_that("findDMRsFromSeeds works with small beta file (in-memory loading)", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    # Run findDMRsFromSeeds with beta_in_mem_threshold_mb=500 (small file loaded in memory)
    options("DMRsegal.beta_in_mem_threshold_mb" = 500)
    expect_no_warning(
        dmrs <- findDMRsFromSeeds(
            .score_dmrs = FALSE,
            beta = beta,
            seeds = dmps,
            pheno = pheno,
            sample_group_col = "Sample_Group",
            min_seeds = 2,
            min_cpgs = 3,
            max_lookup_dist = 1000,
            verbose = 3
        )
    )
    # Assertions
    expect_true(!is.null(dmrs))
    expect_true(inherits(dmrs, "GRanges"))
    expect_true(length(dmrs) > 0)
    expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs))))
})

test_that("findDMRsFromSeeds works with large beta file (tabix indexing)", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    beta_file <- tempfile(fileext = ".tsv")
    withr::defer(unlink(beta_file))
    write.table(as.data.frame(beta), file = beta_file, sep = "\t", col.names = NA, quote = FALSE)
    sorted_beta_file <- sortBetaFileByCoordinates(beta_file, overwrite = TRUE)
    withr::defer(unlink(sorted_beta_file))
    options("DMRsegal.verbose" = 3)
    options("DMRsegal.use_tabix_cache" = FALSE)
    options("DMRsegal.beta_in_mem_threshold_mb" = 1)

    dmrs <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = sorted_beta_file,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        njobs = 1
    )

    expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
    if (!is.null(dmrs) && length(dmrs) > 0) {
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs))))
    }
})

test_that("findDMRsFromSeeds works with BSseq input", {
    set.seed(42)
    n_loci <- 200
    n_cases <- 5
    n_controls <- 5
    n_samples <- n_cases + n_controls
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    base_meth <- runif(n_loci, 0.3, 0.7)
    m_controls <- matrix(rbinom(n_loci * n_controls, size = cov[, 1:n_controls], prob = base_meth), ncol = n_controls)
    dmr_region_idx <- 50:70
    base_meth_dmr <- base_meth
    base_meth_dmr[dmr_region_idx] <- base_meth_dmr[dmr_region_idx] + 0.3
    base_meth_dmr[base_meth_dmr > 1] <- 0.95
    m_cases <- matrix(rbinom(n_loci * n_cases, size = cov[, (n_controls + 1):n_samples], prob = base_meth_dmr), ncol = n_cases)
    met <- cbind(m_controls, m_cases)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(10000, by = 50, length.out = n_loci), width = 1)
    )
    sample_names <- c(paste0("Control", seq_len(n_controls)), paste0("Case", seq_len(n_cases)))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = sample_names
    )
    pheno <- data.frame(
        Sample_ID = sample_names,
        Group = c(rep("Control", n_controls), rep("Case", n_cases)),
        stringsAsFactors = FALSE
    )
    rownames(pheno) <- sample_names
    seeds <- data.frame(
        cpg_id = paste0(seqnames(gr), ":", start(gr))[dmr_region_idx],
        pval = rep(0.001, length(dmr_region_idx))
    )
    expect_warning(
        dmrs <- findDMRsFromSeeds(
            beta = bsseq_obj,
            seeds_id_col = "cpg_id",
            seeds = seeds,
            pheno = pheno,
            sample_group_col = "Group",
            min_cpgs = 5,
            max_lookup_dist = 200,
            njobs = 1,
            annotate_with_genes = FALSE
        ),
        "No DMRs remain after filtering based on min_seeds."
    )
})


test_that("findDMRsFromSeeds works when tabix is not available", {
    skip_if_not_installed("mockery")
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    library(mockery)

    mock_convertBetaToTabix <- mock(NULL) # nolint

    stub(findDMRsFromSeeds, "convertBetaToTabix", mock_convertBetaToTabix)
    options("DMRsegal.use_tabix_cache" = FALSE)
    options("DMRsegal.verbose" = 2)
    options("DMRsegal.beta_in_mem_threshold_mb" = 0.1)
    dmrs <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 2,
        max_lookup_dist = 1000,
        annotate_with_genes = FALSE
    )

    expect_called(mock_convertBetaToTabix, 0)

    expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
    if (!is.null(dmrs) && length(dmrs) > 0) {
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs))))
    }
})


test_that("findDMRsFromSeeds works with minimal bed file", {
    skip_on_ci()
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")
    beta_handler <- getBetaHandler(beta, array = array_type, genome = "hg19")
    beta_mat <- as.matrix(beta_handler$getBeta())
    locs <- beta_handler$getBetaLocs()

    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))
    sample_cols <- rownames(pheno)

    bed_data <- data.frame(
        chrom = as.character(locs$chr),
        start = locs$start,
        stringsAsFactors = FALSE
    )
    for (sample in sample_cols) {
        bed_data[[sample]] <- beta_mat[, sample]
    }

    write.table(bed_data, file = bed_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

    dmps_with_chr_pos <- create_seeds_with_chr_pos(dmps, beta_mat, locs)
    dmrs <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = bed_file,
        seeds = dmps_with_chr_pos,
        seeds_id_col = "ID",
        pheno = pheno,
        sample_group_col = "Sample_Group",
        bed_provided = TRUE,
        bed_chrom_col = "chrom",
        bed_start_col = "start",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        verbose = 2
    )


    expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
    if (!is.null(dmrs) && length(dmrs) > 0) {
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs))))
    }
})
