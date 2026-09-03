options("CMEnt.verbose" = 0)

loadExampleInputDataChr5And11()

test_that("scoreDMRs adds score column to DMRs", {

    example_output_path <- system.file("extdata", "example_outputChr5And11.rds", package = "CMEnt")
    dmrs <- readRDS(example_output_path)
    S4Vectors::mcols(dmrs)$score <- NULL

    expect_false("score" %in% names(S4Vectors::mcols(dmrs)))

    scoring_dmrs <- scoreDMRs(
        dmrs = dmrs,
        beta = beta,
        pheno = pheno,
        array = "450K",
        genome = "hg19",
        sample_group_col = "Sample_Group"
    )

    expect_true("score" %in% names(S4Vectors::mcols(scoring_dmrs)))
    expect_true(all(S4Vectors::mcols(scoring_dmrs)$score >= 0))
    expect_true(all(S4Vectors::mcols(scoring_dmrs)$score <= 1))
    expect_true("score_smoothed" %in% names(S4Vectors::mcols(scoring_dmrs)))
    expect_true("segment_id" %in% names(S4Vectors::mcols(scoring_dmrs)))
    expect_true("segment_slope" %in% names(S4Vectors::mcols(scoring_dmrs)))
    expect_true("block_id" %in% names(S4Vectors::mcols(scoring_dmrs)))
    expect_equal(length(dmrs), length(scoring_dmrs))
})

test_that("scoreDMRs works when called from buildDMRs with score_dmrs=TRUE", {
    skip_if_not_integration_tests()

    dmrs <- buildDMRs(
        score_dmrs = TRUE,
        annotate_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        array = "450K",
        genome = "hg19",
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_sites = 3,
        max_lookup_dist = 1000
    )

    expect_true("score" %in% names(S4Vectors::mcols(dmrs)))
    expect_true(all(S4Vectors::mcols(dmrs)$score >= 0))
    expect_true(all(S4Vectors::mcols(dmrs)$score <= 1))
})

test_that("scoreDMRs accepts a BetaHandler input", {
    skip_if_missing_array_annotation(array = "450K", genome = "hg19")
    dmrs <- readRDS(system.file("extdata", "example_outputChr5And11.rds", package = "CMEnt"))
    beta_handler <- getBetaHandler(beta, array = "450K", genome = "hg19")

    scoring_dmrs <- scoreDMRs(
        dmrs = dmrs,
        beta = beta_handler,
        pheno = pheno,
        sample_group_col = "Sample_Group"
    )

    expect_true("score" %in% names(S4Vectors::mcols(scoring_dmrs)))
    expect_true("cv_accuracy" %in% names(S4Vectors::mcols(scoring_dmrs)))
})

test_that("scoreDMRs can reuse preloaded DMR beta values", {
    dmrs <- readRDS(system.file("extdata", "example_outputChr5And11.rds", package = "CMEnt"))[1:2]
    dmr_sites <- strsplit(as.character(S4Vectors::mcols(dmrs)$sites), ",", fixed = TRUE)
    required_sites <- unique(unlist(dmr_sites, use.names = FALSE))
    preloaded_beta <- beta[required_sites, , drop = FALSE]

    scoring_dmrs <- scoreDMRs(
        dmrs = dmrs,
        beta = beta,
        pheno = pheno,
        array = "450K",
        genome = "hg19",
        sample_group_col = "Sample_Group",
        njobs = 1,
        .dmr_beta = preloaded_beta
    )

    expect_true("score" %in% names(S4Vectors::mcols(scoring_dmrs)))
    expect_error(
        scoreDMRs(
            dmrs = dmrs,
            beta = beta,
            pheno = pheno,
            array = "450K",
            genome = "hg19",
            sample_group_col = "Sample_Group",
            njobs = 1,
            .dmr_beta = preloaded_beta[-1, , drop = FALSE]
        ),
        "missing DMR sites"
    )
})

test_that("scoreDMRs materialization chunk rows are derived from available RAM", {
    chunk_rows <- CMEnt:::.scoringMaterializationChunkRows(
        n_samples = 6,
        njobs = 2,
        n_rows = 5000,
        available_ram_bytes = 1024^2
    )
    expect_equal(chunk_rows, floor(0.9 * 1024^2 / (2 * 6 * 8 * 12)))
    expect_equal(
        CMEnt:::.scoringMaterializationChunkRows(6, 2, 10, available_ram_bytes = 1024^2),
        10L
    )
})

test_that("scoreDMRs materialization chunks respect unique transformed row budget", {
    chunks <- CMEnt:::.scoreDMRIndexChunks(
        list(
            c(1L, 2L),
            c(2L, 3L),
            c(3L, 4L),
            c(10L, 11L, 12L),
            c(12L)
        ),
        max_unique_rows = 3L
    )

    expect_identical(chunks, list(1:2, 3L, 4:5))
    expect_identical(
        CMEnt:::.scoreDMRIndexChunks(list(1:5, 5L), max_unique_rows = 3L),
        list(1L, 2L)
    )
})

test_that("scoreDMRs internal memory job budget tightens materialization chunks", {
    dmrs <- readRDS(system.file("extdata", "example_outputChr5And11.rds", package = "CMEnt"))[1:4]
    dmr_sites <- strsplit(as.character(S4Vectors::mcols(dmrs)$sites), ",", fixed = TRUE)
    required_sites <- unique(unlist(dmr_sites, use.names = FALSE))
    preloaded_beta <- beta[required_sites, , drop = FALSE]

    local_mocked_bindings(.availableRamBytes = function(default_gb = 2) 1024^2)
    expect_message(
        scoreDMRs(
            dmrs = dmrs,
            beta = beta,
            pheno = pheno,
            array = "450K",
            genome = "hg19",
            sample_group_col = "Sample_Group",
            njobs = 1,
            verbose = 3,
            show_progress = FALSE,
            .dmr_beta = preloaded_beta,
            .memory_njobs = 4L
        ),
        "4 budgeted job\\(s\\)"
    )
})

test_that("scoreDMRs infers beta row names for DMRs without sites metadata", {
    old_options <- options(CMEnt.scoring_nfold = 2)
    on.exit(options(old_options), add = TRUE)
    set.seed(1)

    site_ids <- paste0("cg", seq_len(5))
    beta_small <- matrix(
        c(
            0.20, 0.22, 0.18, 0.21,
            0.10, 0.12, 0.88, 0.90,
            0.15, 0.17, 0.85, 0.87,
            0.30, 0.28, 0.31, 0.29,
            0.40, 0.41, 0.39, 0.42
        ),
        nrow = length(site_ids),
        byrow = TRUE,
        dimnames = list(site_ids, paste0("s", seq_len(4)))
    )
    locs <- data.frame(
        chr = "chr1",
        start = seq(10, 50, by = 10),
        end = seq(10, 50, by = 10),
        row.names = site_ids
    )
    pheno_small <- data.frame(
        Sample_Group = c("control", "control", "case", "case"),
        row.names = colnames(beta_small)
    )
    dmrs_no_sites <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 2),
        ranges = IRanges::IRanges(start = c(20, 40), end = c(30, 50))
    )
    GenomeInfoDb::genome(dmrs_no_sites) <- "hg19"

    scoring_dmrs <- scoreDMRs(
        dmrs = dmrs_no_sites,
        beta = beta_small,
        pheno = pheno_small,
        sorted_locs = locs,
        genome = "hg19",
        sample_group_col = "Sample_Group",
        njobs = 1
    )

    expect_setequal(as.character(S4Vectors::mcols(scoring_dmrs)$sites), c("cg2,cg3", "cg4,cg5"))
    expect_true("score" %in% names(S4Vectors::mcols(scoring_dmrs)))
})

test_that("ignored_sample_groups affects detection only, not downstream scoring", {
    skip_if_not_integration_tests()

    dmrs <- buildDMRs(
        score_dmrs = TRUE,
        annotate_with_genes = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        array = "450K",
        genome = "hg19",
        sample_group_col = "Sample_Group",
        ignored_sample_groups = "cancer",
        min_seeds = 2,
        min_sites = 2,
        max_lookup_dist = 2000,
        njobs = 1
    )

    expect_s4_class(dmrs, "GRanges")
    expect_true("score" %in% names(S4Vectors::mcols(dmrs)))
    expect_true("cv_accuracy" %in% names(S4Vectors::mcols(dmrs)))
})

test_that("scoreDMRs score values are meaningful", {

    example_output_path <- system.file("extdata", "example_outputChr5And11.rds", package = "CMEnt")
    dmrs <- readRDS(example_output_path)
    S4Vectors::mcols(dmrs)$score <- NULL

    scoring_dmrs <- scoreDMRs(
        dmrs = dmrs,
        beta = beta,
        pheno = pheno,
        array = "450K",
        genome = "hg19",
        sample_group_col = "Sample_Group"
    )

    expect_true(length(unique(S4Vectors::mcols(scoring_dmrs)$score)) > 1)
    expect_true(any(S4Vectors::mcols(scoring_dmrs)$score > 0.5))
})

test_that("scoreDMRs works with different nfold values", {

    example_output_path <- system.file("extdata", "example_outputChr5And11.rds", package = "CMEnt")
    dmrs <- readRDS(example_output_path)
    S4Vectors::mcols(dmrs)$score <- NULL

    options(CMEnt.scoring_nfold = 3)
    scoring_dmrs_3fold <- scoreDMRs(
        dmrs = dmrs,
        beta = beta,
        pheno = pheno,
        array = "450K",
        genome = "hg19",
        sample_group_col = "Sample_Group"
    )

    options(CMEnt.scoring_nfold = 5)
    scoring_dmrs_5fold <- scoreDMRs(
        dmrs = dmrs,
        beta = beta,
        pheno = pheno,
        array = "450K",
        genome = "hg19",
        sample_group_col = "Sample_Group"
    )

    expect_true("score" %in% names(S4Vectors::mcols(scoring_dmrs_3fold)))
    expect_true("score" %in% names(S4Vectors::mcols(scoring_dmrs_5fold)))
    expect_equal(length(scoring_dmrs_3fold), length(scoring_dmrs_5fold))
})

test_that("scoreDMRs returns DMRs ordered by p-value", {

    example_output_path <- system.file("extdata", "example_outputChr5And11.rds", package = "CMEnt")
    dmrs <- readRDS(example_output_path)
    S4Vectors::mcols(dmrs)$score <- NULL

    scoring_dmrs <- scoreDMRs(
        dmrs = dmrs,
        beta = beta,
        pheno = pheno,
        array = "450K",
        genome = "hg19",
        sample_group_col = "Sample_Group"
    )

    pvals <- S4Vectors::mcols(scoring_dmrs)$pval
    expect_true(all(diff(pvals) >= 0))
})

test_that(".assignDMRBlocksFromScores detects increase-plateau-decrease blocks", {
    set.seed(42)
    x <- seq(100, by = 50, length.out = 30)
    trend <- c(
        seq(0.55, 0.72, length.out = 10),
        rep(0.72, 10),
        seq(0.72, 0.56, length.out = 10)
    )
    y <- trend + stats::rnorm(length(trend), mean = 0, sd = 0.005)

    dmrs <- GenomicRanges::GRanges(
        seqnames = rep("chr1", length(x)),
        ranges = IRanges::IRanges(start = x, width = 30)
    )
    S4Vectors::mcols(dmrs)$score <- y

    with_blocks <- CMEnt:::.assignDMRBlocksFromScores(dmrs)
    block_ids <- S4Vectors::mcols(with_blocks)$block_id

    expect_true(any(!is.na(block_ids)))
    expect_true(length(unique(stats::na.omit(block_ids))) >= 1)
    expect_true(all(is.finite(S4Vectors::mcols(with_blocks)$score_smoothed)))
})

test_that(".assignDMRBlocksFromScores keeps singleton chromosome entries as NA", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(100, 150, 500), width = c(25, 25, 25))
    )
    S4Vectors::mcols(dmrs)$score <- c(0.6, 0.61, 0.9)

    with_blocks <- CMEnt:::.assignDMRBlocksFromScores(dmrs)
    expect_true(is.na(S4Vectors::mcols(with_blocks)$block_id[3]))
})

test_that(".computeChromosomeGapThreshold supports adaptive, fixed, and none modes", {
    dense_x <- seq(1, by = 10000, length.out = 100)
    sparse_x <- seq(1, by = 500000, length.out = 100)

    thr_dense <- CMEnt:::.computeChromosomeGapThreshold(dense_x, mode = "adaptive")
    thr_sparse <- CMEnt:::.computeChromosomeGapThreshold(sparse_x, mode = "adaptive")
    expect_true(thr_sparse > thr_dense)

    thr_fixed <- CMEnt:::.computeChromosomeGapThreshold(dense_x, mode = "fixed", fixed_bp = 123456)
    expect_equal(thr_fixed, 123456)

    thr_none <- CMEnt:::.computeChromosomeGapThreshold(dense_x, mode = "none")
    expect_true(is.infinite(thr_none))
})

test_that("scoreDMRs uses documented adaptive block gap defaults", {
    expect_equal(formals(scoreDMRs)$block_gap_min_bp, 250000)
    expect_equal(formals(scoreDMRs)$block_gap_max_bp, 5000000)
    expect_equal(formals(CMEnt:::.assignDMRBlocksFromScores)$block_gap_min_bp, 250000)
    expect_equal(formals(CMEnt:::.assignDMRBlocksFromScores)$block_gap_max_bp, 5000000)
})

test_that("distance-constrained mode removes over-bridging blocks", {
    set.seed(1)
    x <- c(
        seq(1e6, by = 5e4, length.out = 10),
        seq(2e7, by = 5e4, length.out = 10),
        seq(4e7, by = 5e4, length.out = 10)
    )
    y <- c(
        seq(0.58, 0.74, length.out = 15),
        seq(0.74, 0.58, length.out = 15)
    ) + stats::rnorm(30, mean = 0, sd = 0.002)

    dmrs <- GenomicRanges::GRanges(
        seqnames = rep("chr1", length(x)),
        ranges = IRanges::IRanges(start = x, width = 100)
    )
    S4Vectors::mcols(dmrs)$score <- y

    no_gap <- CMEnt:::.assignDMRBlocksFromScores(dmrs, block_gap_mode = "none")
    fixed_gap <- CMEnt:::.assignDMRBlocksFromScores(
        dmrs,
        block_gap_mode = "fixed",
        block_gap_fixed_bp = 1e6
    )
    expect_true(any(!is.na(S4Vectors::mcols(no_gap)$block_id)))
    expect_true(all(is.na(S4Vectors::mcols(fixed_gap)$block_id)))

    details <- CMEnt:::.computeDMRBlockFormationForChromosome(
        chr = "chr1",
        chr_idx = seq_along(x),
        x_chr = x,
        y_chr = y,
        block_gap_mode = "fixed",
        block_gap_fixed_bp = 1e6
    )
    expect_true(nrow(details$split_events_df) >= 1)
})

test_that("none mode matches effectively unlimited fixed threshold", {
    set.seed(1)
    x <- c(
        seq(1e6, by = 5e4, length.out = 10),
        seq(2e7, by = 5e4, length.out = 10),
        seq(4e7, by = 5e4, length.out = 10)
    )
    y <- c(
        seq(0.58, 0.74, length.out = 15),
        seq(0.74, 0.58, length.out = 15)
    ) + stats::rnorm(30, mean = 0, sd = 0.002)

    dmrs <- GenomicRanges::GRanges(
        seqnames = rep("chr1", length(x)),
        ranges = IRanges::IRanges(start = x, width = 100)
    )
    S4Vectors::mcols(dmrs)$score <- y

    none_mode <- CMEnt:::.assignDMRBlocksFromScores(dmrs, block_gap_mode = "none")
    huge_fixed <- CMEnt:::.assignDMRBlocksFromScores(
        dmrs,
        block_gap_mode = "fixed",
        block_gap_fixed_bp = 1e12
    )
    expect_identical(
        as.character(S4Vectors::mcols(none_mode)$block_id),
        as.character(S4Vectors::mcols(huge_fixed)$block_id)
    )
})

test_that("adaptive mode enforces maximum internal gap per chromosome block", {
    dmrs_path <- system.file("extdata", "example_outputChr5And11.rds", package = "CMEnt")
    if (!nzchar(dmrs_path) || !file.exists(dmrs_path)) {
        skip("Benchmark DMRs not available")
    }
    dmrs <- readRDS(dmrs_path)
    if (length(dmrs) == 0 || !"score" %in% colnames(S4Vectors::mcols(dmrs))) {
        skip("Benchmark DMRs with score not available")
    }

    with_blocks <- CMEnt:::.assignDMRBlocksFromScores(dmrs, block_gap_mode = "adaptive")
    midpoints <- floor((GenomicRanges::start(with_blocks) + GenomicRanges::end(with_blocks)) / 2)
    chrs <- as.character(GenomicRanges::seqnames(with_blocks))
    block_ids <- as.character(S4Vectors::mcols(with_blocks)$block_id)

    for (chr in unique(chrs)) {
        chr_idx <- which(chrs == chr)
        threshold <- CMEnt:::.computeChromosomeGapThreshold(
            midpoints[chr_idx],
            mode = "adaptive",
            quantile = 0.95,
            multiplier = 1.5,
            min_bp = 250000,
            max_bp = 5000000
        )
        chr_blocks <- unique(stats::na.omit(block_ids[chr_idx]))
        for (bid in chr_blocks) {
            idx <- chr_idx[which(block_ids[chr_idx] == bid)]
            block_x <- sort(midpoints[idx])
            max_gap <- if (length(block_x) > 1) max(diff(block_x)) else 0
            expect_lte(max_gap, threshold)
        }
    }
})
