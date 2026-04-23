suppressPackageStartupMessages({
    library(testthat)
    library(CMEnt)
    library(bsseq)
    library(GenomicRanges)
})
options("CMEnt.verbose" = 0)

test_that("BetaHandler can be created from BSseq object", {
    set.seed(123)
    n_loci <- 100
    n_samples <- 10
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 100, length.out = n_loci), width = 1)
    )
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = paste0("Sample", seq_len(n_samples))
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    expect_s3_class(beta_handler, "BetaHandler")
    expect_true(!is.null(beta_handler$sorted_locs))
    expect_equal(nrow(beta_handler$sorted_locs), n_loci)
})

test_that("BetaHandler can extract row names from BSseq object", {
    set.seed(123)
    n_loci <- 50
    n_samples <- 5
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 100, length.out = n_loci), width = 1)
    )
    cpg_names <- paste(seqnames(gr), start(gr), sep = ":")
    names(gr) <- cpg_names
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = paste0("Sample", seq_len(n_samples))
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    row_names <- beta_handler$getBetaRowNames()
    expect_equal(length(row_names), n_loci)
    expect_equal(row_names, cpg_names)
})

test_that("BetaHandler can extract column names from BSseq object", {
    set.seed(123)
    n_loci <- 50
    n_samples <- 5
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 100, length.out = n_loci), width = 1)
    )
    sample_names <- paste0("Sample", seq_len(n_samples))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = sample_names
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    col_names <- beta_handler$getBetaColNames()
    expect_equal(length(col_names), n_samples)
    expect_equal(col_names, sample_names)
})

test_that("BetaHandler can extract beta values from BSseq object", {
    set.seed(123)
    n_loci <- 50
    n_samples <- 5
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 100, length.out = n_loci), width = 1)
    )
    cpg_names <- paste(seqnames(gr), start(gr), sep = ":")
    names(gr) <- cpg_names
    sample_names <- paste0("Sample", seq_len(n_samples))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = sample_names
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    beta_values <- beta_handler$getBeta()
    expect_equal(dim(beta_values), c(n_loci, n_samples))
    expect_equal(rownames(beta_values), cpg_names)
    expect_equal(colnames(beta_values), sample_names)
    expected_beta <- met / cov
    rownames(expected_beta) <- cpg_names
    colnames(expected_beta) <- sample_names
    expect_equal(beta_values, expected_beta, tolerance = 1e-6)
})

test_that("BetaHandler can subset beta values from BSseq object by row names", {
    set.seed(123)
    n_loci <- 50
    n_samples <- 5
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 100, length.out = n_loci), width = 1)
    )
    cpg_names <- paste(seqnames(gr), start(gr), sep = ":")
    names(gr) <- cpg_names
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = paste0("Sample", seq_len(n_samples))
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    subset_cpgs <- cpg_names[c(1, 10, 20)]
    beta_subset <- beta_handler$getBeta(row_names = subset_cpgs)
    expect_equal(nrow(beta_subset), 3)
    expect_equal(rownames(beta_subset), subset_cpgs)
})

test_that("BetaHandler can subset beta values from BSseq object by column names", {
    set.seed(123)
    n_loci <- 50
    n_samples <- 10
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 100, length.out = n_loci), width = 1)
    )
    sample_names <- paste0("Sample", seq_len(n_samples))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = sample_names
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    subset_samples <- c("Sample2", "Sample5", "Sample8")
    beta_subset <- beta_handler$getBeta(col_names = subset_samples)
    expect_equal(ncol(beta_subset), 3)
    expect_equal(colnames(beta_subset), subset_samples)
})

test_that("BetaHandler extracts genomic locations from BSseq object", {
    set.seed(123)
    n_loci <- 50
    n_samples <- 5
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    start_positions <- seq(1000, by = 100, length.out = n_loci)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = start_positions, width = 1)
    )
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = paste0("Sample", seq_len(n_samples))
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    beta_locs <- beta_handler$getBetaLocs()
    expect_equal(nrow(beta_locs), n_loci)
    expect_true(all(c("chr", "start", "end") %in% colnames(beta_locs)))
    expect_equal(beta_locs$chr, rep("chr1", n_loci))
    expect_equal(beta_locs$start, start_positions)
})

test_that("BetaHandler sorts genomic locations from unsorted BSseq input", {
    set.seed(123)
    n_loci <- 8
    n_samples <- 4
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)

    unsorted_starts <- c(300, 100, 800, 200, 700, 400, 600, 500)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = unsorted_starts, width = 1)
    )
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = paste0("Sample", seq_len(n_samples))
    )

    beta_handler <- getBetaHandler(beta = bsseq_obj)
    beta_locs <- beta_handler$getBetaLocs()

    expect_equal(beta_locs$start, sort(unsorted_starts))
    expect_true(all(diff(beta_locs$start) >= 0))
    expect_equal(rownames(beta_locs), paste0("chr1:", sort(unsorted_starts)))
    expect_equal(rownames(beta_locs), beta_handler$getBetaRowNames())
})

test_that("BetaHandler handles missing row names in BSseq gracefully", {
    set.seed(123)
    n_loci <- 50
    n_samples <- 5
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 100, length.out = n_loci), width = 1)
    )
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = paste0("Sample", seq_len(n_samples))
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    row_names <- beta_handler$getBetaRowNames()
    expect_equal(length(row_names), n_loci)
    expect_true(all(grepl("chr1:", row_names)))
})

test_that("BetaHandler throws error when requesting non-existent CpGs from BSseq", {
    set.seed(123)
    n_loci <- 50
    n_samples <- 5
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 100, length.out = n_loci), width = 1)
    )
    names(gr) <- paste0("cg", seq_len(n_loci))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = paste0("Sample", seq_len(n_samples))
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    expect_error(
        beta_handler$getBeta(row_names = c("cg999", "cg1000")),
        "Requested CpG sites not found in BSseq object"
    )
})

test_that("BetaHandler allows missing CpGs from BSseq when allow_missing=TRUE", {
    set.seed(123)
    n_loci <- 50
    n_samples <- 5
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 100, length.out = n_loci), width = 1)
    )
    cpg_names <- paste(seqnames(gr), start(gr), sep = ":")
    names(gr) <- cpg_names
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = paste0("Sample", seq_len(n_samples))
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    beta_subset <- beta_handler$getBeta(
        row_names = c(cpg_names[[1]], cpg_names[[2]], "cg5"),
        allow_missing = TRUE
    )
    expect_equal(nrow(beta_subset), 2)
    expect_equal(rownames(beta_subset), c(cpg_names[[1]], cpg_names[[2]]))
})

test_that("BetaHandler subset returns compact BSseq handler with requested rows and columns", {
    set.seed(123)
    n_loci <- 60
    n_samples <- 6
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 50, length.out = n_loci), width = 1)
    )
    cpg_names <- paste(seqnames(gr), start(gr), sep = ":")
    names(gr) <- cpg_names
    sample_names <- paste0("Sample", seq_len(n_samples))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = sample_names
    )

    beta_handler <- getBetaHandler(beta = bsseq_obj)
    subset_rows <- cpg_names[10:20]
    subset_cols <- sample_names[c(2, 4, 6)]

    subset_handler <- beta_handler$subset(
        row_names = subset_rows,
        col_names = subset_cols
    )

    expect_s3_class(subset_handler, "BetaHandler")
    expect_equal(subset_handler$getBetaRowNames(), subset_rows)
    expect_equal(subset_handler$getBetaColNames(), subset_cols)
    expect_equal(
        subset_handler$getBeta(),
        beta_handler$getBeta(row_names = subset_rows, col_names = subset_cols),
        tolerance = 1e-8
    )
    expect_equal(rownames(subset_handler$getBetaLocs()), subset_rows)
})

test_that("BetaHandler subset supports numeric row indexing", {
    set.seed(99)
    n_loci <- 20
    n_samples <- 4
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr2", n_loci),
        ranges = IRanges(start = seq(5000, by = 25, length.out = n_loci), width = 1)
    )
    sample_names <- paste0("Sample", seq_len(n_samples))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = sample_names
    )

    beta_handler <- getBetaHandler(beta = bsseq_obj)
    subset_handler <- beta_handler$subset(row_names = c(2, 5, 7), col_names = sample_names[1:2])

    expect_equal(
        subset_handler$getBeta(),
        beta_handler$getBeta(row_names = c(2, 5, 7), col_names = sample_names[1:2]),
        tolerance = 1e-8
    )
})
