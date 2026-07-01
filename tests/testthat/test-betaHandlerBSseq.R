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
    expect_null(beta_handler$sorted_locs)
    expect_equal(nrow(beta_handler$getBetaLocs()), n_loci)
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
    site_names <- paste(seqnames(gr), start(gr), sep = ":")
    names(gr) <- site_names
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = paste0("Sample", seq_len(n_samples))
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    row_names <- beta_handler$getBetaRowNames()
    expect_equal(length(row_names), n_loci)
    expect_equal(row_names, site_names)
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
    site_names <- paste(seqnames(gr), start(gr), sep = ":")
    names(gr) <- site_names
    sample_names <- paste0("Sample", seq_len(n_samples))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = sample_names
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    beta_values <- beta_handler$getBeta()
    expect_equal(dim(beta_values), c(n_loci, n_samples))
    expect_equal(rownames(beta_values), site_names)
    expect_equal(colnames(beta_values), sample_names)
    expected_beta <- met / cov
    rownames(expected_beta) <- site_names
    colnames(expected_beta) <- sample_names
    expect_equal(beta_values, expected_beta, tolerance = 1e-6)
})

test_that("BetaHandler returns NA for zero-coverage BSseq beta values", {
    cov <- matrix(c(10L, 0L, 8L, 0L), nrow = 2)
    met <- matrix(c(5L, 0L, 4L, 0L), nrow = 2)
    gr <- GRanges(
        seqnames = rep("chr1", 2),
        ranges = IRanges(start = c(1000L, 1100L), width = 1)
    )
    names(gr) <- paste(seqnames(gr), start(gr), sep = ":")
    bsseq_obj <- BSseq(
        M = met,
        Cov = cov,
        gr = gr,
        sampleNames = c("Sample1", "Sample2")
    )

    beta_values <- getBetaHandler(beta = bsseq_obj)$getBeta()

    expect_equal(beta_values[1, ], c(Sample1 = 0.5, Sample2 = 0.5), tolerance = 1e-8)
    expect_true(all(is.na(beta_values[2, ])))
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
    site_names <- paste(seqnames(gr), start(gr), sep = ":")
    names(gr) <- site_names
    bsseq_obj <- bsseq::BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = paste0("Sample", seq_len(n_samples))
    )
    subset_sites <- site_names[c(1, 10, 20)]

    beta_handler <- getBetaHandler(beta = bsseq_obj)
    beta_subset <- beta_handler$getBeta(row_names = subset_sites)
    expect_equal(nrow(beta_subset), 3)
    expect_equal(rownames(beta_subset), subset_sites)
    bsseq_obj_hdf5 <- HDF5Array::saveHDF5SummarizedExperiment(
        bsseq_obj,
        dir = tempfile("bsseq_hdf5_"),
        replace = TRUE
    )
    beta_handler_hdf5 <- getBetaHandler(beta = bsseq_obj_hdf5)
    beta_subset_hdf5 <- beta_handler_hdf5$getBeta(row_names = subset_sites)
    expect_s4_class(beta_subset_hdf5, "DelayedMatrix")
    expect_equal(nrow(beta_subset_hdf5), 3)
    expect_equal(rownames(beta_subset_hdf5), subset_sites)
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
    expect_equal(colnames(beta_locs), c("chr", "start"))
    expect_equal(as.character(beta_locs$chr), rep("chr1", n_loci))
    expect_equal(beta_locs$start, start_positions)
    # convert the object to disk-backed and check that the genomic locations are still correct
    bsseq_hdf5 <- HDF5Array::saveHDF5SummarizedExperiment(
        bsseq_obj,
        dir = tempfile("bsseq_hdf5_"),
        replace = TRUE
    )
    beta_handler_hdf5 <- getBetaHandler(beta = bsseq_hdf5)
    beta_locs_hdf5 <- beta_handler_hdf5$getBetaLocs()
    expect_equal(as.character(beta_locs_hdf5$chr), rep("chr1", n_loci))
    expect_equal(as.character(beta_locs_hdf5$start), as.character(start_positions))
    expect_equal(colnames(beta_locs_hdf5), c("chr", "start"))
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

test_that("BetaHandler sorts stranded BSseq input by position, not strand", {
    starts <- c(100L, 300L, 101L, 301L)
    gr <- GRanges(
        seqnames = rep("chr1", length(starts)),
        ranges = IRanges(start = starts, width = 1),
        strand = c("+", "+", "-", "-")
    )
    bsseq_obj <- BSseq(
        M = matrix(1L, nrow = length(starts), ncol = 2),
        Cov = matrix(10L, nrow = length(starts), ncol = 2),
        gr = gr,
        sampleNames = c("Sample1", "Sample2")
    )

    beta_locs <- getBetaHandler(beta = bsseq_obj)$getBetaLocs()

    expect_equal(beta_locs$start, sort(starts))
    expect_true(CMEnt:::.bsseqIsSorted(CMEnt:::.prepareBSseqForBetaHandler(bsseq_obj)))
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

test_that("BetaHandler throws error when requesting non-existent sites from BSseq", {
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
        "Requested site sites not found in BSseq object"
    )
})

test_that("BetaHandler allows missing sites from BSseq when allow_missing=TRUE", {
    set.seed(123)
    n_loci <- 50
    n_samples <- 5
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 100, length.out = n_loci), width = 1)
    )
    site_names <- paste(seqnames(gr), start(gr), sep = ":")
    names(gr) <- site_names
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = paste0("Sample", seq_len(n_samples))
    )
    beta_handler <- getBetaHandler(beta = bsseq_obj)
    beta_subset <- beta_handler$getBeta(
        row_names = c(site_names[[1]], site_names[[2]], "cg5"),
        allow_missing = TRUE
    )
    expect_equal(nrow(beta_subset), 2)
    expect_equal(rownames(beta_subset), c(site_names[[1]], site_names[[2]]))
})

test_that("BetaHandler subset materializes compact in-memory BSseq handler", {
    set.seed(123)
    n_loci <- 60
    n_samples <- 6
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 50, length.out = n_loci), width = 1)
    )
    site_names <- paste(seqnames(gr), start(gr), sep = ":")
    names(gr) <- site_names
    sample_names <- paste0("Sample", seq_len(n_samples))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = sample_names
    )

    beta_handler <- getBetaHandler(beta = bsseq_obj)
    subset_rows <- site_names[10:20]
    subset_cols <- sample_names[c(2, 4, 6)]

    subset_handler <- beta_handler$subset(
        row_names = subset_rows,
        col_names = subset_cols
    )

    expect_s3_class(subset_handler, "BetaHandler")
    expect_null(subset_handler$.__enclos_env__$private$.bsseq_object)
    expect_equal(nrow(subset_handler$.__enclos_env__$private$.beta_file_in_memory), length(subset_rows))
    expect_equal(subset_handler$getBetaRowNames(), subset_rows)
    expect_equal(subset_handler$getBetaColNames(), subset_cols)
    expect_equal(
        subset_handler$getBeta(),
        beta_handler$getBeta(row_names = subset_rows, col_names = subset_cols),
        tolerance = 1e-8
    )
    expect_equal(rownames(subset_handler$getBetaLocs()), subset_rows)
})

test_that("BetaHandler subset preserves compact HDF5-backed BSseq views", {
    skip_if_not_installed("HDF5Array")

    set.seed(123)
    n_loci <- 60
    n_samples <- 6
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 50, length.out = n_loci), width = 1)
    )
    site_names <- paste(seqnames(gr), start(gr), sep = ":")
    names(gr) <- site_names
    sample_names <- paste0("Sample", seq_len(n_samples))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = sample_names
    )
    bsseq_hdf5 <- HDF5Array::saveHDF5SummarizedExperiment(
        bsseq_obj,
        dir = tempfile("bsseq_hdf5_"),
        replace = TRUE
    )

    beta_handler <- getBetaHandler(beta = bsseq_hdf5)
    subset_rows <- site_names[10:20]
    subset_cols <- sample_names[c(2, 4, 6)]
    subset_handler <- beta_handler$subset(row_names = subset_rows, col_names = subset_cols)
    subset_bsseq <- subset_handler$.__enclos_env__$private$.bsseq_object

    expect_s4_class(subset_bsseq, "BSseq")
    expect_null(subset_handler$.__enclos_env__$private$.beta_file_in_memory)
    expect_equal(nrow(subset_bsseq), length(subset_rows))
    expect_equal(ncol(subset_bsseq), length(subset_cols))
    expect_equal(
        subset_handler$getBeta(),
        beta_handler$getBeta(row_names = subset_rows, col_names = subset_cols),
        tolerance = 1e-8
    )
})

test_that("disk-backed BSseq window subsetting uses integer row lookup", {
    skip_if_not_installed("HDF5Array")

    set.seed(123)
    n_loci <- 30
    n_samples <- 4
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(1000, by = 50, length.out = n_loci), width = 1)
    )
    sample_names <- paste0("Sample", seq_len(n_samples))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = sample_names
    )
    bsseq_hdf5 <- HDF5Array::saveHDF5SummarizedExperiment(
        bsseq_obj,
        dir = tempfile("bsseq_hdf5_"),
        replace = TRUE
    )
    beta_handler <- getBetaHandler(beta = bsseq_hdf5)
    beta_locs <- beta_handler$getBetaLocs()
    beta_handler$.__enclos_env__$private$.beta_row_names <- paste0("bad", seq_len(n_loci))

    subset <- CMEnt:::.subsetStage2BetaToWindows(
        beta_handler = beta_handler,
        beta_locs = beta_locs,
        col_names = sample_names[1:2],
        expansion_windows = data.frame(chr = "chr1", start = 1250L, end = 1450L)
    )

    expect_s4_class(subset$beta_handler$.__enclos_env__$private$.bsseq_object, "BSseq")
    expect_equal(as.integer(subset$beta_locs$start), seq(1250, 1450, by = 50))
    expect_equal(
        subset$beta_handler$getBeta(row_names = seq_len(nrow(subset$beta_locs))),
        beta_handler$getBeta(row_names = 6:10, col_names = sample_names[1:2]),
        tolerance = 1e-8
    )
})

test_that("HDF5-backed BSseq chromosome tasks carry compact local state", {
    skip_if_not_installed("HDF5Array")

    set.seed(321)
    n_loci <- 20
    n_samples <- 4
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    met <- matrix(rbinom(n_loci * n_samples, size = cov, prob = 0.5), ncol = n_samples)
    gr <- GRanges(
        seqnames = rep(c("chr1", "chr2"), each = n_loci / 2L),
        ranges = IRanges(start = rep(seq(1000, by = 50, length.out = n_loci / 2L), 2L), width = 1)
    )
    site_names <- paste(seqnames(gr), start(gr), sep = ":")
    names(gr) <- site_names
    sample_names <- paste0("Sample", seq_len(n_samples))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = sample_names
    )
    bsseq_hdf5 <- HDF5Array::saveHDF5SummarizedExperiment(
        bsseq_obj,
        dir = tempfile("bsseq_hdf5_"),
        replace = TRUE
    )

    beta_handler <- getBetaHandler(beta = bsseq_hdf5)
    beta_locs <- beta_handler$getBetaLocs()
    beta_chr <- as.character(beta_locs[, "chr"])
    beta_start <- as.numeric(beta_locs[, "start"])
    seed_ids <- site_names[c(3L, 8L, 12L, 18L)]
    seed_beta_index <- CMEnt:::.matchSequencingIdsToBeta(seed_ids, beta_chr, beta_start)
    seeds_locs <- as.data.frame(beta_locs[seed_beta_index, , drop = FALSE])
    rownames(seeds_locs) <- seed_ids

    tasks <- CMEnt:::.buildDMRsChromosomeTasks(
        beta_handler = beta_handler,
        beta_chr = beta_chr,
        beta_locs_rownames = rownames(beta_locs),
        chromosomes = unique(as.character(seeds_locs$chr)),
        seed_ids = seed_ids,
        seed_beta_index = seed_beta_index,
        seed_chr = as.character(seeds_locs$chr),
        beta_col_names = sample_names,
        use_numeric_sequencing_rows = TRUE
    )

    payload_names <- unique(unlist(lapply(tasks, names), use.names = FALSE))
    expect_false(any(c("beta_chr", "beta_start", "seed_chr", "seeds_locs", "beta_row_ids_all", "beta_handler") %in% payload_names))
    expect_equal(tasks$chr1$seed_beta_index, c(3L, 8L))
    expect_equal(tasks$chr2$seed_beta_index, c(2L, 8L))

    for (task in tasks) {
        task_handler <- beta_handler$subset(row_names = task$row_ids, col_names = sample_names)
        task_bsseq <- task_handler$.__enclos_env__$private$.bsseq_object
        task_locs <- as.data.frame(task_handler$getBetaLocs())
        expect_s4_class(task_bsseq, "BSseq")
        expect_null(task_handler$.__enclos_env__$private$.beta_file_in_memory)
        expect_equal(nrow(task_bsseq), nrow(task_locs))
        expect_equal(unique(as.character(task_locs$chr)), task$chr)
    }
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
