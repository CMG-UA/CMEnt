suppressPackageStartupMessages({
    library(testthat)
    library(DMRsegal)
    library(bsseq)
})

makeSyntheticBSseq <- function(n.cpg = 100, n.samples = 20, flat = FALSE,
                               n.clusters = 4, cluster.gap = 2000) {
    cov <- matrix(30L, nrow = n.cpg, ncol = n.samples)
    if (flat) {
        met <- matrix(15L, nrow = n.cpg, ncol = n.samples)
    } else {
        met <- matrix(
            rbinom(n.cpg * n.samples, size = 30, prob = rbeta(n.cpg * n.samples, 8, 8)),
            nrow = n.cpg,
            ncol = n.samples
        )
    }
    sample.names <- paste0("Sample", seq_len(n.samples))
    colnames(met) <- colnames(cov) <- sample.names
    cluster.sizes <- rep(floor(n.cpg / n.clusters), n.clusters)
    cluster.sizes[seq_len(n.cpg %% n.clusters)] <- cluster.sizes[seq_len(n.cpg %% n.clusters)] + 1L
    cluster.starts <- numeric(length(cluster.sizes))
    cluster.starts[[1]] <- 100
    if (length(cluster.sizes) > 1L) {
        for (cluster.idx in 2:length(cluster.sizes)) {
            cluster.starts[[cluster.idx]] <- cluster.starts[[cluster.idx - 1L]] +
                (cluster.sizes[[cluster.idx - 1L]] - 1L) * 100 + cluster.gap
        }
    }
    positions <- unlist(lapply(seq_along(cluster.sizes), function(cluster.idx) {
        seq(cluster.starts[[cluster.idx]], by = 100, length.out = cluster.sizes[[cluster.idx]])
    }))
    bsseq::BSseq(
        chr = rep("chr1", n.cpg),
        pos = positions,
        M = met,
        Cov = cov,
        sampleNames = sample.names
    )
}

centered_neighbor_correlation <- function(sim, bs0, dmr.index = 1) {
    idx <- queryHits(GenomicRanges::findOverlaps(rowRanges(bs0), sim$gr.dmrs[dmr.index]))
    orig.p <- as.matrix(bsseq::getCoverage(bs0, type = "M")[idx, ] /
        bsseq::getCoverage(bs0, type = "Cov")[idx, ])
    sim.p <- as.matrix(bsseq::getCoverage(sim$bs, type = "M")[idx, ] /
        bsseq::getCoverage(sim$bs, type = "Cov")[idx, ])

    half <- ncol(bs0) / 2
    g1 <- seq_len(half)
    g2 <- (half + 1):ncol(bs0)
    s1 <- mean(abs(sim.p[, g1] - orig.p[, g1]))
    s2 <- mean(abs(sim.p[, g2] - orig.p[, g2]))
    affected <- if (s1 > s2) g1 else g2

    delta <- sim.p[, affected, drop = FALSE] - orig.p[, affected, drop = FALSE]
    delta.centered <- delta - rowMeans(delta)
    neigh <- vapply(seq_len(nrow(delta.centered) - 1), function(i) {
        suppressWarnings(cor(delta.centered[i, ], delta.centered[i + 1, ]))
    }, numeric(1))

    mean(neigh, na.rm = TRUE)
}

mean_centered_neighbor_correlation <- function(sim, bs0) {
    mean(vapply(seq_along(sim$gr.dmrs), function(idx) {
        centered_neighbor_correlation(sim, bs0, dmr.index = idx)
    }, numeric(1)))
}

test_that("simDMRsCorrelated uses background calibration by default", {
    set.seed(123)
    bs <- makeSyntheticBSseq()

    set.seed(777)
    sim <- suppressMessages(simDMRsCorrelated(bs, num.dmrs = 3))

    expect_true(all(c(
        "corr_target", "corr_sd_used", "sample_sd_frac_used", "corr_mode_used"
    ) %in% colnames(S4Vectors::mcols(sim$gr.dmrs))))
    expect_true(all(is.finite(S4Vectors::mcols(sim$gr.dmrs)$corr_target)))
    expect_true(all(S4Vectors::mcols(sim$gr.dmrs)$corr_mode_used == "background"))
    expect_true(all(S4Vectors::mcols(sim$gr.dmrs)$corr_sd_used > 0))
    expect_true(all(S4Vectors::mcols(sim$gr.dmrs)$sample_sd_frac_used > 0))
})

test_that("simDMRsCorrelated adds correlated sample-level deviations when enabled", {
    set.seed(123)
    bs <- makeSyntheticBSseq()

    set.seed(777)
    sim.cor <- suppressMessages(simDMRsCorrelated(
        bs,
        num.dmrs = 1,
        use.correlated.effects = TRUE,
        corr.mode = "manual",
        corr.sd = 0.4,
        corr.length = 300,
        sample.sd.frac = 1
    ))
    set.seed(777)
    sim.nocor <- suppressMessages(simDMRsCorrelated(
        bs,
        num.dmrs = 1,
        use.correlated.effects = FALSE,
        corr.mode = "manual",
        corr.sd = 0.4,
        corr.length = 300,
        sample.sd.frac = 1
    ))

    cor.metric <- centered_neighbor_correlation(sim.cor, bs)
    nocor.metric <- centered_neighbor_correlation(sim.nocor, bs)

    expect_gt(cor.metric, nocor.metric + 0.1)
    expect_true(all(S4Vectors::mcols(sim.cor$gr.dmrs)$corr_mode_used == "manual"))
})

test_that("manual settings stay backward compatible when corr.mode is omitted", {
    set.seed(123)
    bs <- makeSyntheticBSseq()

    set.seed(999)
    sim.implicit <- suppressMessages(simDMRsCorrelated(
        bs,
        num.dmrs = 2,
        corr.sd = 0.35,
        sample.sd.frac = 1
    ))
    set.seed(999)
    sim.manual <- suppressMessages(simDMRsCorrelated(
        bs,
        num.dmrs = 2,
        corr.mode = "manual",
        corr.sd = 0.35,
        sample.sd.frac = 1
    ))

    expect_identical(sim.implicit$gr.dmrs, sim.manual$gr.dmrs)
    expect_equal(as.matrix(bsseq::getCoverage(sim.implicit$bs, type = "M")),
        as.matrix(bsseq::getCoverage(sim.manual$bs, type = "M")))
    expect_equal(as.matrix(bsseq::getCoverage(sim.implicit$bs, type = "Cov")),
        as.matrix(bsseq::getCoverage(sim.manual$bs, type = "Cov")))
    expect_true(all(S4Vectors::mcols(sim.implicit$gr.dmrs)$corr_mode_used == "manual"))
})

test_that("corr.rate scales the sampled target and achieved correlation signal", {
    set.seed(123)
    bs <- makeSyntheticBSseq()

    set.seed(2024)
    sim.low <- suppressMessages(simDMRsCorrelated(
        bs,
        num.dmrs = 4,
        corr.rate = 0.5
    ))
    set.seed(2024)
    sim.high <- suppressMessages(simDMRsCorrelated(
        bs,
        num.dmrs = 4,
        corr.rate = 1.5
    ))

    expect_gt(
        mean(S4Vectors::mcols(sim.high$gr.dmrs)$corr_target),
        mean(S4Vectors::mcols(sim.low$gr.dmrs)$corr_target)
    )
    expect_gt(
        mean_centered_neighbor_correlation(sim.high, bs),
        mean_centered_neighbor_correlation(sim.low, bs)
    )
})

test_that("degenerate background fits fall back to tuned manual settings", {
    bs <- makeSyntheticBSseq(flat = TRUE)

    set.seed(321)
    sim <- suppressMessages(simDMRsCorrelated(bs, num.dmrs = 2))

    expect_true(all(S4Vectors::mcols(sim$gr.dmrs)$corr_mode_used == "background"))
    expect_equal(S4Vectors::mcols(sim$gr.dmrs)$corr_sd_used, rep(0.3, 2))
    expect_equal(S4Vectors::mcols(sim$gr.dmrs)$sample_sd_frac_used, rep(0.75, 2))
    expect_true(all(is.finite(S4Vectors::mcols(sim$gr.dmrs)$corr_target)))
})
