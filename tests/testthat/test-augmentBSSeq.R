suppressPackageStartupMessages({
    library(testthat)
    library(DMRsegal)
    library(bsseq)
})

adjacent_sample_correlation <- function(bs_obj, sample_idx) {
    cov <- as.matrix(bsseq::getCoverage(bs_obj, type = "Cov")[, sample_idx, drop = FALSE])
    meth <- as.matrix(bsseq::getCoverage(bs_obj, type = "M")[, sample_idx, drop = FALSE])
    prop <- (meth + 0.5) / (cov + 1)
    signal <- qlogis(pmin(pmax(prop, 1e-6), 1 - 1e-6))

    chr <- as.character(GenomeInfoDb::seqnames(bs_obj))
    pair_idx <- which(chr[-1] == chr[-nrow(signal)])
    corrs <- vapply(pair_idx, function(i) {
        suppressWarnings(stats::cor(signal[i, ], signal[i + 1L, ]))
    }, numeric(1))

    stats::median(corrs, na.rm = TRUE)
}

make_correlated_bsseq <- function(seed = 1) {
    set.seed(seed)

    n_cpg <- 120L
    n_samples <- 16L
    pos <- seq(100, by = 60, length.out = n_cpg)
    base_prob <- plogis(seq(-1.25, 1.25, length.out = n_cpg))

    phi <- exp(-60 / 350)
    latent <- matrix(0, nrow = n_cpg, ncol = n_samples)
    latent[1, ] <- stats::rnorm(n_samples)
    for (i in 2:n_cpg) {
        latent[i, ] <- phi * latent[i - 1L, ] + sqrt(1 - phi^2) * stats::rnorm(n_samples)
    }

    cov <- matrix(35L, nrow = n_cpg, ncol = n_samples)
    prob <- plogis(qlogis(base_prob) + 0.9 * latent)
    meth <- matrix(
        stats::rbinom(n_cpg * n_samples, size = cov, prob = prob),
        nrow = n_cpg,
        ncol = n_samples
    )

    sample_names <- c("synthetic_1", paste0("Sample", 2:n_samples))
    colnames(meth) <- colnames(cov) <- sample_names

    bsseq::BSseq(
        chr = rep("chr1", n_cpg),
        pos = pos,
        M = meth,
        Cov = cov,
        sampleNames = sample_names
    )
}

test_that("augmentBSSeq preserves neighboring CpG correlation", {
    bs <- make_correlated_bsseq(seed = 42)

    aug <- augmentBSSeq(bs, n_new_samples = 30, seed = 99)

    orig_idx <- seq_len(ncol(bs))
    syn_idx <- (ncol(bs) + 1):ncol(aug)

    orig_corr <- adjacent_sample_correlation(aug, orig_idx)
    syn_corr <- adjacent_sample_correlation(aug, syn_idx)

    expect_equal(ncol(aug), ncol(bs) + 30L)
    expect_true(anyDuplicated(colnames(aug)) == 0L)
    expect_gt(syn_corr, 0.25)
    expect_gt(syn_corr, orig_corr - 0.2)
})
