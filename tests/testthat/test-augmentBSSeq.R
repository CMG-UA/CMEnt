suppressPackageStartupMessages({
    library(testthat)
    library(CMEnt)
    library(bsseq)
})
options("CMEnt.verbose" = 0)

#' Create a mock BSseq object for testing
#'
#' @param n_sites Number of site sites
#' @param n_samples Number of samples
#' @param seed Random seed for reproducibility
#' @return A BSseq object
create_mock_bsseq <- function(n_sites = 100, n_samples = 4, seed = 42) {
    set.seed(seed)

    # Create genomic positions
    chr <- rep("chr1", n_sites)
    pos <- sort(sample(1:1000000, n_sites))

    # Create coverage matrix (Poisson distributed)
    cov <- matrix(
        rpois(n_sites * n_samples, lambda = 20),
        nrow = n_sites, ncol = n_samples
    )
    # Ensure some coverage at all sites
    cov[cov == 0] <- 1

    # Create methylation counts (binomial based on coverage)
    # Mix of low and high methylation sites
    meth_prob <- c(rep(0.1, n_sites %/% 2), rep(0.9, n_sites - n_sites %/% 2))
    M <- matrix(0, nrow = n_sites, ncol = n_samples)
    for (i in seq_len(n_sites)) {
        for (j in seq_len(n_samples)) {
            M[i, j] <- rbinom(1, size = cov[i, j], prob = meth_prob[i])
        }
    }

    sample_names <- paste0("sample_", seq_len(n_samples))

    BSseq(
        chr = chr,
        pos = pos,
        M = M,
        Cov = cov,
        sampleNames = sample_names
    )
}

#' Create a BSseq object with some zero coverage sites
#'
#' @param n_sites Number of site sites
#' @param n_samples Number of samples
#' @param zero_frac Fraction of entries to set to zero
#' @param seed Random seed
#' @return A BSseq object with sparse coverage
create_sparse_bsseq <- function(n_sites = 100, n_samples = 4, zero_frac = 0.3, seed = 42) {
    set.seed(seed)

    chr <- rep("chr1", n_sites)
    pos <- sort(sample(1:1000000, n_sites))

    cov <- matrix(
        rpois(n_sites * n_samples, lambda = 15),
        nrow = n_sites, ncol = n_samples
    )

    # Introduce zeros
    n_zeros <- round(n_sites * n_samples * zero_frac)
    zero_idx <- sample(seq_len(n_sites * n_samples), n_zeros)
    cov[zero_idx] <- 0

    # Methylation counts
    meth_prob <- runif(n_sites, 0.1, 0.9)
    M <- matrix(0, nrow = n_sites, ncol = n_samples)
    for (i in seq_len(n_sites)) {
        for (j in seq_len(n_samples)) {
            if (cov[i, j] > 0) {
                M[i, j] <- rbinom(1, size = cov[i, j], prob = meth_prob[i])
            }
        }
    }

    sample_names <- paste0("sample_", seq_len(n_samples))

    BSseq(
        chr = chr,
        pos = pos,
        M = M,
        Cov = cov,
        sampleNames = sample_names
    )
}


test_that("augmentBSSeq returns a BSseq object", {
    bsseq_obj <- create_mock_bsseq(n_sites = 50, n_samples = 3)
    set.seed(123)
    result <- augmentBSSeq(bsseq_obj, n_new_samples = 2)

    expect_s4_class(result, "BSseq")
})

test_that("augmentBSSeq adds correct number of samples", {
    bsseq_obj <- create_mock_bsseq(n_sites = 50, n_samples = 3)
    n_new <- 5
    set.seed(123)
    result <- augmentBSSeq(bsseq_obj, n_new_samples = n_new, min_samples = 2)

    # Original samples + new synthetic samples
    # Note: some sites may be filtered, so we check sample count
    expect_equal(ncol(result), 3 + n_new)
})

test_that("augmentBSSeq synthetic sample names are correct", {
    bsseq_obj <- create_mock_bsseq(n_sites = 50, n_samples = 3)
    set.seed(123)
    result <- augmentBSSeq(bsseq_obj, n_new_samples = 3, min_samples = 2)

    sample_names <- sampleNames(result)
    expect_true(all(paste0("synthetic_", 1:3) %in% sample_names))
})

test_that("augmentBSSeq is reproducible with external seed", {
    bsseq_obj <- create_mock_bsseq(n_sites = 50, n_samples = 3)

    set.seed(42)
    result1 <- augmentBSSeq(bsseq_obj, n_new_samples = 2, min_samples = 2)
    set.seed(42)
    result2 <- augmentBSSeq(bsseq_obj, n_new_samples = 2, min_samples = 2)

    # Coverage should be identical
    expect_equal(getCoverage(result1), getCoverage(result2))
})

test_that("augmentBSSeq produces different results with different external seeds", {
    bsseq_obj <- create_mock_bsseq(n_sites = 50, n_samples = 3)

    set.seed(42)
    result1 <- augmentBSSeq(bsseq_obj, n_new_samples = 2, min_samples = 2)
    set.seed(99)
    result2 <- augmentBSSeq(bsseq_obj, n_new_samples = 2, min_samples = 2)

    # Coverage matrices should differ
    expect_false(identical(getCoverage(result1), getCoverage(result2)))
})

test_that("augmentBSSeq filters sites based on min_samples", {
    bsseq_obj <- create_sparse_bsseq(n_sites = 100, n_samples = 4, zero_frac = 0.5)

    # Count sites with coverage in at least 3 samples
    cov <- getCoverage(bsseq_obj)
    expected_sites <- sum(rowSums(cov > 0) >= 3)

    set.seed(123)
    result <- augmentBSSeq(bsseq_obj, n_new_samples = 1, min_samples = 3)

    expect_equal(nrow(result), expected_sites)
})

test_that("augmentBSSeq coverage values are positive", {
    bsseq_obj <- create_mock_bsseq(n_sites = 50, n_samples = 3)
    set.seed(123)
    result <- augmentBSSeq(bsseq_obj, n_new_samples = 3, min_samples = 2)

    cov <- getCoverage(result)
    # All coverage should be >= 1
    expect_true(all(cov >= 1))
})

test_that("augmentBSSeq methylation counts do not exceed coverage", {
    bsseq_obj <- create_mock_bsseq(n_sites = 50, n_samples = 3)
    set.seed(123)
    result <- augmentBSSeq(bsseq_obj, n_new_samples = 3, min_samples = 2)

    cov <- getCoverage(result)
    M <- getCoverage(result, type = "M")

    # Methylated counts should not exceed coverage
    expect_true(all(M <= cov))
})

test_that("augmentBSSeq methylation counts are non-negative", {
    bsseq_obj <- create_mock_bsseq(n_sites = 50, n_samples = 3)
    set.seed(123)
    result <- augmentBSSeq(bsseq_obj, n_new_samples = 3, min_samples = 2)

    M <- getCoverage(result, type = "M")
    expect_true(all(M >= 0))
})

test_that("augmentBSSeq preserves original samples", {
    bsseq_obj <- create_mock_bsseq(n_sites = 50, n_samples = 3, seed = 42)
    set.seed(123)
    result <- augmentBSSeq(bsseq_obj, n_new_samples = 2, min_samples = 2)

    # Check original sample names are preserved
    orig_names <- sampleNames(bsseq_obj)
    result_names <- sampleNames(result)

    expect_true(all(orig_names %in% result_names))
})

test_that("augmentBSSeq handles single new sample", {
    bsseq_obj <- create_mock_bsseq(n_sites = 50, n_samples = 3)
    set.seed(123)
    result <- augmentBSSeq(bsseq_obj, n_new_samples = 1, min_samples = 2)

    expect_s4_class(result, "BSseq")
    expect_true("synthetic_1" %in% sampleNames(result))
})

test_that("augmentBSSeq handles many new samples", {
    bsseq_obj <- create_mock_bsseq(n_sites = 50, n_samples = 3)
    set.seed(123)
    result <- augmentBSSeq(bsseq_obj, n_new_samples = 20, min_samples = 2)

    expect_s4_class(result, "BSseq")
    expect_equal(ncol(result), 3 + 20)
})

test_that("augmentBSSeq works without an explicit seed argument", {
    bsseq_obj <- create_mock_bsseq(n_sites = 50, n_samples = 3)

    # Should not error
    expect_no_error({
        result <- augmentBSSeq(bsseq_obj, n_new_samples = 2, min_samples = 2)
    })
})

test_that("augmentBSSeq avoids over-dispersion with low-input boundary methylation", {
    n_sites <- 40
    n_samples <- 2

    chr <- rep("chr1", n_sites)
    pos <- seq(100, by = 10, length.out = n_sites)
    cov <- matrix(3L, nrow = n_sites, ncol = n_samples)
    M <- matrix(0L, nrow = n_sites, ncol = n_samples)
    M[(n_sites %/% 2 + 1):n_sites, ] <- 3L

    bsseq_obj <- BSseq(
        chr = chr,
        pos = pos,
        M = M,
        Cov = cov,
        sampleNames = paste0("sample_", seq_len(n_samples))
    )

    set.seed(123)
    result <- augmentBSSeq(bsseq_obj, n_new_samples = 100, min_samples = 2)
    synthetic_idx <- grepl("^synthetic_", sampleNames(result))
    synthetic_meth <- bsseq::getMeth(result, type = "raw")[, synthetic_idx, drop = FALSE]

    low_sites <- seq_len(n_sites %/% 2)
    high_sites <- (n_sites %/% 2 + 1):n_sites
    synthetic_site_means <- rowMeans(synthetic_meth)
    synthetic_site_sd <- apply(synthetic_meth, 1, sd)

    expect_gt(mean(synthetic_site_means[low_sites] < 0.3), 0.9)
    expect_gt(mean(synthetic_site_means[high_sites] > 0.7), 0.9)
    expect_lt(median(synthetic_site_sd), 0.3)
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

lag_correlation_profile <- function(bs_obj, sample_idx, max_lag = 5L) {
    cov <- as.matrix(bsseq::getCoverage(bs_obj, type = "Cov")[, sample_idx, drop = FALSE])
    meth <- as.matrix(bsseq::getCoverage(bs_obj, type = "M")[, sample_idx, drop = FALSE])
    prop <- (meth + 0.5) / (cov + 1)
    signal <- qlogis(pmin(pmax(prop, 1e-6), 1 - 1e-6))

    chr <- as.character(GenomeInfoDb::seqnames(bs_obj))
    profile <- numeric(max_lag)

    for (lag in seq_len(max_lag)) {
        pair_idx <- which(chr[(lag + 1):nrow(signal)] == chr[1:(nrow(signal) - lag)])
        idx_1 <- pair_idx
        idx_2 <- pair_idx + lag

        corrs <- vapply(seq_along(idx_1), function(i) {
            suppressWarnings(stats::cor(signal[idx_1[i], ], signal[idx_2[i], ]))
        }, numeric(1))

        profile[lag] <- stats::median(corrs, na.rm = TRUE)
    }

    profile
}

make_correlated_bsseq <- function(seed = 1) {
    set.seed(seed)

    n_site <- 120L
    n_samples <- 16L
    pos <- seq(100, by = 60, length.out = n_site)
    base_prob <- plogis(seq(-1.25, 1.25, length.out = n_site))

    phi <- exp(-60 / 350)
    latent <- matrix(0, nrow = n_site, ncol = n_samples)
    latent[1, ] <- stats::rnorm(n_samples)
    for (i in 2:n_site) {
        latent[i, ] <- phi * latent[i - 1L, ] + sqrt(1 - phi^2) * stats::rnorm(n_samples)
    }

    cov <- matrix(35L, nrow = n_site, ncol = n_samples)
    prob <- plogis(qlogis(base_prob) + 0.9 * latent)
    meth <- matrix(
        stats::rbinom(n_site * n_samples, size = cov, prob = prob),
        nrow = n_site,
        ncol = n_samples
    )

    sample_names <- c("synthetic_1", paste0("Sample", 2:n_samples))
    colnames(meth) <- colnames(cov) <- sample_names

    bsseq::BSseq(
        chr = rep("chr1", n_site),
        pos = pos,
        M = meth,
        Cov = cov,
        sampleNames = sample_names
    )
}

test_that("augmentBSSeq preserves neighboring site correlation", {
    bs <- make_correlated_bsseq(seed = 42)

    set.seed(99)
    aug <- augmentBSSeq(bs, n_new_samples = 30)

    orig_idx <- seq_len(ncol(bs))
    syn_idx <- (ncol(bs) + 1):ncol(aug)

    orig_corr <- adjacent_sample_correlation(aug, orig_idx)
    syn_corr <- adjacent_sample_correlation(aug, syn_idx)

    expect_equal(ncol(aug), ncol(bs) + 30L)
    expect_true(anyDuplicated(colnames(aug)) == 0L)
    expect_gt(syn_corr, 0.25)
    expect_gt(syn_corr, orig_corr - 0.2)
})

test_that("augmentBSSeq preserves correlation decay with lag", {
    bs <- make_correlated_bsseq(seed = 42)

    set.seed(99)
    aug <- augmentBSSeq(bs, n_new_samples = 30)

    orig_idx <- seq_len(ncol(bs))
    syn_idx <- (ncol(bs) + 1):ncol(aug)

    orig_profile <- lag_correlation_profile(aug, orig_idx, max_lag = 6L)
    syn_profile <- lag_correlation_profile(aug, syn_idx, max_lag = 6L)

    expect_equal(length(orig_profile), 6L)
    expect_equal(length(syn_profile), 6L)
    expect_true(all(diff(orig_profile) <= 1e-8, na.rm = TRUE))
    expect_true(all(diff(syn_profile) <= 1e-8, na.rm = TRUE))
    expect_lt(mean(abs(orig_profile - syn_profile), na.rm = TRUE), 0.25)
})
