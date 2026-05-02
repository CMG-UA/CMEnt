# nolint start: object_name_linter

#' Augment BSseq Object
#'
#' Generate synthetic samples while preserving both per-site coverage and
#' methylation marginals and the local correlation structure of neighboring sites.
#' Coverage and methylation are first fit with per-site Poisson/Beta marginals,
#' then synthetic samples are drawn through chromosome-local Gaussian copula
#' fields whose dependence is estimated from adjacent sites in the observed data.
#' For reproducible synthetic samples, call `set.seed()` before `augmentBSSeq()`.
#'
#' @param bs A BSseq object
#' @param n_new_samples Number of new synthetic samples to generate
#' @param min_samples Minimum number of samples with coverage required per site
#' @param calibrate_correlation Logical. If `TRUE`, iteratively adjusts the
#'     latent Gaussian length scales so adjacent-site correlations are matched
#'     after transforming back through the observed coverage and methylation
#'     sampling layers.
#' @param calibration_iterations Maximum number of bisection iterations used
#'     for each correlation calibration.
#' @param calibration_samples Number of synthetic samples used internally for
#'     correlation calibration. If `NULL`, a capped conservative default is
#'     chosen from the input and requested output sample sizes.
#' @importFrom SummarizedExperiment assays
#' @return A BSseq object with original and synthetic samples
#' @examples
#' \donttest{
#'     # Load example BSseq data
#'     data("BSobj", package = "bsseq")
#'     set.seed(123)
#'     # Augment with 5 synthetic samples
#'     augmented_bs <- augmentBSSeq(BSobj, n_new_samples = 5)
#' }
#' @export
augmentBSSeq <- function(bs, n_new_samples, min_samples = 2,
                         calibrate_correlation = TRUE,
                         calibration_iterations = 8,
                         calibration_samples = NULL) {
    n_new_samples <- as.integer(n_new_samples)
    min_samples <- as.integer(min_samples)
    calibrate_correlation <- as.logical(calibrate_correlation)
    calibration_iterations <- as.integer(calibration_iterations)

    if (length(n_new_samples) != 1L || is.na(n_new_samples) || n_new_samples < 1L) {
        stop("'n_new_samples' must be a positive integer")
    }
    if (length(min_samples) != 1L || is.na(min_samples) || min_samples < 1L) {
        stop("'min_samples' must be a positive integer")
    }
    if (length(calibrate_correlation) != 1L || is.na(calibrate_correlation)) {
        stop("'calibrate_correlation' must be TRUE or FALSE")
    }
    if (
        length(calibration_iterations) != 1L ||
            is.na(calibration_iterations) ||
            calibration_iterations < 0L
    ) {
        stop("'calibration_iterations' must be a non-negative integer")
    }
    if (!is.null(calibration_samples)) {
        calibration_samples <- as.integer(calibration_samples)
        if (
            length(calibration_samples) != 1L ||
                is.na(calibration_samples) ||
                calibration_samples < 3L
        ) {
            stop("'calibration_samples' must be NULL or an integer >= 3")
        }
    }

    # Keep only sites with coverage in at least min_samples
    cov_matrix <- bsseq::getCoverage(bs)
    valid_sites <- rowSums(cov_matrix > 0) >= min_samples
    bsseq_filtered <- bs[valid_sites, ]
    # Keep assay-level sample dimnames aligned with colData before combining objects.
    colnames(SummarizedExperiment::assays(bsseq_filtered)$M) <- colnames(bsseq_filtered)
    colnames(SummarizedExperiment::assays(bsseq_filtered)$Cov) <- colnames(bsseq_filtered)

    if (nrow(bsseq_filtered) == 0L) {
        stop("No site sites have coverage in at least 'min_samples' samples")
    }

    # Extract coverage and methylation counts
    cov <- as.matrix(bsseq::getCoverage(bsseq_filtered))
    M <- as.matrix(bsseq::getCoverage(bsseq_filtered, type = "M"))

    n_sites <- nrow(bsseq_filtered)
    n_orig_samples <- ncol(cov)
    chr <- as.character(seqnames(bsseq_filtered))
    pos <- start(bsseq_filtered)
    if (is.null(calibration_samples)) {
        calibration_samples <- min(max(50L, n_orig_samples, n_new_samples), 200L)
    }

    # Coverage: compute mean per site (excluding zeros)
    cov_nonzero <- cov
    cov_nonzero[cov_nonzero == 0] <- NA_real_
    lambda_per_site <- rowMeans(cov_nonzero, na.rm = TRUE)
    lambda_per_site[is.na(lambda_per_site)] <- 1 # fallback

    # Methylation: use smoothed proportions so boundary values (0/1) are retained.
    valid_cov <- !is.na(cov_nonzero)
    meth_prop <- matrix(NA_real_, nrow = n_sites, ncol = n_orig_samples)
    meth_prop[valid_cov] <- (M[valid_cov] + 0.5) / (cov_nonzero[valid_cov] + 1)

    # Coverage-weighted site means with Jeffreys-style pseudocounts.
    pseudo_total <- rowSums((cov_nonzero + 1) * valid_cov, na.rm = TRUE)
    pseudo_meth <- rowSums((M + 0.5) * valid_cov, na.rm = TRUE)
    site_mean <- pseudo_meth / pmax(pseudo_total, 1)

    # Raw per-site variance across samples (used to estimate Beta concentration).
    meth_obs_n <- rowSums(valid_cov, na.rm = TRUE)
    site_var <- apply(meth_prop, 1, var, na.rm = TRUE)

    # Convert moments to concentration (kappa). Larger kappa => lower variance.
    mean_var_cap <- pmax(site_mean * (1 - site_mean), 1e-6)
    site_var <- pmin(site_var, mean_var_cap * 0.99)
    kappa_obs <- mean_var_cap / pmax(site_var, 1e-6) - 1
    kappa_obs[!is.finite(kappa_obs)] <- NA_real_

    # Empirical prior concentration across informative sites.
    kappa_emp <- median(kappa_obs[meth_obs_n >= 3 & !is.na(kappa_obs)], na.rm = TRUE)
    if (!is.finite(kappa_emp)) kappa_emp <- 20
    kappa_prior <- max(12, kappa_emp)

    # Sample-size-aware shrinkage: low-input sites lean toward conservative prior.
    shrink_w <- pmax(meth_obs_n - 1, 0) / (pmax(meth_obs_n - 1, 0) + 4)
    kappa_site <- ifelse(
        is.na(kappa_obs),
        kappa_prior,
        shrink_w * kappa_obs + (1 - shrink_w) * kappa_prior
    )
    kappa_site <- pmax(pmin(kappa_site, 500), 4)

    alpha_per_site <- pmax(site_mean * kappa_site, 0.1)
    beta_per_site <- pmax((1 - site_mean) * kappa_site, 0.1)

    same_chr <- chr[-1] == chr[-n_sites]
    adjacent_gaps <- pmax(pos[-1] - pos[-n_sites], 1)
    default_gap <- stats::median(adjacent_gaps[same_chr], na.rm = TRUE)
    if (!is.finite(default_gap)) default_gap <- 1
    default_length <- max(5 * default_gap, 50)

    adjacent_pair_correlations <- function(signal) {
        pair_idx <- which(same_chr)
        if (length(pair_idx) == 0L || ncol(signal) < 3L) {
            return(numeric())
        }

        vapply(pair_idx, function(i) {
            x <- signal[i, ]
            y <- signal[i + 1L, ]
            keep <- is.finite(x) & is.finite(y)
            if (sum(keep) < 3L) {
                return(NA_real_)
            }
            if (stats::sd(x[keep]) == 0 || stats::sd(y[keep]) == 0) {
                return(NA_real_)
            }
            suppressWarnings(stats::cor(x[keep], y[keep]))
        }, numeric(1))
    }

    median_pair_correlation <- function(pair_corr, min_pairs = 5L) {
        valid <- is.finite(pair_corr)
        if (sum(valid) < min_pairs) {
            return(NA_real_)
        }
        stats::median(pair_corr[valid], na.rm = TRUE)
    }

    estimate_length_scale <- function(signal, fallback_length) {
        pair_idx <- which(same_chr)
        if (length(pair_idx) == 0L || n_orig_samples < 3L) {
            return(fallback_length)
        }

        pair_corr <- adjacent_pair_correlations(signal)
        valid <- is.finite(pair_corr) & pair_corr > 0
        if (sum(valid) < 5L) {
            return(fallback_length)
        }

        ell <- stats::median(
            -adjacent_gaps[pair_idx][valid] / log(pmin(pair_corr[valid], 0.995)),
            na.rm = TRUE
        )
        if (!is.finite(ell)) {
            return(fallback_length)
        }

        shrink_w <- sum(valid) / (sum(valid) + 25)
        ell <- shrink_w * ell + (1 - shrink_w) * fallback_length
        max(ell, 1)
    }

    clamp_unit_interval <- function(x) {
        pmin(pmax(x, 1e-10), 1 - 1e-10)
    }

    row_correlations <- function(x, y) {
        x_centered <- x - rowMeans(x)
        y_centered <- y - rowMeans(y)
        numerator <- rowSums(x_centered * y_centered)
        denominator <- sqrt(rowSums(x_centered^2) * rowSums(y_centered^2))
        out <- numerator / denominator
        out[!is.finite(out) | denominator <= 0] <- NA_real_
        out
    }

    select_calibration_pairs <- function(pair_idx, max_pairs = 2000L) {
        if (length(pair_idx) <= max_pairs) {
            return(pair_idx)
        }
        pair_idx[unique(round(seq(1, length(pair_idx), length.out = max_pairs)))]
    }

    latent_pair_components <- function(pair_idx, length_scale, base_noise, innovation_noise) {
        if (length(pair_idx) == 0L || !is.finite(length_scale) || length_scale <= 0) {
            return(list(z1 = base_noise, z2 = innovation_noise))
        }

        phi <- exp(-pmax(adjacent_gaps[pair_idx], 1) / length_scale)
        innovation_sd <- sqrt(pmax(1 - phi^2, 1e-8))
        z2 <- sweep(base_noise, 1, phi, `*`) +
            sweep(innovation_noise, 1, innovation_sd, `*`)
        list(z1 = base_noise, z2 = z2)
    }

    simulate_pair_coverage <- function(length_scale, pair_idx, base_noise, innovation_noise) {
        n_pairs <- length(pair_idx)
        n_samples <- ncol(base_noise)
        latent <- latent_pair_components(pair_idx, length_scale, base_noise, innovation_noise)
        u1 <- clamp_unit_interval(stats::pnorm(latent$z1))
        u2 <- clamp_unit_interval(stats::pnorm(latent$z2))

        cov_1 <- matrix(
            stats::qpois(c(u1), lambda = rep(lambda_per_site[pair_idx], n_samples)),
            nrow = n_pairs,
            ncol = n_samples
        )
        cov_2 <- matrix(
            stats::qpois(c(u2), lambda = rep(lambda_per_site[pair_idx + 1L], n_samples)),
            nrow = n_pairs,
            ncol = n_samples
        )

        list(cov_1 = pmax(cov_1, 1L), cov_2 = pmax(cov_2, 1L))
    }

    evaluate_coverage_pair_correlation <- function(length_scale, pair_idx, base_noise, innovation_noise) {
        pair_cov <- simulate_pair_coverage(length_scale, pair_idx, base_noise, innovation_noise)
        median_pair_correlation(row_correlations(log1p(pair_cov$cov_1), log1p(pair_cov$cov_2)))
    }

    evaluate_methylation_pair_correlation <- function(length_scale, pair_idx, base_noise,
                                                      innovation_noise, binom_u_1,
                                                      binom_u_2, pair_cov) {
        n_pairs <- length(pair_idx)
        n_samples <- ncol(base_noise)
        latent <- latent_pair_components(pair_idx, length_scale, base_noise, innovation_noise)
        u1 <- clamp_unit_interval(stats::pnorm(latent$z1))
        u2 <- clamp_unit_interval(stats::pnorm(latent$z2))

        prob_1 <- matrix(
            stats::qbeta(
                c(u1),
                shape1 = rep(alpha_per_site[pair_idx], n_samples),
                shape2 = rep(beta_per_site[pair_idx], n_samples)
            ),
            nrow = n_pairs,
            ncol = n_samples
        )
        prob_2 <- matrix(
            stats::qbeta(
                c(u2),
                shape1 = rep(alpha_per_site[pair_idx + 1L], n_samples),
                shape2 = rep(beta_per_site[pair_idx + 1L], n_samples)
            ),
            nrow = n_pairs,
            ncol = n_samples
        )

        meth_1 <- matrix(
            stats::qbinom(c(binom_u_1), size = c(pair_cov$cov_1), prob = c(prob_1)),
            nrow = n_pairs,
            ncol = n_samples
        )
        meth_2 <- matrix(
            stats::qbinom(c(binom_u_2), size = c(pair_cov$cov_2), prob = c(prob_2)),
            nrow = n_pairs,
            ncol = n_samples
        )

        signal_1 <- qlogis(clamp_unit_interval((meth_1 + 0.5) / (pair_cov$cov_1 + 1)))
        signal_2 <- qlogis(clamp_unit_interval((meth_2 + 0.5) / (pair_cov$cov_2 + 1)))
        median_pair_correlation(row_correlations(signal_1, signal_2))
    }

    calibrate_length_scale <- function(initial_length, target_corr, evaluate_corr,
                                       max_iterations, tolerance = 0.02,
                                       min_target_corr = 0.05) {
        if (
            !is.finite(initial_length) ||
                !is.finite(target_corr) ||
                target_corr <= min_target_corr ||
                max_iterations < 1L
        ) {
            return(initial_length)
        }

        target_corr <- pmin(target_corr, 0.99)
        max_gap <- suppressWarnings(max(adjacent_gaps[same_chr], na.rm = TRUE))
        if (!is.finite(max_gap)) max_gap <- default_gap
        max_length <- max(default_length, initial_length, max_gap, 1) * 1000

        candidate_lengths <- numeric()
        candidate_corrs <- numeric()
        evaluate <- function(length_scale) {
            corr <- evaluate_corr(length_scale)
            candidate_lengths <<- c(candidate_lengths, length_scale)
            candidate_corrs <<- c(candidate_corrs, corr)
            corr
        }
        best_length <- function() {
            valid <- is.finite(candidate_corrs)
            if (!any(valid)) {
                return(initial_length)
            }
            candidate_lengths[valid][which.min(abs(candidate_corrs[valid] - target_corr))]
        }

        initial_corr <- evaluate(initial_length)
        if (!is.finite(initial_corr) || abs(initial_corr - target_corr) <= tolerance) {
            return(best_length())
        }

        if (initial_corr < target_corr) {
            lower <- max(initial_length, 1)
            lower_corr <- initial_corr
            upper <- lower
            upper_corr <- lower_corr
            while (upper_corr < target_corr && upper < max_length) {
                upper <- min(upper * 2, max_length)
                upper_corr <- evaluate(upper)
                if (!is.finite(upper_corr)) {
                    return(best_length())
                }
            }
        } else {
            lower <- 1
            lower_corr <- evaluate(lower)
            upper <- max(initial_length, 1)
            upper_corr <- initial_corr
        }

        if (
            !is.finite(lower_corr) ||
                !is.finite(upper_corr) ||
                lower_corr > target_corr ||
                upper_corr < target_corr
        ) {
            return(best_length())
        }

        for (iteration in seq_len(max_iterations)) {
            mid <- sqrt(lower * upper)
            mid_corr <- evaluate(mid)
            if (!is.finite(mid_corr) || abs(mid_corr - target_corr) <= tolerance) {
                return(best_length())
            }
            if (mid_corr < target_corr) {
                lower <- mid
                lower_corr <- mid_corr
            } else {
                upper <- mid
                upper_corr <- mid_corr
            }
        }

        best_length()
    }

    simulate_latent_field <- function(chr_pos, n_samples, length_scale) {
        n_chr_sites <- length(chr_pos)
        z <- matrix(stats::rnorm(n_chr_sites * n_samples),
            nrow = n_chr_sites,
            ncol = n_samples
        )
        if (n_chr_sites <= 1L || !is.finite(length_scale) || length_scale <= 0) {
            return(z)
        }

        z[1, ] <- stats::rnorm(n_samples)
        phi <- exp(-pmax(diff(chr_pos), 1) / length_scale)
        innovation_sd <- sqrt(pmax(1 - phi^2, 1e-8))
        for (i in 2:n_chr_sites) {
            z[i, ] <- phi[i - 1L] * z[i - 1L, ] + innovation_sd[i - 1L] * z[i, ]
        }
        z
    }

    meth_signal <- matrix(NA_real_, nrow = n_sites, ncol = n_orig_samples)
    meth_signal[valid_cov] <- qlogis(meth_prop[valid_cov])
    cov_signal <- log1p(cov)

    meth_length <- estimate_length_scale(meth_signal, fallback_length = default_length)
    cov_length <- estimate_length_scale(cov_signal, fallback_length = default_length)

    # Marginal Poisson/Beta/Binomial sampling attenuates latent Gaussian
    # dependence, so calibrate on observed adjacent-pair correlations.
    if (calibrate_correlation && calibration_iterations > 0L) {
        calibration_pair_idx <- select_calibration_pairs(which(same_chr))
        if (length(calibration_pair_idx) >= 5L) {
            n_calibration_pairs <- length(calibration_pair_idx)
            cov_target_corr <- median_pair_correlation(adjacent_pair_correlations(cov_signal))
            meth_target_corr <- median_pair_correlation(adjacent_pair_correlations(meth_signal))

            cov_base_noise <- matrix(
                stats::rnorm(n_calibration_pairs * calibration_samples),
                nrow = n_calibration_pairs,
                ncol = calibration_samples
            )
            cov_innovation_noise <- matrix(
                stats::rnorm(n_calibration_pairs * calibration_samples),
                nrow = n_calibration_pairs,
                ncol = calibration_samples
            )
            cov_length <- calibrate_length_scale(
                cov_length,
                cov_target_corr,
                function(length_scale) {
                    evaluate_coverage_pair_correlation(
                        length_scale,
                        calibration_pair_idx,
                        cov_base_noise,
                        cov_innovation_noise
                    )
                },
                calibration_iterations
            )

            calibration_cov <- simulate_pair_coverage(
                cov_length,
                calibration_pair_idx,
                cov_base_noise,
                cov_innovation_noise
            )

            meth_base_noise <- matrix(
                stats::rnorm(n_calibration_pairs * calibration_samples),
                nrow = n_calibration_pairs,
                ncol = calibration_samples
            )
            meth_innovation_noise <- matrix(
                stats::rnorm(n_calibration_pairs * calibration_samples),
                nrow = n_calibration_pairs,
                ncol = calibration_samples
            )
            binom_u_1 <- matrix(
                stats::runif(n_calibration_pairs * calibration_samples),
                nrow = n_calibration_pairs,
                ncol = calibration_samples
            )
            binom_u_2 <- matrix(
                stats::runif(n_calibration_pairs * calibration_samples),
                nrow = n_calibration_pairs,
                ncol = calibration_samples
            )
            meth_length <- calibrate_length_scale(
                meth_length,
                meth_target_corr,
                function(length_scale) {
                    evaluate_methylation_pair_correlation(
                        length_scale,
                        calibration_pair_idx,
                        meth_base_noise,
                        meth_innovation_noise,
                        binom_u_1,
                        binom_u_2,
                        calibration_cov
                    )
                },
                calibration_iterations
            )
        }
    }

    original_names <- colnames(bsseq_filtered)
    if (is.null(original_names)) {
        original_names <- paste0("sample_", seq_len(n_orig_samples))
    }
    new_sample_names <- make.unique(c(
        original_names,
        paste0("synthetic_", seq_len(n_new_samples))
    ))
    new_sample_names <- tail(new_sample_names, n_new_samples)

    new_cov <- matrix(NA_real_,
        nrow = n_sites,
        ncol = n_new_samples,
        dimnames = list(NULL, new_sample_names)
    )
    new_meth <- matrix(NA_real_,
        nrow = n_sites,
        ncol = n_new_samples,
        dimnames = list(NULL, new_sample_names)
    )

    chr_groups <- split(seq_len(n_sites), chr)
    for (idx in chr_groups) {
        chr_pos <- pos[idx]
        cov_latent <- simulate_latent_field(chr_pos, n_new_samples, cov_length)
        meth_latent <- simulate_latent_field(chr_pos, n_new_samples, meth_length)

        cov_u <- pmin(pmax(stats::pnorm(cov_latent), 1e-10), 1 - 1e-10)
        meth_u <- pmin(pmax(stats::pnorm(meth_latent), 1e-10), 1 - 1e-10)

        new_cov[idx, ] <- matrix(
            stats::qpois(c(cov_u), lambda = rep(lambda_per_site[idx], n_new_samples)),
            nrow = length(idx),
            ncol = n_new_samples
        )
        new_cov[idx, ] <- pmax(new_cov[idx, ], 1L)

        new_meth[idx, ] <- matrix(
            stats::qbeta(
                c(meth_u),
                shape1 = rep(alpha_per_site[idx], n_new_samples),
                shape2 = rep(beta_per_site[idx], n_new_samples)
            ),
            nrow = length(idx),
            ncol = n_new_samples
        )
    }

    new_M <- matrix(
        stats::rbinom(
            n_sites * n_new_samples,
            size = c(new_cov),
            prob = c(new_meth)
        ),
        nrow = n_sites,
        ncol = n_new_samples,
        dimnames = list(NULL, new_sample_names)
    )
    storage.mode(new_cov) <- "integer"
    storage.mode(new_M) <- "integer"

    synthetic_bsseq <- bsseq::BSseq(
        chr = chr,
        pos = pos,
        M = new_M,
        Cov = new_cov,
        sampleNames = new_sample_names
    )

    # Combine original and synthetic
    combined <- bsseq::combineList(list(bsseq_filtered, synthetic_bsseq))

    combined
}
# nolint end
