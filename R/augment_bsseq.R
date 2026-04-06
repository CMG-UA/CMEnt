# nolint start: object_name_linter

#' Augment BSseq Object
#'
#' Generate synthetic samples while preserving both per-site coverage and
#' methylation marginals and the local correlation structure of neighboring CpGs.
#' Coverage and methylation are first fit with per-site Poisson/Beta marginals,
#' then synthetic samples are drawn through chromosome-local Gaussian copula
#' fields whose dependence is estimated from adjacent CpGs in the observed data.
#'
#' @param bs A BSseq object
#' @param n_new_samples Number of new synthetic samples to generate
#' @param seed Random seed for reproducibility
#' @param min_samples Minimum number of samples with coverage required per site
#' @importFrom SummarizedExperiment assays
#' @return A BSseq object with original and synthetic samples
#' @export
#' @examples
#' \dontrun{
#' # Load example BSseq data
#' data("BSobj", package = "bsseq")
#' # Augment with 5 synthetic samples
#' augmented_bs <- augmentBSSeq(BSobj, n_new_samples = 5, seed = 123)
#' }
augmentBSSeq <- function(bs, n_new_samples, seed = NULL, min_samples = 2) {
    if (!is.null(seed)) set.seed(seed)
    n_new_samples <- as.integer(n_new_samples)
    min_samples <- as.integer(min_samples)

    if (length(n_new_samples) != 1L || is.na(n_new_samples) || n_new_samples < 1L) {
        stop("'n_new_samples' must be a positive integer")
    }
    if (length(min_samples) != 1L || is.na(min_samples) || min_samples < 1L) {
        stop("'min_samples' must be a positive integer")
    }

    # Keep only sites with coverage in at least min_samples
    cov_matrix <- bsseq::getCoverage(bs)
    valid_sites <- rowSums(cov_matrix > 0) >= min_samples
    bsseq_filtered <- bs[valid_sites, ]
    # Keep assay-level sample dimnames aligned with colData before combining objects.
    colnames(assays(bsseq_filtered)$M) <- colnames(bsseq_filtered)
    colnames(assays(bsseq_filtered)$Cov) <- colnames(bsseq_filtered)

    if (nrow(bsseq_filtered) == 0L) {
        stop("No CpG sites have coverage in at least 'min_samples' samples")
    }

    # Extract coverage and methylation counts
    cov <- as.matrix(bsseq::getCoverage(bsseq_filtered))
    M <- as.matrix(bsseq::getCoverage(bsseq_filtered, type = "M"))

    n_sites <- nrow(bsseq_filtered)
    n_orig_samples <- ncol(cov)
    chr <- as.character(seqnames(bsseq_filtered))
    pos <- start(bsseq_filtered)

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

    estimate_length_scale <- function(signal, fallback_length) {
        pair_idx <- which(same_chr)
        if (length(pair_idx) == 0L || n_orig_samples < 2L) {
            return(fallback_length)
        }

        pair_corr <- vapply(pair_idx, function(i) {
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

    original_names <- colnames(cov)
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
