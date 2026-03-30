# nolint start: object_name_linter

#' Augment BSseq Object (Per-Site)
#'
#' Generate synthetic samples by sampling from per-site coverage and methylation
#' distributions. For each CpG site, coverage is sampled from a Poisson distribution
#' fitted to that site's coverage values, and methylation is sampled from a Beta
#' distribution fit with sample-size-aware shrinkage to avoid over-dispersion when
#' only a few input samples are available.
#'
#' @param bs A BSseq object
#' @param n_new_samples Number of new synthetic samples to generate
#' @param seed Random seed for reproducibility
#' @param min_samples Minimum number of samples with coverage required per site
#' @return A BSseq object with original and synthetic samples
#' @export
#' @examples
#' \dontrun{
#' # Load example BSseq data
#' data("BSobj", package = "bsseq")
#' # Augment with 5 synthetic samples
#' augmented_bs <- augmenBSSeq(BSobj, n_new_samples = 5, seed = 123)
#' }
augmenBSSeq <- function(bs, n_new_samples, seed = NULL, min_samples = 2) {
    if (!is.null(seed)) set.seed(seed)

    # Keep only sites with coverage in at least min_samples
    cov_matrix <- bsseq::getCoverage(bs)
    valid_sites <- rowSums(cov_matrix > 0) >= min_samples
    bsseq_filtered <- bs[valid_sites, ]

    # Extract coverage and methylation counts
    cov <- bsseq::getCoverage(bsseq_filtered)
    M <- bsseq::getCoverage(bsseq_filtered, type = "M")

    n_sites <- nrow(bsseq_filtered)
    n_orig_samples <- ncol(cov)

    # Coverage: compute mean per site (excluding zeros)
    cov[cov == 0] <- NA
    lambda_per_site <- rowMeans(cov, na.rm = TRUE)
    lambda_per_site[is.na(lambda_per_site)] <- 1 # fallback

    # Methylation: use smoothed proportions so boundary values (0/1) are retained.
    valid_cov <- !is.na(cov)
    meth_prop <- matrix(NA_real_, nrow = n_sites, ncol = n_orig_samples)
    meth_prop[valid_cov] <- (M[valid_cov] + 0.5) / (cov[valid_cov] + 1)

    # Coverage-weighted site means with Jeffreys-style pseudocounts.
    pseudo_total <- rowSums((cov + 1) * valid_cov, na.rm = TRUE)
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

    # Generate synthetic samples (vectorized)
    # Coverage: sample from Poisson with per-site lambda
    new_cov <- matrix(
        rpois(n_sites * n_new_samples, lambda = rep(lambda_per_site, n_new_samples)),
        nrow = n_sites, ncol = n_new_samples
    )
    new_cov <- pmax(new_cov, 1L) # Ensure at least 1 coverage

    # Methylation: sample from Beta with per-site parameters
    new_meth <- matrix(
        rbeta(n_sites * n_new_samples,
            shape1 = rep(alpha_per_site, n_new_samples),
            shape2 = rep(beta_per_site, n_new_samples)
        ),
        nrow = n_sites, ncol = n_new_samples
    )

    # Convert methylation to counts
    new_M <- round(new_meth * new_cov)

    # Create new BSseq object with synthetic samples
    new_sample_names <- paste0("synthetic_", seq_len(n_new_samples))

    synthetic_bsseq <- bsseq::BSseq(
        chr = as.character(seqnames(bsseq_filtered)),
        pos = start(bsseq_filtered),
        M = new_M,
        Cov = new_cov,
        sampleNames = new_sample_names
    )

    # Combine original and synthetic
    combined <- bsseq::combineList(list(bsseq_filtered, synthetic_bsseq))

    combined
}
# nolint end