#' Combine spatially correlated p-values with comb-p style SLK correction
#'
#' `combinePvalues()` is an R port of the comb-p Stouffer-Liptak-Kechris
#' (SLK) p-value correction step. For each locus, nearby loci within
#' `max_dist` are collected, an autocorrelation-derived covariance term is
#' computed from their genomic distances, and the local p-values are combined
#' with the SLK z-score formula used by comb-p.
#'
#' If `acf` is not supplied, lagged correlations are estimated from `pvals` and
#' `positions` using comb-p's default convention of Spearman correlations on
#' `-log10(p)`. If `acf` is supplied, it should contain columns named
#' `lag_min`, `lag_max`, and `correlation`, or have those values in its first
#' three columns.
#'
#' @param pvals Numeric vector of p-values.
#' @param positions Numeric vector of genomic positions, one per p-value.
#' @param chr Optional chromosome vector. If `NULL`, all positions are treated
#'   as coming from the same chromosome.
#' @param ends Optional end positions. Defaults to `positions`, which treats
#'   each locus as a single point.
#' @param max_dist Maximum distance used for local SLK neighborhoods and ACF
#'   estimation when `acf` is not supplied. Default is 1000.
#' @param lag_step Width of distance bins used to estimate the ACF when `acf`
#'   is not supplied. Default is 50.
#' @param min_dist Lower bound of the first ACF distance bin. Default is 0.
#' @param lags Optional numeric vector of lag breakpoints. If supplied, it
#'   overrides `min_dist`, `max_dist`, and `lag_step` for ACF estimation.
#' @param acf Optional precomputed ACF table with lag minimum, lag maximum, and
#'   correlation columns.
#' @param method Either `"slk"`/`"none"` to return SLK p-values, or any method
#'   accepted by [stats::p.adjust()] such as `"fdr"`/`"BH"` to adjust the SLK
#'   p-values.
#' @param correlation Correlation method used for ACF estimation. Default is
#'   `"spearman"`, matching comb-p.
#' @param transform Logical. If `TRUE`, estimate correlations on `-log10(p)`.
#'   Default is `TRUE`, matching comb-p.
#' @param partial Logical. If `TRUE`, estimate each lag interval from only
#'   pairs in that interval. If `FALSE`, outer intervals include pairs from
#'   shorter distances as in comb-p's non-partial ACF mode.
#' @param p_floor Minimum p-value used before q-normal transformation.
#' @param na_correlation Correlation value used for ACF bins whose correlation
#'   cannot be estimated. Default is 0.
#'
#' @return A numeric vector in the same order as `pvals`. The estimated or
#'   supplied ACF is attached as the `"acf"` attribute, and the unadjusted SLK
#'   p-values are attached as `"slk_pvalues"` when `method` requests an
#'   additional multiple-testing adjustment.
#'
#' @references
#' Pedersen BS, Schwartz DA, Yang IV, Kechris KJ. Comb-p: software for
#' combining, analyzing, grouping and correcting spatially correlated P-values.
#' Bioinformatics. 2012;28(22):2986-2988.
#'
#' @examples
#' p <- c(0.01, 0.03, 0.4, 0.02)
#' pos <- c(100, 150, 900, 940)
#' combinePvalues(p, pos, max_dist = 100)
#'
#' @export
combinePvalues <- function(
    pvals,
    positions,
    chr = NULL,
    ends = positions,
    max_dist = 1000,
    lag_step = 50,
    min_dist = 0,
    lags = NULL,
    acf = NULL,
    method = c("slk", "none", stats::p.adjust.methods),
    correlation = c("spearman", "pearson", "kendall"),
    transform = TRUE,
    partial = TRUE,
    p_floor = 9e-117,
    na_correlation = 0
) {
    method <- match.arg(method)
    correlation <- match.arg(correlation)

    if (!is.numeric(pvals)) {
        stop("pvals must be a numeric vector.")
    }
    if (!is.numeric(positions)) {
        stop("positions must be a numeric vector.")
    }
    if (length(pvals) != length(positions)) {
        stop("pvals and positions must have the same length.")
    }
    if (length(ends) != length(pvals)) {
        stop("ends must have the same length as pvals.")
    }
    if (is.null(chr)) {
        chr <- rep("chr", length(pvals))
    }
    if (length(chr) != length(pvals)) {
        stop("chr must have the same length as pvals.")
    }

    out <- rep(NA_real_, length(pvals))
    names(out) <- names(pvals)
    if (length(pvals) == 0L) {
        return(out)
    }

    invalid_p <- !is.na(pvals) & (pvals < 0 | pvals > 1)
    if (any(invalid_p)) {
        stop("pvals must be between 0 and 1.")
    }

    dat <- data.frame(
        .index = seq_along(pvals),
        chr = as.character(chr),
        start = as.numeric(positions),
        end = as.numeric(ends),
        p = as.numeric(pvals),
        stringsAsFactors = FALSE
    )
    dat <- dat[
        is.finite(dat$start) &
            is.finite(dat$end) &
            is.finite(dat$p) &
            !is.na(dat$chr),
        ,
        drop = FALSE
    ]
    if (nrow(dat) == 0L) {
        return(out)
    }
    swap_idx <- dat$end < dat$start
    if (any(swap_idx)) {
        tmp <- dat$start[swap_idx]
        dat$start[swap_idx] <- dat$end[swap_idx]
        dat$end[swap_idx] <- tmp
    }
    dat$p <- .clipCombPvalues(dat$p, p_floor)
    dat <- dat[order(dat$chr, dat$start, dat$end, dat$.index), , drop = FALSE]

    if (is.null(acf)) {
        acf_tbl <- .estimateCombPAcf(
            dat = dat,
            max_dist = max_dist,
            lag_step = lag_step,
            min_dist = min_dist,
            lags = lags,
            correlation = correlation,
            transform = transform,
            partial = partial,
            na_correlation = na_correlation
        )
    } else {
        acf_tbl <- .normaliseCombPAcf(acf, na_correlation = na_correlation)
    }
    if (nrow(acf_tbl) == 0L) {
        stop("acf must contain at least one lag interval.")
    }
    acf_tbl <- acf_tbl[order(acf_tbl$lag_min, acf_tbl$lag_max), , drop = FALSE]

    slk <- .slkCombP(dat, acf_tbl)
    out[dat$.index] <- slk
    attr(out, "acf") <- acf_tbl

    if (method %in% c("slk", "none")) {
        return(out)
    }

    adjusted <- stats::p.adjust(out, method = method)
    names(adjusted) <- names(out)
    attr(adjusted, "acf") <- acf_tbl
    attr(adjusted, "slk_pvalues") <- out
    adjusted
}

.clipCombPvalues <- function(pvals, p_floor) {
    p_floor <- as.numeric(p_floor)
    if (!is.finite(p_floor) || p_floor <= 0 || p_floor >= 1) {
        stop("p_floor must be a finite value between 0 and 1.")
    }
    pmax(pmin(pvals, 1 - 9e-16), p_floor)
}

.makeCombPLags <- function(max_dist, lag_step, min_dist, lags) {
    if (!is.null(lags)) {
        lags <- sort(unique(as.numeric(lags)))
        if (length(lags) < 2L || any(!is.finite(lags))) {
            stop("lags must contain at least two finite breakpoints.")
        }
        return(lags)
    }

    max_dist <- as.numeric(max_dist)
    lag_step <- as.numeric(lag_step)
    min_dist <- as.numeric(min_dist)
    if (!is.finite(max_dist) || max_dist <= 0) {
        stop("max_dist must be a positive finite value.")
    }
    if (!is.finite(lag_step) || lag_step <= 0) {
        stop("lag_step must be a positive finite value.")
    }
    if (!is.finite(min_dist) || min_dist < 0) {
        stop("min_dist must be a non-negative finite value.")
    }
    if (min_dist >= max_dist) {
        stop("min_dist must be less than max_dist.")
    }

    lags <- seq(min_dist, max_dist, by = lag_step)
    if (tail(lags, 1L) < max_dist) {
        lags <- c(lags, max_dist)
    }
    unique(lags)
}

.estimateCombPAcf <- function(
    dat,
    max_dist,
    lag_step,
    min_dist,
    lags,
    correlation,
    transform,
    partial,
    na_correlation
) {
    lag_breaks <- .makeCombPLags(max_dist, lag_step, min_dist, lags)
    lag_min <- head(lag_breaks, -1L)
    lag_max <- tail(lag_breaks, -1L)
    max_lag <- max(lag_max)
    by_chr <- split(dat, dat$chr, drop = TRUE)

    n_pairs <- integer(length(lag_min))
    for (chrom_dat in by_chr) {
        n <- nrow(chrom_dat)
        if (n < 2L) {
            next
        }
        for (i in seq_len(n - 1L)) {
            j <- i + 1L
            while (j <= n) {
                dist <- chrom_dat$start[j] - chrom_dat$end[i]
                if (dist > max_lag) {
                    break
                }
                if (dist < 0) {
                    dist <- 0
                }
                bin <- which(lag_min <= dist & dist < lag_max)
                if (length(bin) > 0L) {
                    n_pairs[bin[[1L]]] <- n_pairs[bin[[1L]]] + 1L
                }
                j <- j + 1L
            }
        }
    }

    pair_x <- lapply(n_pairs, numeric)
    pair_y <- lapply(n_pairs, numeric)
    fill_idx <- integer(length(lag_min))
    for (chrom_dat in by_chr) {
        n <- nrow(chrom_dat)
        if (n < 2L) {
            next
        }
        for (i in seq_len(n - 1L)) {
            j <- i + 1L
            while (j <= n) {
                dist <- chrom_dat$start[j] - chrom_dat$end[i]
                if (dist > max_lag) {
                    break
                }
                if (dist < 0) {
                    dist <- 0
                }
                bin <- which(lag_min <= dist & dist < lag_max)
                if (length(bin) > 0L) {
                    bin <- bin[[1L]]
                    fill_idx[bin] <- fill_idx[bin] + 1L
                    pair_x[[bin]][fill_idx[bin]] <- chrom_dat$p[i]
                    pair_y[[bin]][fill_idx[bin]] <- chrom_dat$p[j]
                }
                j <- j + 1L
            }
        }
    }

    corr <- rep(na_correlation, length(pair_x))
    cumulative_x <- numeric()
    cumulative_y <- numeric()
    for (i in seq_along(pair_x)) {
        x <- pair_x[[i]]
        y <- pair_y[[i]]
        if (!partial) {
            cumulative_x <- c(cumulative_x, x)
            cumulative_y <- c(cumulative_y, y)
            x <- cumulative_x
            y <- cumulative_y
        }
        n_pairs[i] <- length(x)
        if (length(x) >= 2L && stats::sd(x) > 0 && stats::sd(y) > 0) {
            if (transform) {
                x <- -log10(pmax(x, 1e-12))
                y <- -log10(pmax(y, 1e-12))
            }
            corr_i <- suppressWarnings(stats::cor(x, y, method = correlation))
            if (is.finite(corr_i)) {
                corr[i] <- corr_i
            }
        }
    }

    data.frame(
        lag_min = lag_min,
        lag_max = lag_max,
        correlation = corr,
        n = n_pairs,
        stringsAsFactors = FALSE
    )
}

.normaliseCombPAcf <- function(acf, na_correlation) {
    acf <- as.data.frame(acf, stringsAsFactors = FALSE)
    if (ncol(acf) < 3L) {
        stop("acf must have at least three columns: lag_min, lag_max, correlation.")
    }
    nms <- names(acf)
    if (all(c("lag_min", "lag_max", "correlation") %in% nms)) {
        out <- acf[, c("lag_min", "lag_max", "correlation"), drop = FALSE]
    } else {
        out <- acf[, seq_len(3L), drop = FALSE]
        names(out) <- c("lag_min", "lag_max", "correlation")
    }
    out$lag_min <- as.numeric(out$lag_min)
    out$lag_max <- as.numeric(out$lag_max)
    out$correlation <- as.numeric(out$correlation)
    if (any(!is.finite(out$lag_min)) || any(!is.finite(out$lag_max))) {
        stop("acf lag_min and lag_max values must be finite.")
    }
    if (any(out$lag_min > out$lag_max)) {
        stop("acf lag_min values must be less than or equal to lag_max values.")
    }
    out$correlation[!is.finite(out$correlation)] <- na_correlation
    out
}

.slkCombP <- function(dat, acf_tbl) {
    lag_max <- max(acf_tbl$lag_max)
    result <- numeric(nrow(dat))
    by_chr_idx <- split(seq_len(nrow(dat)), dat$chr, drop = TRUE)
    for (idx in by_chr_idx) {
        chrom_dat <- dat[idx, , drop = FALSE]
        n <- nrow(chrom_dat)
        imin <- 1L
        imax <- 1L
        for (i in seq_len(n)) {
            while (imin < n && chrom_dat$start[i] - chrom_dat$end[imin] > lag_max) {
                imin <- imin + 1L
            }
            if (imax < i) {
                imax <- i
            }
            while (imax <= n && chrom_dat$start[imax] - chrom_dat$end[i] < lag_max) {
                imax <- imax + 1L
            }
            neighbors <- imin:(imax - 1L)
            result[idx[i]] <- .zScoreCombineCombP(
                pvals = chrom_dat$p[neighbors],
                starts = chrom_dat$start[neighbors],
                ends = chrom_dat$end[neighbors],
                acf_tbl = acf_tbl
            )
        }
    }
    result
}

.zScoreCombineCombP <- function(pvals, starts, ends, acf_tbl) {
    n <- length(pvals)
    if (n == 1L) {
        return(pvals[[1L]])
    }

    z <- mean(stats::qnorm(pvals, lower.tail = FALSE))
    corr_sum <- 0
    for (i in seq_len(n - 1L)) {
        for (j in (i + 1L):n) {
            dist <- starts[j] - ends[i]
            if (dist < 0) {
                dist <- 0
            }
            corr_sum <- corr_sum + .getCombPCorrelation(dist, acf_tbl)
        }
    }

    variance_factor <- n + 2 * corr_sum
    if (!is.finite(variance_factor) || variance_factor <= 0) {
        return(NA_real_)
    }
    scale <- sqrt(variance_factor) / n
    stats::pnorm(z / scale, lower.tail = FALSE)
}

.getCombPCorrelation <- function(dist, acf_tbl) {
    if (dist < acf_tbl$lag_min[[1L]]) {
        return(acf_tbl$correlation[[1L]])
    }
    idx <- which(acf_tbl$lag_min <= dist & dist <= acf_tbl$lag_max)
    if (length(idx) == 0L) {
        return(0)
    }
    acf_tbl$correlation[[idx[[1L]]]]
}
