#' @keywords internal
#' @noRd
.buildStratifiedFolds <- function(groups, nfold = getOption("DMRsegal.ranking_nfold", 5)) {
    groups <- as.factor(groups)
    if (nlevels(groups) < 2) {
        stop("Ranking requires at least two classes in '__casecontrol__'.")
    }
    set.seed(getOption("DMRsegal.random_seed", 42))
    group_folds <- split(seq_along(groups), groups)
    for (g in names(group_folds)) {
        if (length(group_folds[[g]]) < nfold) {
            gsize <- length(group_folds[[g]])
            stop(paste0(
                "Number of samples in group (", gsize, ") '", g, "' is less than nfold = ", nfold,
                ". Cannot perform stratified cross-prediction. Reduce nfold using options(DMRsegal.ranking_nfold=", gsize,
                ") or increase number of samples in this group."
            ))
        }
    }
    folds <- integer(length(groups))
    for (g in names(group_folds)) {
        idx <- group_folds[[g]]
        folds[idx] <- sample(rep(seq_len(nfold), length.out = length(idx)))
    }
    folds
}

#' @keywords internal
#' @noRd
.performCrossPrediction <- function(beta_mat, groups, folds = NULL, nfold = getOption("DMRsegal.ranking_nfold", 5)) {
    groups <- as.factor(groups)
    if (ncol(beta_mat) != length(groups)) {
        stop(
            "Mismatch between beta matrix columns (", ncol(beta_mat),
            ") and group labels (", length(groups), ")."
        )
    }
    if (nlevels(groups) < 2) {
        stop("Ranking requires at least two classes in '__casecontrol__'.")
    }
    if (is.null(folds)) {
        folds <- .buildStratifiedFolds(groups, nfold = nfold)
    } else {
        if (length(folds) != ncol(beta_mat)) {
            stop(
                "Mismatch between folds length (", length(folds),
                ") and number of samples (", ncol(beta_mat), ")."
            )
        }
        folds <- as.integer(folds)
        if (anyNA(folds) || any(folds < 1L)) {
            stop("Fold IDs must be positive integers without NA values.")
        }
        nfold <- max(folds)
    }

    beta_mat_t <- t(beta_mat)
    predictions <- vector("character", ncol(beta_mat))
    decision_values <- rep(NA_real_, ncol(beta_mat))
    groups_chr <- as.character(groups)

    for (fold in seq_len(nfold)) {
        test_indices <- which(folds == fold)
        if (length(test_indices) == 0L) {
            next
        }
        train_indices <- which(folds != fold)
        train_groups <- as.factor(groups[train_indices])
        model <- e1071::svm(
            beta_mat_t[train_indices, , drop = FALSE],
            train_groups,
            kernel = "radial",
            scale = TRUE
        )
        fold_pred <- predict(model, beta_mat_t[test_indices, , drop = FALSE], decision.values = TRUE)
        fold_pred_chr <- as.character(fold_pred)
        predictions[test_indices] <- fold_pred_chr

        fold_decision <- as.numeric(attr(fold_pred, "decision.values"))
        if (length(fold_decision) != length(test_indices)) {
            fold_decision <- rep(NA_real_, length(test_indices))
        } else {
            train_levels <- levels(train_groups)
            pred_from_sign <- ifelse(fold_decision >= 0, train_levels[1], train_levels[2])
            if (mean(pred_from_sign == fold_pred_chr, na.rm = TRUE) < 0.5) {
                fold_decision <- -fold_decision
            }
        }
        decision_values[test_indices] <- fold_decision
    }
    cv_accuracy <- mean(predictions == groups_chr)

    # Margin-sensitive score (in (0, 1]) based on logistic loss of SVM decision values.
    # This breaks ties among perfect-accuracy DMRs by rewarding larger separating margins.
    y <- ifelse(groups_chr == levels(groups)[1], 1, -1)
    finite_mask <- is.finite(decision_values)
    if (!any(finite_mask)) {
        margin_score <- cv_accuracy
    } else {
        logistic_loss <- log1p(exp(-y[finite_mask] * decision_values[finite_mask]))
        margin_score <- exp(-mean(logistic_loss))
    }
    c(score = margin_score, cv_accuracy = cv_accuracy)
}

#' @keywords internal
#' @noRd
.gaussianSmoothAdaptiveKNN <- function(x, y, k = 5L) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    n <- length(x)
    if (n == 0L) {
        return(numeric(0))
    }
    if (n == 1L) {
        return(y)
    }
    if (length(y) != n) {
        stop("x and y must have the same length.")
    }
    if (any(!is.finite(x)) || any(!is.finite(y))) {
        stop("x and y must contain only finite numeric values.")
    }

    k <- as.integer(k)
    if (!is.finite(k) || k < 1L) {
        k <- 1L
    }
    k <- min(k, n)

    unique_x <- sort(unique(x))
    dx <- diff(unique_x)
    min_pos_dx <- suppressWarnings(min(dx[dx > 0], na.rm = TRUE))
    if (!is.finite(min_pos_dx) || min_pos_dx <= 0) {
        min_pos_dx <- 1
    }

    smoothed <- numeric(n)
    for (i in seq_len(n)) {
        d <- abs(x - x[i])
        knearest <- order(d)[seq_len(k)]
        local_d <- d[knearest]
        h <- max(local_d, na.rm = TRUE)
        if (!is.finite(h) || h <= 0) {
            h <- min_pos_dx
        }
        # Use a tighter local kernel on k-nearest neighbors to reduce over-smoothing.
        h <- max(h * 0.75, min_pos_dx)
        w <- exp(-0.5 * (local_d / h)^2)
        w_sum <- sum(w)
        if (!is.finite(w_sum) || w_sum <= 0) {
            smoothed[i] <- y[i]
        } else {
            smoothed[i] <- sum(w * y[knearest]) / w_sum
        }
    }
    smoothed
}

#' @keywords internal
#' @noRd
.buildLinearCostPrefixes <- function(x, y) {
    list(
        sx = c(0, cumsum(x)),
        sy = c(0, cumsum(y)),
        sxx = c(0, cumsum(x * x)),
        sxy = c(0, cumsum(x * y)),
        syy = c(0, cumsum(y * y))
    )
}

#' @keywords internal
#' @noRd
.linearSegmentCost <- function(prefixes, start_idx, end_idx) {
    if (end_idx < start_idx) {
        return(0)
    }
    n <- end_idx - start_idx + 1
    s <- start_idx
    e <- end_idx
    sx <- prefixes$sx[e + 1L] - prefixes$sx[s]
    sy <- prefixes$sy[e + 1L] - prefixes$sy[s]
    sxx <- prefixes$sxx[e + 1L] - prefixes$sxx[s]
    sxy <- prefixes$sxy[e + 1L] - prefixes$sxy[s]
    syy <- prefixes$syy[e + 1L] - prefixes$syy[s]
    if (n <= 1L) {
        return(0)
    }

    tss <- syy - (sy * sy) / n
    if (!is.finite(tss)) {
        tss <- 0
    }
    tss <- max(tss, 0)
    denom <- sxx - (sx * sx) / n
    if (!is.finite(denom) || abs(denom) <= 1e-12) {
        return(tss)
    }
    numer <- sxy - (sx * sy) / n
    rss <- tss - (numer * numer) / denom
    if (!is.finite(rss)) {
        rss <- tss
    }
    max(rss, 0)
}

#' @keywords internal
#' @noRd
.estimateSegmentSlope <- function(x, y) {
    n <- length(x)
    if (n <= 1L) {
        return(0)
    }
    x <- as.numeric(x)
    y <- as.numeric(y)
    x_centered <- x - mean(x)
    denom <- sum(x_centered * x_centered)
    if (!is.finite(denom) || denom <= 1e-12) {
        return(0)
    }
    sum(x_centered * (y - mean(y))) / denom
}

#' @keywords internal
#' @noRd
.segmentLinearPELT <- function(x, y, min_segment_size = 2L, penalty = NULL) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    n <- length(x)
    if (n == 0L) {
        return(list(
            starts = integer(0),
            ends = integer(0),
            slopes = numeric(0),
            segment_id = integer(0)
        ))
    }
    if (length(y) != n) {
        stop("x and y must have the same length.")
    }
    if (n == 1L) {
        return(list(
            starts = 1L,
            ends = 1L,
            slopes = 0,
            segment_id = 1L
        ))
    }

    min_segment_size <- as.integer(min_segment_size)
    if (!is.finite(min_segment_size) || min_segment_size < 1L) {
        min_segment_size <- 1L
    }

    if (is.null(penalty)) {
        y_var <- stats::var(y)
        if (!is.finite(y_var) || y_var <= 0) {
            y_var <- 1e-8
        }
        penalty <- (log(max(n, 2)) / max(n, 2)) * y_var
    }
    penalty <- max(as.numeric(penalty), 0)

    prefixes <- .buildLinearCostPrefixes(x, y)

    # F[t + 1] stores optimal cost up to index t; t = 0..n
    F <- rep(Inf, n + 1L)
    F[1] <- -penalty
    cps <- vector("list", n + 1L)
    R <- 0L
    K <- 0

    for (t in seq_len(n)) {
        eligible <- R[R <= (t - min_segment_size)]
        ineligible <- R[R > (t - min_segment_size)]

        costs <- numeric(0)
        finite_mask <- logical(0)
        if (length(eligible) > 0L) {
            costs <- vapply(eligible, function(tau) {
                prev_cost <- F[tau + 1L]
                if (!is.finite(prev_cost)) {
                    return(Inf)
                }
                prev_cost + .linearSegmentCost(prefixes, tau + 1L, t) + penalty
            }, numeric(1))
            finite_mask <- is.finite(costs)
        }

        if (length(costs) > 0L && any(finite_mask)) {
            best_local <- which.min(costs)
            best_tau <- eligible[best_local]
            F[t + 1L] <- costs[best_local]
            cps[[t + 1L]] <- c(cps[[best_tau + 1L]], best_tau)

            keep_flags <- rep(FALSE, length(eligible))
            keep_flags[finite_mask] <- costs[finite_mask] <= (F[t + 1L] + K)
            keep_eligible <- eligible[keep_flags]
        } else {
            keep_eligible <- eligible
        }

        R <- sort(unique(c(ineligible, keep_eligible, t)))
    }

    best_cps <- cps[[n + 1L]]
    if (length(best_cps) == 0L) {
        ends <- n
    } else {
        best_cps <- best_cps[best_cps > 0L & best_cps < n]
        ends <- c(best_cps, n)
    }
    starts <- c(1L, head(ends, -1L) + 1L)
    nseg <- length(starts)

    slopes <- numeric(nseg)
    segment_id <- integer(n)
    for (s in seq_len(nseg)) {
        idx <- starts[s]:ends[s]
        slopes[s] <- .estimateSegmentSlope(x[idx], y[idx])
        segment_id[idx] <- s
    }

    list(
        starts = starts,
        ends = ends,
        slopes = slopes,
        segment_id = segment_id
    )
}

#' @keywords internal
#' @noRd
.computeChromosomeGapThreshold <- function(x_chr,
                                           mode = c("adaptive", "fixed", "none"),
                                           fixed_bp = NULL,
                                           quantile = 0.95,
                                           multiplier = 1.5,
                                           min_bp = 250000,
                                           max_bp = 5000000) {
    mode <- strex::match_arg(mode, ignore_case = TRUE)
    x_chr <- sort(unique(as.numeric(x_chr)))
    x_chr <- x_chr[is.finite(x_chr)]
    min_bp <- max(as.numeric(min_bp), 0)
    max_bp <- max(as.numeric(max_bp), min_bp)
    if (!is.finite(quantile) || quantile <= 0 || quantile >= 1) {
        stop("quantile must be a numeric value in (0, 1).")
    }
    if (!is.finite(multiplier) || multiplier <= 0) {
        stop("multiplier must be a positive numeric value.")
    }

    if (mode == "none") {
        return(Inf)
    }
    if (mode == "fixed") {
        if (is.null(fixed_bp) || !is.finite(fixed_bp) || fixed_bp <= 0) {
            stop("For block_gap_mode='fixed', block_gap_fixed_bp must be a positive numeric value.")
        }
        return(as.numeric(fixed_bp))
    }

    gaps <- diff(x_chr)
    gaps <- gaps[is.finite(gaps) & gaps > 0]
    if (length(gaps) == 0L) {
        return(max_bp)
    }
    qgap <- as.numeric(stats::quantile(gaps, probs = quantile, names = FALSE, na.rm = TRUE, type = 8))
    if (!is.finite(qgap) || qgap <= 0) {
        qgap <- max(stats::median(gaps), 1)
    }
    threshold <- qgap * multiplier
    threshold <- max(min_bp, min(max_bp, threshold))
    as.numeric(threshold)
}

#' @keywords internal
#' @noRd
.findSlopeCandidates <- function(seg, d) {
    if (length(seg$slopes) < 2L) {
        return(data.frame(candidate_id = integer(0), seg_start = integer(0), seg_end = integer(0)))
    }
    seg_i <- 1L
    nseg <- length(seg$slopes)
    candidates <- list()
    candidate_id <- 0L
    while (seg_i <= (nseg - 1L)) {
        slope_a <- seg$slopes[seg_i]
        if (!is.finite(slope_a) || !(slope_a > d)) {
            seg_i <- seg_i + 1L
            next
        }
        seg_j <- seg_i + 1L
        while (seg_j <= nseg && is.finite(seg$slopes[seg_j]) && abs(seg$slopes[seg_j]) <= d) {
            seg_j <- seg_j + 1L
        }
        if (seg_j <= nseg && is.finite(seg$slopes[seg_j]) && seg$slopes[seg_j] < -d) {
            candidate_id <- candidate_id + 1L
            candidates[[candidate_id]] <- data.frame(
                candidate_id = candidate_id,
                seg_start = seg_i,
                seg_end = seg_j
            )
            seg_i <- seg_j + 1L
        } else {
            seg_i <- seg_i + 1L
        }
    }
    if (length(candidates) == 0L) {
        return(data.frame(candidate_id = integer(0), seg_start = integer(0), seg_end = integer(0)))
    }
    do.call(rbind, candidates)
}

#' @keywords internal
#' @noRd
.splitCandidateByGaps <- function(x_candidate, gap_threshold_bp) {
    x_candidate <- as.numeric(x_candidate)
    n <- length(x_candidate)
    if (n == 0L) {
        return(list(
            ranges = data.frame(local_start = integer(0), local_end = integer(0)),
            split_events = data.frame(split_after_index = integer(0), left_bp = numeric(0), right_bp = numeric(0), gap_bp = numeric(0))
        ))
    }
    if (n == 1L || !is.finite(gap_threshold_bp)) {
        return(list(
            ranges = data.frame(local_start = 1L, local_end = n),
            split_events = data.frame(split_after_index = integer(0), left_bp = numeric(0), right_bp = numeric(0), gap_bp = numeric(0))
        ))
    }

    gaps <- diff(x_candidate)
    split_points <- which(gaps > gap_threshold_bp)
    starts <- c(1L, split_points + 1L)
    ends <- c(split_points, n)
    split_events <- data.frame(
        split_after_index = split_points,
        left_bp = x_candidate[split_points],
        right_bp = x_candidate[split_points + 1L],
        gap_bp = gaps[split_points]
    )
    list(
        ranges = data.frame(local_start = starts, local_end = ends),
        split_events = split_events
    )
}

#' @keywords internal
#' @noRd
.acceptClusterBySlopeOrder <- function(cluster_segment_ids, seg_slopes, d) {
    ids <- sort(unique(as.integer(cluster_segment_ids)))
    ids <- ids[is.finite(ids) & ids >= 1L & ids <= length(seg_slopes)]
    if (length(ids) == 0L) {
        return(FALSE)
    }
    s <- as.numeric(seg_slopes[ids])
    pos_idx <- which(is.finite(s) & s > d)
    neg_idx <- which(is.finite(s) & s < -d)
    if (length(pos_idx) == 0L || length(neg_idx) == 0L) {
        return(FALSE)
    }
    any(vapply(pos_idx, function(p) any(neg_idx > p), logical(1)))
}

#' @keywords internal
#' @noRd
.computeDMRBlockFormationForChromosome <- function(chr,
                                                   chr_idx,
                                                   x_chr,
                                                   y_chr,
                                                   k_neighbors = 5L,
                                                   min_segment_size = 2L,
                                                   block_gap_mode = "adaptive",
                                                   block_gap_fixed_bp = NULL,
                                                   block_gap_quantile = 0.95,
                                                   block_gap_multiplier = 1.5,
                                                   block_gap_min_bp = 2500,
                                                   block_gap_max_bp = 50000) {
    f_chr <- .gaussianSmoothAdaptiveKNN(x_chr, y_chr, k = k_neighbors)
    seg <- .segmentLinearPELT(
        x = x_chr,
        y = f_chr,
        min_segment_size = min_segment_size
    )

    if (length(seg$segment_id) != length(x_chr)) {
        seg$segment_id <- rep(NA_integer_, length(x_chr))
        seg$slopes <- numeric(0)
        seg$starts <- integer(0)
        seg$ends <- integer(0)
    }

    slope_threshold <- stats::mad(seg$slopes, center = stats::median(seg$slopes), na.rm = TRUE)
    if (!is.finite(slope_threshold) || slope_threshold < 0) {
        slope_threshold <- 0
    }

    gap_threshold_bp <- .computeChromosomeGapThreshold(
        x_chr = x_chr,
        mode = block_gap_mode,
        fixed_bp = block_gap_fixed_bp,
        quantile = block_gap_quantile,
        multiplier = block_gap_multiplier,
        min_bp = block_gap_min_bp,
        max_bp = block_gap_max_bp
    )

    candidates_df <- .findSlopeCandidates(seg, d = slope_threshold)
    if (nrow(candidates_df) > 0) {
        cand_ranges <- lapply(seq_len(nrow(candidates_df)), function(i) {
            seg_start <- candidates_df$seg_start[i]
            seg_end <- candidates_df$seg_end[i]
            idx <- which(seg$segment_id >= seg_start & seg$segment_id <= seg_end)
            data.frame(
                dmr_start_idx = if (length(idx) > 0) min(idx) else NA_integer_,
                dmr_end_idx = if (length(idx) > 0) max(idx) else NA_integer_,
                start_bp = if (length(idx) > 0) min(x_chr[idx]) else NA_real_,
                end_bp = if (length(idx) > 0) max(x_chr[idx]) else NA_real_
            )
        })
        cand_ranges <- do.call(rbind, cand_ranges)
        candidates_df <- cbind(candidates_df, cand_ranges)
    }

    split_events_list <- list()
    blocks_list <- list()
    block_local_ids <- rep(NA_integer_, length(x_chr))
    local_block_num <- 0L

    if (nrow(candidates_df) > 0) {
        for (i in seq_len(nrow(candidates_df))) {
            seg_start <- candidates_df$seg_start[i]
            seg_end <- candidates_df$seg_end[i]
            candidate_id <- candidates_df$candidate_id[i]
            candidate_idx <- which(seg$segment_id >= seg_start & seg$segment_id <= seg_end)
            if (length(candidate_idx) == 0L) {
                next
            }

            split_result <- .splitCandidateByGaps(x_chr[candidate_idx], gap_threshold_bp)
            if (nrow(split_result$split_events) > 0) {
                ev <- split_result$split_events
                ev$candidate_id <- candidate_id
                ev$split_after_global_idx <- candidate_idx[ev$split_after_index]
                split_events_list[[length(split_events_list) + 1L]] <- ev
            }

            for (r in seq_len(nrow(split_result$ranges))) {
                cluster_idx <- candidate_idx[split_result$ranges$local_start[r]:split_result$ranges$local_end[r]]
                cluster_segment_ids <- sort(unique(seg$segment_id[cluster_idx]))
                accepted <- .acceptClusterBySlopeOrder(cluster_segment_ids, seg$slopes, d = slope_threshold)
                if (!accepted) {
                    next
                }
                local_block_num <- local_block_num + 1L
                block_local_ids[cluster_idx] <- local_block_num
                cluster_x <- x_chr[cluster_idx]
                cluster_gaps <- diff(cluster_x)
                blocks_list[[length(blocks_list) + 1L]] <- data.frame(
                    block_local_id = local_block_num,
                    candidate_id = candidate_id,
                    dmr_start_idx = min(cluster_idx),
                    dmr_end_idx = max(cluster_idx),
                    start_bp = min(cluster_x),
                    end_bp = max(cluster_x),
                    n_dmrs = length(cluster_idx),
                    span_bp = max(cluster_x) - min(cluster_x),
                    max_gap_bp = if (length(cluster_gaps) > 0) max(cluster_gaps) else 0,
                    segment_ids = paste(cluster_segment_ids, collapse = ","),
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    split_events_df <- if (length(split_events_list) > 0) {
        do.call(rbind, split_events_list)
    } else {
        data.frame(
            split_after_index = integer(0),
            left_bp = numeric(0),
            right_bp = numeric(0),
            gap_bp = numeric(0),
            candidate_id = integer(0),
            split_after_global_idx = integer(0)
        )
    }

    blocks_df <- if (length(blocks_list) > 0) {
        do.call(rbind, blocks_list)
    } else {
        data.frame(
            block_local_id = integer(0),
            candidate_id = integer(0),
            dmr_start_idx = integer(0),
            dmr_end_idx = integer(0),
            start_bp = numeric(0),
            end_bp = numeric(0),
            n_dmrs = integer(0),
            span_bp = numeric(0),
            max_gap_bp = numeric(0),
            segment_ids = character(0),
            stringsAsFactors = FALSE
        )
    }

    segments_df <- if (length(seg$slopes) > 0) {
        data.frame(
            segment_id = seq_along(seg$slopes),
            dmr_start_idx = seg$starts,
            dmr_end_idx = seg$ends,
            start_bp = x_chr[seg$starts],
            end_bp = x_chr[seg$ends],
            slope = seg$slopes
        )
    } else {
        data.frame(
            segment_id = integer(0),
            dmr_start_idx = integer(0),
            dmr_end_idx = integer(0),
            start_bp = numeric(0),
            end_bp = numeric(0),
            slope = numeric(0)
        )
    }

    dmr_df <- data.frame(
        chr = chr,
        dmr_order = seq_along(x_chr),
        chr_index = chr_idx,
        midpoint = x_chr,
        score_raw = y_chr,
        score_smoothed = f_chr,
        segment_id = seg$segment_id,
        segment_slope = if (length(seg$slopes) > 0) seg$slopes[seg$segment_id] else NA_real_,
        block_local_id = block_local_ids
    )

    list(
        dmr_df = dmr_df,
        segments_df = segments_df,
        candidates_df = candidates_df,
        split_events_df = split_events_df,
        blocks_df = blocks_df,
        slope_threshold = slope_threshold,
        gap_threshold_bp = gap_threshold_bp
    )
}

#' @keywords internal
#' @noRd
.assignDMRBlocksFromScores <- function(
    dmrs,
    score_col = "score",
    k_neighbors = 5L,
    min_segment_size = 2L,
    block_gap_mode = c("adaptive", "fixed", "none"),
    block_gap_fixed_bp = NULL,
    block_gap_quantile = 0.95,
    block_gap_multiplier = 1.5,
    block_gap_min_bp = 250000,
    block_gap_max_bp = 5000000,
    return_details = FALSE,
    njobs = getOption("DMRsegal.njobs", min(8, future::availableCores() - 1)),
    verbose = getOption("DMRsegal.verbose", 1L)
) {
    if (!inherits(dmrs, "GRanges")) {
        stop("dmrs must be a GRanges object.")
    }
    if (!(score_col %in% colnames(S4Vectors::mcols(dmrs)))) {
        stop("score_col '", score_col, "' not found in DMR metadata.")
    }
    block_gap_mode <- strex::match_arg(block_gap_mode, ignore_case = TRUE)

    block_id <- rep(NA_character_, length(dmrs))
    score_smooth <- rep(NA_real_, length(dmrs))
    segment_id <- rep(NA_integer_, length(dmrs))
    segment_slope <- rep(NA_real_, length(dmrs))
    details <- list()

    chr_values <- as.character(GenomicRanges::seqnames(dmrs))
    scores <- suppressWarnings(as.numeric(S4Vectors::mcols(dmrs)[[score_col]]))
    midpoints <- floor((GenomicRanges::start(dmrs) + GenomicRanges::end(dmrs)) / 2)
    chromosomes <- unique(chr_values)
    p_con <- NULL
    if (verbose > 0) {
        # check if version of progressr is equal or higher than >= 0.17.0-9002, otherwise p_con will not be used

        if (utils::packageVersion("progressr") >= "0.17.0-9002") {
            p_con <- progressr::progressor(steps = length(dmrs), message = "Assigning DMRs to blocks..")
        }
    }
    fun <- function(chr) {
        chr_idx <- which(chr_values == chr)
        if (length(chr_idx) < 3L) {
            if (!is.null(p_con)) {
                p_con(length(chr_idx))
            }
            return (NULL)
        }
        x_chr <- midpoints[chr_idx]
        y_chr <- scores[chr_idx]
        valid <- is.finite(x_chr) & is.finite(y_chr)
        if (sum(valid) < 3L) {
            if (!is.null(p_con)) {
                p_con(sum(valid))
            }
            return (NULL)
        }

        chr_idx <- chr_idx[valid]
        x_chr <- x_chr[valid]
        y_chr <- y_chr[valid]
        ord <- order(x_chr, seq_along(x_chr))
        chr_idx <- chr_idx[ord]
        x_chr <- x_chr[ord]
        y_chr <- y_chr[ord]

        chr_details <- .computeDMRBlockFormationForChromosome(
            chr = chr,
            chr_idx = chr_idx,
            x_chr = x_chr,
            y_chr = y_chr,
            k_neighbors = k_neighbors,
            min_segment_size = min_segment_size,
            block_gap_mode = block_gap_mode,
            block_gap_fixed_bp = block_gap_fixed_bp,
            block_gap_quantile = block_gap_quantile,
            block_gap_multiplier = block_gap_multiplier,
            block_gap_min_bp = block_gap_min_bp,
            block_gap_max_bp = block_gap_max_bp
        )
        if (!is.null(p_con)) {
            p_con(length(chr_idx))
        }
        list(chr = chr, chr_details = chr_details, chr_idx = chr_idx)
    }


    if (length(chromosomes) == 0L) {
        warning("No valid chromosomes found in DMRs for block assignment.")
    } else {
        if (njobs > 1L && length(chromosomes) > 1L) {
            .setupParallel()
            ret <- future.apply::future_lapply(
                chromosomes,
                fun,
                future.seed = TRUE,
                future.stdout = NA,
                future.globals = c(
                    "chr_values", "scores", "midpoints", ".computeDMRBlockFormationForChromosome",
                    "k_neighbors", "min_segment_size", "block_gap_mode", "block_gap_fixed_bp",
                    "block_gap_quantile", "block_gap_multiplier", "block_gap_min_bp", "block_gap_max_bp",
                    "njobs", "verbose"
                )
            )
        } else {
            ret <- lapply(chromosomes, fun)
        }
    }
    details <- list()
    for (res in ret) {
        if (is.null(res)) {
            next
        }
        chr <- res$chr
        chr_details <- res$chr_details
        details[[chr]] <- chr_details
        chr_idx <- res$chr_idx
        score_smooth[chr_idx] <- chr_details$dmr_df$score_smoothed
        segment_id[chr_idx] <- chr_details$dmr_df$segment_id
        segment_slope[chr_idx] <- chr_details$dmr_df$segment_slope
        local_ids <- chr_details$dmr_df$block_local_id
        has_block <- is.finite(local_ids)
        if (any(has_block)) {
            block_id[chr_idx[has_block]] <- paste0(chr, "_block", as.integer(local_ids[has_block]))
        }
    }


    S4Vectors::mcols(dmrs)$score_smoothed <- score_smooth
    S4Vectors::mcols(dmrs)$segment_id <- segment_id
    S4Vectors::mcols(dmrs)$segment_slope <- segment_slope
    S4Vectors::mcols(dmrs)$block_id <- block_id
    if (isTRUE(return_details)) {
        return(list(dmrs = dmrs, details = details))
    }
    dmrs
}

#' Rank DMRs Based on Classification Score
#'
#' @description Ranks Differentially Methylated Regions (DMRs) based on their ability to
#' discriminate between sample groups using cross-validated Support Vector Machine (SVM)
#' classification. For each DMR, this function performs stratified k-fold cross-prediction
#' using an RBF kernel SVM and computes a margin-sensitive classification score based on
#' decision values, which serves as a measure of the DMR's discriminative power.
#'
#' @param dmrs Data frame or GRanges object containing DMR coordinates and metadata
#' @param beta Character. Path to beta value file, tabix file, beta matrix, BetaHandler object, or bed file
#' @param pheno Data frame. Phenotype data containing sample group information
#' @param genome Character. Genome version (e.g., "hg19", "hg38", "mm10"). Default is "hg19"
#' @param array Character. Array platform type (e.g., "450K", "EPIC", "EPICv2"). Default is "450K"
#' @param sorted_locs Data frame. Optional pre-computed sorted genomic locations. Default is NULL
#' @param sample_group_col Character. Column name in pheno containing sample group information. Default is "Sample_Group"
#' @param block_gap_mode Character. Distance rule for block construction:
#' `"adaptive"` (default), `"fixed"`, or `"none"`.
#' @param block_gap_fixed_bp Numeric. Maximum allowed midpoint gap (bp) when
#' `block_gap_mode = "fixed"`. Ignored otherwise.
#' @param block_gap_quantile Numeric in `(0, 1)`. Quantile of chromosome DMR midpoint
#' gaps used in adaptive thresholding. Default is `0.95`.
#' @param block_gap_multiplier Numeric > 0. Multiplier applied to the adaptive
#' gap quantile. Default is `1.5`.
#' @param block_gap_min_bp Numeric >= 0. Lower clamp for adaptive gap threshold (bp).
#' Default is `250000`.
#' @param block_gap_max_bp Numeric >= `block_gap_min_bp`. Upper clamp for adaptive
#' gap threshold (bp). Default is `5000000`.
#'
#' @return GRanges object with DMRs ordered by score and additional metadata columns:
#' \itemize{
#'   \item score: Margin-sensitive cross-validated classification score for the DMR
#'   \item cv_accuracy: Raw cross-validated classification accuracy for the DMR
#'   \item score_smoothed: Gaussian-kNN smoothed score trajectory per chromosome
#'   \item segment_id: Piecewise-linear segment index estimated with PELT
#'   \item segment_slope: Estimated slope of the segment that each DMR belongs to
#'   \item block_id: Localized DMR block label (NA for DMRs not assigned to a block)
#' }
#'
#' @details
#' The function uses stratified k-fold cross-prediction to ensure balanced representation
#' of sample groups in each fold. The number of folds can be controlled using the
#' option "DMRsegal.ranking_nfold" (default is 5). An RBF (Radial Basis Function) kernel
#' SVM is trained on the beta values of CpG sites within each DMR.
#'
#' The `score` combines classification correctness and margin confidence,
#' making it more sensitive than plain cross-validated accuracy when many DMRs
#' classify perfectly. The `cv_accuracy` column stores the raw cross-validated
#' accuracy for reference. Blocks are detected from smoothed score profiles and
#' split at large midpoint gaps using the selected `block_gap_mode`.
#'
#' @examples
#' # Load example data
#' beta <- loadExampleInputData("beta")
#' pheno <- loadExampleInputData("pheno")
#'
#' # Load pre-computed DMRs
#' dmrs <- readRDS(system.file("extdata", "example_output.rds", package = "DMRsegal"))
#'
#' # Rank DMRs
#' ranked_dmrs <- rankDMRs(
#'     dmrs = dmrs,
#'     beta = beta,
#'     pheno = pheno,
#'     sample_group_col = "Sample_Group"
#' )
#'
#' @export
rankDMRs <- function(
    dmrs, beta, pheno, covariates = NULL,
    genome = "hg19", array = "450K", sorted_locs = NULL,
    sample_group_col = "Sample_Group",
    block_gap_mode = c("adaptive", "fixed", "none"),
    block_gap_fixed_bp = NULL,
    block_gap_quantile = 0.95,
    block_gap_multiplier = 1.5,
    block_gap_min_bp = 2500,
    block_gap_max_bp = 50000,
    njobs = getOption("DMRsegal.njobs", min(8, future::availableCores() - 1)),
    verbose = getOption("DMRsegal.verbose", 1L)
) {
    options("DMRsegal.verbose" = verbose)
    df_provided <- inherits(dmrs, "data.frame") && !inherits(dmrs, "GRanges")
    dmrs <- convertToGRanges(dmrs, genome = genome)
    beta_handler <- getBetaHandler(beta, array = array, genome = genome, sorted_locs = sorted_locs)
    beta_col_names <- beta_handler$getBetaColNames()
    missing_pheno_samples <- setdiff(beta_col_names, rownames(pheno))
    if (length(missing_pheno_samples) > 0) {
        stop(
            "The following beta samples are missing from pheno row names: ",
            paste(head(missing_pheno_samples, 10), collapse = ","),
            if (length(missing_pheno_samples) > 10) " ..." else ""
        )
    }
    pheno <- pheno[beta_col_names, , drop = FALSE]
    if (! "__casecontrol__" %in% colnames(pheno)) {
        pheno[, "__casecontrol__"] <- pheno[, sample_group_col] != pheno[1, sample_group_col]
    }
    class_values <- unique(pheno[, "__casecontrol__"])
    class_values <- class_values[!is.na(class_values)]
    if (length(class_values) < 2) {
        stop("Ranking requires at least two classes in '__casecontrol__'.")
    }
    groups <- pheno[, "__casecontrol__"]
    nfold <- getOption("DMRsegal.ranking_nfold", 5)
    folds <- .buildStratifiedFolds(groups, nfold = nfold)
    dmr_cpgs <- strsplit(as.character(mcols(dmrs)$cpgs), split = ",", fixed = TRUE)
    covariate_model <- .prepareCovariateModel(pheno = pheno, covariates = covariates)
    .log_step("Transforming beta values for DMR ranking", level = 2)
    dmrs_m <- .transformBeta(beta_handler$getBeta(
        row_names = unique(unlist(dmr_cpgs)),
        col_names = beta_col_names
    ), pheno = pheno, covariate_model = covariate_model)
    .log_success("Beta values transformed", level = 2)
    if (njobs > 1L) {
        .setupParallel()
        on.exit(.finalizeParallel(), add = TRUE)
    }

    .log_step("Extracting DMR-specific beta matrices for classification", level = 2)
    dmrs_m_values <- lapply(seq_along(dmrs), function(i) {
        dmr_cpgs_i <- dmr_cpgs[[i]]
        if (length(dmr_cpgs_i) == 0L) {
            return(NULL)
        }
        dmrs_m[dmr_cpgs_i, , drop = FALSE]
    })
    .log_success("DMR-specific beta matrices extracted", level = 2)
    process_args <- list(
        X = dmrs_m_values,
        FUN = function(dmrs_m) {
            cv_results <- .performCrossPrediction(dmrs_m, groups = groups, folds = folds, nfold = nfold)
            if (verbose > 0 && !is.null(p_con)) p_con()
            cv_results
        }
    )
    f <- lapply
    if (njobs > 1L) {
        process_args <- c(
            process_args,
            list(
                future.seed = TRUE,
                future.globals = c(
                    "pheno", "beta_col_names", "p_con",
                    ".performCrossPrediction",
                    "groups", "folds", "nfold", "dmr_cpgs"
                ),
                future.stdout = NA
            )
        )
        f <- future.apply::future_lapply
    }
    p_con <- NULL
    if (verbose > 0) {
        # check if version of progressr is equal or higher than >= 0.17.0-9002, otherwise p_con will not be used

        if (utils::packageVersion("progressr") >= "0.17.0-9002") {
            p_con <- progressr::progressor(steps = length(dmrs), message = "Ranking DMRs..")
        }
    }
    .log_step("Computing cross-validated classification scores for DMRs", level = 2)
    cv_metrics <- do.call(f, process_args)
    cv_metrics <- do.call(rbind, cv_metrics)

    .log_success("Cross-validated classification scores computed", level = 2)
    mcols(dmrs)$score <- as.numeric(cv_metrics[, "score"])
    mcols(dmrs)$cv_accuracy <- as.numeric(cv_metrics[, "cv_accuracy"])
    mcols(dmrs)$rank <- as.numeric(as.factor(rank(-mcols(dmrs)$score, ties.method = "first")))
    .log_step("Assigning DMRs to blocks based on smoothed score profiles", level = 2)

    dmrs <- .assignDMRBlocksFromScores(
        dmrs = dmrs,
        score_col = "score",
        block_gap_mode = block_gap_mode,
        block_gap_fixed_bp = block_gap_fixed_bp,
        block_gap_quantile = block_gap_quantile,
        block_gap_multiplier = block_gap_multiplier,
        block_gap_min_bp = block_gap_min_bp,
        block_gap_max_bp = block_gap_max_bp,
        njobs = njobs,
        verbose = verbose
    )
    .log_success("DMRs assigned to blocks", level = 2)

    dmrs <- dmrs[order(mcols(dmrs)$rank)]
    if (df_provided) {
        dmrs <- convertToDataFrame(dmrs)
    }
    dmrs
}
