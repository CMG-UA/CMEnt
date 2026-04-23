# nolint start: object_name_linter

#' Simulate DMRs in a BSseq object while preserving local structure
#'
#' This is a BSseq DMR simulator in the spirit of `dmrseq::simDMRs()`, but the
#' differential signal is added on the logit methylation scale instead of by
#' independently perturbing CpG-level beta values. Existing sample-to-sample and
#' local CpG structure in the input object is therefore preserved inside the
#' spiked regions, while the truth DMRs remain ordinary regional mean shifts.
#'
#' @param bs A BSseq object.
#' @param num_dmrs Number of DMRs to spike in.
#' @param delta_max0 Baseline maximum beta-scale effect size near the center of
#'   each DMR.
#' @param groups Optional sample group vector. If `NULL`, the first half of
#'   samples are assigned to `Condition1` and the remaining samples to
#'   `Condition2`.
#' @param case_group Group receiving the differential shift. Defaults to the
#'   second group level.
#' @param max_gap Maximum gap, in bp, used to form candidate CpG clusters.
#' @param min_cpgs Minimum number of CpGs per candidate DMR cluster.
#' @param max_cpgs Maximum number of CpGs per candidate DMR cluster.
#' @param truth_min_delta_beta Minimum intended beta-scale perturbation for a
#'   CpG to define the reported truth interval. Set to `0` to report the full
#'   selected cluster.
#' @param delta_jitter Width of the random effect-size jitter around
#'   `delta_max0`.
#' @param profile Shape of the regional effect, either `"triweight"` or
#'   `"flat"`.
#' @param profile_degree Degree used by the triweight profile.
#' @param flank_fraction Fraction of the selected cluster width added on both
#'   sides before evaluating the triweight profile. This mimics
#'   `dmrseq::simDMRs()` and avoids a hard zero at the cluster edges.
#' @param neutral_region_sd Optional logit-scale regional sample effect applied
#'   to both groups before the differential shift. The default is `0`, which
#'   avoids adding DMR-specific correlation beyond what is already present in
#'   `bs`.
#' @param matched_null_regions Number of non-DMR candidate clusters per DMR that
#'   should receive only the neutral regional sample effect when
#'   `neutral_region_sd > 0`. This can be used to avoid making local correlation
#'   unique to true DMRs.
#' @param seed Optional random seed.
#' @param rename_samples If `TRUE`, samples are reordered and renamed using the
#'   `dmrseq::simDMRs()` convention: controls first as `Condition1_Rep1`,
#'   `Condition1_Rep2`, etc., followed by cases as `Condition2_Rep1`,
#'   `Condition2_Rep2`, etc.
#' @param resample_counts If `TRUE`, methylated counts in touched regions are
#'   redrawn from the shifted probabilities and original coverage. If `FALSE`,
#'   counts are rounded deterministically.
#' @param num.dmrs,deltamax0,delta.max0 Compatibility aliases for common
#'   `dmrseq::simDMRs()` argument names.
#'
#' @return A list with `bs`, `gr.dmrs`, `dmr.mncov`, `dmr.L`, `delta`,
#'   `truth`, `selected_regions`, `groups`, and `case_group`.
#' @export
simulateDMRsBSSeq <- function(
    bs,
    num_dmrs = 3000L,
    delta_max0 = 0.3,
    groups = NULL,
    case_group = NULL,
    max_gap = 500L,
    min_cpgs = 5L,
    max_cpgs = 500L,
    truth_min_delta_beta = 0.05,
    delta_jitter = 1 / 3,
    profile = c("triweight", "flat"),
    profile_degree = 4L,
    flank_fraction = 0.2,
    neutral_region_sd = 0,
    matched_null_regions = 0L,
    seed = NULL,
    rename_samples = TRUE,
    resample_counts = TRUE
) {
    if (!inherits(bs, "BSseq")) {
        stop("'bs' must be a BSseq object.")
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }

    profile <- match.arg(profile)
    num_dmrs <- as.integer(num_dmrs)
    max_gap <- as.integer(max_gap)
    min_cpgs <- as.integer(min_cpgs)
    max_cpgs <- as.integer(max_cpgs)
    profile_degree <- as.integer(profile_degree)
    matched_null_regions <- as.integer(matched_null_regions)

    if (length(num_dmrs) != 1L || is.na(num_dmrs) || num_dmrs < 1L) {
        stop("'num_dmrs' must be a positive integer.")
    }
    if (length(delta_max0) != 1L || !is.finite(delta_max0) || delta_max0 <= 0 || delta_max0 >= 1) {
        stop("'delta_max0' must be a numeric scalar in (0, 1).")
    }
    if (length(max_gap) != 1L || is.na(max_gap) || max_gap < 0L) {
        stop("'max_gap' must be a non-negative integer.")
    }
    if (length(min_cpgs) != 1L || is.na(min_cpgs) || min_cpgs < 1L) {
        stop("'min_cpgs' must be a positive integer.")
    }
    if (length(max_cpgs) != 1L || is.na(max_cpgs) || max_cpgs < min_cpgs) {
        stop("'max_cpgs' must be an integer greater than or equal to 'min_cpgs'.")
    }
    if (length(truth_min_delta_beta) != 1L || !is.finite(truth_min_delta_beta) || truth_min_delta_beta < 0) {
        stop("'truth_min_delta_beta' must be a non-negative numeric scalar.")
    }
    if (length(delta_jitter) != 1L || !is.finite(delta_jitter) || delta_jitter < 0) {
        stop("'delta_jitter' must be a non-negative numeric scalar.")
    }
    if (length(profile_degree) != 1L || is.na(profile_degree) || profile_degree < 1L) {
        stop("'profile_degree' must be a positive integer.")
    }
    if (length(flank_fraction) != 1L || !is.finite(flank_fraction) || flank_fraction < 0) {
        stop("'flank_fraction' must be a non-negative numeric scalar.")
    }
    if (length(neutral_region_sd) != 1L || !is.finite(neutral_region_sd) || neutral_region_sd < 0) {
        stop("'neutral_region_sd' must be a non-negative numeric scalar.")
    }
    if (length(matched_null_regions) != 1L || is.na(matched_null_regions) || matched_null_regions < 0L) {
        stop("'matched_null_regions' must be a non-negative integer.")
    }

    bs <- sort(bs)
    n_samples <- ncol(bs)
    if (n_samples < 2L) {
        stop("'bs' must contain at least two samples.")
    }

    groups <- .resolveSimulationGroups(
        groups = groups,
        n_samples = n_samples,
        case_group = case_group
    )
    case_group <- attr(groups, "case_group")
    group_levels <- unique(as.character(groups))
    case_samples <- which(groups == case_group)
    control_samples <- which(groups != case_group)
    if (length(case_samples) == 0L || length(control_samples) == 0L) {
        stop("'case_group' must define at least one case and one non-case sample.")
    }

    sample_names <- colnames(bs)
    if (is.null(sample_names)) {
        sample_names <- paste0("sample_", seq_len(n_samples))
    }
    output_order <- seq_len(n_samples)
    output_groups <- as.character(groups)
    output_case_group <- case_group
    if (isTRUE(rename_samples)) {
        output_order <- c(control_samples, case_samples)
        output_groups <- c(
            rep("Condition1", length(control_samples)),
            rep("Condition2", length(case_samples))
        )
        output_case_group <- "Condition2"
        sample_names <- .makeSimulationSampleNames(output_groups)
    }

    meth_mat <- as.matrix(bsseq::getCoverage(bs, type = "M"))
    cov_mat <- as.matrix(bsseq::getCoverage(bs, type = "Cov"))
    gr <- GenomicRanges::granges(bs)
    chr <- as.character(GenomeInfoDb::seqnames(gr))
    pos <- GenomicRanges::start(gr)

    meth_mat[is.na(meth_mat)] <- 0
    cov_mat[is.na(cov_mat)] <- 0
    meth_mat <- round(meth_mat)
    cov_mat <- round(cov_mat)
    cov_mat[cov_mat < 0] <- 0
    meth_mat[meth_mat < 0] <- 0
    meth_over_cov <- meth_mat > cov_mat
    meth_mat[meth_over_cov] <- cov_mat[meth_over_cov]

    collapsed_input <- .collapseSimulationDuplicateLoci(
        meth = meth_mat,
        cov = cov_mat,
        chr = chr,
        pos = pos
    )
    meth_mat <- collapsed_input$meth
    cov_mat <- collapsed_input$cov
    chr <- collapsed_input$chr
    pos <- collapsed_input$pos

    clusters <- .makeSimulationClusters(chr = chr, pos = pos, max_gap = max_gap)
    index_by_cluster <- split(seq_along(clusters), clusters)
    cluster_lengths <- lengths(index_by_cluster)
    eligible <- which(cluster_lengths >= min_cpgs & cluster_lengths <= max_cpgs)
    if (length(eligible) < num_dmrs) {
        stop(
            "Only ", length(eligible), " eligible candidate clusters found, but ",
            num_dmrs, " DMRs were requested. Decrease 'num_dmrs' or relax ",
            "'min_cpgs'/'max_cpgs'/'max_gap'."
        )
    }

    cluster_weights <- vapply(index_by_cluster[eligible], function(idx) {
        beta_region <- .simulationBetaFromCounts(
            meth_mat[idx, , drop = FALSE],
            cov_mat[idx, , drop = FALSE]
        )
        p <- stats::median(rowMeans(beta_region, na.rm = TRUE), na.rm = TRUE)
        if (!is.finite(p)) {
            return(0)
        }
        pmax(1 - sqrt(2) * abs(0.5 - p)^0.5, 0)
    }, numeric(1))
    if (!any(cluster_weights > 0)) {
        cluster_weights <- rep(1, length(eligible))
    }

    dmr_cluster_pos <- sample(seq_along(eligible), num_dmrs, replace = FALSE, prob = cluster_weights)
    dmr_cluster_ids <- eligible[dmr_cluster_pos]
    dmr_indices <- index_by_cluster[dmr_cluster_ids]

    null_indices <- list()
    if (neutral_region_sd > 0 && matched_null_regions > 0L) {
        null_pool <- setdiff(eligible, dmr_cluster_ids)
        n_null <- min(length(null_pool), num_dmrs * matched_null_regions)
        if (n_null > 0L) {
            null_indices <- index_by_cluster[sample(null_pool, n_null, replace = FALSE)]
        }
    }

    touched_indices <- c(null_indices, dmr_indices)
    if (neutral_region_sd > 0 && length(touched_indices) > 0L) {
        for (idx in touched_indices) {
            eta <- .simulationLogitFromCounts(meth_mat[idx, , drop = FALSE], cov_mat[idx, , drop = FALSE])
            prof <- .simulationEffectProfile(pos[idx], profile, profile_degree, flank_fraction)
            eta <- eta + neutral_region_sd * outer(prof, stats::rnorm(n_samples))
            meth_mat[idx, ] <- .simulationCountsFromEta(
                eta = eta,
                cov = cov_mat[idx, , drop = FALSE],
                resample_counts = resample_counts
            )
        }
    }

    truth_rows <- vector("list", num_dmrs)
    truth_gr <- vector("list", num_dmrs)
    selected_gr <- vector("list", num_dmrs)
    deltas <- numeric(num_dmrs)
    dmr_mncov <- numeric(num_dmrs)
    dmr_lengths <- integer(num_dmrs)

    for (dmr_i in seq_along(dmr_indices)) {
        idx <- dmr_indices[[dmr_i]]
        eta0 <- .simulationLogitFromCounts(meth_mat[idx, , drop = FALSE], cov_mat[idx, , drop = FALSE])
        p0 <- plogis(eta0)
        prof <- .simulationEffectProfile(pos[idx], profile, profile_degree, flank_fraction)

        delta_max <- delta_max0 + (stats::rbeta(1L, 2, 2) - 0.5) * delta_jitter
        delta_max <- pmin(pmax(delta_max, 0.01), 0.95)

        reference_beta <- stats::median(p0[, case_samples, drop = FALSE], na.rm = TRUE)
        if (!is.finite(reference_beta)) {
            reference_beta <- stats::median(p0, na.rm = TRUE)
        }
        if (!is.finite(reference_beta)) {
            reference_beta <- 0.5
        }
        reference_beta <- .clampSimulationBeta(reference_beta)

        direction <- sample(c(-1, 1), 1L)
        if (direction > 0 && reference_beta > 1 - delta_max) {
            direction <- -1
        } else if (direction < 0 && reference_beta < delta_max) {
            direction <- 1
        }
        target_beta <- .clampSimulationBeta(reference_beta + direction * delta_max)
        effect_logit <- stats::qlogis(target_beta) - stats::qlogis(reference_beta)

        eta1 <- eta0
        eta1[, case_samples] <- eta1[, case_samples, drop = FALSE] + outer(prof, rep(effect_logit, length(case_samples)))

        intended_delta <- rowMeans(plogis(eta1[, case_samples, drop = FALSE]), na.rm = TRUE) -
            rowMeans(plogis(eta0[, case_samples, drop = FALSE]), na.rm = TRUE)
        intended_delta[!is.finite(intended_delta)] <- 0
        truth_local <- which(abs(intended_delta) >= truth_min_delta_beta)
        if (length(truth_local) < min(min_cpgs, length(idx))) {
            truth_local <- order(abs(intended_delta), decreasing = TRUE)[seq_len(min(min_cpgs, length(idx)))]
        }
        truth_local <- sort(unique(truth_local))
        truth_idx <- idx[truth_local]

        meth_mat[idx, ] <- .simulationCountsFromEta(
            eta = eta1,
            cov = cov_mat[idx, , drop = FALSE],
            resample_counts = resample_counts
        )

        case_beta <- rowMeans(.simulationBetaFromCounts(meth_mat[truth_idx, case_samples, drop = FALSE], cov_mat[truth_idx, case_samples, drop = FALSE]), na.rm = TRUE)
        control_beta <- rowMeans(.simulationBetaFromCounts(meth_mat[truth_idx, control_samples, drop = FALSE], cov_mat[truth_idx, control_samples, drop = FALSE]), na.rm = TRUE)
        observed_delta <- stats::median(case_beta - control_beta, na.rm = TRUE)
        if (!is.finite(observed_delta)) {
            observed_delta <- direction * delta_max
        }

        dmr_mncov[[dmr_i]] <- mean(rowMeans(cov_mat[truth_idx, , drop = FALSE], na.rm = TRUE), na.rm = TRUE)
        dmr_lengths[[dmr_i]] <- length(truth_idx)
        deltas[[dmr_i]] <- observed_delta

        selected_gr[[dmr_i]] <- GenomicRanges::GRanges(
            seqnames = chr[idx[[1L]]],
            ranges = IRanges::IRanges(start = min(pos[idx]), end = max(pos[idx]))
        )
        truth_gr[[dmr_i]] <- GenomicRanges::GRanges(
            seqnames = chr[truth_idx[[1L]]],
            ranges = IRanges::IRanges(start = min(pos[truth_idx]), end = max(pos[truth_idx]))
        )
        truth_rows[[dmr_i]] <- data.frame(
            seqnames = chr[truth_idx[[1L]]],
            start = min(pos[truth_idx]),
            end = max(pos[truth_idx]),
            width = max(pos[truth_idx]) - min(pos[truth_idx]) + 1L,
            num_cpgs = length(truth_idx),
            cov = dmr_mncov[[dmr_i]],
            delta_beta = observed_delta,
            delta_beta_abs = abs(observed_delta),
            intended_delta_beta = stats::median(intended_delta[truth_local], na.rm = TRUE),
            intended_delta_beta_abs = stats::median(abs(intended_delta[truth_local]), na.rm = TRUE),
            direction = ifelse(observed_delta >= 0, "hyper", "hypo"),
            case_group = output_case_group,
            stringsAsFactors = FALSE
        )
    }

    gr_dmrs <- suppressWarnings(do.call(c, truth_gr))
    selected_regions <- suppressWarnings(do.call(c, selected_gr))
    truth <- do.call(rbind, truth_rows)
    rownames(truth) <- NULL

    GenomicRanges::mcols(gr_dmrs)$cov <- dmr_mncov
    GenomicRanges::mcols(gr_dmrs)$num_cpgs <- dmr_lengths
    GenomicRanges::mcols(gr_dmrs)$delta_beta <- deltas
    GenomicRanges::mcols(gr_dmrs)$delta_beta_abs <- abs(deltas)
    GenomicRanges::mcols(gr_dmrs)$case_group <- output_case_group

    meth_mat <- meth_mat[, output_order, drop = FALSE]
    cov_mat <- cov_mat[, output_order, drop = FALSE]
    colnames(meth_mat) <- colnames(cov_mat) <- sample_names
    bs_new <- bsseq::BSseq(
        chr = chr,
        pos = pos,
        M = meth_mat,
        Cov = cov_mat,
        sampleNames = sample_names
    )
    GenomeInfoDb::seqinfo(bs_new) <- GenomeInfoDb::seqinfo(bs)
    SummarizedExperiment::colData(bs_new) <- S4Vectors::DataFrame(
        Sample_ID = sample_names,
        Sample_Group = output_groups,
        casecontrol = output_groups == output_case_group,
        row.names = sample_names
    )

    list(
        gr.dmrs = gr_dmrs,
        dmr.mncov = dmr_mncov,
        dmr.L = dmr_lengths,
        bs = bs_new,
        delta = deltas,
        truth = truth,
        selected_regions = selected_regions,
        groups = stats::setNames(output_groups, sample_names),
        case_group = output_case_group,
        input_groups = stats::setNames(as.character(groups)[output_order], sample_names),
        input_case_group = case_group,
        duplicate_loci_collapsed = collapsed_input$n_collapsed,
        null_regions = if (length(null_indices) == 0L) GenomicRanges::GRanges() else {
            suppressWarnings(do.call(c, lapply(null_indices, function(idx) {
                GenomicRanges::GRanges(
                    seqnames = chr[idx[[1L]]],
                    ranges = IRanges::IRanges(start = min(pos[idx]), end = max(pos[idx]))
                )
            })))
        }
    )
}

.collapseSimulationDuplicateLoci <- function(meth, cov, chr, pos) {
    key <- paste(chr, pos, sep = "\t")
    duplicate_mask <- duplicated(key)
    if (!any(duplicate_mask)) {
        return(list(
            meth = meth,
            cov = cov,
            chr = chr,
            pos = pos,
            n_collapsed = 0L
        ))
    }

    first_idx <- match(unique(key), key)
    meth_collapsed <- rowsum(meth, group = key, reorder = FALSE)
    cov_collapsed <- rowsum(cov, group = key, reorder = FALSE)
    colnames(meth_collapsed) <- colnames(meth)
    colnames(cov_collapsed) <- colnames(cov)

    list(
        meth = as.matrix(meth_collapsed),
        cov = as.matrix(cov_collapsed),
        chr = chr[first_idx],
        pos = pos[first_idx],
        n_collapsed = as.integer(sum(duplicate_mask))
    )
}

.resolveSimulationGroups <- function(groups, n_samples, case_group = NULL) {
    if (is.null(groups)) {
        n_control <- floor(n_samples / 2)
        groups <- c(
            rep("Condition1", n_control),
            rep("Condition2", n_samples - n_control)
        )
    }
    groups <- as.character(groups)
    if (length(groups) != n_samples) {
        stop("'groups' must have length equal to ncol(bs).")
    }
    if (anyNA(groups) || any(!nzchar(groups))) {
        stop("'groups' must not contain missing or empty values.")
    }
    group_levels <- unique(groups)
    if (length(group_levels) < 2L) {
        stop("'groups' must contain at least two groups.")
    }
    if (is.null(case_group)) {
        case_group <- group_levels[[2L]]
    }
    case_group <- as.character(case_group)[[1L]]
    if (!case_group %in% group_levels) {
        stop("'case_group' must be one of: ", paste(group_levels, collapse = ", "))
    }
    attr(groups, "case_group") <- case_group
    groups
}

.makeSimulationSampleNames <- function(groups) {
    reps <- ave(seq_along(groups), groups, FUN = seq_along)
    paste0(groups, "_Rep", as.integer(reps))
}

.makeSimulationClusters <- function(chr, pos, max_gap) {
    if (length(chr) != length(pos)) {
        stop("'chr' and 'pos' must have the same length.")
    }
    if (length(chr) == 0L) {
        return(integer(0))
    }
    starts <- c(TRUE, chr[-1L] != chr[-length(chr)] | diff(pos) > max_gap)
    cumsum(starts)
}

.simulationEffectProfile <- function(pos, profile, degree, flank_fraction) {
    if (length(pos) == 1L || identical(profile, "flat")) {
        return(rep(1, length(pos)))
    }
    first <- min(pos)
    last <- max(pos)
    width <- max(last - first, 1)
    first <- first - flank_fraction * width
    last <- last + flank_fraction * width
    width <- max(last - first, 1)
    center <- first + width / 2
    scaled <- abs(pos - center) / (width / 2)
    out <- (1 - pmin(scaled, 1)^degree)^degree
    out[!is.finite(out)] <- 0
    pmax(out, 0)
}

.simulationBetaFromCounts <- function(meth, cov) {
    beta <- (meth + 0.5) / (cov + 1)
    beta[cov <= 0] <- NA_real_
    beta
}

.simulationLogitFromCounts <- function(meth, cov) {
    stats::qlogis(.clampSimulationBeta(.simulationBetaFromCounts(meth, cov)))
}

.clampSimulationBeta <- function(x) {
    x_dim <- dim(x)
    x_dimnames <- dimnames(x)
    out <- pmin(pmax(x, 1e-6), 1 - 1e-6)
    if (!is.null(x_dim)) {
        dim(out) <- x_dim
        dimnames(out) <- x_dimnames
    }
    out
}

.simulationCountsFromEta <- function(eta, cov, resample_counts) {
    prob <- plogis(eta)
    prob[!is.finite(prob)] <- 0.5
    cov <- round(cov)
    cov[cov < 0] <- 0
    if (isTRUE(resample_counts)) {
        out <- matrix(
            stats::rbinom(length(cov), size = c(cov), prob = c(prob)),
            nrow = nrow(cov),
            ncol = ncol(cov)
        )
    } else {
        out <- round(cov * prob)
    }
    out[out < 0] <- 0
    out_over_cov <- out > cov
    out[out_over_cov] <- cov[out_over_cov]
    out
}

# nolint end: object_name_linter
