# nolint start: object_name_linter

#' Simulate DMRs while preserving local methylation structure
#'
#' This simulator is inspired by `dmrseq::simDMRs()`, but adds differential
#' signal on the logit methylation scale. The same regional perturbation engine
#' is used across methylation assays by operating on methylated-count and
#' coverage-count representations. For microarray-style beta values, pseudo
#' counts are constructed first, then the same simulation mechanism is applied.
#'
#' @param beta A BSseq object, a [BetaHandler] object, a beta matrix/data
#'   frame, a path to a beta/tabix file, or a path to an `.rds` file
#'   containing a BSseq object.
#' @param num_dmrs Number of DMRs to spike in.
#' @param delta_max0 Baseline maximum beta-scale effect size near the center of
#'   each DMR.
#' @param groups Optional sample group vector. If `NULL`, the first half of
#'   samples are assigned to `Condition1` and the remaining samples to
#'   `Condition2`.
#' @param case_group Group receiving the differential shift. Defaults to the
#'   second group level.
#' @param max_gap Maximum gap, in bp, used to form candidate site clusters.
#' @param min_sites Minimum number of sites per candidate DMR cluster.
#' @param max_sites Maximum number of sites per candidate DMR cluster.
#' @param truth_min_delta_beta Minimum intended beta-scale perturbation for a
#'   site to define the reported truth interval. Set to `0` to report the full
#'   selected cluster.
#' @param delta_jitter Width of the random effect-size jitter around
#'   `delta_max0`.
#' @param profile Shape of the regional effect, either `"triweight"` or
#'   `"flat"`.
#' @param profile_degree Degree used by the triweight profile.
#' @param flank_fraction Fraction of the selected cluster width added on both
#'   sides before evaluating the triweight profile.
#' @param seed Optional random seed.
#' @param rename_samples If `TRUE`, samples are reordered and renamed using the
#'   `dmrseq::simDMRs()` convention: controls first as `Condition1_Rep1`,
#'   `Condition1_Rep2`, etc., followed by cases as `Condition2_Rep1`,
#'   `Condition2_Rep2`, etc.
#' @param resample_counts If `TRUE`, methylated counts in touched regions are
#'   redrawn from shifted probabilities and coverage. If `FALSE`, counts are
#'   rounded deterministically.
#' @param array Array platform type. Used only when `beta` is not a BSseq object
#'   or a self-contained [BetaHandler].
#' @param genome Reference genome. Used only when `beta` is not a BSseq object
#'   or a self-contained [BetaHandler].
#' @param sorted_locs Optional genomic locations with `chr`, `start`, and
#'   optionally `end` columns. Used only for non-BSseq inputs.
#' @param beta_row_names_file Optional file with beta row names. Used only for
#'   non-BSseq file inputs.
#' @param chrom_col Chromosome column name for tabix inputs.
#' @param start_col Start column name for tabix inputs.
#' @param njobs Number of parallel jobs for loading non-BSseq inputs.
#'
#' @return A list with simulated output (`simulated`), optional genomic
#'   locations for non-BSseq inputs (`beta_locs`), and dmrseq-like metadata:
#'   `gr.dmrs`, `dmr.mncov`, `dmr.L`, `delta`, `truth`, `selected_regions`,
#'   `groups`, and `case_group`.
#' @importFrom stats qlogis plogis
#' @export
simulateDMRs <- function(
    beta,
    num_dmrs = 3000L,
    delta_max0 = 0.3,
    groups = NULL,
    case_group = NULL,
    max_gap = 500L,
    min_sites = 5L,
    max_sites = 500L,
    truth_min_delta_beta = 0.05,
    delta_jitter = 1 / 3,
    profile = c("triweight", "flat"),
    profile_degree = 4L,
    flank_fraction = 0.2,
    seed = NULL,
    rename_samples = TRUE,
    resample_counts = TRUE,
    array = c("450K", "27K", "EPIC", "EPICv2"),
    genome = c("hg38", "hg19", "hs1", "mm10", "mm39"),
    sorted_locs = NULL,
    beta_row_names_file = NULL,
    chrom_col = "#chrom",
    start_col = "start",
    njobs = getOption("CMEnt.njobs", min(8, future::availableCores() - 1))
) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    profile <- match.arg(profile)
    num_dmrs <- as.integer(num_dmrs)
    max_gap <- as.integer(max_gap)
    min_sites <- as.integer(min_sites)
    max_sites <- as.integer(max_sites)
    profile_degree <- as.integer(profile_degree)

    if (length(num_dmrs) != 1L || is.na(num_dmrs) || num_dmrs < 1L) {
        stop("'num_dmrs' must be a positive integer.")
    }
    if (length(delta_max0) != 1L || !is.finite(delta_max0) || delta_max0 <= 0 || delta_max0 >= 1) {
        stop("'delta_max0' must be a numeric scalar in (0, 1).")
    }
    if (length(max_gap) != 1L || is.na(max_gap) || max_gap < 0L) {
        stop("'max_gap' must be a non-negative integer.")
    }
    if (length(min_sites) != 1L || is.na(min_sites) || min_sites < 1L) {
        stop("'min_sites' must be a positive integer.")
    }
    if (length(max_sites) != 1L || is.na(max_sites) || max_sites < min_sites) {
        stop("'max_sites' must be an integer greater than or equal to 'min_sites'.")
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

    input <- .prepareSimulationInput(
        beta = beta,
        array = array,
        genome = genome,
        sorted_locs = sorted_locs,
        beta_row_names_file = beta_row_names_file,
        chrom_col = chrom_col,
        start_col = start_col,
        njobs = njobs
    )

    meth_mat <- input$meth
    cov_mat <- input$cov
    chr <- input$chr
    pos <- input$pos
    end_pos <- input$end
    row_ids <- input$row_ids

    n_samples <- ncol(meth_mat)
    if (n_samples < 2L) {
        stop("Input methylation data must contain at least two samples.")
    }

    groups <- .resolveSimulationGroups(
        groups = groups,
        n_samples = n_samples,
        case_group = case_group
    )
    case_group <- attr(groups, "case_group")
    case_samples <- which(groups == case_group)
    control_samples <- which(groups != case_group)
    if (length(case_samples) == 0L || length(control_samples) == 0L) {
        stop("'case_group' must define at least one case and one non-case sample.")
    }

    sample_names <- input$sample_names
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
        pos = pos,
        end = end_pos,
        row_ids = row_ids
    )
    meth_mat <- collapsed_input$meth
    cov_mat <- collapsed_input$cov
    chr <- collapsed_input$chr
    pos <- collapsed_input$pos
    end_pos <- collapsed_input$end
    row_ids <- collapsed_input$row_ids

    clusters <- .makeSimulationClusters(chr = chr, pos = pos, max_gap = max_gap)
    index_by_cluster <- split(seq_along(clusters), clusters)
    cluster_lengths <- lengths(index_by_cluster)
    eligible <- which(cluster_lengths >= min_sites & cluster_lengths <= max_sites)
    if (length(eligible) < num_dmrs) {
        stop(
            "Only ", length(eligible), " eligible candidate clusters found, but ",
            num_dmrs, " DMRs were requested. Decrease 'num_dmrs' or relax ",
            "'min_sites'/'max_sites'/'max_gap'."
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
        if (length(truth_local) < min(min_sites, length(idx))) {
            truth_local <- order(abs(intended_delta), decreasing = TRUE)[seq_len(min(min_sites, length(idx)))]
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
            ranges = IRanges::IRanges(start = min(pos[idx]), end = max(end_pos[idx]))
        )
        truth_gr[[dmr_i]] <- GenomicRanges::GRanges(
            seqnames = chr[truth_idx[[1L]]],
            ranges = IRanges::IRanges(start = min(pos[truth_idx]), end = max(end_pos[truth_idx]))
        )
        truth_rows[[dmr_i]] <- data.frame(
            seqnames = chr[truth_idx[[1L]]],
            start = min(pos[truth_idx]),
            end = max(end_pos[truth_idx]),
            width = max(end_pos[truth_idx]) - min(pos[truth_idx]) + 1L,
            num_sites = length(truth_idx),
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
    GenomicRanges::mcols(gr_dmrs)$num_sites <- dmr_lengths
    GenomicRanges::mcols(gr_dmrs)$delta_beta <- deltas
    GenomicRanges::mcols(gr_dmrs)$delta_beta_abs <- abs(deltas)
    GenomicRanges::mcols(gr_dmrs)$case_group <- output_case_group

    meth_mat <- meth_mat[, output_order, drop = FALSE]
    cov_mat <- cov_mat[, output_order, drop = FALSE]
    colnames(meth_mat) <- colnames(cov_mat) <- sample_names

    output <- .buildSimulationOutputObject(
        input = input,
        meth = meth_mat,
        cov = cov_mat,
        chr = chr,
        pos = pos,
        end_pos = end_pos,
        row_ids = row_ids,
        sample_names = sample_names,
        output_groups = output_groups,
        output_case_group = output_case_group
    )

    result <- list(
        gr.dmrs = gr_dmrs,
        dmr.mncov = dmr_mncov,
        dmr.L = dmr_lengths,
        simulated = output$simulated,
        assay = input$assay,
        delta = deltas,
        truth = truth,
        selected_regions = selected_regions,
        groups = stats::setNames(output_groups, sample_names),
        case_group = output_case_group,
        input_groups = stats::setNames(as.character(groups)[output_order], sample_names),
        input_case_group = case_group,
        duplicate_loci_collapsed = collapsed_input$n_collapsed
    )

    if (!is.null(output$beta_locs)) {
        result$beta_locs <- output$beta_locs
    }

    result
}

.prepareSimulationInput <- function(beta,
                                    array,
                                    genome,
                                    sorted_locs,
                                    beta_row_names_file,
                                    chrom_col,
                                    start_col,
                                    njobs) {
    if (
        is.character(beta) &&
            length(beta) == 1L &&
            file.exists(beta) &&
            tolower(tools::file_ext(beta)) == "rds"
    ) {
        beta_rds <- readRDS(beta)
        if (inherits(beta_rds, "BSseq")) {
            beta <- beta_rds
        }
    }

    if (inherits(beta, "BSseq")) {
        bs <- sort(beta)
        meth_mat <- as.matrix(bsseq::getCoverage(bs, type = "M"))
        cov_mat <- as.matrix(bsseq::getCoverage(bs, type = "Cov"))
        gr <- GenomicRanges::granges(bs)
        chr <- as.character(GenomeInfoDb::seqnames(gr))
        pos <- GenomicRanges::start(gr)
        end_pos <- GenomicRanges::end(gr)
        row_ids <- paste0(chr, ":", pos)
        sample_names <- colnames(bs)
        if (is.null(sample_names)) {
            sample_names <- paste0("sample_", seq_len(ncol(meth_mat)))
        }
        return(list(
            assay = "BSseq",
            source = bs,
            seqinfo = GenomeInfoDb::seqinfo(bs),
            sample_names = sample_names,
            row_ids = row_ids,
            chr = chr,
            pos = pos,
            end = end_pos,
            meth = meth_mat,
            cov = cov_mat
        ))
    }

    beta_handler <- getBetaHandler(
        beta = beta,
        array = array,
        genome = genome,
        beta_row_names_file = beta_row_names_file,
        sorted_locs = sorted_locs,
        chrom_col = chrom_col,
        start_col = start_col,
        njobs = njobs
    )

    beta_mat <- as.matrix(beta_handler$getBeta())
    beta_locs <- as.data.frame(beta_handler$getBetaLocs())
    if (!all(c("chr", "start") %in% colnames(beta_locs))) {
        stop("Non-BSseq simulation requires genomic locations with 'chr' and 'start' columns.")
    }
    if (!("end" %in% colnames(beta_locs))) {
        beta_locs$end <- beta_locs$start
    }

    if (is.null(rownames(beta_locs)) || any(!nzchar(rownames(beta_locs)))) {
        beta_locs_ids <- paste0(as.character(beta_locs$chr), ":", as.integer(beta_locs$start))
        rownames(beta_locs) <- make.unique(beta_locs_ids)
    }

    if (is.null(rownames(beta_mat)) || any(!nzchar(rownames(beta_mat)))) {
        if (nrow(beta_mat) != nrow(beta_locs)) {
            stop("Unable to align non-BSseq beta values with genomic locations.")
        }
        rownames(beta_mat) <- rownames(beta_locs)
    }

    common_ids <- rownames(beta_locs)[rownames(beta_locs) %in% rownames(beta_mat)]
    if (length(common_ids) == 0L) {
        stop("No overlapping row identifiers between beta values and genomic locations.")
    }
    beta_locs <- beta_locs[common_ids, , drop = FALSE]
    beta_mat <- beta_mat[common_ids, , drop = FALSE]

    suppressWarnings(storage.mode(beta_mat) <- "double")
    invalid_mask <- !is.finite(beta_mat)
    beta_mat[invalid_mask] <- 0
    beta_mat <- .clampSimulationBeta(beta_mat)

    cov_mat <- matrix(
        100L,
        nrow = nrow(beta_mat),
        ncol = ncol(beta_mat),
        dimnames = dimnames(beta_mat)
    )
    cov_mat[invalid_mask] <- 0
    meth_mat <- round(beta_mat * cov_mat)
    meth_mat[invalid_mask] <- 0

    sample_names <- colnames(beta_mat)
    if (is.null(sample_names)) {
        sample_names <- paste0("sample_", seq_len(ncol(beta_mat)))
    }

    list(
        assay = "microarray",
        source = beta_handler,
        sample_names = sample_names,
        row_ids = rownames(beta_mat),
        chr = as.character(beta_locs$chr),
        pos = as.integer(beta_locs$start),
        end = as.integer(beta_locs$end),
        meth = meth_mat,
        cov = cov_mat,
        beta_locs = beta_locs
    )
}

.buildSimulationOutputObject <- function(input,
                                         meth,
                                         cov,
                                         chr,
                                         pos,
                                         end_pos,
                                         row_ids,
                                         sample_names,
                                         output_groups,
                                         output_case_group) {
    if (identical(input$assay, "BSseq")) {
        bs_new <- bsseq::BSseq(
            chr = chr,
            pos = pos,
            M = meth,
            Cov = cov,
            sampleNames = sample_names
        )
        if (!is.null(input$seqinfo)) {
            GenomeInfoDb::seqinfo(bs_new) <- input$seqinfo
        }
        SummarizedExperiment::colData(bs_new) <- S4Vectors::DataFrame(
            Sample_ID = sample_names,
            Sample_Group = output_groups,
            casecontrol = output_groups == output_case_group,
            row.names = sample_names
        )
        return(list(
            simulated = bs_new,
            beta_locs = NULL
        ))
    }

    beta_new <- .simulationBetaFromCounts(meth = meth, cov = cov)
    rownames(beta_new) <- row_ids
    colnames(beta_new) <- sample_names
    beta_locs <- data.frame(
        chr = chr,
        start = as.integer(pos),
        end = as.integer(end_pos),
        row.names = row_ids,
        stringsAsFactors = FALSE
    )

    list(
        simulated = beta_new,
        beta_locs = beta_locs
    )
}

.collapseSimulationDuplicateLoci <- function(meth, cov, chr, pos, end = NULL, row_ids = NULL) {
    key <- paste(chr, pos, sep = "\t")
    duplicate_mask <- duplicated(key)
    if (is.null(end)) {
        end <- pos
    }
    if (is.null(row_ids)) {
        row_ids <- paste0(chr, ":", pos)
    }

    if (!any(duplicate_mask)) {
        return(list(
            meth = meth,
            cov = cov,
            chr = chr,
            pos = pos,
            end = end,
            row_ids = row_ids,
            n_collapsed = 0L
        ))
    }

    first_idx <- match(unique(key), key)
    meth_collapsed <- rowsum(meth, group = key, reorder = FALSE)
    cov_collapsed <- rowsum(cov, group = key, reorder = FALSE)
    colnames(meth_collapsed) <- colnames(meth)
    colnames(cov_collapsed) <- colnames(cov)

    row_ids_collapsed <- as.character(row_ids[first_idx])
    if (anyDuplicated(row_ids_collapsed) > 0L) {
        row_ids_collapsed <- make.unique(row_ids_collapsed)
    }

    list(
        meth = as.matrix(meth_collapsed),
        cov = as.matrix(cov_collapsed),
        chr = chr[first_idx],
        pos = pos[first_idx],
        end = end[first_idx],
        row_ids = row_ids_collapsed,
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
        stop("'groups' must have length equal to the number of samples.")
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

#' @importFrom stats ave
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
