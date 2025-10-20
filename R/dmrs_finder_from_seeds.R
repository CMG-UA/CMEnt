#' Find Differentially Methylated Regions (DMRs) from Pre-computed DMPs
#'
#' @name findDMRsFromSeeds
#' @description This function identifies Differentially Methylated Regions (DMRs) from pre-computed
#' Differentially Methylated Positions (DMPs) using a correlation-based approach. It expands
#' significant DMPs into regions, considering both statistical significance and biological
#' relevance of methylation changes.
#'
#' @section Important Note on Input Data:
#' Do not apply heavy filtering to your DMPs prior to using this function, particularly based on
#' beta values or effect sizes. The function works by expanding regions around significant DMPs
#' and connecting nearby CpGs into larger regions. Filtering out DMPs with smaller effect sizes
#' may remove important CpGs that could serve as "bridges" to connect more significant DMPs into
#' larger, biologically meaningful DMRs. For optimal results, include all statistically
#' significant DMPs (e.g., adjusted p-value < 0.05) and let the function handle region expansion
#' and filtering internally using the min_cpg_delta_beta parameter if needed.
#'
#' @param beta_file Path to the methylation beta values file or a data matrix with beta values
#' @param dmps_file Path to the pre-computed DMPs file or a data frame with DMPs
#' @param pheno Data frame. Phenotype data.
#' @param dmps_tsv_id_col Character. Column name for DMP identifiers in the DMPs TSV file. Default is NULL.
#' @param dmp_groups_info Named list. Required when `dmps_tsv_id_col` is given. List of DMP group information, where names are group identifiers, found in dmps_tsv_id_col column, and values are the samples names, found in the beta values columns. Default is NULL.
#' @param pval_col Column name in DMPs file containing p-values (default: "pval_adj")
#' @param sample_group_col Column in pheno for sample grouping (default: "Sample_Group")
#' @param dmp_group_col Column in DMPs file for grouping DMPs (default: NULL)
#' @param casecontrol_col Column in pheno for case/control status (default: "casecontrol")
#' @param min_cpg_delta_beta Minimum delta beta threshold for CpGs (default: 0)
#' @param expansion_step Distance in bp to expand regions during search (default: 500)
#' @param expansion_relaxation Relaxation parameter for region expansion (default: 0)
#' @param array Array platform, either "450K" or "EPIC" (default: c("450K", "EPIC"))
#' @param genome Reference genome, "hg19" or "hg38" (default: c("hg19", "hg38", "mm10"))
#' @param max_pval Maximum p-value threshold for DMPs (default: 0.05)
#' @param max_lookup_dist Maximum distance for region expansion in bp (default: 10000)
#' @param min_dmps Minimum number of DMPs required per region (default: 1)
#' @param min_adj_dmps Minimum number of adjacent DMPs required per region (default: 1)
#' @param min_cpgs Minimum number of CpGs required per region, existing in the array (default: 50)
#' @param ignored_sample_groups Sample groups to ignore during analysis (default: NULL)
#' @param output_prefix Optional identifier prefix for output files (default: NULL)
#' @param aggfun Aggregation function for computing group means, either "median" or "mean" (default: c("median", "mean"))
#' @param njobs Number of parallel jobs (default: detectCores())
#' @param memory_threshold_mb Memory threshold in MB for switching between in-memory and file-based processing (default: 500)
#' @param verbose Numeric. Level of verbosity for logging messages, from 0 (not verbose) to 3 (very verbose). Default is 1.
#' @param beta_row_names_file Optional file with beta value row names (default: NULL)
#' @param tabix_file Path to tabix-indexed beta values file (alternative to beta_file, default: NULL)
#'
#' @return A GRanges object containing identified DMRs with metadata columns:
#' \itemize{
#'   \item cpgs_num: Number of CpGs in the region
#'   \item dmps_num: Number of DMPs in the region
#'   \item delta_beta: Mean methylation difference (aggregated using the specified aggregation function)
#'   \item delta_beta_min: Minimum methylation difference
#'   \item delta_beta_max: Maximum methylation difference
#'   \item delta_beta_sd: Standard deviation of methylation differences
#'   \item pval_adj: Adjusted p-value (if available in DMPs file)
#'   \item cases_beta: Mean beta value in cases
#'   \item controls_beta: Mean beta value in controls
#'   \item start_dmp: ID of the first DMP in the region
#'   \item end_dmp: ID of the last DMP in the region
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(example_beta)
#' data(example_dmps)
#' data(example_pheno)
#'
#' # Write beta values to file
#' beta_file <- tempfile(fileext = ".txt")
#' write.table(cbind(ID = rownames(example_beta), example_beta),
#'     file = beta_file, sep = "\t", quote = FALSE, row.names = FALSE
#' )
#'
#' # Find DMRs
#' dmrs <- findDMRsFromSeeds(
#'     beta_file = beta_file,
#'     dmps_file = example_dmps,
#'     pheno = example_pheno
#' )
#' }
#'
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom future.apply future_lapply
#' @importFrom future availableCores
#' @importFrom progressr progressor
#' @importFrom stringr str_count str_order
#' @importFrom readr read_tsv
#' @importFrom data.table fread fwrite
#' @importFrom psych corr.test
#' @importFrom dplyr %>% filter select mutate
#' @importFrom bedr tabix
#' @importFrom BSgenome getSeq
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rtracklayer import.chain
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom utils write.table read.table
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom GenomicFeatures genes promoters
#' @importFrom AnnotationDbi mapIds
#' @importFrom S4Vectors queryHits subjectHits

#' @keywords internal
#' @noRd
.expandDMRs <- function(dmr,
                        beta_file_handler,
                        beta_row_names,
                        beta_col_names,
                        sorted_locs,
                        max_pval,
                        min_cpg_delta_beta = 0,
                        casecontrol = NULL,
                        min_cpgs = 3,
                        expansion_step = 500,
                        expansion_relaxation = 0,
                        max_lookup_dist = 1e7,
                        group_inds = NULL,
                        extreme_verbosity = FALSE,
                        aggfun = mean,
                        pval_mode = c("parametric", "empirical"),
                        empirical_strategy = c("auto", "montecarlo"),
                        nperm = 200L,
                        mid_p = FALSE,
                        perm_seed = NULL) {
    .log_step("Expanding DMR..", level=4)
    pval_mode <- match.arg(pval_mode)
    empirical_strategy <- match.arg(empirical_strategy)
    sorted_locs <- sorted_locs[beta_row_names, ]

    dmr_start <- dmr["start_dmp"]
    dmr_end <- dmr["end_dmp"]
    dmr_start_ind <- which(beta_row_names == dmr_start)
    dmr_end_ind <- which(beta_row_names == dmr_end)
    if (length(dmr_start_ind) == 0) {
        stop("Could not find the start CpG ", dmr_start, " in the beta file row names.")
    }
    if (length(dmr_end_ind) == 0) {
        stop("Could not find the end CpG ", dmr_end, " in the beta file row names.")
    }


    .check_upstream <- function(ustream_end_lookup_site_ind, exp_step) {
        .log_step("Setting up..", level = 6)
        ustream_stop_reason <- NULL
        if (ustream_end_lookup_site_ind < 0) {
            ustream_stop_reason <- "end-of-input"
            ustream_exp <- 1
            return(list(
                ustream_end_lookup_site_ind = ustream_end_lookup_site_ind,
                ustream_stop_reason = ustream_stop_reason, ustream_exp = ustream_exp
            ))
        }
        ustream_start_lookup_site_ind <- max(0, ustream_end_lookup_site_ind - exp_step)
        x <- which(sorted_locs[ustream_start_lookup_site_ind:ustream_end_lookup_site_ind, "chr"] == dmr["chr"])
        if (length(x) == 0) {
            ustream_stop_reason <- "end-of-input"
            ustream_exp <- 1
            return(list(
                ustream_end_lookup_site_ind = ustream_end_lookup_site_ind,
                ustream_stop_reason = ustream_stop_reason, ustream_exp = ustream_exp
            ))
        }
        x <- x[[1]]
        ustream_start_lookup_site_ind <- ustream_start_lookup_site_ind + x - 1
        .log_success("Setup done.", level=6)
        .log_step("Getting beta values...", level=6)
        ustream_betas <- beta_file_handler$getBeta(row_names = rownames(sorted_locs)[ustream_start_lookup_site_ind:ustream_end_lookup_site_ind], col_names = beta_col_names)
        .log_success("Beta values retrieved.", level=6)
        if (nrow(ustream_betas) == 1) {
            ustream_stop_reason <- "end-of-input"
            return(list(
                ustream_end_lookup_site_ind = ustream_end_lookup_site_ind,
                ustream_stop_reason = ustream_stop_reason, ustream_exp = ustream_exp
            ))
        }
        .log_step("Preparing beta values...", level=6)
        ustream_betas_reversed <- ustream_betas[rev(seq_len(nrow(ustream_betas))), , drop = FALSE]
        .log_success("Beta values prepared.", level=6)
        .log_step("Testing connectivity...", level=6)
        corr_ret <- .testConnectivityBatch(
            sites_beta = ustream_betas_reversed,
            group_inds = group_inds,
            casecontrol = casecontrol,
            max_pval = max_pval,
            min_delta_beta = min_cpg_delta_beta,
            max_lookup_dist = max_lookup_dist,
            sites_locs = sorted_locs[ustream_start_lookup_site_ind:ustream_end_lookup_site_ind, ],
            aggfun = aggfun,
            pval_mode = pval_mode,
            empirical_strategy = empirical_strategy,
            nperm = nperm,
            mid_p = mid_p,
            perm_seed = if (is.null(perm_seed)) NULL else as.integer(perm_seed)
        )
        .log_success("Connectivity tested.", level=6)
        .log_step("Analyzing connectivity results...", level=6)
        # find consecutive FALSE counts in connected
        consecutive_fails <- rle(!corr_ret$connected)
        # pick the first run of FALSE with length > expansion_relaxation
        fail_runs <- which(consecutive_fails$values & consecutive_fails$lengths > expansion_relaxation)
        if (length(fail_runs) > 0) {
            first_fail_run <- fail_runs[[1]]
            # Calculate the index in the original vector
            fail_start_idx <- sum(consecutive_fails$lengths[seq_len(first_fail_run - 1)]) + 1
            ustream_exp <- ustream_end_lookup_site_ind - fail_start_idx + 1
            ustream_stop_reason <- corr_ret$reason[fail_start_idx]
            .log_success("Connectivity results analyzed.", level=6)
            .log_success("DMR expanded.", level=4)
            return(list(
                ustream_end_lookup_site_ind = ustream_end_lookup_site_ind,
                ustream_stop_reason = ustream_stop_reason, ustream_exp = ustream_exp
            ))
        }
        .log_success("DMR expanded.", level = 4)
        .log_success("Connectivity results analyzed.", level=6)
        ustream_end_lookup_site_ind <- ustream_end_lookup_site_ind - exp_step
        ret <- list(
            ustream_end_lookup_site_ind = ustream_end_lookup_site_ind,
            ustream_stop_reason = ustream_stop_reason, ustream_exp = ustream_exp
        )
        ret
    }

    .check_downstream <- function(dstream_start_lookup_site_ind, exp_step) {
        dstream_stop_reason <- NULL
        if (dstream_start_lookup_site_ind > length(beta_row_names)) {
            dstream_stop_reason <- "end-of-input"
            dstream_exp <- dstream_start_lookup_site_ind - 1
            return(list(
                dstream_start_lookup_site_ind = dstream_start_lookup_site_ind,
                dstream_stop_reason = dstream_stop_reason, dstream_exp = dstream_exp
            ))
        }
        dstream_end_lookup_site_ind <- min(dstream_start_lookup_site_ind + exp_step, length(beta_row_names))
        x <- which(rev(sorted_locs[dstream_start_lookup_site_ind:dstream_end_lookup_site_ind, "chr"]) == dmr["chr"])
        if (length(x) == 0) {
            dstream_stop_reason <- "end-of-input"
            dstream_exp <- dstream_start_lookup_site_ind - 1
            return(list(
                dstream_start_lookup_site_ind = dstream_start_lookup_site_ind,
                dstream_stop_reason = dstream_stop_reason, dstream_exp = dstream_exp
            ))
        }
        x <- x[[1]]
        dstream_end_lookup_site_ind <- dstream_end_lookup_site_ind - x + 1
        if (dstream_end_lookup_site_ind <= dstream_start_lookup_site_ind + 1) {
            dstream_stop_reason <- "end-of-input"
            dstream_exp <- dstream_start_lookup_site_ind - 1
            return(list(
                dstream_start_lookup_site_ind = dstream_start_lookup_site_ind,
                dstream_stop_reason = dstream_stop_reason, dstream_exp = dstream_exp
            ))
        }
        dstream_betas <- beta_file_handler$getBeta(row_names = rownames(sorted_locs)[dstream_start_lookup_site_ind:dstream_end_lookup_site_ind], col_names = beta_col_names)
        corr_ret <- .testConnectivityBatch(
            sites_beta = dstream_betas,
            group_inds = group_inds,
            casecontrol = casecontrol,
            max_pval = max_pval,
            min_delta_beta = min_cpg_delta_beta,
            max_lookup_dist = max_lookup_dist,
            sites_locs = sorted_locs[dstream_start_lookup_site_ind:dstream_end_lookup_site_ind, ],
            aggfun = aggfun,
            pval_mode = pval_mode,
            empirical_strategy = empirical_strategy,
            nperm = nperm,
            mid_p = mid_p,
            perm_seed = if (is.null(perm_seed)) NULL else as.integer(perm_seed)
        )
        # find consecutive FALSE counts in connected
        consecutive_fails <- rle(!corr_ret$connected)
        # pick the first run of FALSE with length > expansion_relaxation
        fail_runs <- which(consecutive_fails$values & consecutive_fails$lengths > expansion_relaxation)
        if (length(fail_runs) > 0) {
            first_fail_run <- fail_runs[[1]]
            # Calculate the index in the original vector
            fail_start_idx <- sum(consecutive_fails$lengths[seq_len(first_fail_run - 1)]) + 1
            dstream_exp <- dstream_start_lookup_site_ind + fail_start_idx - 1
            dstream_stop_reason <- corr_ret$reason[fail_start_idx]
            return(list(
                dstream_start_lookup_site_ind = dstream_start_lookup_site_ind,
                dstream_stop_reason = dstream_stop_reason, dstream_exp = dstream_exp
            ))
        }
        prev_start_lookup_site_ind <- dstream_start_lookup_site_ind
        dstream_start_lookup_site_ind <- dstream_start_lookup_site_ind + exp_step - x + 1
        if (dstream_start_lookup_site_ind < prev_start_lookup_site_ind) {
            stop("BUG: dstream_start_lookup_site_ind < prev_start_lookup_site_ind during downstream expansion. dstream_start_lookup_site_ind: ", dstream_start_lookup_site_ind, " prev_start_lookup_site_ind: ", prev_start_lookup_site_ind, " expansion_step: ", exp_step, " x: ", x)
        }
        return(list(
            dstream_start_lookup_site_ind = dstream_start_lookup_site_ind,
            dstream_stop_reason = dstream_stop_reason, dstream_exp = dstream_exp
        ))
    }

    ustream_end_lookup_site_ind <- dmr_start_ind[[1]]
    ustream_exp <- ustream_end_lookup_site_ind
    ustream_stop_reason <- NULL
    dstream_start_lookup_site_ind <- dmr_end_ind[[1]]
    dstream_exp <- dstream_start_lookup_site_ind
    dstream_stop_reason <- NULL
    
    t <- 0
    while (TRUE) {
        exp_step <- expansion_step
        if (t == 0) { # first iteration, use min_cpgs + expansion_relaxation and remove the DMRs that are not long enough
            ccpgs <- dstream_exp - ustream_exp + 1
            if (ccpgs < (min_cpgs)) {
                .log_info("DMR ", dmr["dmr_id"], " too short (", ccpgs, " CpGs). Expanding to reach min_cpgs=", min_cpgs, ".", level = 3)
                exp_step <- min_cpgs - ccpgs + expansion_relaxation
            }
            .log_info("Number of CpGs in DMR: ", ccpgs, level = 5)

        }
        .log_info("Expansion step size: ", exp_step, " bp.", level = 5)
        .log_step("Checking upstream expansion...", level = 5)
        if (is.null(ustream_stop_reason)) {
            res <- .check_upstream(ustream_end_lookup_site_ind, exp_step)
            ustream_end_lookup_site_ind <- res$ustream_end_lookup_site_ind
            ustream_stop_reason <- res$ustream_stop_reason
            ustream_exp <- res$ustream_exp
        }
        .log_success("Upstream expansion checked.", level = 5)
        .log_step("Checking downstream expansion...", level = 5)
        if (is.null(dstream_stop_reason)) {
            res <- .check_downstream(dstream_start_lookup_site_ind, exp_step)
            dstream_start_lookup_site_ind <- res$dstream_start_lookup_site_ind
            dstream_stop_reason <- res$dstream_stop_reason
            dstream_exp <- res$dstream_exp
        }

        if (t == 0) {
            new_ccpgs <- dstream_exp - ustream_exp + 1
            .log_info("Number of CpGs in expanded DMR: ", new_ccpgs, " from ", ccpgs,  level = 5)
            if (new_ccpgs < min_cpgs) {
                ustream_stop_reason <- "min_cpgs_reached"
                dstream_stop_reason <- "min_cpgs_reached"
            }

            t <- 1
        }
        .log_success("Downstream expansion checked.", level = 5)
        if (!is.null(ustream_stop_reason) && !is.null(dstream_stop_reason)) {
            break
        }
    }
    .log_step("Finalizing expanded DMR ", dmr["dmr_id"], ".", level = 5)
    dmr["start_cpg"] <- beta_row_names[ustream_exp]
    dmr["end_cpg"] <- beta_row_names[dstream_exp]
    dmr["start"] <- sorted_locs[dmr["start_cpg"], "pos"]
    dmr["end"] <- sorted_locs[dmr["end_cpg"], "pos"]
    dmr["downstream_cpg_expansion_stop_reason"] <- dstream_stop_reason
    dmr["upstream_cpg_expansion_stop_reason"] <- ustream_stop_reason
    .log_success("Expanded DMR finalized: ", dmr["dmr_id"], " (start_cpg: ", dmr["start_cpg"], ", end_cpg: ", dmr["end_cpg"], ").", level = 5)
    dmr
}

#' Vectorized connectivity testing for consecutive site pairs, given their beta values
#'
#' @param sites_beta Matrix where each row is a site's beta values across samples
#' @param group_inds List of sample indices for each group
#' @param max_pval Maximum p-value threshold
#' @param casecontrol Optional case/control indicator vector
#' @param min_delta_beta Minimum delta beta threshold
#' @param max_lookup_dist Maximum distance between consecutive sites
#' @param sites_locs Data frame with chr and pos columns for each site
#' @param aggfun Aggregation function for computing group means (default: mean)
#' @param pval_mode Character. "parametric" (default) for t-based p-values, or "empirical" for permutation-based p-values.
#' @param empirical_strategy Character. When pval_mode = "empirical": "auto" (default) uses Monte Carlo for rows with >=4 valid samples and falls back to parametric for rows with <=3; "montecarlo" forces Monte Carlo for all eligible rows.
#' @param nperm Integer. Number of permutations when using empirical Monte Carlo. Default is 200.
#' @param mid_p Logical. Use mid-p adjustment for empirical p-values (reduces tie conservatism). Default is FALSE.
#' @param perm_seed Integer or NULL. RNG seed for reproducible permutations when pval_mode = "empirical". Default is NULL.
#'
#' @return Data frame with columns: connected, pval, delta_beta, reason, first_failing_group, stop_reason
#' @keywords internal
#' @noRd
.testConnectivityBatch <- function(sites_beta, group_inds, max_pval,
                                   casecontrol = NULL, min_delta_beta = 0,
                                   max_lookup_dist = NULL, sites_locs = NULL, aggfun = mean,
                                   pval_mode = c("parametric", "empirical"),
                                   empirical_strategy = c("auto", "montecarlo"),
                                   nperm = 200L, mid_p = FALSE, perm_seed = NULL) {
    pval_mode <- match.arg(pval_mode)
    empirical_strategy <- match.arg(empirical_strategy)
    n_sites <- nrow(sites_beta)
    if (n_sites < 2) {
        return(data.frame(
            connected = logical(0),
            pval = numeric(0),
            delta_beta = numeric(0),
            reason = character(0),
            first_failing_group = character(0)
        ))
    }

    n_pairs <- n_sites - 1
    n_groups <- length(group_inds)
    max_pval_corrected <- max_pval / n_groups

    # Initialize result vectors
    connected <- rep(TRUE, n_pairs)
    pvals <- rep(0, n_pairs)
    delta_betas <- rep(NA_real_, n_pairs)
    reasons <- rep("", n_pairs)
    failing_groups <- rep("", n_pairs)

    # Check distance condition if provided (vectorized)
    if (!is.null(max_lookup_dist) && !is.null(sites_locs)) {
        dists <- sites_locs$pos[2:n_sites] - sites_locs$pos[1:(n_sites - 1)]
        exceeded_dist <- dists > max_lookup_dist
        connected[exceeded_dist] <- FALSE
        reasons[exceeded_dist] <- "exceeded max distance"
    }

    # Fully vectorized correlation testing for each group
    for (g in names(group_inds)) {
        idx <- group_inds[[g]]
        if (length(idx) < 3) next

        # Get data for this group - subset columns
        group_beta <- sites_beta[, idx, drop = FALSE]
        group_m <- log(group_beta / (1 - group_beta + 1e-6) + 1e-6) # M-values transformation

        # Extract consecutive pairs matrices
        x_mat <- group_m[1:(n_sites - 1), , drop = FALSE] # Sites i
        y_mat <- group_m[2:n_sites, , drop = FALSE] # Sites i+1

        # Compute means for each pair (vectorized)
        x_means <- rowMeans(x_mat, na.rm = TRUE)
        y_means <- rowMeans(y_mat, na.rm = TRUE)

        # Center the data
        x_centered <- x_mat - x_means
        y_centered <- y_mat - y_means

        # Compute sum of products (numerator of correlation)
        sum_xy <- rowSums(x_centered * y_centered, na.rm = TRUE)
        sum_x2 <- rowSums(x_centered^2, na.rm = TRUE)
        sum_y2 <- rowSums(y_centered^2, na.rm = TRUE)

        # Compute correlations (fully vectorized)
        cors <- sum_xy / sqrt(sum_x2 * sum_y2)

        # Compute degrees of freedom (vectorized)
        # Count non-NA pairs for each row
        valid_pairs <- !is.na(x_mat) & !is.na(y_mat)
        n_valid <- rowSums(valid_pairs)
        dfs <- n_valid - 2L

        # Identify failures (vectorized logical operations)
        na_r <- is.na(cors) & connected
        connected[na_r] <- FALSE
        reasons[na_r] <- "na r"
        failing_groups[na_r] <- g

        low_df <- (dfs < 1) & connected
        connected[low_df] <- FALSE
        reasons[low_df] <- "df<1"
        failing_groups[low_df] <- g

        ps <- rep(NA_real_, n_pairs)
        # Precompute parametric p-values as fallback when empirical is not feasible/resolved
        tstats <- cors * sqrt(dfs / pmax(1e-12, 1 - cors * cors))
        ps_parametric <- rep(NA_real_, n_pairs)
        finite_tstat_mask <- !is.na(tstats)
        valid_df_mask <- dfs >= 1L
        finite_and_valid_df <- finite_tstat_mask & valid_df_mask
        if (any(finite_and_valid_df)) {
            ps_parametric[finite_and_valid_df] <- -2 * expm1(pt(abs(tstats[finite_and_valid_df]), df = dfs[finite_and_valid_df], log.p = TRUE))
        }

        if (pval_mode == "parametric") {
            # Compute t-statistics (vectorized)
            na_tstat <- is.na(tstats) & connected
            connected[na_tstat] <- FALSE
            reasons[na_tstat] <- "na tstat"
            failing_groups[na_tstat] <- g

            # Compute p-values (vectorized), only for the non-na and connected tstats
            finite_tstat_and_connected <- finite_tstat_mask & connected
            ps[finite_tstat_and_connected] <- ps_parametric[finite_tstat_and_connected]
        } else {
            # Empirical p-values via permutations of sample labels within group
            # Only compute for rows that are still connected and have finite cors
            mask <- is.finite(cors) & connected
            if (any(mask)) {
                # RNG management
                old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                if (old_seed_exists) {
                    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
                }
                if (!is.null(perm_seed)) {
                    set.seed(perm_seed)
                }
                on.exit({
                    if (old_seed_exists) {
                        assign(".Random.seed", old_seed, envir = .GlobalEnv)
                    }
                }, add = TRUE)
                counts_ge <- integer(n_pairs)
                counts_eq <- integer(n_pairs)
                # Number of samples in this group
                m <- ncol(y_mat)
                # Empirical strategy: auto -> fallback to parametric when n_valid <= 3 (impossible to get p < 0.05)
                do_empirical <- mask
                if (empirical_strategy == "auto") {
                    do_empirical <- do_empirical & (n_valid >= 4L)
                }
                if (nperm > 0L && m > 2L && any(do_empirical)) {
                    for (b in seq_len(nperm)) {
                        # Permute sample labels (columns) only for y; x remains fixed
                        perm <- sample.int(m, size = m, replace = FALSE)
                        yp <- y_mat[, perm, drop = FALSE]
                        ym <- rowMeans(yp, na.rm = TRUE)
                        yc <- sweep(yp, 1L, ym, FUN = "-")
                        sxy <- rowSums(x_centered * yc, na.rm = TRUE)
                        sy2 <- rowSums(yc^2, na.rm = TRUE)
                        rperm <- sxy / sqrt(sum_x2 * sy2)
                        comp_mask <- do_empirical & is.finite(rperm)
                        if (any(comp_mask)) {
                            ap <- abs(rperm[comp_mask])
                            ao <- abs(cors[comp_mask])
                            counts_ge[comp_mask] <- counts_ge[comp_mask] + (ap > ao)
                            counts_eq[comp_mask] <- counts_eq[comp_mask] + (ap == ao)
                        }
                    }
                }
                if (mid_p) {
                    ps[do_empirical] <- (counts_ge[do_empirical] + 0.5 * counts_eq[do_empirical] + 1) / (nperm + 1)
                } else {
                    ps[do_empirical] <- (counts_ge[do_empirical] + counts_eq[do_empirical] + 1) / (nperm + 1)
                }
                # Fallback to parametric for rows where empirical was not applied (auto) or unavailable
                fallback_mask <- mask & !do_empirical
                if (any(fallback_mask)) {
                    ps[fallback_mask] <- ps_parametric[fallback_mask]
                }
            }
        }

        na_p <- is.na(ps) & connected
        connected[na_p] <- FALSE
        reasons[na_p] <- "na pval"
        failing_groups[na_p] <- g

        # Update maximum p-values across groups (vectorized)
        pvals <- pmax(pvals, ps, na.rm = TRUE)

        # Check p-value threshold (vectorized)
        exceed_pval <- (pvals > max_pval_corrected) & connected
        connected[exceed_pval] <- FALSE
        reasons[exceed_pval] <- "pval>max_pval (corrected)"
        failing_groups[exceed_pval] <- g
    }
    ret <- data.frame(
        connected = connected,
        pval = pvals,
        reason = reasons,
        first_failing_group = failing_groups,
        stringsAsFactors = FALSE
    )
    # Vectorized delta beta check if needed
    if (!is.null(casecontrol) && min_delta_beta > 0) {
        # Extract site2 beta values for all pairs
        site2_beta_mat <- sites_beta[2:n_sites, , drop = FALSE]

        # Vectorized mean computation across case/control
        case_betas <- apply(site2_beta_mat[, casecontrol == 1, drop = FALSE], 1, aggfun, na.rm = TRUE)
        control_betas <- apply(site2_beta_mat[, casecontrol == 0, drop = FALSE], 1, aggfun, na.rm = TRUE)
        delta_betas <- case_betas - control_betas

        # Check threshold (vectorized)
        low_delta <- (is.na(delta_betas) | abs(delta_betas) < min_delta_beta) & connected
        connected[low_delta] <- FALSE
        reasons[low_delta] <- "delta_beta<min_delta_beta"
        ret[["delta_beta"]] <- delta_betas
    }

    ret
}

#' Extract CpG Information from DMR Results
#'
#' @description This function extracts detailed CpG methylation information for
#' identified DMRs, providing beta values and additional metadata for each CpG
#' site within the DMR regions.
#'
#' @param dmrs GRanges object containing DMR results from findDMRsFromSeeds
#' @param beta_file_handler BetaFileHandler object for the methylation beta values file
#' @param beta_row_names Character vector. Row names for beta values
#' @param beta_col_names Character vector. Column names for beta values
#' @param sorted_locs GRanges or data frame with sorted genomic locations
#' @param output_prefix Character. Optional prefix for output files (default: NULL)
#' @param njobs Integer. Number of cores for parallel processing (default: 1)
#'
#' @return A data frame containing detailed CpG information for each DMR
#'
#' @examples
#' \dontrun{
#' # Extract CpG info from DMR results
#' cpg_info <- extractCpgInfoFromResultDMRs(
#'     dmrs = dmr_results,
#'     beta_file_handler = beta_file_handler,
#'     beta_row_names = row_names,
#'     beta_col_names = col_names,
#'     sorted_locs = genomic_locations
#' )
#' }
#'
#' @export
extractCpgInfoFromResultDMRs <- function(dmrs,
                                         beta_file_handler,
                                         beta_row_names,
                                         beta_col_names,
                                         sorted_locs,
                                         output_prefix = NULL,
                                         njobs = 1) {
    if (!is.null(output_prefix)) {
        if (!dir.exists(dirname(output_prefix))) {
            dir.create(dirname(output_prefix), recursive = TRUE)
        }
    }
    dmrs_cpgs <- data.frame()
    dmrs_sites <- c()
    .log_step("Finding constituent CpGs of DMRs present in beta file..")
    for (i in seq_len(nrow(dmrs))) {
        start_ind <- dmrs$start_ind[i]
        end_ind <- dmrs$end_ind[i]
        cpgs <- rownames(sorted_locs)[start_ind:end_ind]
        cpgs <- cpgs[cpgs %in% beta_row_names]
        dmrs_cpgs <- rbind(dmrs_cpgs, data.frame(dmr = dmrs$id[i], cpgs = cpgs))
        dmrs_sites <- c(dmrs_sites, cpgs)
    }

    if (!is.null(output_prefix)) {
        dmrs_cpgs_file <- paste0(output_prefix, ".cpgs.tsv.gz")
        .log_step("Saving DMR-CpG mapping to ", dmrs_cpgs_file, " ...")
        gz <- gzfile(dmrs_cpgs_file, "w")
        write.table(
            dmrs_cpgs,
            gz,
            sep = "\t",
            quote = FALSE,
            qmethod = "double",
            col.names = TRUE,
            row.names = FALSE
        )
        close(gz)
    }

    dmrs_sites <- unique(dmrs_sites)
    dmrs_sites <- rownames(sorted_locs)[rownames(sorted_locs) %in% dmrs_sites]

    dmrs_beta_file <- paste0(output_prefix, "dmrs_beta.tsv.gz")
    .log_step("Saving constituent CpG betas", dmrs_beta_file, " ...")
    dmrs_beta <- beta_file_handler$getBeta(
        row_names = dmrs_sites,
        col_names = beta_col_names
    )
    rownames(dmrs_beta) <- dmrs_sites
    if (!is.null(output_prefix)) {
        gz <- gzfile(dmrs_beta_file, "w")
        write.table(
            dmrs_beta,
            gz,
            sep = "\t",
            quote = FALSE,
            qmethod = "double",
            col.names = NA,
            row.names = TRUE
        )
        close(gz)
    }
    ret <- list(
        dmrs_cpgs = dmrs_cpgs,
        dmrs_beta = dmrs_beta
    )
    invisible(ret)
}



#' Find Differentially Methylated Regions (DMRs) from Differentially Methylated Positions (DMPs)
#'
#' This function identifies DMRs from a given set of DMPs and a beta value file.
#'
#' @param beta_file Character. Path to the beta value file. Either this or tabix_file must be provided.
#' @param dmps_file Character. Path to the DMPs TSV file.
#' @param pheno Data frame. Phenotype data.
#' @param dmps_tsv_id_col Character. Column name for DMP identifiers in the DMPs TSV file. Default is NULL.
#' @param pval_col Character. Column name for p-values in the DMPs file. Default is "pval_adj".
#' @param sample_group_col Character. Column name for sample group information in the phenotype data. Default is NULL.
#' @param dmp_group_col Character. Column name for DMP group information in the DMPs TSV file. Default is NULL.
#' @param dmp_groups_info Named list. Required when `dmp_group_col` is given. List of DMP group information, where names are group identifiers, found in dmps_tsv_id_col column, and values are the samples names, found in the beta values columns. Default is NULL.
#' @param casecontrol_col Character. Column name for case-control information in the phenotype data. Default is "casecontrol".
#' @param min_cpg_delta_beta Numeric. Minimum delta beta value for CpGs. Default is 0.
#' @param expansion_step Numeric. Step size for expanding DMRs. Increasing it means higher memory usage and faster computation. Default is 500.
#' @param expansion_relaxation Numeric. Maximum number of intermittent CpGs allowed to not be significantly correlated, to increase the extended DMR size. Default is 0.
#' @param array Character. Type of array used (e.g., "450K", "EPIC", "EPICv2", "27K"). Ignored if using mm10 genome.
#' @param genome Character. Genome version (e.g., "hg19", "hg38", "mm10"). Default is "hg19".
#' @param max_pval Numeric. Maximum p-value to assume DMPs correlation is significant. Default is 0.05.
#' @param pval_mode Character. "parametric" (default) to use t-based correlation p-values during connectivity testing, or "empirical" to use permutation-based p-values.
#' @param nperm Integer. Number of permutations when pval_mode = "empirical". Default is 0 (disabled).
#' @param perm_seed Integer or NULL. RNG seed for reproducibility when pval_mode = "empirical". Default is NULL.
#' @param max_lookup_dist Numeric. Maximum distance to look up for adjacent DMPs belonging to the same DMR. Default is 10000.
#' @param min_dmps Numeric. Minimum number of connected DMPs in a DMR. Default is 1.
#' @param min_adj_dmps Numeric. Minimum number of DMPs, adjusted by CpG density, in a DMR after extension. Default is 1.
#' @param min_cpgs Numeric. Minimum number of CpGs in a DMR after extension. Default is 50.
#' @param aggfun Character. Aggregation function to use ("median" or "mean") when calculating delta beta values and p-values of DMRs. Default is "median".
#' @param ignored_sample_groups Character vector. Sample groups to ignore, separated by commas. Default is NULL.
#' @param output_prefix Character. Identifier for the output files. If not provided, no output will be saved. Default is NULL.
#' @param njobs Numeric. Number of parallel jobs to use. Default is the number of available cores.
#' @param verbose Numeric. Level of verbosity for logging messages, from 0 (not verbose) to 3 (very verbose). Default is 1.
#' @param memory_threshold_mb Numeric. Memory threshold in MB for loading beta files. Default is 500.
#' @param beta_row_names_file Character. Path to a file containing row names for the beta values. If not provided, row names will be read from the beta file. Default is NULL.
#' @param tabix_file Character. Path to a tabix-indexed beta file. Either this or beta_file must be provided. Default is NULL.
#'
#' @return Data frame of identified DMRs.
#' @export
findDMRsFromSeeds <- function(beta_file = NULL,
                              dmps_file = NULL,
                              pheno = NULL,
                              dmps_tsv_id_col = NULL,
                              pval_col = "pval_adj",
                              sample_group_col = "Sample_Group",
                              dmp_group_col = NULL,
                              dmp_groups_info = NULL,
                              casecontrol_col = "casecontrol",
                              min_cpg_delta_beta = 0,
                              expansion_step = 500,
                              expansion_relaxation = 0,
                              array = c("450K", "27K", "EPIC", "EPICv2"),
                              genome = c("hg19", "hg38", "mm10", "mm39"),
                              max_pval = 0.05,
                              pval_mode = c("parametric", "empirical"),
                              empirical_strategy = c("auto", "montecarlo"),
                              nperm = 200L,
                              mid_p = FALSE,
                              perm_seed = NULL,
                              max_lookup_dist = 10000,
                              min_dmps = 1,
                              min_adj_dmps = 1,
                              min_cpgs = 50,
                              aggfun = c("median", "mean"),
                              ignored_sample_groups = NULL,
                              output_prefix = NULL,
                              njobs = future::availableCores(),
                              memory_threshold_mb = 500,
                              verbose = 1,
                              beta_row_names_file = NULL,
                              tabix_file = NULL) {
    pval_mode <- match.arg(pval_mode)
    empirical_strategy <- match.arg(empirical_strategy)
    # Clean up any zombie processes on exit
    includes <- "#include <sys/wait.h>"
    code <- "int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};"
    wait <- inline::cfunction(body = code, includes = includes, convention = ".C")
    on.exit(wait(), add = TRUE)

    # Bridge verbose to logging option for consistent styled logs
    old_opt <- options(DMRSegal.verbose = verbose)
    on.exit(options(old_opt), add = TRUE)
    if (is.null(dmps_file) || is.null(pheno)) {
        stop("dmps_file and pheno parameters are required")
    }
    beta_file_handler <- BetaFileHandler$new(
        beta_file = beta_file,
        tabix_file = tabix_file,
        beta_row_names_file = beta_row_names_file,
        njobs = njobs,
        memory_threshold_mb = memory_threshold_mb,
        verbose = verbose
    )

    aggfun <- match.arg(aggfun)
    aggfun <- switch(aggfun,
        median = stats::median,
        mean = mean
    )
    stopifnot(!is.null(max_pval))
    stopifnot(!is.null(min_dmps))
    stopifnot(!is.null(min_cpgs))
    stopifnot(!is.null(expansion_step))
    stopifnot(!is.null(min_cpg_delta_beta))
    stopifnot(!is.null(max_lookup_dist))

    array <- match.arg(array)
    genome <- match.arg(genome)

    stopifnot(file.exists(dmps_file))
    stopifnot(casecontrol_col %in% colnames(pheno))
    stopifnot(sample_group_col %in% colnames(pheno))
    if (is.null(ignored_sample_groups)) {
        ignored_sample_groups <- c()
    } else {
        ignored_sample_groups <- unlist(strsplit(ignored_sample_groups, ","))
    }
    if (!is.null(output_prefix)) {
        output_dir <- dirname(output_prefix)
        dir.create(output_dir, showWarnings = FALSE)
        output_prefix <- paste0(output_prefix, ".")
    } else {
        output_dir <- NULL
        output_prefix <- NULL
    }

    .log_step("Preparing inputs...")
    .log_step("Reading DMP tsv..", level = 2)
    dmps_tsv <- try(read.table(
        dmps_file,
        header = TRUE,
        sep = "\t",
        check.names = FALSE,
        quote = "",
        comment.char = "",
        row.names = NULL,
    ))
    if (inherits(dmps_tsv, "try-error")) {
        .log_warn("Provided DMPs file is empty or does not exist. Not proceeding.")
        if (!is.null(output_prefix)) {
            for (f in c(".methylation.tsv.gz", ".tsv.gz")) {
                gzfile <- gzfile(paste0(output_prefix, f), "w", compression = 2)
                close(gzfile)
            }
        }
        return(NULL)
    }
    
    # Check if dmps_tsv has any rows
    if (nrow(dmps_tsv) == 0) {
        .log_warn("Provided DMPs file has no data rows. Not proceeding.")
        if (!is.null(output_prefix)) {
            for (f in c(".methylation.tsv.gz", ".tsv.gz")) {
                gzfile <- gzfile(paste0(output_prefix, f), "w", compression = 2)
                close(gzfile)
            }
        }
        return(GenomicRanges::GRanges())
    }
    
    if (!is.null(dmp_group_col)) {
        if (!dmp_group_col %in% colnames(dmps_tsv)) {
            stop(
                "DMP group column '", dmp_group_col,
                "' does not reside in the DMPs file columns: ",
                paste(colnames(dmps_tsv), collapse = ",")
            )
        }
        if (is.null(dmp_groups_info)) {
            stop(
                "dmp_groups_info parameter is required when dmp_group_col is provided"
            )
        }
        unique_dmp_groups <- unique(dmps_tsv[, dmp_group_col])
        if (!all(unique_dmp_groups %in% names(dmp_groups_info))) {
            missing <- unique_dmp_groups[!(unique_dmp_groups %in% names(dmp_groups_info))]
            stop("The following DMP groups do not exist in the names of the dmp_groups_info list'", dmp_group_col, "' column: ", paste(missing, collapse = ","))
        }
    } else {
        dmp_group_col <- "_DUMMY_DMP_GROUP_COL_"
        dmps_tsv[, dmp_group_col] <- "all"
    }

    if (!pval_col %in% colnames(dmps_tsv)) {
        stop(
            "P-value column '", pval_col, "' does not reside in the DMPs file columns: ",
            paste(colnames(dmps_tsv), collapse = ",")
        )
    }
    stopifnot(pval_col %in% colnames(dmps_tsv))

    .log_step("Reading beta file characteristics..", level = 2)

    beta_row_names <- beta_file_handler$getBetaRowNames()
    beta_col_names <- beta_file_handler$getBetaColNames()


    samples_selection_mask <- !(pheno[, sample_group_col] %in% ignored_sample_groups)
    beta_col_names <- beta_col_names[samples_selection_mask]
    pheno <- pheno[beta_col_names, ]
    .log_info("Samples to process: ", length(beta_col_names), level = 1)

    pheno <- pheno[beta_col_names, ]
    sample_groups <- factor(pheno[beta_col_names, sample_group_col])
    group_inds <- split(seq_along(sample_groups), sample_groups)
    case_mask <- pheno[beta_col_names, casecontrol_col] == 1

    .log_step("Reordering DMPs by genomic location...", level = 2)

    if (is.null(dmps_tsv_id_col)) {
        dmps_tsv_id_col <- "row.names"
    }
    if (!dmps_tsv_id_col %in% colnames(dmps_tsv)) {
        stop(
            "DMP id column '", dmps_tsv_id_col,
            "' does not reside in the DMPs file columns: ",
            paste(colnames(dmps_tsv), collapse = ",")
        )
    }
    sorted_locs <- getSortedGenomicLocs(array = array, genome = genome)
    dmps_tsv <- dmps_tsv[orderByLoc(dmps_tsv[, dmps_tsv_id_col], genome = genome, genomic_locs = sorted_locs), , drop = FALSE]

    dmps <- unique(dmps_tsv[, dmps_tsv_id_col])
    # Filter DMPs not present in array annotation first (prevents NA logical indices later)
    missing_in_annotation <- setdiff(dmps, rownames(sorted_locs))
    if (length(missing_in_annotation) > 0) {
        .log_warn(
            "Dropping ", length(missing_in_annotation), " DMP(s) not found in the array annotation: ",
            paste(head(missing_in_annotation, 10), collapse = ","),
            if (length(missing_in_annotation) > 10) " ..." else ""
        )
        dmps_tsv <- dmps_tsv[!(dmps_tsv[, dmps_tsv_id_col] %in% missing_in_annotation), , drop = FALSE]
        dmps <- setdiff(dmps, missing_in_annotation)
    }
    if (length(dmps) == 0) {
        stop("No DMPs remain after filtering against array annotation.")
    }
    if (!all(dmps %in% beta_row_names)) {
        missing_in_beta <- dmps[!(dmps %in% beta_row_names)]
        stop("Some of the DMPs are not present in the beta file. DMPs: ", paste(missing_in_beta, collapse = ","))
    }
    dmps <- dmps[orderByLoc(dmps, genome = genome, genomic_locs = sorted_locs)]
    .log_step("Validating beta file sorting by position...", level = 2)


    .log_step("Subsetting beta matrix for DMPs...", level = 2)
    dmps_locs <- sorted_locs[dmps, , drop = FALSE]
    dmps_tsv[, "chr"] <- dmps_locs[dmps_tsv[, dmps_tsv_id_col], "chr"]
    dmps_tsv[, "pos"] <- dmps_locs[dmps_tsv[, dmps_tsv_id_col], "pos"]
    if (!is.null(output_prefix)) {
        dmps_beta_output_file <- paste0(output_prefix, "dmps_beta.tsv.gz")
    }
    dmps_beta <- beta_file_handler$getBeta(row_names = dmps, col_names = beta_col_names)
    .log_step("Calculating delta_beta related columns in the DMPs table...", level = 2)
    if (dmp_group_col == "_DUMMY_DMP_GROUP_COL_") {
        dmp_groups_info <- list(all = NULL)
    }
    for (dmp_group in unique(dmps_tsv[, dmp_group_col])) {
        group_mask <- dmps_tsv[, dmp_group_col] == dmp_group
        if (!is.null(dmp_groups_info[[dmp_group]])) {
            beta_mask <- rownames(dmps_beta) %in% dmp_groups_info[[dmp_group]]
        } else {
            beta_mask <- rep(TRUE, nrow(dmps_beta))
        }
        group_beta_subset <- dmps_beta[beta_mask, , drop = FALSE]
        cases_beta <- apply(group_beta_subset[, pheno[beta_col_names, casecontrol_col] == 1, drop = FALSE], 1, aggfun, na.rm = TRUE)
        controls_beta <- apply(group_beta_subset[, pheno[beta_col_names, casecontrol_col] == 0, drop = FALSE], 1, aggfun, na.rm = TRUE)
        case_sd <- apply(group_beta_subset[, pheno[beta_col_names, casecontrol_col] == 1, drop = FALSE], 1, sd, na.rm = TRUE)
        control_sd <- apply(group_beta_subset[, pheno[beta_col_names, casecontrol_col] == 0, drop = FALSE], 1, sd, na.rm = TRUE)
        cases_num <- rowSums(!is.na(group_beta_subset[, pheno[beta_col_names, casecontrol_col] == 1, drop = FALSE]))
        controls_num <- rowSums(!is.na(group_beta_subset[, pheno[beta_col_names, casecontrol_col] == 0, drop = FALSE]))
        dmps_tsv[group_mask, "cases_num"] <- cases_num[dmps_tsv[group_mask, dmps_tsv_id_col]]
        dmps_tsv[group_mask, "controls_num"] <- controls_num[dmps_tsv[group_mask, dmps_tsv_id_col]]
        dmps_tsv[group_mask, "cases_beta"] <- cases_beta[dmps_tsv[group_mask, dmps_tsv_id_col]]
        dmps_tsv[group_mask, "controls_beta"] <- controls_beta[dmps_tsv[group_mask, dmps_tsv_id_col]]
        dmps_tsv[group_mask, "cases_beta_sd"] <- case_sd[dmps_tsv[group_mask, dmps_tsv_id_col]]
        dmps_tsv[group_mask, "controls_beta_sd"] <- control_sd[dmps_tsv[group_mask, dmps_tsv_id_col]]
        dmps_tsv[group_mask, "delta_beta"] <- dmps_tsv[group_mask, "cases_beta"] - dmps_tsv[group_mask, "controls_beta"]
    }
    .log_success("Delta beta columns calculated", level = 2)
    if (!is.null(output_prefix)) {
        .log_step("Saving dmps beta to file: ", dmps_beta_output_file, " ...", level = 2)
        gz <- gzfile(dmps_beta_output_file, "w")
        write.table(dmps_beta, gz, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
        close(gz)
    }

    if (nrow(dmps_locs) != nrow(dmps_beta)) {
        stop(
            "Number of rows in the queried dmps beta file does not match the number of DMPs. Number of rows in beta file: ",
            nrow(dmps_beta),
            " Number of rows in DMPs: ",
            nrow(dmps_locs)
        )
    }
    # Diagnostic: identify any rows with all NA betas (should not happen with synthetic test data)
    all.na.rows <- apply(dmps_beta, 1, function(r) all(is.na(r)))
    if (any(all.na.rows)) {
        stop(
            "Beta extraction failure: the following DMP rows have all NA beta values: ",
            paste(rownames(dmps_beta)[all.na.rows], collapse = ","),
            ". This indicates a mismatch between requested CpG IDs and beta file columns or a parsing issue."
        )
    }

    .log_info("Subset size: ", paste(dim(dmps_beta), collapse = ","), level = 2)
    .log_info("Number of provided DMPs: ", length(dmps), level = 2)

    .log_success("Input preparation complete.", level = 1)
    .log_step("Connecting DMPs to form initial DMRs..", level = 1)

    # Set up progress tracking for DMP connection
    chromosomes <- unique(dmps_locs$chr)

    # Set up future plan for parallel processing
    if (njobs < 0) {
        njobs <- future::availableCores() + njobs
    }
    if (njobs > 1) {
        if (future::availableCores("multicore") > 1L) {
            .log_info("Using multicore parallelization with ", njobs, " workers")
            future::plan(future::multicore, workers = njobs)
        } else {
            .log_info("Using multisession parallelization with ", njobs, " workers")
            future::plan(future::multisession, workers = njobs)
        }
        on.exit(future::plan(future::sequential), add = TRUE)
    } else {
        .log_info("Using sequential process ing (njobs=1)")
        future::plan(future::sequential)
    }
    # Use progressr for cross-platform progress reporting
    .log_info("Processing ", length(chromosomes), " chromosomes...", level = 2)

    if (verbose > 0) {
        p_con <- progressr::progressor(steps = length(chromosomes))
    }
    ret <- future.apply::future_lapply(
        X = chromosomes,
        future.seed=TRUE,
        FUN = function(chr) {
            op <- options(warn = 2)$warn
            .log_step("Subsetting DMPs for chromosome ", chr, " ...", level = 3)

            m <- dmps_locs$chr == chr
            cdmps_tsv <- dmps_tsv[(dmps_tsv[, dmps_tsv_id_col] %in% rownames(dmps_locs)), , drop = FALSE]
            unique_groups <- unique(cdmps_tsv[, dmp_group_col])
            seed_group_inds <- split(seq_len(nrow(cdmps_tsv)), cdmps_tsv[, dmp_group_col])
            gdmps_tsv_list <- list()
            for (dmp_group in unique_groups) {
                gdmps_tsv_list[[dmp_group]] <- cdmps_tsv[
                    seed_group_inds[[dmp_group]], ,
                    drop = FALSE
                ]
                rownames(gdmps_tsv_list[[dmp_group]]) <- gdmps_tsv_list[[dmp_group]][, dmps_tsv_id_col]
            }
            cdmps <- dmps[m]
            cdmps_beta <- dmps_beta[m, , drop = FALSE]
            cdmps_locs <- dmps_locs[m, , drop = FALSE]
            .log_step("Transforming beta subset to matrix...", level = 3)
            # Ensure plain numeric matrix to avoid per-iteration unlist()
            if (!is.matrix(cdmps_beta)) cdmps_beta <- as.matrix(cdmps_beta)
            storage.mode(cdmps_beta) <- "double"
            .log_success("Beta subset transformed to matrix", level = 3)
            # Collect rows and combine once (avoid repeated data.frame rbind)
            dmr_list <- vector("list", length = 128L)
            dmr_n <- 0L

            start_ind <- 1L
            corr_pval <- 1
            dmr_dmps_inds <- integer(0)
            .log_step("Testing DMP connectivity on chromosome ", chr, " ...", level = 3)
            corr_ret <- .testConnectivityBatch(
                sites_beta = cdmps_beta,
                group_inds = group_inds,
                casecontrol = case_mask,
                max_lookup_dist = max_lookup_dist,
                sites_locs = cdmps_locs,
                max_pval = max_pval,
                aggfun = aggfun,
                pval_mode = pval_mode,
                empirical_strategy = empirical_strategy,
                nperm = nperm,
                mid_p = mid_p,
                perm_seed = if (is.null(perm_seed)) NULL else as.integer(perm_seed)
            )
            stopifnot(nrow(corr_ret) == nrow(cdmps_beta) - 1)
            if (verbose >= 3) {
                dir.create("debug", showWarnings = FALSE)
                write.table(corr_ret,
                    file = paste0("debug/dmp_connectivity_chr", chr, ".tsv"),
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE,
                    quote = FALSE
                )
            }
            .log_success("DMP connectivity tested.", level = 3)
            breakpoints <- c(which(!corr_ret$connected), nrow(cdmps_beta))
            start_ind <- 1
            for (bp_ind in seq_len(length(breakpoints))) {
                end_ind <- breakpoints[bp_ind]
                if (end_ind - start_ind + 1 < min_dmps) {
                    .log_info("Skipping DMR from DMP ", start_ind, " to DMP ", end_ind, " (id: ", chr, ":", cdmps_locs[start_ind, "pos"], "-", cdmps_locs[end_ind, "pos"], ") due to insufficient number of connected DMPs (", end_ind - start_ind + 1, " < ", min_dmps, ").", level = 4)
                    start_ind <- breakpoints[bp_ind] + 1
                    next
                }
                .log_step("Registering ", bp_ind, "/", (length(breakpoints) - 1), " DMR from DMP ", start_ind, " to DMP ", end_ind, " (id: ", chr, ":", cdmps_locs[start_ind, "pos"], "-", cdmps_locs[end_ind, "pos"], ")", level = 3)
                dmr_dmps_inds <- seq.int(start_ind, end_ind)
                if (end_ind == start_ind) {
                    corr_pval <- NA
                } else {
                    corr_pval <- aggfun(corr_ret$pval[dmr_dmps_inds[-length(dmr_dmps_inds)]], na.rm = TRUE)
                }
                if (end_ind < nrow(cdmps_beta)) {
                    stop_reason <- corr_ret$reason[[end_ind]]
                } else {
                    stop_reason <- "end of chromosome"
                }
                for (dmp_group in unique_groups) {
                    gdmps_tsv <- gdmps_tsv_list[[dmp_group]]
                    dmr_dmps_tsv <- gdmps_tsv[cdmps[dmr_dmps_inds], , drop = FALSE]
                    for (num.col in c("cases_num", "controls_num")) {
                        if (!num.col %in% colnames(dmr_dmps_tsv)) {
                            dmr_dmps_tsv[, num.col] <- sum(!is.na(dmr_dmps_tsv[, "cases_beta"]))
                        }
                    }
                    for (opt.col in c("cases_beta_sd", "controls_beta_sd")) {
                        if (!opt.col %in% colnames(dmr_dmps_tsv)) {
                            dmr_dmps_tsv[, opt.col] <- NA
                        }
                    }
                    new_dmr <- list(
                        chr = chr,
                        start_dmp = cdmps[[start_ind]],
                        end_dmp = cdmps[[end_ind]],
                        start_dmp_pos = cdmps_locs[start_ind, "pos"],
                        end_dmp_pos = cdmps_locs[end_ind, "pos"],
                        dmps_num = length(dmr_dmps_inds),
                        delta_beta = aggfun(abs(dmr_dmps_tsv[, "delta_beta"])) * sign(sum(sign(dmr_dmps_tsv[, "delta_beta"]))),
                        delta_beta_sd = stats::sd(dmr_dmps_tsv[, "delta_beta"]),
                        delta_beta_min = min(dmr_dmps_tsv[, "delta_beta"]),
                        delta_beta_max = max(dmr_dmps_tsv[, "delta_beta"]),
                        cases_beta = aggfun(abs(dmr_dmps_tsv[, "cases_beta"])) * sign(sum(sign(dmr_dmps_tsv[, "cases_beta"]))),
                        controls_beta = aggfun(abs(dmr_dmps_tsv[, "controls_beta"])) * sign(sum(sign(dmr_dmps_tsv[, "controls_beta"]))),
                        cases_num = min(dmr_dmps_tsv[, "cases_num"]),
                        controls_num = min(dmr_dmps_tsv[, "controls_num"]),
                        corr_pval = corr_pval,
                        stop_connection_reason = stop_reason,
                        dmps = paste(cdmps[dmr_dmps_inds], collapse = ",")
                    )
                    new_dmr[[pval_col]] <- aggfun(dmr_dmps_tsv[, pval_col])
                    if (dmp_group_col != "_DUMMY_DMP_GROUP_COL_") {
                        new_dmr[[dmp_group_col]] <- dmp_group
                    }
                    dmr_n <- dmr_n + 1L
                    if (dmr_n > length(dmr_list)) length(dmr_list) <- length(dmr_list) * 2L
                    dmr_list[[dmr_n]] <- new_dmr
                    start_ind <- breakpoints[bp_ind] + 1

                }
                .log_success("DMR registered.",
                    level = 3
                )
            }
            if (verbose > 0 && exists("p_con")) p_con()
            dmrs <- if (dmr_n > 0L) data.table::rbindlist(dmr_list[seq_len(dmr_n)], fill = TRUE) else data.frame()
            if (nrow(dmrs)) {
                rownames(dmrs) <- seq_len(nrow(dmrs))
                dmrs[, "chr"] <- chr
            }
            options(warn = op)
            dmrs
        }
    )


    if (inherits(ret[[1]], "try-error")) {
        stop(ret)
    }
    dmrs <- as.data.frame(do.call(rbind, ret))
    .log_info("Summary:\n\t", paste(capture.output(summary(dmrs)), collapse = "\n\t"), level = 1)
    if (nrow(dmrs) == 0) {
        .log_warn("No DMRs remain after filtering based on connected DMP number.")
        if (!is.null(output_prefix)) {
            for (f in c(".methylation.tsv.gz", ".tsv.gz")) {
                gzfile <- gzfile(paste0(output_prefix, f), "w", compression = 2)
                close(gzfile)
            }
        }
        return(NULL)
    }
    if (verbose >= 2) {
        dir.create("debug", showWarnings = FALSE)
        write.table(dmrs,
            file = file.path("debug", "01_dmrs_from_connected_dmps.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
    }
    cases_num <- dmrs$cases_num
    controls_num <- dmrs$controls_num
    if (anyNA(c(cases_num, controls_num))) {
        .log_warn("NAs introduced while coercing cases_num / controls_num to numeric; replacing NAs with 1 to avoid division errors.")
        cases_num[is.na(cases_num)] <- 1
        controls_num[is.na(controls_num)] <- 1
    }

    if (dmp_group_col %in% colnames(dmrs)) {
        ungrouped_dmrs <- dmrs[
            dmrs[, dmp_group_col] == dmrs[1, dmp_group_col],
            c("chr", "start_dmp", "end_dmp", "start_dmp_pos", "end_dmp_pos", "dmps_num", "corr_pval")
        ]
    } else {
        ungrouped_dmrs <- dmrs[, c("chr", "start_dmp", "end_dmp", "start_dmp_pos", "end_dmp_pos", "dmps_num", "corr_pval")]
    }
    .log_success("Initial DMRs formed: ", nrow(dmrs), level = 1)
    .log_info("Summary:\n\t", paste(capture.output(summary(ungrouped_dmrs)), collapse = "\n\t"), level = 1)
    .log_step("Expanding DMRs on neighborhood CpGs..", level = 1)
    # Set up progress tracking for DMR expansion
    n_dmrs <- nrow(ungrouped_dmrs)

    if (verbose > 0) {
        p_ext <- progressr::progressor(steps = n_dmrs)
    }
    .log_step("Expanding ", n_dmrs, " DMRs using up to ", njobs, " parallel jobs...", level = 2)
    ret <- future.apply::future_apply(
        X = ungrouped_dmrs,
        MARGIN = 1,
        simplify = FALSE,
        future.seed=TRUE,
        future.scheduling = ceiling(n_dmrs / njobs),
        FUN = function(dmr) {
            op <- options(warn = 2)$warn
            x <- .expandDMRs(
                dmr = dmr,
                group_inds = group_inds,
                max_pval = max_pval,
                beta_file_handler = beta_file_handler,
                beta_row_names = beta_row_names,
                beta_col_names = beta_col_names,
                casecontrol = case_mask,
                expansion_step = expansion_step,
                expansion_relaxation = expansion_relaxation,
                min_cpgs = min_cpgs,
                min_cpg_delta_beta = min_cpg_delta_beta,
                sorted_locs = sorted_locs,
                max_lookup_dist = max_lookup_dist,
                aggfun = aggfun,
                pval_mode = pval_mode,
                empirical_strategy = empirical_strategy,
                nperm = nperm,
                mid_p = mid_p,
                perm_seed = if (is.null(perm_seed)) NULL else as.integer(perm_seed)
            )
            options(warn = op)
            if (verbose > 0 && exists("p_ext")) p_ext()
            x
        }
    )
    .log_success("DMR expansion complete.", level = 2)

    .log_step("Post-processing extended DMRs..", level = 2)
    if (inherits(ret, "try-error")) {
        stop(ret)
    }
    extended_dmrs <- as.data.frame(do.call(rbind, ret))
    extended_dmrs$end <- as.numeric(extended_dmrs$end)
    extended_dmrs$start <- as.numeric(extended_dmrs$start)
    extended_dmrs$start_dmp_pos <- as.numeric(extended_dmrs$start_dmp_pos)
    extended_dmrs$end_dmp_pos <- as.numeric(extended_dmrs$end_dmp_pos)
    extended_dmrs$dmps_num <- as.numeric(extended_dmrs$dmps_num)
    extended_dmrs$corr_pval <- as.numeric(extended_dmrs$corr_pval)

    end_less_than_start <- extended_dmrs$end - extended_dmrs$start < 0

    if (any(end_less_than_start)) {
        .log_warn(
            paste(
                sum(end_less_than_start),
                "DMRs have been assigned an end larger than start ! (CODE BUG TO BE REPORTED)"
            )
        )
        .log_warn(
            "Those are: \n\t",
            paste0(capture.output(print(extended_dmrs[end_less_than_start, ])), collapse = "\n\t")
        )

        .log_warn("Removing them..")
        extended_dmrs <- extended_dmrs[!end_less_than_start, ]
        .log_warn("Remaining: ", nrow(extended_dmrs))
    }
    all_locs_inds <- rownames(sorted_locs)
    names(all_locs_inds) <- all_locs_inds
    all_locs_inds[seq_along(all_locs_inds)] <- seq_along(all_locs_inds)
    extended_dmrs$start_ind <- as.numeric(all_locs_inds[extended_dmrs$start_cpg])
    extended_dmrs$end_ind <- as.numeric(all_locs_inds[extended_dmrs$end_cpg])
    extended_dmrs$sup_cpgs_num <- extended_dmrs$end_ind - extended_dmrs$start_ind + 1

    extended_dmrs$id <- paste0(extended_dmrs$chr, ":", extended_dmrs$start, "-", extended_dmrs$end)
    
    .log_success("Post-processing complete.", level = 2)
    .log_success("Extended DMRs formed: ", nrow(extended_dmrs), level = 1)

    .log_step("Finding GC content of DMRs..", level = 1)

    extended_dmrs_ranges <- GenomicRanges::makeGRangesFromDataFrame(
        extended_dmrs,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqnames.field = "chr",
        seqinfo = GenomeInfoDb::Seqinfo(genome = genome),
        na.rm = TRUE
    )
    sequences <- getDMRSequences(extended_dmrs_ranges, genome)
    extended_dmrs$cpgs_num <- stringr::str_count(sequences, "(CG)|(GC)")
    colnames(extended_dmrs)[colnames(extended_dmrs) == "seqnames"] <- "chr"

    extended_dmrs[extended_dmrs$cpgs_num == 0, "cpgs_num"] <- 1

    extended_dmrs$dmps_num_adj <- ceiling(extended_dmrs$cpgs_num / extended_dmrs$sup_cpgs_num * extended_dmrs$dmps_num)
    .log_success("CpG content calculated.", level = 1)
    .log_info("Summary of extended DMRs before filtering based on supporting CpGs and adjusted DMPs number:\n\t", paste(capture.output(summary(extended_dmrs)), collapse = "\n\t"), level = 1)
    if (verbose >= 2) {
        dir.create("debug", showWarnings = FALSE)
        write.table(extended_dmrs,
            file = file.path("debug", "02_extended_dmrs_prior_filtering.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
    }
    filtered_dmrs <- extended_dmrs[
        extended_dmrs$dmps_num_adj >= min_adj_dmps &
            extended_dmrs$sup_cpgs_num >= min_cpgs, ,
        drop = FALSE
    ]
    .log_info(
        "Keeping ",
        nrow(filtered_dmrs),
        " out of ",
        nrow(extended_dmrs),
        " with at least ",
        min_adj_dmps,
        " (adjusted) supporting DMPs and ",
        min_cpgs,
        " supporting CpGs.",
        level = 1
    )
    extended_dmrs <- filtered_dmrs
    if (nrow(extended_dmrs) == 0) {
        .log_warn("No DMRs passed the filtering based on min_cpgs and min_adj_dmps criteria.")
        if (!is.null(output_prefix)) {
            for (f in c(".methylation.tsv.gz", ".tsv.gz")) {
                gzfile <- gzfile(file.path(output_dir, paste0(output_prefix, f)), "w", compression = 2)
                close(gzfile)
            }
        }
        return(NULL)
    }
    .log_step("Merging extended DMRs with original DMRs table to fill in missing information..", level = 2)
    ne_dmrs_cols_to_keep <- colnames(dmrs)
    ne_dmrs_cols_to_keep <- c("start_dmp", "end_dmp", ne_dmrs_cols_to_keep[!ne_dmrs_cols_to_keep %in% colnames(extended_dmrs)])
    dmrs <- merge(extended_dmrs, dmrs[, ne_dmrs_cols_to_keep], by = c("start_dmp", "end_dmp"))
    dmrs <- dmrs[stringr::str_order(dmrs[, "id"], numeric = TRUE), ]
    .log_success("Merging done.", level = 2)
    .log_info("Summary of final DMRs:\n\t", paste(capture.output(summary(dmrs)), collapse = "\n\t"), level = 1)
    dmrs_granges <- GenomicRanges::makeGRangesFromDataFrame(
        dmrs,
        keep.extra.columns = TRUE,
        seqnames.field = "chr",
        seqinfo = GenomeInfoDb::Seqinfo(genome = genome),
        na.rm = TRUE
    )

    if (!is.null(output_prefix)) {
        dmrs_file <- paste0(output_prefix, "dmrs.tsv.gz")
        .log_step("Saving DMRs to ", dmrs_file, "..")
        gz <- gzfile(dmrs_file, "w")
        write.table(
            dmrs,
            gz,
            sep = "\t",
            quote = FALSE,
            qmethod = "double",
            col.names = TRUE,
            row.names = FALSE
        )
        close(gz)
        .log_success("DMRs saved.")
    }
    invisible(dmrs_granges)
}
