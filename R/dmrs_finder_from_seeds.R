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
#' @param beta Character. Path to the beta value file, or a tabix file, or a beta matrix, or a BetaHandler object, or a bed file. If a bed file is provided, it must at least contain bed_chrom_col and bed_chrom_start, followed by samples names in the given pheno, with corresponging beta values, and it will be converted to a tabix-indexed beta file internall, with the locations separately saved and queried as a bigmatrix.
#' @param dmps Character. Path to the DMPs TSV file or the dmps dataframe, in a format like the one produced by dmpFinder.
#' @param pheno Data frame. Phenotype data.
#' @param dmps_tsv_id_col Character. Column name or index for DMP identifiers in the DMPs TSV file. Default is 1.
#' @param sample_group_col Character. Column name for sample group information in the phenotype data. Default is NULL.
#' @param casecontrol_col Boolean Column in pheno for case (TRUE/1) / control (FALSE/0) status . If NULL, controls will be assumed to be the first level of sample_group_col. Default is NULL.
#' @param min_cpg_delta_beta Numeric. Minimum delta beta value for CpGs. Default is 0.
#' @param expansion_step Numeric. Step size for expanding DMRs. Increasing it means higher memory usage and faster computation. Default is 500.
#' @param array Character. Type of array used (e.g., "450K", "EPIC", "EPICv2", "27K"). Ignored if using mm10 genome.
#' @param genome Character. Genome version (e.g., "hg19", "hg38", "mm10"). Default is "hg19".
#' @param max_pval Numeric. Maximum p-value to assume DMPs correlation is significant. Default is 0.05.
#' @param pval_mode Character. "parametric" (default) to use t-based correlation p-values during connectivity testing, or "empirical" to use permutation-based p-values.
#' @param empirical_strategy Character. When pval_mode = "empirical": "auto" (default) uses Monte Carlo for groups with <6 samples and permutations for groups with >=6 samples; "montecarlo" always uses Monte Carlo; "permutations" always uses permutations.
#' @param ntries Integer. Number of permutations when pval_mode = "empirical". Default is 0 (disabled).
#' @param tries_seed Integer or NULL. RNG seed for reproducibility when pval_mode = "empirical". Default is NULL.
#' @param max_lookup_dist Numeric. Maximum distance to look up for adjacent DMPs belonging to the same DMR. Default is 10000.
#' @param min_dmps Numeric. Minimum number of connected DMPs in a DMR. Default is 1.
#' @param min_adj_dmps Numeric. Minimum number of DMPs, adjusted by CpG density, in a DMR after extension. Default is 1.
#' @param min_cpgs Numeric. Minimum number of CpGs in a DMR after extension. Default is 50.
#' @param aggfun Function or character. Aggregation function to use when calculating delta beta values and p-values of DMRs. Can be "median", "mean", or a function (e.g., median, mean). Default is "median".
#' @param ignored_sample_groups Character vector. Sample groups to ignore, separated by commas. Default is NULL.
#' @param output_prefix Character. Identifier for the output files. If not provided, no output will be saved. Default is NULL.
#' @param njobs Numeric. Number of parallel jobs to use. Default is the number of available cores.
#' @param verbose Numeric. Level of verbosity for logging messages, from 0 (not verbose) to 5 (very very verbose). Default is 1.
#' @param memory_threshold_mb Numeric. Memory threshold in MB for loading beta files. Default is 500.
#' @param beta_row_names_file Character. Path to a file containing row names for the beta values. If not provided, row names will be read from the beta file. Default is NULL.
#' @param annotate_with_genes Logical. Whether to annotate DMRs with overlapping genes. Default is TRUE.
#' @param bed_provided Logical. Whether the beta file is provided as a BED file. Default is FALSE.
#' @param bed_chrom_col Character. Column name for chromosome in the BED file. Default is "chrom".
#' @param bed_start_col Character. Column name for start position in the BED file. Default is "start".
#' @param bed_end_col Character or NULL. Column name for end position in the BED file. If NULL, start_col + 1 will be used. Default is NULL, which means that the column is assumed missing.
#' @param bed_id_col Character or NULL. Column name for CpG identifier in the BED file. Default is NULL, which means that the column is assumed missing.
#' @param bed_score_col Character or NULL. Column name for score in the BED file. Default is NULL, which means that the column is assumed missing.
#' @param bed_strand_col Character or NULL. Column name for strand in the BED file. Default is NULL, which means that the column is assumed missing.
#' @param .load_debug Logical. If TRUE, enables debug mode for loading beta files. Default is FALSE.
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
#' if (!requireNamespace("DMRsegaldata", quietly = TRUE)) {
#'     remotes::install_github("CMG-UA/DMRsegaldata")
#' }
#' library(DMRsegaldata)
#' beta <- DMRsegaldata::beta
#' dmps <- DMRsegaldata::dmps
#' pheno <- DMRsegaldata::pheno
#' # Find DMRs
#' dmrs <- findDMRsFromSeeds(
#'     beta = beta,
#'     dmps = dmps,
#'     pheno = pheno
#' )
#' }
#'
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom future.apply future_lapply future_mapply
#' @importFrom future availableCores
#' @importFrom progressr progressor
#' @importFrom stringr str_count str_order
#' @importFrom readr read_tsv
#' @importFrom data.table fread fwrite
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


.buildConnectivityArray <- function(
    beta_handler, group_inds, sorted_locs, max_pval = 0.05, min_delta_beta = 0, casecontrol = NULL, max_lookup_dist = 1000,
    chunk_size = 1000, aggfun = median, empirical_strategy = "auto",
    pval_mode = "empirical", ntries = 500, mid_p = TRUE, tries_seed = 42, njobs = 1, bed_provided = FALSE
) {
    # split sorted_locs into chunks of chunk_size, respecting chromosome boundaries
    splits <- c()
    chr_ends <- c()
    for (chr in unique(sorted_locs$chr)) {
        chr_inds <- which(sorted_locs$chr == chr)
        chr_start <- min(chr_inds)
        chr_end <- max(chr_inds)
        chr_ends <- c(chr_ends, chr_end)
        for (chunk_start in seq(chr_start, chr_end, by = chunk_size)) {
            chunk_end <- min(chunk_start + chunk_size + 1, chr_end)
            splits <- rbind(splits, c(chunk_start, chunk_end))
        }
    }
    verbose <- getOption("DMRsegal.verbose", 0)
    if (verbose > 0) {
        p_ext <- progressr::progressor(steps = nrow(splits))
    }
    ret <- future.apply::future_lapply(
        X = seq_len(nrow(splits)),
        future.seed = TRUE,
        future.stdout = NA,
        future.globals = c(
            "beta_handler",
            "bed_provided",
            "group_inds",
            "casecontrol",
            "max_pval",
            "min_delta_beta",
            "aggfun",
            "pval_mode",
            "empirical_strategy",
            "ntries",
            "mid_p",
            "tries_seed",
            "chr_ends",
            "verbose",
            "p_ext"
        ),
        FUN = function(split_ind) {
            split <- splits[split_ind, ]
            beta_locs <- beta_handler$getBetaLocs()
            if (bed_provided) {
                chunk_beta <- beta_handler$getBeta(row_names = rownames(beta_locs)[split[1]:split[2]])
            } else {
                chunk_beta <- beta_handler$getBeta(
                    row_names = split[1]:split[2]
                )
            }
            x <- .testConnectivityBatch(
                sites_beta = chunk_beta,
                group_inds = group_inds,
                casecontrol = casecontrol,
                max_pval = max_pval,
                min_delta_beta = min_delta_beta,
                sites_locs = beta_locs[split[1]:split[2], ],
                aggfun = aggfun,
                pval_mode = pval_mode,
                empirical_strategy = empirical_strategy,
                ntries = ntries,
                mid_p = mid_p,
                tries_seed = if (is.null(tries_seed)) NULL else as.integer(tries_seed)
            )
            if (split[2] %in% chr_ends) {
                add <- data.frame(
                    connected = FALSE,
                    pval = NA,
                    reason = "end-of-input",
                    first_failing_group = NA
                )

                if (min_delta_beta > 0 && !is.null(casecontrol)) {
                    add["delta_beta"] <- NA
                }
                x <- rbind(x, add)
            }
            if (verbose > 0) {
                p_ext()
            }
            x
        }
    )
    do.call(rbind, ret)
}


#' @keywords internal
#' @noRd
.expandDMRs <- function(dmr,
                        chr_array,
                        chr_locs,
                        min_cpg_delta_beta = 0,
                        min_cpgs = 3,
                        expansion_step = 500,
                        bed_provided = FALSE) {
    .log_step("Expanding DMR..", level = 4)
    dmr_start <- dmr["start_dmp"]
    dmr_end <- dmr["end_dmp"]
    if (bed_provided) {
        dmr_start_ind <- as.integer(dmr_start)
        dmr_end_ind <- as.integer(dmr_end)
    } else {
        dmr_start_ind <- which(rownames(chr_locs) == dmr_start)
        dmr_end_ind <- which(rownames(chr_locs) == dmr_end)
    }
    if (length(dmr_start_ind) == 0) {
        stop("Could not find the start CpG ", dmr_start, " in the beta file row names.")
    }
    if (length(dmr_end_ind) == 0) {
        stop("Could not find the end CpG ", dmr_end, " in the beta file row names.")
    }


    .check_upstream <- function(ustream_exp, exp_step) {
        ustream_stop_reason <- NULL
        ustream_end_lookup_site_ind <- ustream_exp - 1 # corr_ret[i] corresponds to the connection between i and i+1, so corr_ret[ustream_exp - 1] reports the connection between ustream_exp - 1 and ustream_exp
        if (ustream_end_lookup_site_ind < 0) {
            ustream_stop_reason <- "end-of-input"
            return(list(
                ustream_stop_reason = ustream_stop_reason, ustream_exp = ustream_exp
            ))
        }
        ustream_start_lookup_site_ind <- max(1, ustream_end_lookup_site_ind - exp_step)
        corr_ret <- chr_array[ustream_start_lookup_site_ind:ustream_end_lookup_site_ind, , drop = FALSE]
        corr_ret <- corr_ret[rev(seq_len(nrow(corr_ret))), , drop = FALSE]
        if (nrow(corr_ret) == 1) {
            ustream_stop_reason <- "end-of-input"
            return(list(
                ustream_end_lookup_site_ind = ustream_end_lookup_site_ind,
                ustream_stop_reason = ustream_stop_reason,
                ustream_exp = ustream_exp
            ))
        }
        # pick the first run of FALSE
        fail_runs <- which(!corr_ret$connected)
        # [----+---0] where - is connected, + is not connected, 0 stands for the current DMR start (ustream_exp), then fail_start_idx is the first + from the right, 4 in this case. That means that the 4th from the right failed to connect to the 3rd from the right.
        if (length(fail_runs) > 0) {
            fail_start_idx <- fail_runs[[1]] # first failing index in the reversed corr_ret
            ustream_exp <- ustream_end_lookup_site_ind - fail_start_idx + 2
            ustream_stop_reason <- corr_ret$reason[fail_start_idx]
            return(list(
                ustream_stop_reason = ustream_stop_reason, ustream_exp = ustream_exp
            ))
        }
        ustream_exp <- ustream_start_lookup_site_ind
        list(
            ustream_stop_reason = ustream_stop_reason, ustream_exp = ustream_exp
        )
    }

    .check_downstream <- function(dstream_exp, exp_step) {
        dstream_stop_reason <- NULL
        dstream_start_lookup_site_ind <- dstream_exp
        dstream_end_lookup_site_ind <- min(dstream_start_lookup_site_ind + exp_step, nrow(chr_locs))
        if (dstream_end_lookup_site_ind == dstream_start_lookup_site_ind) {
            dstream_stop_reason <- "end-of-input"
            return(list(
                dstream_stop_reason = dstream_stop_reason, dstream_exp = dstream_exp
            ))
        }
        corr_ret <- chr_array[dstream_start_lookup_site_ind:dstream_end_lookup_site_ind, , drop = FALSE]
        # pick the first run of FALSE
        fail_runs <- which(!corr_ret$connected)
        # [0----+---] where - is connected, + is not connected, 0 stands for the current DMR start (dstream_exp), then fail_start_idx is the first + from the left, 5 in this case. That means that the 5th from the left failed to connect to the 6th from the left.
        # That means we need to expand the DMR to include the 5th from the left, which is dstream_exp + fail_start_idx - 1
        if (length(fail_runs) > 0) {
            fail_start_idx <- fail_runs[[1]]
            dstream_exp <- dstream_start_lookup_site_ind + fail_start_idx - 1
            dstream_stop_reason <- corr_ret$reason[fail_start_idx]
            return(list(
                dstream_stop_reason = dstream_stop_reason, dstream_exp = dstream_exp
            ))
        }
        dstream_exp <- dstream_exp + exp_step
        list(
            dstream_stop_reason = dstream_stop_reason, dstream_exp = dstream_exp
        )
    }

    ustream_exp <- dmr_start_ind[[1]]
    ustream_stop_reason <- NULL
    dstream_exp <- dmr_end_ind[[1]]
    dstream_stop_reason <- NULL

    t <- 0
    while (TRUE) {
        exp_step <- expansion_step
        if (t == 0) { # first iteration, use min_cpgs and remove the DMRs that are not long enough
            ccpgs <- dstream_exp - ustream_exp + 1
            if (ccpgs < (min_cpgs)) {
                .log_info("DMR  too short (", ccpgs, " CpGs). Expanding to reach min_cpgs=", min_cpgs, ".", level = 3)
                exp_step <- min_cpgs - ccpgs
            }
            .log_info("Number of CpGs in DMR: ", ccpgs, level = 5)
        }
        .log_info("Expansion step size: ", exp_step, " bp.", level = 5)
        .log_step("Checking upstream expansion...", level = 5)
        if (is.null(ustream_stop_reason)) {
            res <- .check_upstream(ustream_exp, exp_step)
            ustream_stop_reason <- res$ustream_stop_reason
            if (res$ustream_exp > ustream_exp) {
                .log_info("Upstream expanded by ", ustream_exp - res$ustream_exp, " CpGs.", level = 5)
            }
            ustream_exp <- res$ustream_exp
        }
        .log_success("Upstream expansion checked.", level = 5)
        .log_step("Checking downstream expansion...", level = 5)
        if (is.null(dstream_stop_reason)) {
            res <- .check_downstream(dstream_exp, exp_step)
            dstream_stop_reason <- res$dstream_stop_reason
            if (res$dstream_exp > dstream_exp) {
                .log_info("Downstream expanded by ", res$dstream_exp - dstream_exp, " CpGs.", level = 5)
            }
            dstream_exp <- res$dstream_exp
        }
        .log_success("Downstream expansion checked.", level = 5)
        if (t == 0) {
            new_ccpgs <- dstream_exp - ustream_exp + 1
            .log_info("Number of CpGs in expanded DMR: ", new_ccpgs, " from ", ccpgs, level = 3)
            if (new_ccpgs < min_cpgs) {
                ustream_stop_reason <- "min_cpgs_reached"
                dstream_stop_reason <- "min_cpgs_reached"
                .log_info("DMR could not reach min_cpgs=", min_cpgs, " after expansion (", new_ccpgs, "). Stopping expansion.", level = 3)
            }

            t <- 1
        }
        if (!is.null(ustream_stop_reason) && !is.null(dstream_stop_reason)) {
            break
        }
    }
    .log_step("Finalizing expanded DMR.", level = 5)
    if (!bed_provided) {
        dmr["start_cpg"] <- rownames(chr_locs)[ustream_exp]
        dmr["end_cpg"] <- rownames(chr_locs)[dstream_exp]
    } else {
        dmr["start_cpg"] <- ustream_exp
        dmr["end_cpg"] <- dstream_exp
    }
    dmr["start"] <- chr_locs[dmr["start_cpg"], "start"]
    dmr["end"] <- chr_locs[dmr["end_cpg"], "start"]
    dmr["downstream_cpg_expansion_stop_reason"] <- dstream_stop_reason
    dmr["upstream_cpg_expansion_stop_reason"] <- ustream_stop_reason
    .log_success("Expanded DMR finalized: (start_cpg: ", dmr["start_cpg"], ", end_cpg: ", dmr["end_cpg"], ").", level = 5)
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
#' @param sites_locs Data frame with chr and start columns for each site
#' @param aggfun Aggregation function for computing group means (default: mean)
#' @param pval_mode Character. "parametric" (default) for t-based p-values, or "empirical" for permutation-based p-values.
#' @param empirical_strategy Character. When pval_mode = "empirical": "auto" (default) uses Monte Carlo for groups with <6 samples and permutations for groups with >=6 samples
#' @param ntries Integer. Number of tries when pval_mode = "empirical". Default is 500. If set to 0, uses 500 for Monte Carlo and min(2m!,500) for permutations, where m is the number of samples.
#' @param mid_p Logical. Use mid-p adjustment for empirical p-values (reduces tie conservatism). Default is FALSE.
#' @param tries_seed Integer or NULL. RNG seed for reproducible permutations when pval_mode = "empirical". Default is NULL.
#'
#' @return Data frame with columns: connected, pval, delta_beta, reason, first_failing_group, stop_reason
#' @keywords internal
#' @noRd
.testConnectivityBatch <- function(sites_beta, group_inds, max_pval,
                                   casecontrol = NULL, min_delta_beta = 0,
                                   max_lookup_dist = NULL, sites_locs = NULL, aggfun = mean,
                                   pval_mode = c("parametric", "empirical"),
                                   empirical_strategy = c("auto", "montecarlo", "permutations"),
                                   ntries = 0, mid_p = FALSE, tries_seed = NULL) {
    pval_mode <- strex::match_arg(pval_mode, ignore_case = TRUE)
    empirical_strategy <- strex::match_arg(empirical_strategy, ignore_case = TRUE)
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
    pvals <- rep(NA_real_, n_pairs)
    reasons <- rep("", n_pairs)
    failing_groups <- rep(NA_character_, n_pairs)

    # Check distance condition if provided (vectorized)
    if (!is.null(max_lookup_dist) && !is.null(sites_locs)) {
        dists <- sites_locs$start[2:n_sites] - sites_locs$start[1:(n_sites - 1)]
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


        if (pval_mode == "parametric") {
            tstats <- cors * sqrt(dfs / pmax(1e-12, 1 - cors * cors))
            # Compute t-statistics (vectorized)
            na_tstat <- is.na(tstats) & connected
            connected[na_tstat] <- FALSE
            reasons[na_tstat] <- "na tstat"
            failing_groups[na_tstat] <- g

            ps[connected] <- -2 * expm1(pt(abs(tstats[connected]), df = dfs[connected], log.p = TRUE))
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
                if (!is.null(tries_seed)) {
                    set.seed(tries_seed)
                }
                on.exit(
                    {
                        if (old_seed_exists) {
                            assign(".Random.seed", old_seed, envir = .GlobalEnv)
                        }
                    },
                    add = TRUE
                )
                counts_ge <- integer(n_pairs)
                counts_eq <- integer(n_pairs)
                # Number of samples in this group
                m <- ncol(y_mat)
                # Empirical strategy: auto -> MonteCarlo when n_valid < 6  or permutations when n_valid >= 6
                if (empirical_strategy == "auto") {
                    do_permutations <- ncol(x_mat) >= 6
                } else if (empirical_strategy == "montecarlo") {
                    do_permutations <- FALSE
                } else if (empirical_strategy == "permutations") {
                    do_permutations <- TRUE
                } else {
                    stop("Unknown empirical_strategy: ", empirical_strategy)
                }
                if (do_permutations) {
                    min_possible_pval <- 1 / (1 + factorial(m))
                    if (min_possible_pval > max_pval_corrected) {
                        .log_warn("Skipping empirical p-value computation for group '", g, "' since minimum possible p-value (", min_possible_pval, ") exceeds max_pval_corrected (", max_pval_corrected, ").")
                        next
                    }
                }
                .log_info("Computing empirical p-values for group '", g, "' using ", if (do_permutations) "permutations" else "Monte Carlo", " with ", ntries, " tries.", level = 4)
                if (ntries == 0) {
                    if (do_permutations) {
                        ntries <- min(500L, factorial(m))
                    } else {
                        ntries <- 500
                    }
                }
                for (b in seq_len(ntries)) {
                    # Permute sample labels (columns) only for y; x remains fixed
                    if (do_permutations) {
                        perm <- sample.int(m, size = m, replace = FALSE)
                        yp <- y_mat[, perm, drop = FALSE]
                    } else {
                        yp <- matrix(runif(n = nrow(y_mat) * m, min = min(y_mat, na.rm = TRUE), max = max(y_mat, na.rm = TRUE)), nrow = nrow(y_mat), ncol = m)
                    }
                    ym <- rowMeans(yp, na.rm = TRUE)
                    yc <- sweep(yp, 1L, ym, FUN = "-")
                    sxy <- rowSums(x_centered * yc, na.rm = TRUE)
                    sy2 <- rowSums(yc^2, na.rm = TRUE)
                    rperm <- sxy / sqrt(sum_x2 * sy2)
                    comp_mask <- is.finite(rperm)
                    if (any(comp_mask)) {
                        ap <- abs(rperm[comp_mask])
                        ao <- abs(cors[comp_mask])
                        counts_ge[comp_mask] <- counts_ge[comp_mask] + (ap > ao)
                        counts_eq[comp_mask] <- counts_eq[comp_mask] + (ap == ao)
                    }
                }
                if (mid_p) {
                    ps <- (counts_ge + 0.5 * counts_eq + 1) / (ntries + 1)
                } else {
                    ps <- (counts_ge + counts_eq + 1) / (ntries + 1)
                }
            }
        }


        na_p <- is.na(ps) & connected
        connected[na_p] <- FALSE
        reasons[na_p] <- "na pval"
        failing_groups[na_p] <- g

        # Update maximum p-values across groups (vectorized)
        pvals <- pmax(pvals, ps, na.rm = TRUE)
        pvals[is.na(ps)] <- NA_real_

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
#' @param beta_handler BetaHandler object for the methylation beta values file
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
#'     beta_handler = beta_handler,
#'     beta_row_names = row_names,
#'     beta_col_names = col_names,
#'     sorted_locs = genomic_locations
#' )
#' }
#'
#' @export
extractCpgInfoFromResultDMRs <- function(dmrs,
                                         beta_handler,
                                         beta_row_names,
                                         beta_col_names,
                                         sorted_locs,
                                         output_prefix = NULL,
                                         bed_provided = FALSE,
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
        start_cpg_ind <- dmrs$start_cpg_ind[i]
        end_cpg_ind <- dmrs$end_cpg_ind[i]
        if (!bed_provided) {
            cpgs <- rownames(sorted_locs)[start_cpg_ind:end_cpg_ind]
        } else {
            cpgs <- start_cpg_ind:end_cpg_ind
        }
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
    if (!bed_provided) {
        dmrs_sites <- rownames(sorted_locs)[rownames(sorted_locs) %in% dmrs_sites]
    }

    dmrs_beta_file <- paste0(output_prefix, "dmrs_beta.tsv.gz")
    .log_step("Saving constituent CpG betas", dmrs_beta_file, " ...")
    dmrs_beta <- beta_handler$getBeta(
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

.calculateBetaStats <- function(beta_values, beta_col_names, pheno,
                                casecontrol_col,
                                aggfun) {
        cases_beta <- apply(beta_values[, pheno[beta_col_names, casecontrol_col] == 1, drop = FALSE], 1, aggfun, na.rm = TRUE)
        controls_beta <- apply(beta_values[, pheno[beta_col_names, casecontrol_col] == 0, drop = FALSE], 1, aggfun, na.rm = TRUE)
        case_sd <- apply(beta_values[, pheno[beta_col_names, casecontrol_col] == 1, drop = FALSE], 1, sd, na.rm = TRUE)
        control_sd <- apply(beta_values[, pheno[beta_col_names, casecontrol_col] == 0, drop = FALSE], 1, sd, na.rm = TRUE)
        cases_num <- rowSums(!is.na(beta_values[, pheno[beta_col_names, casecontrol_col] == 1, drop = FALSE]))
        controls_num <- rowSums(!is.na(beta_values[, pheno[beta_col_names, casecontrol_col] == 0, drop = FALSE]))
        list(
            cases_beta = cases_beta,
            controls_beta = controls_beta,
            cases_beta_sd = case_sd,
            controls_beta_sd = control_sd,
            cases_num = cases_num,
            controls_num = controls_num
        )
}


#' Find Differentially Methylated Regions (DMRs) from Differentially Methylated Positions (DMPs)
#'
#' This function identifies DMRs from a given set of DMPs and a beta value file.
#'
#' @param beta Character. Path to the beta value file, or a tabix file, or a beta matrix, or a BetaHandler object, or a bed file. If a bed file is provided, it must at least contain bed_chrom_col and bed_chrom_start, followed by samples names in the given pheno, with corresponging beta values, and it will be converted to a tabix-indexed beta file internall, with the locations separately saved and queried as a bigmatrix.
#' @param dmps Character. Path to the DMPs TSV file or the dmps dataframe, in a format like the one produced by dmpFinder.
#' @param pheno Data frame. Phenotype data.
#' @param dmps_tsv_id_col Character. Column name or index for DMP identifiers in the DMPs TSV file. Default is NULL, which corresponds to the rows names if existing, or the first column if not.
#' @param sample_group_col Character. Column name for sample group information in the phenotype data. Default is NULL.
#' @param casecontrol_col Boolean Column in pheno for case (TRUE/1) / control (FALSE/0) status . If NULL, controls will be assumed to be the first level of sample_group_col. Default is NULL.
#' @param min_cpg_delta_beta Numeric. Minimum delta beta value for CpGs. Default is 0.
#' @param expansion_step Numeric. Step size for expanding DMRs. Increasing it means higher memory usage and faster computation. Default is 500.
#' @param array Character. Type of array used (e.g., "450K", "EPIC", "EPICv2", "27K"). Ignored if using mm10 genome.
#' @param genome Character. Genome version (e.g., "hg19", "hg38", "mm10"). Default is "hg19".
#' @param max_pval Numeric. Maximum p-value to assume DMPs correlation is significant. Default is 0.05.
#' @param pval_mode Character. "parametric" (default) to use t-based correlation p-values during connectivity testing, or "empirical" to use permutation-based p-values.
#' @param empirical_strategy Character. When pval_mode = "empirical": "auto" (default) uses Monte Carlo for groups with <6 samples and permutations for groups with >=6 samples; "montecarlo" always uses Monte Carlo; "permutations" always uses permutations.
#' @param ntries Integer. Number of permutations when pval_mode = "empirical". Default is 0 (disabled).
#' @param tries_seed Integer or NULL. RNG seed for reproducibility when pval_mode = "empirical". Default is NULL.
#' @param max_lookup_dist Numeric. Maximum distance to look up for adjacent DMPs belonging to the same DMR. Default is 10000.
#' @param min_dmps Numeric. Minimum number of connected DMPs in a DMR. Default is 1.
#' @param min_adj_dmps Numeric. Minimum number of DMPs, adjusted by CpG density, in a DMR after extension. Default is 1.
#' @param min_cpgs Numeric. Minimum number of CpGs in a DMR after extension. Default is 50.
#' @param aggfun Character. Aggregation function to use ("median" or "mean") when calculating delta beta values and p-values of DMRs. Default is "median".
#' @param ignored_sample_groups Character vector. Sample groups to ignore, separated by commas. Default is NULL.
#' @param output_prefix Character. Identifier for the output files. If not provided, no output will be saved. Default is NULL.
#' @param njobs Numeric. Number of parallel jobs to use. Default is the number of available cores.
#' @param verbose Numeric. Level of verbosity for logging messages, from 0 (not verbose) to 5 (very very verbose). Default is 1.
#' @param memory_threshold_mb Numeric. Memory threshold in MB for loading beta files. Default is 500.
#' @param beta_row_names_file Character. Path to a file containing row names for the beta values. If not provided, row names will be read from the beta file. Default is NULL.
#' @param annotate_with_genes Logical. Whether to annotate DMRs with overlapping genes. Default is TRUE.
#' @param bed_provided Logical. Whether the beta file is provided as a BED file. Default is FALSE.
#' @param bed_chrom_col Character. Column name for chromosome in the BED file. Default is "chrom".
#' @param bed_start_col Character. Column name for start position in the BED file. Default is "start".
#' @param bed_end_col Character or NULL. Column name for end position in the BED file. If NULL, start_col + 1 will be used. Default is NULL, which means that the column is assumed missing.
#' @param bed_id_col Character or NULL. Column name for CpG identifier in the BED file. Default is NULL, which means that the column is assumed missing.
#' @param bed_score_col Character or NULL. Column name for score in the BED file. Default is NULL, which means that the column is assumed missing.
#' @param bed_strand_col Character or NULL. Column name for strand in the BED file. Default is NULL, which means that the column is assumed missing.
#' @param .load_debug Logical. If TRUE, enables debug mode for loading beta files. Default is FALSE.

#'
#' @return Data frame of identified DMRs.
#' @export
findDMRsFromSeeds <- function(beta = NULL,
                              dmps = NULL,
                              pheno = NULL,
                              dmps_tsv_id_col = NULL,
                              sample_group_col = "Sample_Group",
                              casecontrol_col = NULL,
                              min_cpg_delta_beta = 0,
                              expansion_step = 500,
                              array = c("450K", "27K", "EPIC", "EPICv2"),
                              genome = c("hg19", "hg38", "mm10", "mm39"),
                              max_pval = 0.05,
                              pval_mode = c("parametric", "empirical"),
                              empirical_strategy = c("auto", "montecarlo", "permutations"),
                              ntries = 200L,
                              mid_p = FALSE,
                              tries_seed = NULL,
                              max_lookup_dist = 10000,
                              min_dmps = 1,
                              min_adj_dmps = 1,
                              min_cpgs = 50,
                              aggfun = c("median", "mean"),
                              ignored_sample_groups = NULL,
                              output_prefix = NULL,
                              njobs = future::availableCores(),
                              memory_threshold_mb = 500,
                              beta_row_names_file = NULL,
                              annotate_with_genes = TRUE,
                              bed_provided = FALSE,
                              bed_chrom_col = "chrom",
                              bed_start_col = "start",
                              bed_end_col = NULL,
                              bed_id_col = NULL,
                              bed_score_col = NULL,
                              bed_strand_col = NULL,
                              .load_debug = FALSE,
                              verbose = 1) {
    pval_mode <- strex::match_arg(pval_mode, ignore_case = TRUE)
    empirical_strategy <- strex::match_arg(empirical_strategy, ignore_case = TRUE)
    # Clean up any zombie processes on exit
    if (Sys.info()[["sysname"]] != "Windows") {
        includes <- "#include <sys/wait.h>"
        code <- "int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};"
        wait <- inline::cfunction(body = code, includes = includes, convention = ".C")
        on.exit(wait(), add = TRUE)
    }
    options(DMRsegal.verbose = verbose)
    # Set up future plan for parallel processing
    if (njobs < 0) {
        njobs <- future::availableCores() + njobs
    }
    old_globals_maxsize <- getOption("future.globals.maxSize", 0)
    if (njobs > 1) {
        if (future::availableCores("multicore") > 1L) {
            .log_info("Using multicore parallelization with ", njobs, " workers", level = 1)
            future::plan(future::multicore, workers = njobs)
        } else {
            .log_info("Using multisession parallelization with ", njobs, " workers", level = 1)
            future::plan(future::multisession, workers = njobs)
        }
        on.exit(future::plan(future::sequential), add = TRUE)
        globals_maxsize <- max(max(memory_threshold_mb * 10 * 1024^2, old_globals_maxsize, na.rm = TRUE), 1024^2 * 500) # at least 500 MB
        .log_info("Setting future.globals.maxSize to ", globals_maxsize / 1024^2, " MB", level = 2)
        options(future.globals.maxSize = globals_maxsize)
    } else {
        .log_info("Using sequential processing (njobs=1)", level = 1)
        future::plan(future::sequential)
        options(future.globals.maxSize = Inf)
    }
    on.exit(options(future.globals.maxSize = old_globals_maxsize), add = TRUE)

    if (is.null(dmps) || is.null(pheno)) {
        stop("dmps and pheno parameters are required")
    }

    .log_step("Preparing inputs...")
    .log_step("Reading DMP tsv..", level = 2)
    if (is.character(dmps) && length(dmps) == 1) {
        dmps_tsv <- try(read.table(
            dmps,
            header = TRUE,
            sep = "\t",
            check.names = FALSE,
            quote = "",
            comment.char = "",
            row.names = NULL
        ))
    } else if (is.data.frame(dmps)) {
        dmps_tsv <- dmps
    } else {
        stop("dmps must be either a file path or a data frame")
    }
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

    if (is.null(dmps_tsv_id_col)) {
        dmps_tsv_id_col <- "_DMP_ROW_NAMES_"
        if (!is.null(rownames(dmps_tsv))) {
            dmps_tsv[, dmps_tsv_id_col] <- rownames(dmps_tsv)
        } else {
            dmps_tsv_id_col <- colnames(dmps_tsv)[1]
        }
        rownames(dmps_tsv) <- NULL
    }
    if (is.numeric(dmps_tsv_id_col)) {
        dmps_tsv_id_col <- colnames(dmps_tsv)[dmps_tsv_id_col]
    }
    if (!dmps_tsv_id_col %in% colnames(dmps_tsv)) {
        stop(
            "DMP id column '", dmps_tsv_id_col,
            "' does not reside in the DMPs file columns: ",
            paste(colnames(dmps_tsv), collapse = ",")
        )
    }

    if (inherits(beta, "BetaHandler")) {
        beta_handler <- beta
        sorted_locs <- beta_handler$getGenomicLocs()
    } else {
        sorted_locs <- NULL
        if (is.character(beta) && length(beta) == 1 && file.exists(beta)) {
            beta_file_ext <- tools::file_ext(beta)
            if (beta_file_ext == "bed" || bed_provided) {
                # Make sure that the dmp tsv_id col has entries that are chr:pos format
                bed_provided <- TRUE
                dmp_ids <- dmps[, dmps_tsv_id_col]
                if (!all(grepl("^(chr)?[0-9XYM]+:[0-9]+$", dmp_ids))) {
                    stop("When providing a bed file as beta input, the DMP IDs in the dmps file/dataframe (using dmps_tsv_id_col: ", dmps_tsv_id_col, ") must be in 'chr:pos' format (e.g., chr1:123456).")
                }
                .log_step("Converting bed beta file to tabix-indexed beta file...")
                ret <- processMethylationBedData(
                    bed_file = beta,
                    pheno = pheno,
                    chrom_col = bed_chrom_col,
                    start_col = bed_start_col,
                    end_col = bed_end_col,
                    id_col = bed_id_col,
                    score_col = bed_score_col,
                    strand_col = bed_strand_col
                )
                beta <- ret$tabix_file
                sorted_locs <- readRDS(ret$locations_file)
                .log_success("Conversion to tabix-indexed beta file completed.")
                # Converting dmps ids to sorted_locs indices, based on their position
                converted_dmp_ids <- c()
                for (dmp_id in dmp_ids) {
                    parts <- unlist(strsplit(dmp_id, ":"))
                    chrom <- as.integer(factor(parts[1], levels = CHROMOSOMES))
                    pos <- as.numeric(parts[2])
                    loc_index <- which(sorted_locs$chr == chrom & sorted_locs$start == pos)

                    if (length(loc_index) == 1) {
                        converted_dmp_ids <- c(converted_dmp_ids, loc_index)
                    } else {
                        .log_warn("DMP ID '", dmp_id, "' could not be found in the bed beta file. It will be skipped.")
                        converted_dmp_ids <- c(converted_dmp_ids, NA)
                    }
                }
                dmps_tsv[, dmps_tsv_id_col] <- converted_dmp_ids
                dmps_tsv <- dmps_tsv[!is.na(dmps_tsv[, dmps_tsv_id_col]), ]
                if (nrow(dmps_tsv) == 0) {
                    stop("No valid DMPs found after converting IDs to bed positions.")
                }
            }
        }
        beta_handler <- getBetaHandler(
            beta = beta,
            beta_row_names_file = beta_row_names_file,
            njobs = njobs,
            memory_threshold_mb = memory_threshold_mb,
            sorted_locs = sorted_locs,
            array = array,
            genome = genome
        )
    }

    if (!is.function(aggfun)) {
        aggfun_choice <- strex::match_arg(aggfun, ignore_case = TRUE)
        aggfun <- switch(aggfun_choice,
            median = stats::median,
            mean = mean
        )
    }
    stopifnot(!is.null(max_pval))
    stopifnot(!is.null(min_dmps))
    stopifnot(!is.null(min_cpgs))
    stopifnot(!is.null(expansion_step))
    stopifnot(!is.null(min_cpg_delta_beta))
    stopifnot(!is.null(max_lookup_dist))

    array <- strex::match_arg(array, ignore_case = TRUE)
    genome <- strex::match_arg(genome, ignore_case = TRUE)
    if (is.character(dmps) && length(dmps) == 1) {
        stopifnot(file.exists(dmps))
    }
    stopifnot(sample_group_col %in% colnames(pheno))
    if (is.null(casecontrol_col)) {
        casecontrol_col <- "_CASECONTROL_INFERRED_"
        pheno[, casecontrol_col] <- ifelse(pheno[, sample_group_col] == levels(as.factor(pheno[, sample_group_col]))[1], 0, 1)
    }
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

    .log_step("Reading beta file characteristics..", level = 2)

    beta_row_names <- beta_handler$getBetaRowNames()
    beta_col_names <- beta_handler$getBetaColNames()


    samples_selection_mask <- !(pheno[, sample_group_col] %in% ignored_sample_groups)
    beta_col_names <- beta_col_names[samples_selection_mask]
    pheno <- pheno[beta_col_names, ]
    .log_info("Samples to process: ", length(beta_col_names), level = 1)

    pheno <- pheno[beta_col_names, ]
    sample_groups <- factor(pheno[beta_col_names, sample_group_col])
    group_inds <- split(seq_along(sample_groups), sample_groups)
    case_mask <- pheno[beta_col_names, casecontrol_col] == 1

    .log_step("Reordering DMPs by genomic location...", level = 2)


    if (!all(dmps_tsv[, dmps_tsv_id_col] %in% beta_row_names)) {
        if (!any(dmps_tsv[, dmps_tsv_id_col] %in% beta_row_names)) {
            # incorrect dmps_tsv_id_col, figure out the correct one
            dmps_tsv_id_col_found <- NULL
            for (col in colnames(dmps_tsv)) {
                if (all(dmps_tsv[, col] %in% beta_row_names)) {
                    dmps_tsv_id_col_found <- col
                    break
                }
            }
            if (is.null(dmps_tsv_id_col_found)) {
                stop("None of the IDs in the DMPs dmps_tsv_id_col match the beta file row names.")
            } else {
                .log_warn("The provided dmps_tsv_id_col '", dmps_tsv_id_col, "' is incorrect. The correct column is '", dmps_tsv_id_col_found, "', switching to that one.")
                dmps_tsv_id_col <- dmps_tsv_id_col_found
            }
        } else {
            missing_in_beta <- dmps_tsv[, dmps_tsv_id_col][!(dmps_tsv[, dmps_tsv_id_col] %in% beta_row_names)]
            .log_warn(length(missing_in_beta), " DMPs are not present in the beta file or the array annotation (maybe due to liftOver?). DMPs: ", paste(missing_in_beta, collapse = ","))
            .log_warn("Ignoring them..")
            dmps_tsv <- dmps_tsv[dmps_tsv[, dmps_tsv_id_col] %in% beta_row_names, , drop = FALSE]
        }
    }
    if (is.null(sorted_locs)) {
        sorted_locs <- getSortedGenomicLocs(array = array, genome = genome)
    }
    if (!bed_provided) {
        sorted_locs <- sorted_locs[beta_row_names, ]
        dmps_tsv <- dmps_tsv[orderByLoc(dmps_tsv[, dmps_tsv_id_col], genomic_locs = sorted_locs), , drop = FALSE]
    } else {
        dmps_tsv <- dmps_tsv[order(dmps_tsv[, dmps_tsv_id_col]), , drop = FALSE]
    }

    # Filter DMPs not present in array annotation first (prevents NA logical indices later)
    if (!bed_provided) {
        dmps <- unique(dmps_tsv[, dmps_tsv_id_col])
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
        dmps <- dmps[orderByLoc(dmps, genome = genome, genomic_locs = sorted_locs)]
    } else {
        dmps <- dmps_tsv[, dmps_tsv_id_col]
    }
    .log_step("Validating beta file sorting by position...", level = 2)


    .log_step("Subsetting beta matrix for DMPs...", level = 2)
    dmps_locs <- sorted_locs[dmps, , drop = FALSE]
    dmps_beta <- beta_handler$getBeta(row_names = dmps, col_names = beta_col_names)
    if (!is.null(output_prefix)) {
        dmps_beta_output_file <- paste0(output_prefix, "dmps_beta.tsv.gz")
    }

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

    if (verbose > 1 && .load_debug && file.exists(file.path("debug", "01_dmrs_from_connected_dmps.tsv"))) {
        .log_info("Loading debug DMRs from file...", level = 2)
        dmrs <- read.table(
            file.path("debug", "01_dmrs_from_connected_dmps.tsv"),
            header = TRUE,
            sep = "\t",
            check.names = FALSE,
            quote = "",
            comment.char = "",
            row.names = NULL
        )
    } else {
        # Use progressr for cross-platform progress reporting
        .log_info("Processing ", length(chromosomes), " chromosomes...", level = 2)

        if (verbose > 0) {
            p_con <- progressr::progressor(steps = length(chromosomes))
        }
        # Split by chromosome for parallel processing
        dmps_list <- split(dmps, dmps_locs$chr)
        dmps_list <- dmps_list[chromosomes]
        dmps_locs_list <- split(dmps_locs, dmps_locs$chr)
        dmps_locs_list <- dmps_locs_list[chromosomes]

        if (!is.matrix(dmps_beta)) dmps_beta <- as.matrix(dmps_beta)
        storage.mode(dmps_beta) <- "double"
        dmps_beta_list <- lapply(split(dmps_beta, dmps_locs$chr), matrix, ncol = ncol(dmps_beta))
        dmps_beta_list <- dmps_beta_list[chromosomes]
        ret <- future.apply::future_mapply(
            chromosomes,
            dmps_list,
            dmps_beta_list,
            dmps_locs_list,
            SIMPLIFY = FALSE,
            future.seed = TRUE,
            future.stdout = NA,
            future.globals = c(
                "group_inds",
                "case_mask",
                "max_lookup_dist",
                "max_pval",
                "aggfun",
                "pval_mode",
                "empirical_strategy",
                "ntries",
                "mid_p",
                "tries_seed",
                "min_dmps",
                ".testConnectivityBatch",
                ".log_step",
                ".log_success",
                ".log_info"
            ),
            FUN = function(chr, cdmps, cdmps_beta, cdmps_locs) {
                op <- options(warn = 2)$warn
                dmr_list <- vector("list", length = 128L)
                dmr_n <- 0L
                connection_corr_pval <- NA
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
                    ntries = ntries,
                    mid_p = mid_p,
                    tries_seed = if (is.null(tries_seed)) NULL else as.integer(tries_seed)
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
                start_dmp_ind <- 1
                for (bp_ind in seq_len(length(breakpoints))) {
                    end_dmp_ind <- breakpoints[bp_ind]
                    if (end_dmp_ind - start_dmp_ind + 1 < min_dmps) {
                        .log_info("Skipping DMR from DMP ", start_dmp_ind, " to DMP ", end_dmp_ind, " (id: ", chr, ":", cdmps_locs[start_dmp_ind, "start"], "-", cdmps_locs[end_dmp_ind, "start"], ") due to insufficient number of connected DMPs (", end_dmp_ind - start_dmp_ind + 1, " < ", min_dmps, ").", level = 4)
                        start_dmp_ind <- breakpoints[bp_ind] + 1
                        next
                    }
                    .log_step("Registering ", bp_ind, "/", (length(breakpoints) - 1), " DMR from DMP ", start_dmp_ind, " to DMP ", end_dmp_ind, " (id: ", chr, ":", cdmps_locs[start_dmp_ind, "start"], "-", cdmps_locs[end_dmp_ind, "start"], ")", level = 3)
                    dmr_dmps_inds <- seq.int(start_dmp_ind, end_dmp_ind)
                    if (end_dmp_ind == start_dmp_ind) {
                        connection_corr_pval <- NA
                    } else {
                        connection_corr_pval <- aggfun(corr_ret$pval[dmr_dmps_inds[-length(dmr_dmps_inds)]], na.rm = TRUE)
                    }
                    if (end_dmp_ind < nrow(cdmps_beta)) {
                        stop_reason <- corr_ret$reason[[end_dmp_ind]]
                    } else {
                        stop_reason <- "end of chromosome"
                    }
                    new_dmr <- list(
                        chr = chr,
                        start_dmp = cdmps[[start_dmp_ind]],
                        end_dmp = cdmps[[end_dmp_ind]],
                        start_dmp_pos = cdmps_locs[start_dmp_ind, "start"],
                        end_dmp_pos = cdmps_locs[end_dmp_ind, "start"],
                        dmps_num = length(dmr_dmps_inds),
                        connection_corr_pval = connection_corr_pval,
                        stop_connection_reason = stop_reason,
                        dmps = paste(cdmps[dmr_dmps_inds], collapse = ",")
                    )
                    dmr_n <- dmr_n + 1L
                    if (dmr_n > length(dmr_list)) length(dmr_list) <- length(dmr_list) * 2L
                    dmr_list[[dmr_n]] <- new_dmr
                    start_dmp_ind <- breakpoints[bp_ind] + 1
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
    }
    cases_num <- dmrs$cases_num
    controls_num <- dmrs$controls_num
    if (anyNA(c(cases_num, controls_num))) {
        .log_warn("NAs introduced while coercing cases_num / controls_num to numeric; replacing NAs with 1 to avoid division errors.")
        cases_num[is.na(cases_num)] <- 1
        controls_num[is.na(controls_num)] <- 1
    }

    .log_success("Initial DMRs formed: ", nrow(dmrs), level = 1)
    .log_info("Summary:\n\t", paste(capture.output(summary(dmrs)), collapse = "\n\t"), level = 1)
    .log_step("Expanding DMRs on neighborhood CpGs..", level = 1)
    # Set up progress tracking for DMR expansion
    n_dmrs <- nrow(dmrs)
    if (verbose > 1 && .load_debug && file.exists(file.path("debug", "connectivity_array.rds"))) {
        .log_info("Loading debug connectivity array from file...", level = 2)
        connectivity_array <- readRDS(file.path("debug", "connectivity_array.rds"))
    } else {
        .log_step("Building connectivity array..", level = 1)
        connectivity_array <- .buildConnectivityArray(
            beta_handler = beta_handler,
            group_inds = group_inds,
            casecontrol = case_mask,
            sorted_locs = sorted_locs,
            max_lookup_dist = max_lookup_dist,
            max_pval = max_pval,
            min_delta_beta = min_cpg_delta_beta,
            aggfun = aggfun,
            pval_mode = pval_mode,
            empirical_strategy = empirical_strategy,
            ntries = ntries,
            mid_p = mid_p,
            tries_seed = if (is.null(tries_seed)) NULL else as.integer(tries_seed),
            njobs = njobs,
            bed_provided = bed_provided
        )
        if (verbose > 0) {
            p_ext <- progressr::progressor(steps = n_dmrs)
        }
        .log_success("Connectivity array built.", level = 2)
    }
    .log_info("Number of connected CpGs found: ", sum(connectivity_array$connected), level = 2)
    if (verbose > 1) {
        dir.create("debug", showWarnings = FALSE)
        saveRDS(connectivity_array, file = file.path("debug", "connectivity_array.rds"))
    }
    stopifnot(nrow(connectivity_array) != nrow(sorted_locs))
    .log_step("Expanding ", n_dmrs, " DMRs using up to ", njobs, " parallel jobs...", level = 2)
    ret <- list()
    for (chr in unique(dmrs$chr)) {
        chr_mask <- sorted_locs$chr == chr
        chr_dmrs <- dmrs[dmrs$chr == chr, ]
        chr_locs <- sorted_locs[chr_mask, , drop = FALSE]
        chr_array <- connectivity_array[chr_mask, , drop = FALSE]

        chr_ret <- future.apply::future_apply(
            X = chr_dmrs,
            MARGIN = 1,
            simplify = FALSE,
            future.seed = TRUE,
            future.scheduling = ceiling(n_dmrs / njobs),
            FUN = function(dmr) {
                op <- options(warn = 2)$warn
                x <- .expandDMRs(
                    dmr = dmr,
                    chr_array = chr_array,
                    expansion_step = expansion_step,
                    min_cpgs = min_cpgs,
                    min_cpg_delta_beta = min_cpg_delta_beta,
                    chr_locs = chr_locs,
                    bed_provided = bed_provided
                )
                options(warn = op)
                if (verbose > 0 && exists("p_ext")) p_ext()
                x
            }
        )
        ret <- c(ret, chr_ret)
    }
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
    extended_dmrs$connection_corr_pval <- as.numeric(extended_dmrs$connection_corr_pval)

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
    if (!bed_provided) {
        all_locs_inds <- rownames(sorted_locs)
        names(all_locs_inds) <- all_locs_inds
        all_locs_inds[seq_along(all_locs_inds)] <- seq_along(all_locs_inds)
        extended_dmrs$start_cpg_ind <- as.numeric(all_locs_inds[extended_dmrs$start_cpg])
        extended_dmrs$end_cpg_ind <- as.numeric(all_locs_inds[extended_dmrs$end_cpg])
    } else {
        extended_dmrs$start_cpg_ind <- extended_dmrs$start_cpg
        extended_dmrs$end_cpg_ind <- extended_dmrs$end_cpg
    }


    .log_success("Post-processing complete.", level = 2)


    extended_dmrs_ranges <- GenomicRanges::makeGRangesFromDataFrame(
        extended_dmrs,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqnames.field = "chr",
        seqinfo = GenomeInfoDb::Seqinfo(genome = genome),
        na.rm = TRUE
    )

    .log_step("Merging overlapping extended DMRs..", level = 2)
    extended_dmrs_ranges_reduced <- GenomicRanges::reduce(extended_dmrs_ranges, ignore.strand = TRUE)
    hits <- GenomicRanges::findOverlaps(extended_dmrs_ranges_reduced, extended_dmrs_ranges, ignore.strand = TRUE)
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    orig_mcols <- as.data.frame(GenomicRanges::mcols(extended_dmrs_ranges))
    mcol_names <- colnames(orig_mcols)
    agg_df <- data.frame(matrix(NA, nrow = length(extended_dmrs_ranges_reduced), ncol = length(mcol_names)))
    colnames(agg_df) <- mcol_names

    for (i in seq_len(length(extended_dmrs_ranges_reduced))) {
        inds <- sh[qh == i]
        cols_vals <- orig_mcols[inds, ]
        agg_df[i, "start_cpg_ind"] <- min(cols_vals$start_cpg_ind)
        agg_df[i, "end_cpg_ind"] <- max(cols_vals$end_cpg_ind)
        agg_df[i, "sup_cpgs_num"] <- agg_df[i, "end_cpg_ind"] - agg_df[i, "start_cpg_ind"] + 1
        agg_df[i, "start_dmp"] <- cols_vals$start_dmp[[1]]
        agg_df[i, "end_dmp"] <- cols_vals$end_dmp[[length(inds)]]
        agg_df[i, "start_dmp_pos"] <- cols_vals$start_dmp_pos[[1]]
        agg_df[i, "end_dmp_pos"] <- cols_vals$end_dmp_pos[[length(inds)]]
        agg_dmps <- unique(unlist(strsplit(cols_vals$dmps, ",")))
        agg_df[i, "dmps"] <- paste(agg_dmps, collapse = ",")
        agg_df[i, "dmps_num"] <- length(agg_dmps)
        agg_df[i, "connection_corr_pval"] <- aggfun(cols_vals$connection_corr_pval, na.rm = TRUE)
        agg_df[i, "stop_connection_reason"] <- paste(cols_vals$stop_connection_reason, collapse = ",")
        agg_df[i, "start_cpg"] <- cols_vals$start_cpg[[1]]
        agg_df[i, "end_cpg"] <- cols_vals$end_cpg[[length(inds)]]
        agg_df[i, "downstream_cpg_expansion_stop_reason"] <- paste(cols_vals$downstream_cpg_expansion_stop_reason, collapse = ",")
        agg_df[i, "upstream_cpg_expansion_stop_reason"] <- paste(cols_vals$upstream_cpg_expansion_stop_reason, collapse = ",")
        agg_df[i, "merged_dmrs_num"] <- length(inds)
    }

    agg_df[, "id"] <- paste0(seqnames(extended_dmrs_ranges_reduced), ":", agg_df$start_cpg_ind, "-", agg_df$end_cpg_ind)

    GenomicRanges::mcols(extended_dmrs_ranges_reduced) <- agg_df
    extended_dmrs_ranges <- extended_dmrs_ranges_reduced

    .log_step("Finding GC content of DMRs..", level = 1)
    # increase end by 1 to include last base in getDMRSequences, in case there is a C, belonging to a CpG, at the end
    GenomicRanges::end(extended_dmrs_ranges) <- GenomicRanges::end(extended_dmrs_ranges) + 1
    sequences <- getDMRSequences(extended_dmrs_ranges, genome)
    extended_dmrs <- as.data.frame(extended_dmrs_ranges)
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

    .log_step("Adding DMR delta-beta information..", level = 2)
    all_selected_cpgs <- unique(unlist(lapply(seq_len(nrow(extended_dmrs)), function(i) {
        seq(extended_dmrs$start_cpg_ind[i], extended_dmrs$end_cpg_ind[i])
    })))
    if (!bed_provided) {
        all_selected_cpgs <- rownames(sorted_locs)[all_selected_cpgs]
    }
    all_selected_cpgs_beta <- beta_handler$getBeta(row_names = all_selected_cpgs, col_names = beta_col_names)

    beta_stats <- .calculateBetaStats(
        beta = all_selected_cpgs_beta,
        beta_col_names = beta_col_names,
        pheno = pheno,
        casecontrol_col = casecontrol_col,
        aggfun = aggfun
    )

    beta_stats <- as.data.frame(beta_stats)
    rownames(beta_stats) <- all_selected_cpgs
    dmrs_with_beta_stats <- list()
    for (dmr_ind in seq_len(nrow(extended_dmrs))) {
        dmr <- extended_dmrs[dmr_ind, ]
        dmr_dmps <- strsplit(dmr$dmps, ",")[[1]]
        dmr$cases_beta <- aggfun(abs(beta_stats[dmr_dmps, "cases_beta"])) * sign(sum(sign(beta_stats[dmr_dmps, "cases_beta"])))
        dmr$controls_beta <- aggfun(abs(beta_stats[dmr_dmps, "controls_beta"])) * sign(sum(sign(beta_stats[dmr_dmps, "controls_beta"])))
        dmr$delta_beta <- dmr$cases_beta - dmr$controls_beta
        dmr$cases_beta_sd <- aggfun(beta_stats[dmr_dmps, "cases_beta_sd"], na.rm = TRUE)
        dmr$controls_beta_sd <- aggfun(beta_stats[dmr_dmps, "controls_beta_sd"], na.rm = TRUE)
        dmr$cases_beta_min <- min(beta_stats[dmr_dmps, "cases_beta"], na.rm = TRUE)
        dmr$cases_beta_max <- max(beta_stats[dmr_dmps, "cases_beta"], na.rm = TRUE)
        dmr$controls_beta_min <- min(beta_stats[dmr_dmps, "controls_beta"], na.rm = TRUE)
        dmr$controls_beta_max <- max(beta_stats[dmr_dmps, "controls_beta"], na.rm = TRUE)
        if (!bed_provided) {
            dmr_cpgs <- rownames(sorted_locs)[seq.int(dmr$start_cpg_ind, dmr$end_cpg_ind)]
        } else {
            dmr_cpgs <- seq.int(dmr$start_cpg_ind, dmr$end_cpg_ind)
        }
        dmr$cpgs_cases_beta <- aggfun(abs(beta_stats[dmr_cpgs, "cases_beta"])) * sign(sum(sign(beta_stats[dmr_cpgs, "cases_beta"])))
        dmr$cpgs_controls_beta <- aggfun(abs(beta_stats[dmr_cpgs, "controls_beta"])) * sign(sum(sign(beta_stats[dmr_cpgs, "controls_beta"])))
        dmr$cpgs_delta_beta <- dmr$cpgs_cases_beta - dmr$cpgs_controls_beta
        dmr$cpgs_cases_beta_sd <- aggfun(beta_stats[dmr_cpgs, "cases_beta_sd"], na.rm = TRUE)
        dmr$cpgs_controls_beta_sd <- aggfun(beta_stats[dmr_cpgs, "controls_beta_sd"], na.rm = TRUE)
        dmr$cpgs_cases_beta_min <- min(beta_stats[dmr_cpgs, "cases_beta"], na.rm = TRUE)
        dmr$cpgs_cases_beta_max <- max(beta_stats[dmr_cpgs, "cases_beta"], na.rm = TRUE)
        dmr$cpgs_controls_beta_min <- min(beta_stats[dmr_cpgs, "controls_beta"], na.rm = TRUE)
        dmr$cpgs_controls_beta_max <- max(beta_stats[dmr_cpgs, "controls_beta"], na.rm = TRUE)

        dmrs_with_beta_stats[[dmr_ind]] <- dmr
    }
    dmrs <- do.call(rbind, dmrs_with_beta_stats)
    .log_success("DMR delta-beta information added.", level = 2)

    dmrs_granges <- GenomicRanges::makeGRangesFromDataFrame(
        dmrs,
        keep.extra.columns = TRUE,
        seqnames.field = "chr",
        seqinfo = GenomeInfoDb::Seqinfo(genome = genome),
        na.rm = TRUE
    )

    if (annotate_with_genes) {
        .log_step("Annotating DMRs with gene information...", level = 2)
        dmrs_granges <- annotateDMRsWithGenes(dmrs_granges, genome = genome)
        .log_success("DMR annotation completed.", level = 2)
    }
    dmrs <- as.data.frame(dmrs_granges)
    colnames(dmrs)[colnames(dmrs) == "seqnames"] <- "chr"
    if (bed_provided) {
        dmrs$start_dmp <-paste0(sorted_locs[dmrs$start_dmp, "chr"], ":", sorted_locs[dmrs$start_dmp, "start"])
        dmrs$end_dmp <-paste0(sorted_locs[dmrs$end_dmp, "chr"], ":", sorted_locs[dmrs$end_dmp, "start"])
        dmrs$start_cpg <- paste0(sorted_locs[dmrs$start_cpg, "chr"], ":", dmrs$start_cpg)
        dmrs$end_cpg <- paste0(sorted_locs[dmrs$end_cpg, "chr"], ":", dmrs$end_cpg)
        dmrs$dmps <- sapply(dmrs$dmps, function(dmp_ids) {
            dmp_ids_split <- unlist(strsplit(dmp_ids, ","))
            dmp_ids_annotated <- paste0(sorted_locs[dmp_ids_split, "chr"], ":", sorted_locs[dmp_ids_split, "start"])
            paste(dmp_ids_annotated, collapse = ",")
        })

    }
    .log_info("Summary of final DMRs:\n\t", paste(capture.output(summary(dmrs)), collapse = "\n\t"), level = 1)
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
