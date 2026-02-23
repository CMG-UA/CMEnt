#' Find Differentially Methylated Regions (DMRs) from Pre-computed seeds
#'
#' @name findDMRsFromSeeds
#' @description This function identifies Differentially Methylated Regions (DMRs) from pre-computed
#' Differentially Methylated Positions (seeds) using a correlation-based approach. It expands
#' significant seeds into regions, considering both statistical significance and biological
#' relevance of methylation changes.
#'
#' @section Important Note on Input Data:
#' Do not apply heavy filtering to your seeds prior to using this function, particularly based on
#' beta values or effect sizes. The function works by expanding regions around significant seeds
#' and connecting nearby CpGs into larger regions. Filtering out seeds with smaller effect sizes
#' may remove important CpGs that could serve as "bridges" to connect more significant seeds into
#' larger, biologically meaningful DMRs. For optimal results, include all statistically
#' significant seeds (e.g., adjusted p-value < 0.05) and let the function handle region expansion
#' and filtering internally using the min_cpg_delta_beta parameter if needed.
#'
#' @param beta Character. Path to the beta value file, or a tabix file, or a beta matrix, or a BetaHandler object, or a bed file, or a BSseq object. If a bed file is provided, it must at least contain bed_chrom_col and bed_chrom_start, followed by samples names in the given pheno, with corresponging beta values, and it will be converted to a tabix-indexed beta file internall, with the locations separately saved and queried as a bigmemmory::big.matrix. If a BSseq object is provided, genomic locations and methylation values will be extracted using bsseq methods.
#' @param seeds Character. Path to the seeds TSV file or the seeds dataframe, in a format like the one produced by dmpFinder.
#' @param pheno Data frame. Phenotype data.
#' @param seeds_id_col Character. Column name or index for Seed identifiers in the seeds TSV file. Default is 1.
#' @param sample_group_col Character. Column name for sample group information in the phenotype data. Default is NULL.
#' @param casecontrol_col Boolean Column in pheno for case (TRUE/1) / control (FALSE/0) status . If NULL, controls will be assumed to be the first level of sample_group_col. Default is NULL.
#' @param covariates Character vector of column names in pheno to adjust for (e.g. "age", "sex"). When provided, correlations are computed on residuals after regressing M-values on these covariates within each group
#' @param min_cpg_delta_beta Numeric. Minimum delta beta value for CpGs. Default is 0.1.
#' @param adaptive_min_cpg_delta_beta Logical. Whether to adaptively increase min_cpg_delta_beta from seed-level delta-beta distribution (never below min_cpg_delta_beta). Default is TRUE.
#' @param expansion_step Numeric. Index-specific step size for expanding DMRs. Increasing it means higher memory usage and faster computation. Default is 500.
#' @param array Character. Type of array used (e.g., "450K", "EPIC", "EPICv2", "27K"). Ignored if using mm10 genome.
#' @param genome Character. Genome version (e.g., "hg19", "hg38", "mm10"). Default is "hg19".
#' @param max_pval Numeric. Maximum p-value to assume seeds correlation is significant. Default is 0.05.
#' @param pval_mode Character. "auto" (default) selects between t-based correlation p-values and empirical p-values per sample group using data diagnostics. You can also force "parametric" for t-based correlation p-values or "empirical" for permutation-based p-values.
#' @param empirical_strategy Character. When pval_mode = "empirical": "auto" (default) uses Monte Carlo for groups with <6 samples and permutations for groups with >=6 samples; "montecarlo" always uses Monte Carlo; "permutations" always uses permutations.
#' @param ntries Integer. Number of permutations when pval_mode = "empirical". Default is 0 (disabled).
#' @param max_lookup_dist Numeric. Maximum distance to look up for adjacent seeds belonging to the same DMR during Stage 1. Default is 10000 (10 kb).
#' @param connectivity_window_bp Numeric. Stage 2 connectivity is computed only in windows centered on seed-derived Stage 1 DMR neighborhoods, with this total window width in bp. Set <=0 for genome-wide connectivity. Default is -1 for microarrays and 10000 (10 kb) for NGS datasets.
#' @param max_bridge_seeds_gaps Integer. Maximum number of consecutive failed seed-to-seed edges to bridge during Stage 1 when both flanking edges are connected and failures are p-value driven. Set to 0 to disable. Default is 1.
#' @param max_bridge_extension_gaps Integer. Maximum gap size to consider during Stage 2 extension. Default is 1 (i.e., at most 1 consecutive failing CpG to bridge).
#' @param min_seeds Numeric. Minimum number of connected seeds in a DMR. Minimum is 2. Default is 2.
#' @param min_adj_seeds Numeric. Minimum number of seeds, adjusted by array CpG density, in a DMR after extension. Minimum is 2. Default is 2.
#' @param min_cpgs Numeric. Minimum number of CpGs in a DMR after extension, including the seeds. Minimum is 2. Default is 50.
#' @param aggfun Function or character. Aggregation function to use when calculating delta beta values and p-values of DMRs. Can be "median", "mean", or a function (e.g., median, mean). Default is "median".
#' @param ignored_sample_groups Character vector. Sample groups to ignore during connection and expansion, separated by commas. Can also be "case" or "control". Default is NULL.
#' @param output_prefix Character. Identifier for the output files. If not provided, no output will be saved. Default is NULL.
#' @param njobs Numeric. Number of parallel jobs to use. Default is the number of available cores.
#' @param verbose Numeric. Level of verbosity for logging messages, from 0 (not verbose) to 5 (very very verbose). Default is 1.
#' @param beta_row_names_file Character. Path to a file containing row names for the beta values. If not provided, row names will be read from the beta file. Default is NULL.
#' @param annotate_with_genes Logical. Whether to annotate DMRs with overlapping genes. Default is TRUE.
#' @param bed_provided Logical. Whether the beta file is provided as a BED file. Default is FALSE.
#' @param bed_chrom_col Character. Column name for chromosome in the BED file. Default is "chrom".
#' @param bed_start_col Character. Column name for start position in the BED file. Default is "start".
#' @param .load_debug Logical. If TRUE, enables debug mode for loading beta files. Default is FALSE.
#'
#' @return A GRanges object containing identified DMRs with metadata columns:
#' \itemize{
#'   \item cpgs_num: Number of CpGs in the region
#'   \item seeds_num: Number of seeds in the region
#'   \item delta_beta: Aggregated methylation difference (using the specified aggregation function)
#'   \item delta_beta_min: Minimum methylation difference
#'   \item delta_beta_max: Maximum methylation difference
#'   \item delta_beta_sd: Standard deviation of methylation differences
#'   \item cases_beta: Mean beta value in cases
#'   \item controls_beta: Mean beta value in controls
#'   \item start_seed: ID of the first Seed in the region
#'   \item end_seed: ID of the last Seed in the region
#'   \item start_seed_pos: Genomic position of the first Seed in the region
#'   \item end_seed_pos: Genomic position of the last Seed in the region
#'   \item seeds: Comma-separated list of Seed IDs in the region
#'   \item connection_corr_pval: Aggregated correlation p-value of connected seeds
#'   \item stop_connection_reason: Reasons for stopping the connection of seeds
#'   \item start_cpg: ID of the starting CpG after expansion
#'   \item end_cpg: ID of the ending CpG after expansion
#'   \item upstream_cpg_expansion: Number of CpGs expanded upstream
#'   \item upstream_cpg_expansion_stop_reason: Reason for stopping upstream expansion
#'   \item downstream_cpg_expansion: Number of CpGs expanded downstream
#'   \item downstream_cpg_expansion_stop_reason: Reason for stopping downstream expansion
#'   \item merged_dmrs_num: Number of merged DMRs
#' }
#'
#' @examples
#' # Load example data
#' beta <- loadExampleInputData("beta")
#' dmps <- loadExampleInputData("dmps")
#' pheno <- loadExampleInputData("pheno")
#' # Find DMRs
#' dmrs <- findDMRsFromSeeds(
#'     beta = beta,
#'     seeds = seeds,
#'     pheno = pheno
#' )
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

#' @keywords internal
#' @noRd
.mergeGenomicWindows <- function(windows) {
    if (is.null(windows) || nrow(windows) == 0L) {
        return(data.frame(chr = character(0), start = numeric(0), end = numeric(0)))
    }
    windows <- as.data.frame(windows)
    windows$chr <- as.character(windows$chr)
    windows$start <- as.integer(round(as.numeric(windows$start)))
    windows$end <- as.integer(round(as.numeric(windows$end)))
    windows <- windows[!is.na(windows$chr) & !is.na(windows$start) & !is.na(windows$end), , drop = FALSE]
    if (nrow(windows) == 0L) {
        return(data.frame(chr = character(0), start = numeric(0), end = numeric(0)))
    }
    bad <- windows$end < windows$start
    if (any(bad)) {
        swap <- windows$start[bad]
        windows$start[bad] <- windows$end[bad]
        windows$end[bad] <- swap
    }
    gr <- GenomicRanges::GRanges(
        seqnames = windows$chr,
        ranges = IRanges::IRanges(start = windows$start, end = windows$end)
    )
    gr <- GenomicRanges::reduce(gr, ignore.strand = TRUE)
    data.frame(
        chr = as.character(GenomicRanges::seqnames(gr)),
        start = GenomicRanges::start(gr),
        end = GenomicRanges::end(gr)
    )
}

#' @keywords internal
#' @noRd
.buildConnectivityWindowsFromDMRs <- function(dmrs, connectivity_window_bp) {
    if (is.null(dmrs) || nrow(dmrs) == 0L || !is.finite(connectivity_window_bp) || connectivity_window_bp <= 0) {
        return(data.frame(chr = character(0), start = numeric(0), end = numeric(0)))
    }
    flank <- as.numeric(connectivity_window_bp) / 2
    windows <- data.frame(
        chr = as.character(dmrs$chr),
        start = pmax(1, as.numeric(dmrs$start_seed_pos) - flank),
        end = as.numeric(dmrs$end_seed_pos) + flank
    )
    .mergeGenomicWindows(windows)
}

#' @keywords internal
#' @noRd
.bridgeConnectivityGaps <- function(corr_ret, max_bridge_seeds_gaps = 0L, bridgeable_reasons = c("pval>max_pval", "na pval")) {
    if (is.null(corr_ret) || nrow(corr_ret) == 0L || is.na(max_bridge_seeds_gaps) || max_bridge_seeds_gaps <= 0) {
        attr(corr_ret, "bridged_edges") <- 0L
        return(corr_ret)
    }
    connected <- as.logical(corr_ret$connected)
    reasons <- as.character(corr_ret$reason)
    n <- length(connected)
    if (n < 3L) {
        attr(corr_ret, "bridged_edges") <- 0L
        return(corr_ret)
    }
    runs <- rle(connected)
    run_ends <- cumsum(runs$lengths)
    run_starts <- run_ends - runs$lengths + 1L
    bridged <- 0L
    for (r in which(!runs$values)) {
        gap_start <- run_starts[r]
        gap_end <- run_ends[r]
        gap_len <- runs$lengths[r]
        if (gap_len > max_bridge_seeds_gaps || gap_start == 1L || gap_end == n) {
            next
        }
        if (!connected[gap_start - 1L] || !connected[gap_end + 1L]) {
            next
        }
        run_reasons <- reasons[gap_start:gap_end]
        reason_ok <- all(vapply(run_reasons, function(rr) {
            if (is.na(rr) || rr == "") {
                return(FALSE)
            }
            parts <- trimws(strsplit(rr, ";", fixed = TRUE)[[1]])
            all(parts %in% bridgeable_reasons)
        }, logical(1)))
        if (!reason_ok) {
            next
        }
        connected[gap_start:gap_end] <- TRUE
        reasons[gap_start:gap_end] <- "bridged"
        bridged <- bridged + gap_len
    }
    corr_ret$connected <- connected
    corr_ret$reason <- reasons
    attr(corr_ret, "bridged_edges") <- bridged
    corr_ret
}

#' @keywords internal
#' @noRd
.buildConnectivityArray <- function(
    beta_handler, pheno, group_inds, col_names = NULL, max_pval = 0.05,
    min_delta_beta = 0, covariates = NULL, max_lookup_dist = 1000,
    chunk_size = getOption("DMRsegal.chunk_size", 1000), entanglement = "strong",
    aggfun = median, empirical_strategy = "auto",
    pval_mode = "auto", ntries = 500, mid_p = TRUE, njobs = 1,
    connectivity_windows = NULL, connectivity_array = NULL, gap = 1L, splits = NULL
) {
    beta_locs <- beta_handler$getBetaLocs()
    n_sites <- nrow(beta_locs)
    if (n_sites < 2L) {
        ret <- data.frame(
            connected = rep(FALSE, n_sites),
            pval = rep(NA_real_, n_sites),
            reason = rep("end-of-input", n_sites),
            stringsAsFactors = FALSE
        )
        if (identical(entanglement, "strong")) {
            ret$first_failing_group <- rep("", n_sites)
        } else {
            ret$failing_groups <- rep("", n_sites)
        }
        if (min_delta_beta > 0) {
            ret$delta_beta <- rep(NA_real_, n_sites)
        }
        return(list(connectivity_array = ret, splits = NULL))
    }
    beta_chr <- as.character(beta_locs[, "chr"])
    chr_ends <- as.integer(vapply(split(seq_len(n_sites), beta_chr), max, integer(1)))
    window_mode <- !is.null(connectivity_windows) && nrow(connectivity_windows) > 0L
    default_reason <- if (window_mode) "outside_connectivity_window" else ""

    .make_template <- function(nrows, reason_default) {
        ret <- data.frame(
            connected = rep(FALSE, nrows),
            pval = rep(NA_real_, nrows),
            reason = rep(reason_default, nrows),
            stringsAsFactors = FALSE
        )
        if (identical(entanglement, "strong")) {
            ret$first_failing_group <- rep("", nrows)
        } else {
            ret$failing_groups <- rep("", nrows)
        }
        if (min_delta_beta > 0) {
            ret$delta_beta <- rep(NA_real_, nrows)
        }
        ret
    }

    .get_chunk_beta <- function(site_start, site_end, mask = NULL) {
        if (!is.null(mask)) {
            site_indices <- if (!bigmemory::is.big.matrix(beta_locs)) {
                rownames(beta_locs)[site_start:site_end][mask]
            } else {
                (site_start:site_end)[mask]
            }
        } else {
            site_indices <- if (!bigmemory::is.big.matrix(beta_locs)) {
                rownames(beta_locs)[site_start:site_end]
            } else {
                site_start:site_end
            }
        }
        beta_handler$getBeta(
            row_names = site_indices,
            col_names = col_names,
            check_mem = TRUE
        )
    }

    .chunk_pair_ranges <- function(pair_ranges_df) {
        if (is.null(pair_ranges_df) || nrow(pair_ranges_df) == 0L) {
            return(matrix(numeric(0), ncol = 2))
        }
        out <- vector("list", nrow(pair_ranges_df) * 2L)
        out_n <- 0L
        for (i in seq_len(nrow(pair_ranges_df))) {
            ps <- as.integer(pair_ranges_df$start_pair[i])
            pe <- as.integer(pair_ranges_df$end_pair[i])
            if (pe < ps) {
                next
            }
            for (chunk_ps in seq(ps, pe, by = chunk_size)) {
                chunk_pe <- min(chunk_ps + chunk_size - 1L, pe)
                out_n <- out_n + 1L
                out[[out_n]] <- c(chunk_ps, chunk_pe)
            }
        }
        do.call(rbind, out[seq_len(out_n)])
    }

    .build_all_pair_ranges <- function() {
        out <- vector("list", length(unique(beta_chr)))
        out_n <- 0L
        for (chr in unique(beta_chr)) {
            chr_inds <- which(beta_chr == chr)
            if (length(chr_inds) < 2L) {
                next
            }
            out_n <- out_n + 1L
            out[[out_n]] <- data.frame(
                start_pair = min(chr_inds),
                end_pair = max(chr_inds) - 1L
            )
        }
        if (out_n == 0L) {
            return(data.frame(start_pair = integer(0), end_pair = integer(0)))
        }
        do.call(rbind, out[seq_len(out_n)])
    }

    .build_window_pair_ranges <- function() {
        wins <- .mergeGenomicWindows(connectivity_windows)
        if (nrow(wins) == 0L) {
            return(data.frame(start_pair = integer(0), end_pair = integer(0)))
        }
        pair_ranges <- vector("list", length(unique(beta_chr)))
        pair_n <- 0L
        beta_start <- as.integer(round(as.numeric(beta_locs[, "start"])))
        for (chr in intersect(unique(wins$chr), unique(beta_chr))) {
            chr_inds <- which(beta_chr == chr)
            if (length(chr_inds) < 2L) {
                next
            }
            chr_wins <- wins[wins$chr == chr, , drop = FALSE]
            if (nrow(chr_wins) == 0L) {
                next
            }
            win_ir <- IRanges::IRanges(
                start = as.integer(round(as.numeric(chr_wins$start))),
                end = as.integer(round(as.numeric(chr_wins$end)))
            )
            site_ir <- IRanges::IRanges(
                start = beta_start[chr_inds],
                width = 1L
            )
            ov <- IRanges::findOverlaps(win_ir, site_ir)
            if (length(ov) == 0L) {
                next
            }
            qh <- S4Vectors::queryHits(ov)
            sh <- S4Vectors::subjectHits(ov)
            min_rel <- as.integer(tapply(sh, qh, min))
            max_rel <- as.integer(tapply(sh, qh, max))
            pair_start <- chr_inds[min_rel]
            pair_end <- chr_inds[max_rel] - 1L
            keep <- !is.na(pair_start) & !is.na(pair_end) & pair_end >= pair_start
            if (!any(keep)) {
                next
            }
            pair_n <- pair_n + 1L
            pair_ranges[[pair_n]] <- data.frame(
                start_pair = pair_start[keep],
                end_pair = pair_end[keep]
            )
        }
        if (pair_n == 0L) {
            return(data.frame(start_pair = integer(0), end_pair = integer(0)))
        }
        pair_ranges <- do.call(rbind, pair_ranges[seq_len(pair_n)])
        merged <- IRanges::reduce(IRanges::IRanges(
            start = pair_ranges$start_pair,
            end = pair_ranges$end_pair
        ))
        data.frame(
            start_pair = IRanges::start(merged),
            end_pair = IRanges::end(merged)
        )
    }
    if (is.null(splits)) {
        pair_ranges <- if (window_mode) .build_window_pair_ranges() else .build_all_pair_ranges()
        splits <- .chunk_pair_ranges(pair_ranges)

    }
    if (nrow(splits) == 0L) {
        connectivity_array <- .make_template(n_sites, default_reason)
        connectivity_array[chr_ends, "connected"] <- FALSE
        connectivity_array[chr_ends, "reason"] <- "end-of-input"
        return(list(connectivity_array = connectivity_array, splits = splits))
    }

    if (is.null(connectivity_array)) {
        connectivity_array <- .make_template(n_sites, default_reason)
        connectivity_array[chr_ends, "connected"] <- FALSE
        connectivity_array[chr_ends, "reason"] <- "end-of-input"
        revisited_mask <- NULL
    } else {
        # If connectivity array is provided, we assume it has already been filled for all sites up to the chromosome ends.
        # Instead, we re-assess the disconnected sites on the edges of the connected regions,
        # comparing i with i + 2 instead of i with i + 1, to see if we can connect them by bridging the gap of one site.
        # The connectivity array is then updated with the bridged connections.
        connected <- connectivity_array[, "connected"]
        runs <- rle(connected)
        run_ends <- cumsum(runs$lengths)
        run_starts <- run_ends - runs$lengths + 1L
        m <- connectivity_array[run_ends, "reason"] != "end-of-input"
        run_ends <- run_ends[m]
        run_starts <- run_starts[m]
        values <- runs$values[m]
        lengths <- runs$lengths[m]
        run_mask <- lengths > gap & values == 1
        run_ends <- run_ends[run_mask]
        run_starts <- run_starts[run_mask]
        revisited_mask <- rep(FALSE, n_sites)
        # The following indices will be re-checked
        checked_inds <- sort(c(run_starts - gap, run_starts, run_ends, run_ends + gap))
        checked_inds <- checked_inds[checked_inds > 0 & checked_inds <= n_sites]
        revisited_mask[checked_inds] <- TRUE
        splits <- splits[apply(splits, 1, function(r) any(revisited_mask[r[1]:r[2]])), , drop = FALSE]
        .log_info(
            "Re-assessing connectivity for ", sum(revisited_mask), " sites at the edges of existing connected regions to see if we can bridge small gaps.",
            level = 3
        )
    }
    .log_info(
        "Connectivity computation mode: ",
        if (window_mode) "windowed" else "genome-wide",
        "; chunks to evaluate: ", nrow(splits), ".",
        level = 3
    )
    verbose <- getOption("DMRsegal.verbose", 1)
    if (verbose > 0) {
        p_ext <- progressr::progressor(steps = nrow(splits), message = "Computing connectivity array...")
    }
    gc()
    fun <- function(split_ind) {
        split <- splits[split_ind, ]
        pair_start <- as.integer(split[1])
        pair_end <- as.integer(split[2])
        site_start <- pair_start
        site_end <- pair_end + 1L
        locs <- beta_locs[site_start:site_end, , drop = FALSE]
        if (!is.null(revisited_mask)) {
            mask_seg <- revisited_mask[site_start:site_end]
            locs <- locs[mask_seg, , drop = FALSE]
            chunk_beta <- .get_chunk_beta(site_start, site_end, mask_seg)
        } else {
            chunk_beta <- .get_chunk_beta(site_start, site_end)
        }
        if (nrow(chunk_beta) != nrow(locs)) {
            stop("Mismatch between number of rows in beta chunk (", nrow(chunk_beta), ") and number of locations (", nrow(locs), ").")
        }
        x <- .testConnectivityBatch(
            sites_beta = chunk_beta,
            group_inds = group_inds,
            pheno = pheno,
            covariates = covariates,
            max_pval = max_pval,
            min_delta_beta = min_delta_beta,
            max_lookup_dist = max_lookup_dist,
            sites_locs = locs,
            entanglement = entanglement,
            aggfun = aggfun,
            pval_mode = pval_mode,
            empirical_strategy = empirical_strategy,
            ntries = ntries,
            mid_p = mid_p
        )

        if (!is.null(revisited_mask)) {
            # global site indices in this chunk are site_start:site_end
            sel_sites_rel <- which(mask_seg)  # indices relative to site_start
            if (length(sel_sites_rel) >= 2L) {
                # pair indices correspond to the first site of each consecutive selected pair
                recomputed_pairs <- site_start + sel_sites_rel[-length(sel_sites_rel)] - 1L
            } else {
                recomputed_pairs <- integer(0)
            }
            # attach to result so outer loop can map back exactly
            attr(x, "recomputed_pairs") <- recomputed_pairs
        }
        if (verbose > 0) {
            p_ext()
        }
        list(pair_start = pair_start, pair_end = pair_end, result = x)
    }
    if (njobs == 1) {
        ret <- lapply(seq_len(nrow(splits)), fun)
    } else {
        .setupParallel()
        ret <- future.apply::future_lapply(
            X = seq_len(nrow(splits)),
            future.seed = TRUE,
            future.stdout = NA,
            future.globals = c(
                "splits",
                "beta_locs",
                "group_inds",
                "pheno",
                "covariates",
                "max_pval",
                "min_delta_beta",
                "max_lookup_dist",
                "entanglement",
                "aggfun",
                "pval_mode",
                "empirical_strategy",
                "ntries",
                "revisited_mask",
                "mid_p",
                "verbose",
                "p_ext",
                ".get_chunk_beta",
                ".testConnectivityBatch"
            ),
            FUN = fun
        )
        .finalizeParallel()
    }
    for (item in ret) {
        idx <- item$pair_start:item$pair_end
        x <- item$result
        missing_cols <- setdiff(names(x), names(connectivity_array))
        if (length(missing_cols) > 0L) {
            for (col in missing_cols) {
                if (is.character(x[[col]])) {
                    connectivity_array[[col]] <- rep("", n_sites)
                } else {
                    connectivity_array[[col]] <- rep(NA, n_sites)
                }
            }
        }
        if (is.null(revisited_mask)) {
            # First pass: full overwrite
            for (col in names(x)) {
                connectivity_array[idx, col] <- x[[col]]
            }
        } else {
            # Map result rows to the exact global pair indices returned by the worker.
            recomputed_pairs <- attr(x, "recomputed_pairs")
            if (is.null(recomputed_pairs)) recomputed_pairs <- integer(0)
            masked_idx <- recomputed_pairs
            if (length(masked_idx) != nrow(x)) {
                stop("Mismatch between recomputed pair indices (", length(masked_idx), ") and recomputed connectivity rows (", nrow(x), ").")
            }
            update_m <- x[["connected"]]
            if (any(update_m) && length(masked_idx) >= 1L) {
                for (col in names(x)) {
                    connectivity_array[masked_idx[update_m], col] <- x[[col]][update_m]
                }
            }
        }
    }
    if (!is.null(connectivity_array) && gap > 0L) {
        connected <- connectivity_array[, "connected"]
        runs <- rle(connected)

        run_lengths <- runs$lengths
        run_values  <- runs$values

        run_starts <- cumsum(c(1L, head(run_lengths, -1)))
        run_ends   <- cumsum(run_lengths)

        fill_indices <- integer(0)

        for (i in seq_along(run_values)) {

            if (!run_values[i] && run_lengths[i] <= gap) {

                left_run  <- i - 1L
                right_run <- i + 1L

                if (left_run >= 1L && right_run <= length(run_values) &&
                        run_values[left_run] && run_values[right_run]) {

                    fill_indices <- c(
                        fill_indices,
                        run_starts[i]:run_ends[i]
                    )
                }
            }
        }

        if (length(fill_indices) > 0) {
            connectivity_array[fill_indices, "connected"] <- TRUE
            connectivity_array[fill_indices, "reason"]    <- "bridged"
        }
    } else if (!is.null(revisited_mask)) {
        connectivity_array[, "contingently_connected"] <- connectivity_array[, "connected"]
    }

    list(connectivity_array = connectivity_array, splits = splits)
}

#' @keywords internal
#' @noRd
.expandDMR <- function(dmr,
                       chr_array,
                       chr_locs,
                       min_cpg_delta_beta = 0,
                       min_cpgs = 3,
                       expansion_step = 500,
                       chr_start_base = 0) {
    .log_step("Expanding DMR..", level = 4)
    dmr_start <- dmr["start_seed"]
    dmr_end <- dmr["end_seed"]
    if (!is.data.frame(chr_locs) && bigmemory::is.sub.big.matrix(chr_locs)) {
        dmr_start_ind <- as.integer(dmr_start) - chr_start_base
        dmr_end_ind <- as.integer(dmr_end) - chr_start_base
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
        # [----+---0] where - is connected, + is not connected, 0 stands for the current DMR start (ustream_exp),
        # then fail_start_idx is the first + from the right, 4 in this case.
        # That means that the 4th from the right failed to connect to the 3rd from the right.
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
        # [0----+---] where - is connected, + is not connected, 0 stands for the current DMR start (dstream_exp),
        # then fail_start_idx is the first + from the left, 5 in this case.
        # That means that the 5th from the left failed to connect to the 6th from the left.
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
                .log_info("DMR  too short (", ccpgs, " CpGs). Expanding to reach min_cpgs=", min_cpgs, ".", level = 4)
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
            .log_info("Number of CpGs in expanded DMR after first iteration: ", new_ccpgs, " from ", ccpgs, level = 4)
            if (new_ccpgs < min_cpgs) {
                ustream_stop_reason <- "min-cpgs-not-reached"
                dstream_stop_reason <- "min-cpgs-not-reached"
                .log_info("DMR could not reach min_cpgs=", min_cpgs, " after expansion (", new_ccpgs, "). Stopping expansion.", level = 4)
            }
            t <- 1
        }
        if (!is.null(ustream_stop_reason) && !is.null(dstream_stop_reason)) {
            break
        }
    }
    .log_step("Finalizing expanded DMR.", level = 4)
    if (!is.data.frame(chr_locs) && bigmemory::is.sub.big.matrix(chr_locs)) {
        dmr["start_cpg"] <- ustream_exp + chr_start_base
        dmr["end_cpg"] <- dstream_exp + chr_start_base
        dmr["start_cpg_ind"] <- ustream_exp + chr_start_base
        dmr["end_cpg_ind"] <- dstream_exp + chr_start_base
        dmr["start"] <- chr_locs[ustream_exp, "start"]
        dmr["end"] <- chr_locs[dstream_exp, "start"]
    } else {
        dmr["start_cpg"] <- rownames(chr_locs)[ustream_exp]
        dmr["end_cpg"] <- rownames(chr_locs)[dstream_exp]
        dmr["start_cpg_ind"] <- ustream_exp + chr_start_base
        dmr["end_cpg_ind"] <- dstream_exp + chr_start_base
        dmr["start"] <- chr_locs[dmr["start_cpg"], "start"]
        dmr["end"] <- chr_locs[dmr["end_cpg"], "start"]
    }
    to_cpg_ids <- function(local_inds) {
        if (length(local_inds) == 0) {
            return(character(0))
        }
        if (!is.data.frame(chr_locs) && bigmemory::is.sub.big.matrix(chr_locs)) {
            return(as.character(local_inds + chr_start_base))
        }
        rownames(chr_locs)[local_inds]
    }

    dmr["upstream_cpg_expansion_stop_reason"] <- ustream_stop_reason
    upstream_candidate <- if (ustream_exp <= (dmr_start_ind - 1L)) {
        seq.int(ustream_exp, dmr_start_ind - 1L)
    } else {
        integer(0)
    }
    if (length(upstream_candidate) > 0) {
        bridged_upstream_m <- chr_array[upstream_candidate, "connection_reason"] == "bridged"
        upstream_kept <- upstream_candidate[!bridged_upstream_m]
    } else {
        upstream_kept <- integer(0)
    }
    dmr["upstream_cpgs_inds"] <- paste(upstream_kept + chr_start_base, collapse = ",")
    dmr["upstream_cpgs"] <- paste(to_cpg_ids(upstream_kept), collapse = ",")
    dmr["upstream_cpg_expansion"] <- length(upstream_kept)

    dmr["downstream_cpg_expansion_stop_reason"] <- dstream_stop_reason
    downstream_candidate <- if ((dmr_end_ind + 1L) <= dstream_exp) {
        seq.int(dmr_end_ind + 1L, dstream_exp)
    } else {
        integer(0)
    }
    if (length(downstream_candidate) > 0) {
        bridged_downstream_m <- chr_array[downstream_candidate, "connection_reason"] == "bridged"
        downstream_kept <- downstream_candidate[!bridged_downstream_m]
    } else {
        downstream_kept <- integer(0)
    }
    dmr["downstream_cpgs_inds"] <- paste(downstream_kept + chr_start_base, collapse = ",")
    dmr["downstream_cpgs"] <- paste(to_cpg_ids(downstream_kept), collapse = ",")
    dmr["downstream_cpg_expansion"] <- length(downstream_kept)

    .log_success("Expanded DMR finalized: (start_cpg: ", dmr["start_cpg"], ", end_cpg: ", dmr["end_cpg"], ").", level = 4)
    dmr
}

#' Vectorized connectivity testing for consecutive site pairs, given their beta values
#'
#' @return Data frame with columns: connected, pval, delta_beta, reason, first_failing_group, stop_reason
#' @keywords internal
#' @noRd
.testConnectivityBatch <- function(sites_beta, group_inds, pheno,
                                   max_pval, covariates = NULL,
                                   min_delta_beta = 0,
                                   max_lookup_dist = NULL, sites_locs = NULL,
                                   entanglement = "strong",
                                   aggfun = mean,
                                   pval_mode = c("auto", "parametric", "empirical"),
                                   empirical_strategy = c("auto", "montecarlo", "permutations"),
                                   ntries = 0, mid_p = FALSE) {
    pval_mode <- strex::match_arg(pval_mode, ignore_case = TRUE)
    empirical_strategy <- strex::match_arg(empirical_strategy, ignore_case = TRUE)
    n_sites <- nrow(sites_beta)
    strict_mode <- identical(entanglement, "strong")
    if (n_sites < 2) {
        ret <- data.frame(
            connected = logical(0),
            pval = numeric(0),
            delta_beta = numeric(0),
            reason = character(0)
        )
        if (strict_mode) {
            ret$first_failing_group <- character(0)
        } else {
            ret$failing_groups <- character(0)
        }
        return(ret)
    }

    n_pairs <- n_sites - 1
    n_groups <- length(group_inds)
    if (strict_mode) {
        # null hypothesis: all groups must be significant -> bonferroni correction
        max_pval_corrected <- max_pval / n_groups
        .log_info("Using max p-value of ", max_pval_corrected, "(group multi-testing corrected) for connectivity testing.", level = 4)
    } else {
        # null hypothesis: at least one group must be significant -> independent testing
        max_pval_corrected <- max_pval
        .log_info("Using max p-value of ", max_pval_corrected, " for connectivity testing.", level = 4)
    }
    # Initialize result vectors
    connected <- rep(TRUE, n_pairs)
    pvals <- rep(NA_real_, n_pairs)
    reasons <- rep("", n_pairs)
    failing_groups <- rep("", n_pairs)

    # Check distance condition if provided (vectorized)
    if (!is.null(max_lookup_dist) && !is.null(sites_locs)) {
        dists <- sites_locs[2:n_sites, "start"] - sites_locs[1:(n_sites - 1), "start"]
        exceeded_dist <- dists > max_lookup_dist
        connected[exceeded_dist] <- FALSE
        reasons[exceeded_dist] <- "exceeded max distance"
    } else {
        exceeded_dist <- rep(FALSE, n_pairs)
    }
    nexdist_mask <- !exceeded_dist
    .log_info(sum(exceeded_dist), " out of ", n_pairs, " site pairs exceeded the maximum lookup distance and will be marked as not connected.", level = 4)
    if (!strict_mode && n_groups > 0) {
        per_group_reasons <- matrix("", nrow = n_groups, ncol = n_pairs)
        per_group_reasons[, exceeded_dist] <- "exceeded max distance"
        per_group_p <- matrix(NA_real_, nrow = n_groups, ncol = n_pairs)
        rownames(per_group_reasons) <- names(group_inds)
        rownames(per_group_p) <- names(group_inds)
    }
    g_index <- 0
    # Fully vectorized correlation testing for each group
    for (g in names(group_inds)) {
        .log_step("Processing group '", g, "' (", g_index, "/", n_groups, ") for chunk connectivity testing...", level = 4)
        g_index <- g_index + 1
        idx <- group_inds[[g]]
        if (length(idx) < 3) next

        # Get data for this group - subset columns
        group_beta <- sites_beta[, idx, drop = FALSE]
        group_m <- log(group_beta / (1 - group_beta + 1e-6) + 1e-6) # M-values transformation

        if (!is.null(covariates) && length(covariates) > 0L && !is.null(pheno)) {
            if (!all(covariates %in% colnames(pheno))) {
                stop("Not all covariates are present in pheno.")
            }
            xc <- as.data.frame(pheno[idx, covariates, drop = FALSE])
            xc <- data.frame(`(Intercept)` = 1, xc, check.names = FALSE)
            # convert any string columns to factors
            for (col in colnames(xc)) {
                if (is.character(xc[[col]])) {
                    xc[[col]] <- as.numeric(as.factor(xc[[col]]))
                }
            }
            xc <- as.matrix(xc)
            storage.mode(xc) <- "double"
            group_m <- .remove_confounder_effect(group_m, xc)
        }

        # Extract consecutive pairs matrices
        x_mat_full <- group_m[1:(n_sites - 1), , drop = FALSE] # Sites i
        y_mat_full <- group_m[2:n_sites, , drop = FALSE] # Sites i+1

        # Apply distance mask
        x_mat <- x_mat_full[nexdist_mask, , drop = FALSE]
        y_mat <- y_mat_full[nexdist_mask, , drop = FALSE]

        sn_pairs <- nrow(x_mat)
        if (sn_pairs == 0L) {
            next
        }
        g_reasons <- rep("", sn_pairs)
        g_mask <- rep(TRUE, sn_pairs)
        if (strict_mode) {
            g_mask <- connected[nexdist_mask]
        } else {
            g_mask <- rep(TRUE, sn_pairs)
        }

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

        low_df <- (dfs < 1) & g_mask
        g_reasons[low_df] <- ifelse(g_reasons[low_df] == "", "df<1", g_reasons[low_df])
        g_mask[low_df] <- FALSE

        na_r <- is.na(cors) & g_mask
        g_reasons[na_r] <- ifelse(g_reasons[na_r] == "", "na r", g_reasons[na_r])
        g_mask[na_r] <- FALSE

        ps <- rep(NA_real_, sn_pairs)
        effective_pval_mode <- pval_mode
        if (pval_mode == "auto") {
            auto_diag <- .summarizeCorrelationAssumptions(
                x_mat = x_mat,
                y_mat = y_mat,
                n_valid = n_valid
            )
            effective_pval_mode <- if (auto_diag$use_parametric) "parametric" else "empirical"
            .log_info(
                "Auto-selected pval_mode='", effective_pval_mode, "' for group '", g,
                "' (q10_n_valid=", signif(auto_diag$q10_n_valid, 3),
                ", median|pearson-spearman|=", signif(auto_diag$median_abs_delta_spearman, 3),
                ", median|pearson-winsorized|=", signif(auto_diag$median_abs_delta_winsorized, 3),
                ", median|skew|=", signif(auto_diag$median_abs_skew, 3),
                ", median_excess_kurtosis=", signif(auto_diag$median_excess_kurtosis, 3),
                ", pilot_pairs=", auto_diag$n_pairs_used,
                ", reason=", auto_diag$reason, ").",
                level = 4
            )
        }

        # Precompute parametric p-values as fallback when empirical is not feasible/resolved
        if (effective_pval_mode == "parametric") {
            tstats <- cors * sqrt(dfs / pmax(1e-12, 1 - cors * cors))
            # Compute t-statistics (vectorized)
            na_tstat <- is.na(tstats) & g_mask
            g_reasons[na_tstat] <- "na tstat"
            g_mask[na_tstat] <- FALSE
            ps[g_mask] <- -2 * expm1(pt(abs(tstats[g_mask]), df = dfs[g_mask], log.p = TRUE))
        } else {
            # Empirical p-values via permutations of sample labels within group
            # Only compute for rows that are still connected and have finite cors
            mask <- is.finite(cors) & g_mask
            if (any(mask)) {
                set.seed(getOption("DMRsegal.random_seed", 42))
                counts_ge <- integer(sn_pairs)
                counts_eq <- integer(sn_pairs)
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
                skip_empirical <- FALSE
                if (do_permutations) {
                    min_possible_pval <- 1 / (1 + factorial(m))
                    if (min_possible_pval > max_pval_corrected) {
                        .log_warn(
                            "Cannot compute sufficiently small empirical p-values for group '", g,
                            "' because minimum possible p-value (", min_possible_pval,
                            ") exceeds max_pval_corrected (", max_pval_corrected,
                            "). Marking currently eligible pairs as not connected for this group."
                        )
                        ps[g_mask] <- 1
                        skip_empirical <- TRUE
                    }
                }
                if (!skip_empirical) {
                    .log_info("Computing empirical p-values for group '", g, "' using ", if (do_permutations) "permutations" else "Monte Carlo", " with ", ntries, " tries.", level = 4)
                    if (ntries == 0) {
                        if (do_permutations) {
                            ntries <- min(500L, factorial(m))
                        } else {
                            ntries <- 500
                        }
                    }
                    maxval <- max(y_mat, na.rm = TRUE)
                    minval <- min(y_mat, na.rm = TRUE)
                    abs_cors <- abs(cors)
                    for (b in seq_len(ntries)) {
                        # Permute sample labels (columns) only for y; x remains fixed
                        if (do_permutations) {
                            perm <- sample.int(m, size = m, replace = FALSE)
                            yp <- y_mat[, perm, drop = FALSE]
                        } else {
                            yp <- matrix(stats::runif(n = nrow(y_mat) * m, min = minval, max = maxval), nrow = nrow(y_mat), ncol = m)
                        }
                        yc <- yp - rowMeans(yp, na.rm = TRUE)
                        sxy <- rowSums(x_centered * yc, na.rm = TRUE)
                        sy2 <- rowSums(yc^2, na.rm = TRUE)
                        rperm <- sxy / sqrt(sum_x2 * sy2)
                        comp_mask <- is.finite(rperm)
                        if (any(comp_mask)) {
                            ap <- abs(rperm[comp_mask])
                            ao <- abs_cors[comp_mask]
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
        }


        na_p <- is.na(ps) & g_mask
        g_reasons[na_p] <- "na pval"
        g_mask[na_p] <- FALSE
        exceed_pval <- g_mask & (ps > max_pval_corrected)
        g_reasons[exceed_pval] <- "pval>max_pval"
        g_mask[exceed_pval] <- FALSE

        # Update results back to main vectors
        if (strict_mode) {
            sconnected <- connected[nexdist_mask]
            broad_mask <- nexdist_mask & connected
            prev_p <- pvals[broad_mask]
            next_p <- ps[sconnected]
            pvals[broad_mask] <- ifelse(
                is.na(prev_p),
                next_p,
                pmax(prev_p, next_p)
            )
            reasons[broad_mask] <- g_reasons[sconnected]
            failing_groups[broad_mask] <- ifelse(g_mask[sconnected], "", g)
            connected[broad_mask] <- connected[broad_mask] & g_mask[sconnected]
        } else {
            per_group_p[g_index, nexdist_mask] <- ps
            per_group_reasons[g_index, nexdist_mask] <- g_reasons
        }
        .log_success("Finished processing chunk for group '", g, "'.", level = 4)
    }
    if (!strict_mode) {
        not_failed <- per_group_reasons == ""
        connected <- colSums(per_group_reasons == "") > 0
        pvals[connected] <- as.vector(apply(
            per_group_p[, connected, drop = FALSE], 2, function(v) {
                if (all(is.na(v))) {
                    return(NA_real_)
                }
                max(v, na.rm = TRUE)
            }
        ))
        reasons[!connected] <- apply(
            per_group_reasons[, !connected, drop = FALSE], 2, function(v) paste(v, collapse = ";")
        )
        failing_groups[!connected] <- apply(
            !not_failed[, !connected, drop = FALSE], 2, function(v) paste(names(group_inds)[v], collapse = ";")
        )
    }
    ret <- data.frame(
        connected = connected,
        pval = pvals,
        reason = reasons
    )
    if (strict_mode) {
        ret[, "first_failing_group"] <- failing_groups
    } else {
        ret[, "failing_groups"] <- failing_groups
    }
    # Vectorized delta beta check if needed
    if (min_delta_beta > 0 && length(unique(pheno[, "__casecontrol__"])) > 1) {
        # Extract site2 beta values for all pairs
        site2_beta_mat <- sites_beta[2:n_sites, , drop = FALSE]

        # Vectorized mean computation across case/control
        case_betas <- apply(site2_beta_mat[, pheno[, "__casecontrol__"] == 1, drop = FALSE], 1, aggfun, na.rm = TRUE)
        control_betas <- apply(site2_beta_mat[, pheno[, "__casecontrol__"] == 0, drop = FALSE], 1, aggfun, na.rm = TRUE)
        delta_betas <- case_betas - control_betas

        # Check threshold (vectorized)
        low_delta <- (is.na(delta_betas) | abs(delta_betas) < min_delta_beta) & connected
        connected[low_delta] <- FALSE
        reasons[low_delta] <- "delta_beta<min_delta_beta"
        ret[["delta_beta"]] <- delta_betas
    }
    ret$connected <- connected
    ret$reason <- reasons
    if (strict_mode) {
        ret$first_failing_group <- failing_groups
    } else {
        ret$failing_groups <- failing_groups
    }


    ret
}



#' Find Differentially Methylated Regions (DMRs) from Differentially Methylated Positions (seeds)
#'
#' This function identifies DMRs from a given set of seeds and a beta value file.
#'
#' @param beta Character. Path to the beta value file, or a tabix file, or a beta matrix, or a BetaHandler object, or a bed file. If a bed file is provided, it must at least contain bed_chrom_col and bed_chrom_start, followed by samples names in the given pheno, with corresponging beta values, and it will be converted to a tabix-indexed beta file internall, with the locations separately saved and queried as a bigmemmory::big.matrix.
#' @param seeds Character. Path to the seeds (seeds, etc.) TSV file or the seeds dataframe, in a format like the one produced by dmpFinder.
#' @param pheno Data frame. Phenotype data.
#' @param seeds_id_col Character. Column name or index for Seed identifiers in the seeds TSV file. Default is NULL, which corresponds to the rows names if existing, or the first column if not.
#' @param sample_group_col Character. Column name for sample group information in the phenotype data. Default is NULL.
#' @param casecontrol_col Boolean Column in pheno for case (TRUE/1) / control (FALSE/0) status . If NULL, controls will be assumed to be the first level of sample_group_col. Default is NULL.
#' @param covariates Character vector of column names in pheno to adjust for (e.g. "age", "sex"). When provided, correlations are computed on residuals after regressing M-values on these covariates within each group
#' @param min_cpg_delta_beta Numeric. Minimum delta beta value for CpGs. Default is 0.1.
#' @param adaptive_min_cpg_delta_beta Logical. Whether to adaptively increase min_cpg_delta_beta from seed-level delta-beta distribution (never below min_cpg_delta_beta). Default is TRUE.
#' @param expansion_step Numeric. Index-specific step size for expanding DMRs. Increasing it means higher memory usage and faster computation. Default is 500.
#' @param array Character. Type of array used (e.g., "450K", "EPIC", "EPICv2", "27K"). Ignored if using a mouse genome. Also ignored if the beta file is provided as a beta values BED file. Default is "450K".
#' @param genome Character. Genome version. Default is "hg19".
#' @param max_pval Numeric. Maximum p-value to assume seeds correlation is significant. Default is 0.05.
#' @param entanglement Character. "strong" (default) requires all groups to show significant correlation for connectivity; "weak" requires at least one group to show significant correlation.
#' @param pval_mode Character. "auto" (default) selects between t-based correlation p-values and empirical p-values per sample group using data diagnostics. You can also force "parametric" for t-based correlation p-values or "empirical" for permutation-based p-values.
#' @param empirical_strategy Character. When pval_mode = "empirical": "auto" (default) uses Monte Carlo for groups with <6 samples and permutations for groups with >=6 samples; "montecarlo" always uses Monte Carlo; "permutations" always uses permutations.
#' @param ntries Integer. Number of permutations when pval_mode = "empirical". Default is 0 (disabled).
#' @param max_lookup_dist Numeric. Maximum distance to look up for adjacent seeds belonging to the same DMR during Stage 1. Default is 10000 (10 kb).
#' @param connectivity_window_bp Numeric. Stage 2 connectivity is computed only in windows centered on seed-derived Stage 1 DMR neighborhoods, with this total window width in bp. This value sets a maximum effective size of a DMR after stage 2. Set <=0 for genome-wide connectivity. Default is -1 for microarrays and 10000 (10 kb) for NGS datasets.
#' @param max_bridge_seeds_gaps Integer. Maximum number of consecutive failed seed-to-seed edges to bridge during Stage 1 when both flanking edges are connected and failures are p-value driven. Set to 0 to disable. Default is 1.
#' @param max_bridge_extension_gaps Integer. Maximum gap size to consider during Stage 2 extension. Default is 1 (i.e., at most 1 consecutive failing CpG to bridge).
#' @param min_seeds Numeric. Minimum number of connected seeds in a DMR. Minimum is 2. Default is 2.
#' @param min_adj_seeds Numeric. Minimum number of seeds, adjusted by array CpG density, in a DMR after extension. Minimum is 2. Default is 2.
#' @param min_cpgs Numeric. Minimum number of CpGs in a DMR after extension, including the seeds. Minimum is 2. Default is 50.
#' @param aggfun Function or character. Aggregation function to use when calculating delta beta values and p-values of DMRs. Can be "median", "mean", or a function (e.g., median, mean). Default is "median".
#' @param ignored_sample_groups Character vector. Sample groups to ignore during connection and expansion, separated by commas. Can also be "case" or "control". Default is NULL.
#' @param output_prefix Character. Identifier for the output files. If not provided, no output will be saved. Default is NULL.
#' @param njobs Numeric. Number of parallel jobs to use. Default is the number of available cores.
#' @param beta_row_names_file Character. Path to a file containing row names for the beta values. If not provided, row names will be read from the beta file. Default is NULL.
#' @param annotate_with_genes Logical. Whether to annotate DMRs with overlapping genes. Default is TRUE.
#' @param rank_dmrs Logical. Whether to rank DMRs based on significance and effect size. Default is TRUE.
#' @param bed_provided Logical. Whether the beta file is provided as a BED file. Default is FALSE. In case the input has a .bed extension, this will be set to TRUE automatically.
#' @param bed_chrom_col Character. Column name for chromosome in the BED file. Default is "chrom".
#' @param bed_start_col Character. Column name for start position in the BED file. Default is "start".
#' @param chr_levels Character vector. Custom chromosome levels to use. Default is NULL.
#' @param verbose Numeric. Level of verbosity for logging messages, from 0 (not verbose) to 5 (very very verbose). Default is retrieved from option "DMRsegal.verbose".
#' @param .load_debug Logical. If TRUE, enables debug mode for loading beta files. Default is FALSE.
#' @param chunk_size Numeric. Number of CpGs to process in each chunk. Default is retrieved from option "DMRsegal.chunk_size".
#'
#' @return Data frame of identified DMRs.
#' @export
findDMRsFromSeeds <- function(
    beta,
    seeds,
    pheno,
    seeds_id_col = NULL,
    sample_group_col = "Sample_Group",
    casecontrol_col = NULL,
    covariates = NULL,
    min_cpg_delta_beta = 0.1,
    adaptive_min_cpg_delta_beta = TRUE,
    expansion_step = 500,
    array = c("450K", "27K", "EPIC", "EPICv2", "NULL"),
    genome = "hg19",
    max_pval = 0.05,
    entanglement = c("strong", "weak"),
    pval_mode = c("auto", "parametric", "empirical"),
    empirical_strategy = c("auto", "montecarlo", "permutations"),
    ntries = 200L,
    mid_p = FALSE,
    max_lookup_dist = 10000,
    connectivity_window_bp = "auto",
    max_bridge_seeds_gaps = 1L,
    max_bridge_extension_gaps = 1L,
    min_seeds = 2,
    min_adj_seeds = 2,
    min_cpgs = 50,
    aggfun = c("median", "mean"),
    ignored_sample_groups = NULL,
    output_prefix = NULL,
    njobs = getOption("DMRsegal.njobs", min(8, future::availableCores() - 1)),
    chunk_size = getOption("DMRsegal.chunk_size", 10000),
    beta_row_names_file = NULL,
    annotate_with_genes = TRUE,
    rank_dmrs = TRUE,
    bed_provided = FALSE,
    bed_chrom_col = "chrom",
    bed_start_col = "start",
    chr_levels = NULL,
    verbose = getOption("DMRsegal.verbose", 1),
    .load_debug = FALSE
) {
    pval_mode <- strex::match_arg(pval_mode, ignore_case = TRUE)
    empirical_strategy <- strex::match_arg(empirical_strategy, ignore_case = TRUE)

    # Clean up any zombie processes on exit
    if (Sys.info()[["sysname"]] != "Windows") {
        includes <- "#include <sys/wait.h>"
        code <- "int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};"
        wait <- inline::cfunction(body = code, includes = includes, convention = ".C")
        withr::defer(wait())
    }
    options(DMRsegal.verbose = verbose)
    options(cli.num_colors = cli::num_ansi_colors())
    options(future.globals.maxSize = Inf)

    # Set up future plan for parallel processing

    options("DMRsegal.njobs" = njobs)


    .log_step("Preparing inputs...")
    .log_step("Reading Seed tsv..", level = 2)
    if (is.character(seeds) && length(seeds) == 1) {
        seeds_tsv <- try(as.data.frame(read.table(
            seeds,
            header = TRUE,
            sep = "\t",
            check.names = FALSE,
            quote = "",
            comment.char = "",
            row.names = NULL
        )))
    } else if (is.data.frame(seeds)) {
        seeds_tsv <- as.data.frame(seeds)
    } else {
        stop("seeds must be either a file path or a data frame")
    }
    if (inherits(seeds_tsv, "try-error")) {
        .log_warn("Provided seeds file is empty or does not exist. Not proceeding.")
        if (!is.null(output_prefix)) {
            for (f in c(".methylation.tsv.gz", ".tsv.gz")) {
                gzfile <- gzfile(paste0(output_prefix, f), "w", compression = 2)
                close(gzfile)
            }
        }
        return(NULL)
    }

    if (is.null(seeds_id_col)) {
        seeds_id_col <- "_SEED_ROW_NAMES_"
        if (!is.null(rownames(seeds_tsv))) {
            seeds_tsv[, seeds_id_col] <- rownames(seeds_tsv)
        } else {
            seeds_id_col <- colnames(seeds_tsv)[1]
        }
        rownames(seeds_tsv) <- NULL
    }
    if (is.numeric(seeds_id_col)) {
        seeds_id_col <- colnames(seeds_tsv)[seeds_id_col]
    }
    if (!seeds_id_col %in% colnames(seeds_tsv)) {
        stop(
            "Seed id column '", seeds_id_col,
            "' does not reside in the seeds file columns: ",
            paste(colnames(seeds_tsv), collapse = ",")
        )
    }
    if (!is.null(covariates)) {
        missing_covars <- covariates[!covariates %in% colnames(pheno)]
        if (length(missing_covars) > 0) {
            stop("The following covariates are not present in pheno: ", paste(missing_covars, collapse = ", "))
        }
    }
    array_based <- !is.null(array)
    if (array_based) {
        array <- strex::match_arg(array, ignore_case = TRUE)
    }
    if (is.null(chr_levels)) {
        chr_levels <- GenomeInfoDb::getChromInfoFromUCSC(genome)[, 1]
    }
    if (inherits(beta, "BetaHandler")) {
        beta_handler <- beta
    } else {
        beta_locs <- NULL
        if (is.character(beta) && length(beta) == 1 && file.exists(beta)) {
            beta_file_ext <- tools::file_ext(beta)
            if (beta_file_ext == "bed" || bed_provided) {
                # Make sure that the seeds_id_col has entries that are chr:pos format
                bed_provided <- TRUE
                seed_ids <- seeds_tsv[, seeds_id_col]
                if (!all(grepl("^(chr)?[0-9XYM]+:[0-9]+$", seed_ids))) {
                    stop("When providing a bed file as beta input, the seed IDs in the seeds file/dataframe (using seeds_id_col: ", seeds_id_col, ") must be in 'chr:pos' format (e.g., chr1:123456).")
                }
                .log_step("Converting bed beta file to tabix-indexed beta file...")
                ret <- readCustomMethylationBedData(
                    bed_file = beta,
                    pheno = pheno,
                    genome = genome,
                    chrom_col = bed_chrom_col,
                    start_col = bed_start_col,
                    chr_levels = chr_levels
                )
                beta <- ret$tabix_file
                beta_locs <- readRDS(ret$locations_file)
                try(beta_locs <- bigmemory::attach.big.matrix(beta_locs), silent = TRUE)
                .log_success("Conversion to tabix-indexed beta file completed.")
                # Converting seeds ids to beta_locs indices, based on their position
                .log_step("Converting Seed IDs to bed positions...", level = 2)
                converted_seed_ids <- c()
                chroms_and_pos <- strsplit(seed_ids, ":")
                chroms <- sapply(chroms_and_pos, function(x) x[1])
                positions <- as.numeric(sapply(chroms_and_pos, function(x) x[2]))
                # Convert chroms to UCSC style if needed
                if (!all(chroms %in% chr_levels)) {
                    if (all(paste0("chr", chroms) %in% chr_levels)) {
                        chroms <- paste0("chr", chroms)
                    } else if (all(sub("^chr", "", chroms) %in% chr_levels)) {
                        chroms <- sub("^chr", "", chroms)
                    } else {
                        stop("Seed IDs chromosomes in the seeds file/dataframe (using seeds_id_col: ", seeds_id_col, ") do not match the chromosome style of the provided bed beta file.")
                    }
                }
                chroms <- as.integer(factor(chroms, levels = chr_levels))
                seeds_tsv <- seeds_tsv[order(chroms, positions), ]
                converted_seed_ids <- match(
                    paste0(chroms, ":", positions),
                    rownames(beta_locs)
                )
                if (all(is.na(converted_seed_ids))) {
                    stop("None of the Seed IDs could be converted to bed positions based on the provided bed beta file.")
                }
                if (any(is.na(converted_seed_ids))) {
                    .log_warn(sum(is.na(converted_seed_ids)), "out of", length(converted_seed_ids), " Seed IDs could not be matched to the provided bed beta file. They will be ignored.", level = 2)
                    seeds_tsv <- seeds_tsv[!is.na(converted_seed_ids), ]
                    converted_seed_ids <- converted_seed_ids[!is.na(converted_seed_ids)]
                }
                .log_success("Seed IDs successfully converted to bed indices.", level = 2)
                seeds_tsv[, seeds_id_col] <- converted_seed_ids
                seeds_tsv <- seeds_tsv[!is.na(seeds_tsv[, seeds_id_col]), ]
                if (nrow(seeds_tsv) == 0) {
                    stop("No valid seeds found after converting IDs to bed indices.")
                }
            }
        }
        beta_handler <- getBetaHandler(
            beta = beta,
            beta_row_names_file = beta_row_names_file,
            njobs = njobs,
            sorted_locs = beta_locs,
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
    stopifnot(!is.null(min_seeds))
    stopifnot(!is.null(min_cpgs))
    if (min_seeds < 2 && min_cpgs < 2) {
        stop("min_seeds or min_cpgs must be at least 2, to define a DMR")
    }
    stopifnot(!is.null(min_adj_seeds))
    if (min_adj_seeds < 2) {
        stop("min_adj_seeds must be at least 2, to define a DMR")
    }
    stopifnot(!is.null(expansion_step))
    stopifnot(!is.null(min_cpg_delta_beta))
    stopifnot(!is.null(max_lookup_dist))
    stopifnot(!is.null(entanglement))
    if (!is.logical(adaptive_min_cpg_delta_beta) || length(adaptive_min_cpg_delta_beta) != 1 || is.na(adaptive_min_cpg_delta_beta)) {
        stop("adaptive_min_cpg_delta_beta must be TRUE or FALSE.")
    }
    if (connectivity_window_bp == "auto") {
        connectivity_window_bp <- if (array_based) -1 else 10000
    }
    if (!is.numeric(connectivity_window_bp) || length(connectivity_window_bp) != 1 || is.na(connectivity_window_bp)) {
        stop("connectivity_window_bp must be a numeric scalar.")
    }
    if (!is.numeric(max_bridge_seeds_gaps) || length(max_bridge_seeds_gaps) != 1 || is.na(max_bridge_seeds_gaps) || max_bridge_seeds_gaps < 0) {
        stop("max_bridge_seeds_gaps must be a non-negative integer.")
    }
    max_bridge_seeds_gaps <- as.integer(max_bridge_seeds_gaps)
    if (!is.numeric(max_bridge_extension_gaps) || length(max_bridge_extension_gaps) != 1 || is.na(max_bridge_extension_gaps) || max_bridge_extension_gaps < 0) {
        stop("max_bridge_extension_gaps must be a non-negative integer.")
    }
    max_bridge_extension_gaps <- as.integer(max_bridge_extension_gaps)

    entanglement <- strex::match_arg(entanglement, ignore_case = TRUE)

    if (is.character(seeds) && length(seeds) == 1) {
        stopifnot(file.exists(seeds))
    }
    stopifnot(sample_group_col %in% colnames(pheno))
    if (is.null(casecontrol_col)) {
        pheno[, "__casecontrol__"] <- ifelse(pheno[, sample_group_col] == levels(as.factor(pheno[, sample_group_col]))[1], 0, 1)
    } else {
        pheno[, "__casecontrol__"] <- pheno[, casecontrol_col]
    }
    if (is.null(ignored_sample_groups)) {
        ignored_sample_groups <- c()
    } else {
        ignored_sample_groups <- trimws(unlist(strsplit(ignored_sample_groups, ",")))
        ignored_sample_groups <- ignored_sample_groups[nzchar(ignored_sample_groups)]
    }
    if (!is.null(output_prefix)) {
        output_dir <- dirname(output_prefix)
        dir.create(output_dir, showWarnings = FALSE)
        output_prefix <- paste0(output_prefix, ".")
    } else {
        output_dir <- NULL
        output_prefix <- NULL
    }


    # Check if seeds_tsv has any rows
    if (nrow(seeds_tsv) == 0) {
        .log_warn("Provided seeds file has no data rows. Not proceeding.")
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
    missing_pheno_samples <- setdiff(beta_col_names, rownames(pheno))
    if (length(missing_pheno_samples) > 0) {
        stop(
            "The following beta samples are missing from pheno row names: ",
            paste(head(missing_pheno_samples, 10), collapse = ","),
            if (length(missing_pheno_samples) > 10) " ..." else ""
        )
    }
    pheno_all <- pheno[beta_col_names, , drop = FALSE]
    beta_locs <- beta_handler$getBetaLocs()


    samples_selection_mask <- !(pheno_all[, sample_group_col] %in% ignored_sample_groups)
    if ("case" %in% ignored_sample_groups) {
        samples_selection_mask <- samples_selection_mask & (pheno_all[, "__casecontrol__"] != 1)
    }
    if ("control" %in% ignored_sample_groups) {
        samples_selection_mask <- samples_selection_mask & (pheno_all[, "__casecontrol__"] != 0)
    }
    beta_col_names_detection <- beta_col_names[samples_selection_mask]
    if (length(beta_col_names_detection) < 2) {
        stop("At least two samples are required after applying ignored_sample_groups.")
    }
    pheno_detection <- pheno_all[beta_col_names_detection, , drop = FALSE]
    .log_info("Samples to process during DMR detection: ", length(beta_col_names_detection), level = 1)

    sample_groups <- factor(pheno_detection[, sample_group_col])
    group_inds <- split(seq_along(sample_groups), sample_groups)

    .log_step("Reordering seeds by genomic location...", level = 2)


    if (!all(seeds_tsv[, seeds_id_col] %in% beta_row_names)) {
        if (!any(seeds_tsv[, seeds_id_col] %in% beta_row_names)) {
            # incorrect seeds_id_col, figure out the correct one
            seeds_id_col_found <- NULL
            for (col in colnames(seeds_tsv)) {
                if (all(seeds_tsv[, col] %in% beta_row_names)) {
                    seeds_id_col_found <- col
                    break
                }
            }
            if (is.null(seeds_id_col_found)) {
                stop("None of the IDs in the seeds seeds_id_col match the beta file row names.")
            } else {
                .log_warn("The provided seeds_id_col '", seeds_id_col, "' is incorrect. The correct column is '", seeds_id_col_found, "', switching to that one.")
                seeds_id_col <- seeds_id_col_found
            }
        } else {
            missing_in_beta <- seeds_tsv[, seeds_id_col][!(seeds_tsv[, seeds_id_col] %in% beta_row_names)]
            .log_warn(length(missing_in_beta), " seeds are not present in the beta file or the array annotation (maybe due to liftOver?). seeds: ", paste(missing_in_beta, collapse = ","))
            .log_warn("Ignoring them..")
            seeds_tsv <- seeds_tsv[seeds_tsv[, seeds_id_col] %in% beta_row_names, , drop = FALSE]
        }
    }
    if (!bed_provided) {
        seeds_tsv <- seeds_tsv[orderByLoc(seeds_tsv[, seeds_id_col], genomic_locs = beta_handler$getBetaLocs()), , drop = FALSE]
    } else {
        seeds_tsv <- seeds_tsv[order(seeds_tsv[, seeds_id_col]), , drop = FALSE]
    }

    # Filter seeds not present in array annotation first (prevents NA logical indices later)
    if (!bed_provided) {
        seeds <- unique(seeds_tsv[, seeds_id_col])
        missing_in_annotation <- setdiff(seeds, rownames(beta_handler$getGenomicLocs()))
        if (length(missing_in_annotation) > 0) {
            .log_warn(
                "Dropping ", length(missing_in_annotation), " seed(s) not found in the array annotation: ",
                paste(head(missing_in_annotation, 10), collapse = ","),
                if (length(missing_in_annotation) > 10) " ..." else ""
            )
            seeds_tsv <- seeds_tsv[!(seeds_tsv[, seeds_id_col] %in% missing_in_annotation), , drop = FALSE]
            seeds <- setdiff(seeds, missing_in_annotation)
        }
        missing_in_beta <- setdiff(seeds, beta_row_names)
        if (length(missing_in_beta) > 0) {
            .log_warn(
                "Dropping ", length(missing_in_beta), " seed(s) not found in the beta file: ",
                paste(head(missing_in_beta, 10), collapse = ","),
                if (length(missing_in_beta) > 10) " ..." else ""
            )
            seeds_tsv <- seeds_tsv[!(seeds_tsv[, seeds_id_col] %in% missing_in_beta), , drop = FALSE]
            seeds <- setdiff(seeds, missing_in_beta)
        }
        if (length(seeds) == 0) {
            stop("No seeds remain after filtering against array annotation.")
        }
        seeds <- seeds[orderByLoc(seeds, genome = genome, genomic_locs = beta_locs)]
        seeds_inds <- match(seeds, rownames(beta_locs))
    } else {
        seeds <- seeds_tsv[, seeds_id_col]
        seeds_inds <- seeds
    }

    .log_step("Validating beta file sorting by position...", level = 2)


    .log_step("Subsetting beta matrix for seeds...", level = 2)
    seeds_locs <- as.data.frame(beta_locs[seeds, , drop = FALSE])
    seeds_beta <- beta_handler$getBeta(row_names = seeds, col_names = beta_col_names_detection, check_mem = TRUE)
    if (!is.null(output_prefix)) {
        seeds_beta_output_file <- paste0(output_prefix, "seeds_beta.tsv.gz")
    }

    if (!is.null(output_prefix)) {
        .log_step("Saving seeds beta to file: ", seeds_beta_output_file, " ...", level = 2)
        gz <- gzfile(seeds_beta_output_file, "w")
        write.table(seeds_beta, gz, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
        close(gz)
    }

    if (nrow(seeds_locs) != nrow(seeds_beta)) {
        stop(
            "Number of rows in the queried seeds beta file does not match the number of seeds. Number of rows in beta file: ",
            nrow(seeds_beta),
            " Number of rows in seeds: ",
            nrow(seeds_locs)
        )
    }
    .log_info("Checking for seeds with all NA beta values...", level = 3)
    if (bigmemory::is.big.matrix(seeds_beta)) {
        beta_chunk_size <- beta_handler$getBetaChunkSize()
        all.na.rows <- rep(FALSE, nrow(seeds_beta))
        for (start_row in seq(1, nrow(seeds_beta), by = beta_chunk_size)) {
            end_row <- min(start_row + beta_chunk_size - 1, nrow(seeds_beta))
            chunk <- bigmemory::sub.big.matrix(seeds_beta, firstRow = start_row, lastRow = end_row)
            chunk <- as.matrix(chunk)
            na_rows_chunk <- matrixStats::rowAlls(is.na(as.matrix(chunk)))
            if (any(na_rows_chunk)) {
                all.na.rows[start_row:end_row][na_rows_chunk] <- TRUE
                break
            }
        }
    } else {
        all.na.rows <- matrixStats::rowAlls(is.na(as.matrix(seeds_beta)))
    }
    if (any(all.na.rows)) {
        stop(
            "Beta extraction failure: the following Seed rows have all NA beta values: ",
            paste(rownames(seeds_beta)[all.na.rows], collapse = ","),
            ". This indicates a mismatch between requested CpG IDs and beta file columns or a parsing issue."
        )
    }

    .log_info("Subset size: ", paste(dim(seeds_beta), collapse = ","), level = 2)
    .log_info("Number of provided seeds: ", length(seeds), level = 2)
    resolved_min_cpg_delta_beta <- as.numeric(min_cpg_delta_beta)
    if (adaptive_min_cpg_delta_beta) {
        resolved_min_cpg_delta_beta <- .resolveAdaptiveMinCpgDeltaBeta(
            seeds_beta = seeds_beta,
            pheno = pheno_detection,
            aggfun = aggfun,
            base_threshold = resolved_min_cpg_delta_beta
        )
    }
    .log_info(
        "Using min_cpg_delta_beta threshold: ", signif(resolved_min_cpg_delta_beta, 4),
        if (adaptive_min_cpg_delta_beta) " (adaptive)." else " (fixed).",
        level = 2
    )

    .log_success("Input preparation complete.", level = 1)
    .log_step("Stage 1: Connecting seeds to form initial DMRs..", level = 1)

    # Set up progress tracking for Seed connection
    chromosomes <- unique(seeds_locs[, "chr"])

    if (verbose > 1 && .load_debug && file.exists(file.path("debug", "01_dmrs_from_connected_seeds.tsv"))) {
        .log_info("Loading debug DMRs from file...", level = 2)
        dmrs <- read.table(
            file.path("debug", "01_dmrs_from_connected_seeds.tsv"),
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
        p_con <- NULL
        if (verbose > 0) {
            # check if version of progressr is equal or higher than >= 0.17.0-9002, otherwise p_con will not be used

            if (utils::packageVersion("progressr") >= "0.17.0-9002") {
                p_con <- progressr::progressor(steps = length(chromosomes), message = "Connecting seeds to form DMRs..")
            }
        }
        # Split by chromosome for parallel processing
        seeds_list <- split(seeds, seeds_locs[, "chr"])
        seeds_list <- seeds_list[as.character(chromosomes)]

        seeds_locs_list <- split(seeds_locs, seeds_locs[, "chr"])
        seeds_locs_list <- seeds_locs_list[as.character(chromosomes)]
        seeds_inds_list <- split(seeds_inds, seeds_locs[, "chr"])
        seeds_inds_list <- seeds_inds_list[as.character(chromosomes)]
        if (bigmemory::is.big.matrix(seeds_beta)) {
            seeds_beta_list <- lapply(seeds_locs_list, function(loc_df) {
                bigmemory::sub.big.matrix(
                    seeds_beta,
                    firstRow = min(match(rownames(seeds_beta), rownames(loc_df)), na.rm = TRUE),
                    lastRow = max(match(rownames(seeds_beta), rownames(loc_df)), na.rm = TRUE)
                )
            })
        } else {
            if (!is.matrix(seeds_beta)) seeds_beta <- as.matrix(seeds_beta)
            storage.mode(seeds_beta) <- "double"
            seeds_beta_list <- lapply(split(seeds_beta, seeds_locs[, "chr"]), matrix, ncol = ncol(seeds_beta))
            seeds_beta_list <- seeds_beta_list[as.character(chromosomes)]
        }
        fun <- function(chr, cseeds, cseeds_inds, cseeds_beta, cseeds_locs) {
            op <- options(warn = 2)$warn
            dmr_list <- vector("list", length = 128L)
            dmr_n <- 0L
            connection_corr_pval <- NA
            dmr_seeds_inds <- integer(0)
            .log_step("Testing Seed connectivity on chromosome ", chr, " ...", level = 3)
            corr_ret <- .testConnectivityBatch(
                sites_beta = cseeds_beta,
                group_inds = group_inds,
                pheno = pheno_detection,
                max_lookup_dist = max_lookup_dist,
                sites_locs = cseeds_locs,
                max_pval = max_pval,
                min_delta_beta = resolved_min_cpg_delta_beta,
                aggfun = aggfun,
                entanglement = entanglement,
                pval_mode = pval_mode,
                empirical_strategy = empirical_strategy,
                ntries = ntries,
                mid_p = mid_p
            )
            corr_ret <- .bridgeConnectivityGaps(
                corr_ret = corr_ret,
                max_bridge_seeds_gaps = max_bridge_seeds_gaps
            )
            bridged_edges <- attr(corr_ret, "bridged_edges")
            if (!is.null(bridged_edges) && bridged_edges > 0L) {
                .log_info("Bridged ", bridged_edges, " seed-connection gap edge(s) on chromosome ", chr, ".", level = 3)
            }
            stopifnot(nrow(corr_ret) == nrow(cseeds_beta) - 1)
            if (getOption("DMRsegal.make_debug_dir", FALSE)) {
                .log_info("Saving seed connectivity results to debug/seeds_connectivity_chr", chr, ".tsv", level = 3)
                dir.create("debug", showWarnings = FALSE)
                write.table(corr_ret,
                    file = paste0("debug/seeds_connectivity_chr", chr, ".tsv"),
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE,
                    quote = FALSE
                )
            }
            .log_success("Seed connectivity tested.", level = 3)
            breakpoints <- c(which(!corr_ret$connected), nrow(cseeds_beta))
            start_seed_ind <- 1
            single_seeds <- 0
            for (bp_ind in seq_along(breakpoints)) {
                end_seed_ind <- breakpoints[bp_ind]
                if (end_seed_ind - start_seed_ind + 1 < min_seeds) {
                    .log_info("Skipping DMR from Seed ", start_seed_ind, " to Seed ", end_seed_ind, " (id: ", chr, ":", cseeds_locs[start_seed_ind, "start"], "-", cseeds_locs[end_seed_ind, "start"], ") due to insufficient number of connected seeds (", end_seed_ind - start_seed_ind + 1, " < ", min_seeds, ").", level = 4)
                    start_seed_ind <- breakpoints[bp_ind] + 1
                    next
                }
                if (end_seed_ind > start_seed_ind) {
                    .log_step("Registering ", bp_ind, "/", (length(breakpoints) - 1), " DMR from Seed ", start_seed_ind, " to Seed ", end_seed_ind, " (id: ", chr, ":", cseeds_locs[start_seed_ind, "start"], "-", cseeds_locs[end_seed_ind, "start"], ")", level = 3)
                } else {
                    single_seeds <- single_seeds + 1
                }
                dmr_seeds_inds <- seq.int(start_seed_ind, end_seed_ind)
                bridged_inds <- which(corr_ret$reason[dmr_seeds_inds] == "bridged")
                bridged_inds <- dmr_seeds_inds[bridged_inds]
                dmr_seeds_inds <- setdiff(dmr_seeds_inds, bridged_inds)
                if (end_seed_ind == start_seed_ind) {
                    connection_corr_pval <- NA
                } else {
                    connection_corr_pval <- aggfun(corr_ret$pval[dmr_seeds_inds[-length(dmr_seeds_inds)]], na.rm = TRUE)
                }
                if (end_seed_ind < nrow(cseeds_beta)) {
                    stop_reason <- corr_ret$reason[[end_seed_ind]]
                } else {
                    stop_reason <- "end of chromosome"
                }
                new_dmr <- list(
                    chr = chr,
                    start_seed = cseeds[[start_seed_ind]],
                    end_seed = cseeds[[end_seed_ind]],
                    start_seed_pos = cseeds_locs[start_seed_ind, "start"],
                    end_seed_pos = cseeds_locs[end_seed_ind, "start"],
                    seeds_num = length(dmr_seeds_inds),
                    connection_corr_pval = connection_corr_pval,
                    stop_connection_reason = stop_reason,
                    seeds_inds = paste(cseeds_inds[dmr_seeds_inds], collapse = ","),
                    seeds = paste(cseeds[dmr_seeds_inds], collapse = ",")
                )
                dmr_n <- dmr_n + 1L
                if (dmr_n > length(dmr_list)) length(dmr_list) <- length(dmr_list) * 2L
                dmr_list[[dmr_n]] <- new_dmr
                start_seed_ind <- breakpoints[bp_ind] + 1
                if (end_seed_ind > start_seed_ind) {
                    .log_success("DMR registered.",
                        level = 3
                    )
                } else {
                    single_seeds <- single_seeds + 1
                }
            }
            if (single_seeds > 0L) {
                .log_info("Also registered ", single_seeds, " DMR(s) consisting of a single seed on chromosome ", chr, ".", level = 3)
            }
            if (verbose > 0 && !is.null(p_con)) p_con()
            dmrs <- if (dmr_n > 0L) data.table::rbindlist(dmr_list[seq_len(dmr_n)], fill = TRUE) else data.frame()
            if (nrow(dmrs)) {
                rownames(dmrs) <- seq_len(nrow(dmrs))
                dmrs[, "chr"] <- chr
            }
            options(warn = op)
            dmrs
        }
        if (njobs == 1) {
            ret <- lapply(seq_along(chromosomes), function(i) {
                chr <- chromosomes[i]
                cseeds <- seeds_list[[i]]
                cseeds_inds <- seeds_inds_list[[i]]
                cseeds_beta <- seeds_beta_list[[i]]
                cseeds_locs <- seeds_locs_list[[i]]
                fun(
                    chr = chr,
                    cseeds = cseeds,
                    cseeds_inds = cseeds_inds,
                    cseeds_beta = cseeds_beta,
                    cseeds_locs = cseeds_locs
                )
            })
        } else {
            .setupParallel()
            ret <- future.apply::future_mapply(
                chromosomes,
                seeds_list,
                seeds_inds_list,
                seeds_beta_list,
                seeds_locs_list,
                SIMPLIFY = FALSE,
                future.seed = TRUE,
                future.stdout = NA,
                future.globals = c(
                    "group_inds",
                    "case_mask",
                    "max_lookup_dist",
                    "max_pval",
                    "resolved_min_cpg_delta_beta",
                    "aggfun",
                    "entanglement",
                    "pval_mode",
                    "empirical_strategy",
                    "ntries",
                    "mid_p",
                    "max_bridge_seeds_gaps",
                    "min_seeds",
                    ".testConnectivityBatch",
                    ".bridgeConnectivityGaps",
                    ".log_step",
                    ".log_success",
                    ".log_info",
                    "p_con",
                    "verbose"
                ),
                FUN = fun
            )
            .finalizeParallel()
        }

        if (inherits(ret[[1]], "try-error")) {
            stop(ret)
        }
        dmrs <- as.data.frame(do.call(rbind, ret))
        if (nrow(dmrs) == 0) {
            .log_warn("No DMRs remain after filtering based on min_seeds.")
            if (!is.null(output_prefix)) {
                for (f in c(".methylation.tsv.gz", ".tsv.gz")) {
                    gzfile <- gzfile(paste0(output_prefix, f), "w", compression = 2)
                    close(gzfile)
                }
            }
            return(NULL)
        }
        if (getOption("DMRsegal.make_debug_dir", FALSE)) {
            .log_info("Saving initial DMRs from connected seeds to debug/01_dmrs_from_connected_seeds.tsv", level = 1)
            dir.create("debug", showWarnings = FALSE)
            write.table(dmrs,
                file = file.path("debug", "01_dmrs_from_connected_seeds.tsv"),
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE
            )
        }
    }
    .log_success("Initial DMRs formed: ", nrow(dmrs), level = 1)
    .log_info("Summary of connected seeds DMRs:\n\t", paste(capture.output(summary(dmrs)), collapse = "\n\t"), level = 3)
    .log_step("Stage 2: Expanding DMRs on neighborhood CpGs..", level = 1)
    # Set up progress tracking for DMR expansion
    n_dmrs <- nrow(dmrs)
    if (verbose > 1 && .load_debug && file.exists(file.path("debug", "connectivity_array.rds"))) {
        .log_info("Loading debug connectivity array from file...", level = 2)
        connectivity_array <- readRDS(file.path("debug", "connectivity_array.rds"))
    } else {
        connectivity_windows <- NULL
        if (is.finite(connectivity_window_bp) && connectivity_window_bp > 0) {
            connectivity_windows <- .buildConnectivityWindowsFromDMRs(
                dmrs = dmrs,
                connectivity_window_bp = connectivity_window_bp
            )
            if (nrow(connectivity_windows) > 0L) {
                .log_info(
                    "Stage 2 connectivity restricted to ", nrow(connectivity_windows),
                    " seed-derived windows (total span: ",
                    format(sum(connectivity_windows$end - connectivity_windows$start + 1), big.mark = ","),
                    " bp).",
                    level = 2
                )
            } else {
                .log_info("No connectivity windows were generated; Stage 2 will return disconnected CpGs outside chromosome termini.", level = 2)
            }
        } else {
            .log_info("Stage 2 connectivity computed genome-wide (connectivity_window_bp <= 0).", level = 2)
        }
        .log_step("Building connectivity array..", level = 2)
        connectivity_array <- NULL
        splits <- NULL
        for (gap in seq(0L, max_bridge_extension_gaps)) {
            if (gap == 0L) {
                .log_info("Building initial connectivity array with no gap bridging.", level = 2)
            } else {
                .log_info("Building bridged connectivity array allowing up to ", gap, " gap(s) between connected seeds.", level = 2)
            }
            build_ret <- .buildConnectivityArray(
                beta_handler = beta_handler,
                group_inds = group_inds,
                col_names = beta_col_names_detection,
                pheno = pheno_detection,
                covariates = covariates,
                max_lookup_dist = max_lookup_dist,
                max_pval = max_pval,
                min_delta_beta = resolved_min_cpg_delta_beta,
                entanglement = entanglement,
                aggfun = aggfun,
                pval_mode = pval_mode,
                empirical_strategy = empirical_strategy,
                chunk_size = chunk_size,
                ntries = ntries,
                mid_p = mid_p,
                njobs = njobs,
                connectivity_windows = connectivity_windows,
                connectivity_array = connectivity_array,
                splits = splits,
                gap = gap
            )
            connectivity_array <- build_ret$connectivity_array
            splits <- build_ret$splits
            .log_info("Connectivity array built with gap allowance of ", gap, " (", sum(connectivity_array$connected), " connected CpGs).", level = 2)
        }
    }
    .log_success("Connectivity array built.", level = 2)
    .log_info("Number of underlying connected CpGs found: ", sum(connectivity_array$connected), level = 1)
    if (getOption("DMRsegal.make_debug_dir", FALSE)) {
        .log_info("Saving connectivity array to debug/connectivity_array.rds", level = 1)
        dir.create("debug", showWarnings = FALSE)
        saveRDS(connectivity_array, file = file.path("debug", "connectivity_array.rds"))
    }
    stopifnot(nrow(connectivity_array) == nrow(beta_locs))
    .log_info("Expanding ", n_dmrs, " DMRs using ", njobs, " jobs...", level = 2)
    if (verbose > 0) {
        p_ext <- NULL
        if (verbose > 0) {
            # check if version of progressr is equal or higher than >= 0.17.0-9002, otherwise p_ext will not be used

            if (utils::packageVersion("progressr") >= "0.17.0-9002") {
                p_ext <- progressr::progressor(steps = n_dmrs, message = "Expanding DMRs to proximal CpGs..")
            }
        }
    }
    ret <- list()
    for (chr in unique(dmrs$chr)) {
        .log_info("Processing ", chr, level = 2)
        chr_dmrs <- dmrs[dmrs$chr == chr, ]
        chr_mask <- beta_locs[, "chr"] == chr
        first_row <- which(chr_mask)[1]
        last_row <- which(chr_mask)[sum(chr_mask)]
        chr_start_base <- first_row - 1
        if (bigmemory::is.big.matrix(beta_locs)) {
            chr_locs <- bigmemory::sub.big.matrix(
                beta_locs,
                firstRow = first_row,
                lastRow = last_row
            )
        } else {
            chr_locs <- beta_locs[chr_mask, , drop = FALSE]
        }
        chr_array <- connectivity_array[chr_mask, , drop = FALSE]
        .setupParallel()
        chr_ret <- future.apply::future_apply(
            X = chr_dmrs,
            MARGIN = 1,
            FUN = function(dmr) {
                op <- options(warn = 2)$warn
                x <- .expandDMR(
                    dmr = dmr,
                    chr_array = chr_array,
                    expansion_step = expansion_step,
                    min_cpgs = min_cpgs,
                    min_cpg_delta_beta = min_cpg_delta_beta,
                    chr_locs = chr_locs,
                    chr_start_base = chr_start_base
                )
                options(warn = op)
                if (verbose > 0 && !is.null(p_ext)) p_ext()
                x
            },
            simplify = FALSE,
            future.seed = TRUE,
            future.globals = c(
                ".expandDMR",
                "chr_array",
                "expansion_step",
                "min_cpgs",
                "min_cpg_delta_beta",
                "chr_locs",
                "chr_start_base",
                "verbose",
                "p_ext"
            )
        )
        .finalizeParallel()
        .log_info("Chromosome ", chr, ": Number of DMRs processed: ", length(chr_ret), level = 2)
        ret <- c(ret, chr_ret)
    }
    if (inherits(ret, "try-error")) {
        stop(ret)
    }
    upstream_cpg_expansion_table <- table(sapply(ret, function(x) x[["upstream_cpg_expansion"]]))
    downstream_cpg_expansion_table <- table(sapply(ret, function(x) x[["downstream_cpg_expansion"]]))
    # sort tables by names (expansion sizes)
    upstream_cpg_expansion_table <- upstream_cpg_expansion_table[order(as.integer(names(upstream_cpg_expansion_table)))]
    downstream_cpg_expansion_table <- downstream_cpg_expansion_table[order(as.integer(names(downstream_cpg_expansion_table)))]
    .log_info("Table of upstream_cpg_expansion:\n\t", paste(capture.output(upstream_cpg_expansion_table), collapse = "\n\t"), level = 2)
    .log_info("Table of downstream_cpg_expansion:\n\t", paste(capture.output(downstream_cpg_expansion_table), collapse = "\n\t"), level = 2)

    .log_step("Post-processing extended DMRs..", level = 2)

    extended_dmrs <- as.data.frame(do.call(rbind, ret))
    extended_dmrs$end <- as.numeric(extended_dmrs$end)
    extended_dmrs$start <- as.numeric(extended_dmrs$start)
    extended_dmrs$start_seed_pos <- as.numeric(extended_dmrs$start_seed_pos)
    extended_dmrs$end_seed_pos <- as.numeric(extended_dmrs$end_seed_pos)
    extended_dmrs$seeds_num <- as.numeric(extended_dmrs$seeds_num)
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
        extended_dmrs$start_cpg_ind <- match(extended_dmrs$start_cpg, rownames(beta_locs))
        extended_dmrs$end_cpg_ind <- match(extended_dmrs$end_cpg, rownames(beta_locs))
    } else {
        extended_dmrs$start_cpg_ind <- as.integer(extended_dmrs$start_cpg)
        extended_dmrs$end_cpg_ind <- as.integer(extended_dmrs$end_cpg)
    }


    .log_success("Post-processing complete.", level = 2)
    .log_success("DMR expansion complete.", level = 1)
    .log_info("Summary of extended DMRs:\n\t", paste(capture.output(summary(extended_dmrs)), collapse = "\n\t"), level = 3)
    if (getOption("DMRsegal.make_debug_dir", FALSE)) {
        .log_info("Saving extended DMRs prior to filtering to debug/03_extended_dmrs.tsv", level = 1)
        dir.create("debug", showWarnings = FALSE)
        write.table(extended_dmrs,
            file = file.path("debug", "02_extended_dmrs.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
    }
    .log_step("Stage 3: Merging overlapping extended DMRs..", level = 1)

    if (bed_provided) {
        # Change seqnames from factor to character
        extended_dmrs[, "chr"] <- chr_levels[as.integer(extended_dmrs[, "chr"])]
    }
    extended_dmrs_ranges <- GenomicRanges::makeGRangesFromDataFrame(
        extended_dmrs,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqnames.field = "chr",
        seqinfo = GenomeInfoDb::Seqinfo(genome = genome),
        na.rm = TRUE
    )

    merged_dmrs_ranges <- GenomicRanges::reduce(extended_dmrs_ranges, ignore.strand = TRUE)
    hits <- GenomicRanges::findOverlaps(merged_dmrs_ranges, extended_dmrs_ranges, ignore.strand = TRUE)
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    orig_mcols <- as.data.frame(GenomicRanges::mcols(extended_dmrs_ranges))
    mcol_names <- colnames(orig_mcols)
    agg_df <- data.frame(matrix(NA, nrow = length(merged_dmrs_ranges), ncol = length(mcol_names) + 1))
    colnames(agg_df) <- c(mcol_names, "merged_dmrs_num")
    # ranges that do not need to be merged will have only one hit
    tqh <- table(qh)
    .log_info("Frequency of N-DMRs overlap:\n", paste(capture.output(print(table(tqh))), collapse = "\n"), level = 2)
    single_hits <- which(tqh == 1)
    if (length(single_hits) > 0) {
        # get the corresponding indices in the original extended_dmrs_ranges
        .log_info("Copying over ", length(single_hits), " non-overlapping extended DMRs...", level = 2)
        # do it vectorizedly
        agg_df[single_hits, ] <- orig_mcols[sh[qh %in% single_hits], ]
        agg_df[single_hits, "merged_dmrs_num"] <- 1
    }

    multiple_hits <- which(tqh > 1)
    .log_info("Merging ", length(multiple_hits), " overlapping extended DMRs...", level = 2)
    for (i in multiple_hits) {
        inds <- sh[qh == i]
        cols_vals <- orig_mcols[inds, ]
        agg_df[i, "start_cpg_ind"] <- min(as.integer(cols_vals$start_cpg_ind))
        agg_df[i, "end_cpg_ind"] <- max(as.integer(cols_vals$end_cpg_ind))
        agg_df[i, "start_seed"] <- cols_vals$start_seed[[1]]
        agg_df[i, "end_seed"] <- cols_vals$end_seed[[length(inds)]]
        agg_df[i, "start_seed_pos"] <- cols_vals$start_seed_pos[[1]]
        agg_df[i, "end_seed_pos"] <- cols_vals$end_seed_pos[[length(inds)]]
        agg_seeds <- unique(unlist(lapply(cols_vals$seeds, .splitCsvValues), use.names = FALSE))
        agg_df[i, "seeds"] <- paste(agg_seeds, collapse = ",")
        agg_df[i, "seeds_num"] <- length(agg_seeds)
        agg_seeds_inds <- suppressWarnings(as.integer(unlist(lapply(cols_vals$seeds_inds, .splitCsvValues), use.names = FALSE)))
        agg_seeds_inds <- sort(unique(agg_seeds_inds[!is.na(agg_seeds_inds)]))
        agg_df[i, "seeds_inds"] <- paste(agg_seeds_inds, collapse = ",")
        agg_df[i, "connection_corr_pval"] <- aggfun(as.double(cols_vals$connection_corr_pval), na.rm = TRUE)
        agg_df[i, "stop_connection_reason"] <- paste(cols_vals$stop_connection_reason, collapse = ",")
        agg_df[i, "start_cpg"] <- cols_vals$start_cpg[[1]]
        agg_df[i, "end_cpg"] <- cols_vals$end_cpg[[length(inds)]]
        agg_upstream_cpgs <- unique(unlist(lapply(cols_vals$upstream_cpgs, .splitCsvValues), use.names = FALSE))
        agg_df[i, "upstream_cpg_expansion"] <- paste(cols_vals$upstream_cpg_expansion, collapse = ",")
        agg_df[i, "upstream_cpgs"] <- paste(agg_upstream_cpgs, collapse = ",")
        agg_upstream_inds <- suppressWarnings(as.integer(unlist(lapply(cols_vals$upstream_cpgs_inds, .splitCsvValues), use.names = FALSE)))
        agg_upstream_inds <- sort(unique(agg_upstream_inds[!is.na(agg_upstream_inds)]))
        agg_df[i, "upstream_cpgs_inds"] <- paste(agg_upstream_inds, collapse = ",")
        agg_df[i, "upstream_cpg_expansion_stop_reason"] <- paste(cols_vals$upstream_cpg_expansion_stop_reason, collapse = ",")
        agg_df[i, "downstream_cpg_expansion"] <- paste(cols_vals$downstream_cpg_expansion, collapse = ",")
        agg_downsteam_cpgs <- unique(unlist(lapply(cols_vals$downstream_cpgs, .splitCsvValues), use.names = FALSE))
        agg_df[i, "downstream_cpgs"] <- paste(agg_downsteam_cpgs, collapse = ",")
        agg_downstream_inds <- suppressWarnings(as.integer(unlist(lapply(cols_vals$downstream_cpgs_inds, .splitCsvValues), use.names = FALSE)))
        agg_downstream_inds <- sort(unique(agg_downstream_inds[!is.na(agg_downstream_inds)]))
        agg_df[i, "downstream_cpgs_inds"] <- paste(agg_downstream_inds, collapse = ",")
        agg_df[i, "downstream_cpg_expansion_stop_reason"] <- paste(cols_vals$downstream_cpg_expansion_stop_reason, collapse = ",")
        agg_df[i, "merged_dmrs_num"] <- length(inds)
    }
    agg_df[, "cpgs"] <- apply(agg_df[, c("seeds", "upstream_cpgs", "downstream_cpgs")], 1, function(x) {
        vals <- unique(unlist(lapply(x, .splitCsvValues), use.names = FALSE))
        paste(vals, collapse = ",")
    })
    agg_df[, "supporting_cpgs_num"] <- vapply(agg_df$cpgs, function(x) {
        length(.splitCsvValues(x))
    }, integer(1))
    agg_df[, "cpgs_num"] <- agg_df[, "end_cpg_ind"] - agg_df[, "start_cpg_ind"] + 1
    agg_df[, "cpgs_inds"] <- apply(agg_df[, c("seeds_inds", "upstream_cpgs_inds", "downstream_cpgs_inds")], 1, function(x) {
        vals <- suppressWarnings(as.integer(unlist(lapply(x, .splitCsvValues), use.names = FALSE)))
        vals <- sort(unique(vals[!is.na(vals)]))
        paste(vals, collapse = ",")
    })
    agg_df[, "id"] <- paste0(seqnames(merged_dmrs_ranges), ":", agg_df$start_cpg, "-", agg_df$end_cpg)

    GenomicRanges::mcols(merged_dmrs_ranges) <- agg_df
    .log_success("Overlapping extended DMRs merged: ", length(merged_dmrs_ranges), " resulting DMRs.", level = 1)
    merged_dmrs <- as.data.frame(merged_dmrs_ranges)
    colnames(merged_dmrs)[colnames(merged_dmrs) == "seqnames"] <- "chr"
    .log_info("Summary of merged DMRs:\n\t", paste(capture.output(summary(merged_dmrs)), collapse = "\n\t"), level = 3)

    if (getOption("DMRsegal.make_debug_dir", FALSE)) {
        .log_info("Saving merged DMRs prior to filtering to debug/03_merged_dmrs.tsv", level = 1)
        dir.create("debug", showWarnings = FALSE)
        write.table(merged_dmrs,
            file = file.path("debug", "03_merged_dmrs.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
    }

    .log_step("Stage 4: Filtering and annotating resulting DMRs..", level = 1)

    if (min_cpgs > 1 || min_seeds > 1) {
        filtered_dmrs_ranges <- merged_dmrs_ranges[
            GenomicRanges::mcols(merged_dmrs_ranges)$seeds_num >= min_seeds &
                GenomicRanges::mcols(merged_dmrs_ranges)$cpgs_num >= min_cpgs
        ]
        .log_info(
            "Keeping ",
            length(filtered_dmrs_ranges),
            " out of ",
            length(merged_dmrs_ranges),
            " with at least ",
            min_seeds,
            " supporting seeds and ",
            min_cpgs,
            " CpGs.",
            level = 1
        )
    } else {
        filtered_dmrs_ranges <- merged_dmrs_ranges
    }
    filtered_dmrs <- as.data.frame(filtered_dmrs_ranges)
    colnames(filtered_dmrs)[colnames(filtered_dmrs) == "seqnames"] <- "chr"

    array_based <- !is.null(array)
    underlying_cpgs_num <- filtered_dmrs$end_cpg_ind - filtered_dmrs$start_cpg_ind + 1

    if (array_based && min_adj_seeds > min_seeds) {
        .log_step("Calculating CpG content and adjusted seeds number..", level = 2)
        filtered_dmrs$cpgs_num_bg <- getCpGBackgroundCounts(filtered_dmrs_ranges, genome)

        underlying_cpgs_num[underlying_cpgs_num == 0] <- 1

        filtered_dmrs$seeds_num_adj <- ceiling(underlying_cpgs_num / filtered_dmrs[, "cpgs_num_bg"] * filtered_dmrs[, "seeds_num"])

        .log_success("CpG content calculated.", level = 2)
        adj_filtered_dmrs <- filtered_dmrs[
            filtered_dmrs$seeds_num_adj >= min_adj_seeds, ,
            drop = FALSE
        ]
        .log_info(
            "Keeping ",
            nrow(adj_filtered_dmrs),
            " out of ",
            nrow(filtered_dmrs),
            " with at least ",
            min_adj_seeds,
            " (adjusted) supporting seeds.",
            level = 1
        )
        filtered_dmrs <- adj_filtered_dmrs
    } else {
        .log_info("Skipping adjusted seeds number calculation as min_adj_seeds <= min_seeds.", level = 2)
        filtered_dmrs$cpgs_num_bg <- NA
        filtered_dmrs$seeds_num_adj <- NA
    }
    if (nrow(filtered_dmrs) == 0) {
        .log_warn("No DMRs passed the filtering step.")
        if (!is.null(output_prefix)) {
            for (f in c(".methylation.tsv.gz", ".tsv.gz")) {
                gzfile <- gzfile(file.path(output_dir, paste0(output_prefix, f)), "w", compression = 2)
                close(gzfile)
            }
        }
        return(NULL)
    }

    annotated_dmrs <- filtered_dmrs

    all_selected_cpgs <- unique(unlist(lapply(seq_len(nrow(annotated_dmrs)), function(i) {
        seq(annotated_dmrs$start_cpg_ind[i], annotated_dmrs$end_cpg_ind[i])
    })))
    all_selected_cpgs_beta <- beta_handler$getBeta(row_names = all_selected_cpgs, col_names = beta_col_names)
    .log_step("Calculating per-CpG beta statistics..", level = 2)
    beta_stats <- .calculateBetaStats(
        beta_values = all_selected_cpgs_beta,
        pheno = pheno_all,
        aggfun = aggfun
    )
    .log_success("Per-CpG beta statistics calculated.", level = 2)

    beta_stats <- as.data.frame(beta_stats)
    rownames(beta_stats) <- all_selected_cpgs
    dmrs_seeds_inds <- strsplit(annotated_dmrs$seeds_inds, ",")
    .log_step("Adding DMR delta-beta information..", level = 2)

    dmr_seeds_list <- lapply(dmrs_seeds_inds, as.character)
    dmr_seeds_indices <- unlist(dmr_seeds_list, use.names = FALSE)
    dmr_seeds_groups <- rep(seq_along(dmr_seeds_list), lengths(dmr_seeds_list))

    beta_stats_seeds <- beta_stats[dmr_seeds_indices, , drop = FALSE]
    beta_stats_seeds$dmr_id <- dmr_seeds_groups

    seeds_agg <- data.table::as.data.table(beta_stats_seeds)[, .(
        cases_beta = aggfun(abs(cases_beta)) * sign(sum(sign(cases_beta))),
        controls_beta = aggfun(abs(controls_beta)) * sign(sum(sign(controls_beta))),
        cases_beta_sd = aggfun(cases_beta_sd, na.rm = TRUE),
        controls_beta_sd = aggfun(controls_beta_sd, na.rm = TRUE),
        cases_beta_min = min(cases_beta, na.rm = TRUE),
        cases_beta_max = max(cases_beta, na.rm = TRUE),
        controls_beta_min = min(controls_beta, na.rm = TRUE),
        controls_beta_max = max(controls_beta, na.rm = TRUE)
    ), by = dmr_id]

    annotated_dmrs$cases_beta <- seeds_agg$cases_beta
    annotated_dmrs$controls_beta <- seeds_agg$controls_beta
    annotated_dmrs$delta_beta <- annotated_dmrs$cases_beta - annotated_dmrs$controls_beta
    annotated_dmrs$cases_beta_sd <- seeds_agg$cases_beta_sd
    annotated_dmrs$controls_beta_sd <- seeds_agg$controls_beta_sd
    annotated_dmrs$cases_beta_min <- seeds_agg$cases_beta_min
    annotated_dmrs$cases_beta_max <- seeds_agg$cases_beta_max
    annotated_dmrs$controls_beta_min <- seeds_agg$controls_beta_min
    annotated_dmrs$controls_beta_max <- seeds_agg$controls_beta_max

    dmr_cpgs_list <- lapply(seq_len(nrow(annotated_dmrs)), function(i) {
        as.character(c(
            seq.int(annotated_dmrs$start_cpg_ind[i], as.integer(dmrs_seeds_inds[[i]][1])),
            seq.int(as.integer(dmrs_seeds_inds[[i]][length(dmrs_seeds_inds[[i]])]), annotated_dmrs$end_cpg_ind[i])
        ))
    })
    dmr_cpgs_indices <- unlist(dmr_cpgs_list, use.names = FALSE)
    dmr_cpgs_groups <- rep(seq_along(dmr_cpgs_list), lengths(dmr_cpgs_list))

    beta_stats_cpgs <- beta_stats[dmr_cpgs_indices, , drop = FALSE]
    beta_stats_cpgs$dmr_id <- dmr_cpgs_groups

    cpgs_agg <- data.table::as.data.table(beta_stats_cpgs)[, .(
        cpgs_cases_beta = aggfun(abs(cases_beta)) * sign(sum(sign(cases_beta))),
        cpgs_controls_beta = aggfun(abs(controls_beta)) * sign(sum(sign(controls_beta))),
        cpgs_cases_beta_sd = aggfun(cases_beta_sd, na.rm = TRUE),
        cpgs_controls_beta_sd = aggfun(controls_beta_sd, na.rm = TRUE),
        cpgs_cases_beta_min = min(cases_beta, na.rm = TRUE),
        cpgs_cases_beta_max = max(cases_beta, na.rm = TRUE),
        cpgs_controls_beta_min = min(controls_beta, na.rm = TRUE),
        cpgs_controls_beta_max = max(controls_beta, na.rm = TRUE)
    ), by = dmr_id]

    annotated_dmrs$cpgs_cases_beta <- cpgs_agg$cpgs_cases_beta
    annotated_dmrs$cpgs_controls_beta <- cpgs_agg$cpgs_controls_beta
    annotated_dmrs$cpgs_delta_beta <- annotated_dmrs$cpgs_cases_beta - annotated_dmrs$cpgs_controls_beta
    annotated_dmrs$cpgs_cases_beta_sd <- cpgs_agg$cpgs_cases_beta_sd
    annotated_dmrs$cpgs_controls_beta_sd <- cpgs_agg$cpgs_controls_beta_sd
    annotated_dmrs$cpgs_cases_beta_min <- cpgs_agg$cpgs_cases_beta_min
    annotated_dmrs$cpgs_cases_beta_max <- cpgs_agg$cpgs_cases_beta_max
    annotated_dmrs$cpgs_controls_beta_min <- cpgs_agg$cpgs_controls_beta_min
    annotated_dmrs$cpgs_controls_beta_max <- cpgs_agg$cpgs_controls_beta_max

    .log_success("DMR delta-beta information added.", level = 2)


    if (bed_provided) {
        annotated_dmrs$chr <- chr_levels[annotated_dmrs$chr]
        start_seeds <- as.integer(annotated_dmrs$start_seed)
        end_seeds <- as.integer(annotated_dmrs$end_seed)
        annotated_dmrs$start_seed <- paste0(chr_levels[beta_locs[start_seeds, "chr"]], ":", beta_locs[start_seeds, "start"])
        annotated_dmrs$end_seed <- paste0(chr_levels[beta_locs[end_seeds, "chr"]], ":", beta_locs[end_seeds, "start"])
        start_cpgs <- as.integer(annotated_dmrs$start_cpg)
        end_cpgs <- as.integer(annotated_dmrs$end_cpg)
        annotated_dmrs$start_cpg <- paste0(chr_levels[beta_locs[start_cpgs, "chr"]], ":", beta_locs[start_cpgs, "start"])
        annotated_dmrs$end_cpg <- paste0(chr_levels[beta_locs[end_cpgs, "chr"]], ":", beta_locs[end_cpgs, "start"])
        annotated_dmrs$upstream_cpgs <- sapply(seq_len(nrow(annotated_dmrs)), function(i) {
            cpgs <- annotated_dmrs$upstream_cpgs[i]
            cpg_inds <- suppressWarnings(as.integer(.splitCsvValues(cpgs)))
            cpg_inds <- cpg_inds[!is.na(cpg_inds)]
            if (length(cpg_inds) == 0) {
                return("")
            }
            paste0(chr_levels[beta_locs[cpg_inds, "chr"]], ":", beta_locs[cpg_inds, "start"], collapse = ",")
        })
        annotated_dmrs$downstream_cpgs <- sapply(seq_len(nrow(annotated_dmrs)), function(i) {
            cpgs <- annotated_dmrs$downstream_cpgs[i]
            cpg_inds <- suppressWarnings(as.integer(.splitCsvValues(cpgs)))
            cpg_inds <- cpg_inds[!is.na(cpg_inds)]
            if (length(cpg_inds) == 0) {
                return("")
            }
            paste0(chr_levels[beta_locs[cpg_inds, "chr"]], ":", beta_locs[cpg_inds, "start"], collapse = ",")
        })
        annotated_dmrs$cpgs <- sapply(seq_len(nrow(annotated_dmrs)), function(i) {
            cpgs <- annotated_dmrs$cpgs[i]
            cpg_inds <- suppressWarnings(as.integer(.splitCsvValues(cpgs)))
            cpg_inds <- cpg_inds[!is.na(cpg_inds)]
            if (length(cpg_inds) == 0) {
                return("")
            }
            paste0(chr_levels[beta_locs[cpg_inds, "chr"]], ":", beta_locs[cpg_inds, "start"], collapse = ",")
        })
        annotated_dmrs$tmp_id <- seq_len(nrow(annotated_dmrs))
        annotated_dmrs$seeds_unlisted <- sapply(annotated_dmrs$seeds, function(seed_ids) {
            as.integer(unlist(strsplit(seed_ids, ",")))
        })
        # explode on seeds column to annotate
        seeds_exploded <- annotated_dmrs[, "tmp_id", drop = FALSE]
        seeds_exploded <- seeds_exploded[rep(seq_len(nrow(seeds_exploded)), times = sapply(annotated_dmrs$seeds_unlisted, length)), , drop = FALSE]
        seeds_exploded$seeds_unlisted <- as.integer(unlist(annotated_dmrs$seeds_unlisted))
        seeds_exploded$chr <- chr_levels[beta_locs[seeds_exploded$seeds_unlisted, "chr"]]
        seeds_exploded$pos <- beta_locs[seeds_exploded$seeds_unlisted, "start"]
        seeds_exploded$annotated_seed <- paste0(seeds_exploded$chr, ":", seeds_exploded$pos)
        seeds_annotated_list <- split(seeds_exploded$annotated_seed, seeds_exploded$tmp_id)
        annotated_dmrs$seeds <- sapply(seeds_annotated_list, function(x) {
            paste(x, collapse = ",")
        })
        annotated_dmrs$seeds_unlisted <- NULL
        annotated_dmrs$tmp_id <- NULL
        annotated_dmrs$id <- paste0(annotated_dmrs$chr, ":", beta_locs[start_cpgs, "start"], "-", beta_locs[end_cpgs, "start"])
    }
    annotated_dmrs_granges <- GenomicRanges::makeGRangesFromDataFrame(
        annotated_dmrs,
        keep.extra.columns = TRUE,
        seqnames.field = "chr",
        seqinfo = GenomeInfoDb::Seqinfo(genome = genome),
        na.rm = TRUE
    )

    if (annotate_with_genes) {
        .log_step("Annotating DMRs with gene information...", level = 1)
        annotated_dmrs_granges <- annotateDMRsWithGenes(annotated_dmrs_granges, genome = genome)
        .log_success("DMR annotation completed.", level = 1)
    }

    if (!bed_provided) {
        .log_step("Converting DMR start/end CpG indices to be relative to all CpGs...", level = 2)
        sorted_locs <- beta_handler$getGenomicLocs()
        start_cpg_inds <- mcols(annotated_dmrs_granges)$start_cpg_ind
        end_cpg_inds <- mcols(annotated_dmrs_granges)$end_cpg_ind

        mcols(annotated_dmrs_granges)$start_cpg_ind_abs <- match(rownames(beta_locs)[start_cpg_inds], rownames(sorted_locs))
        mcols(annotated_dmrs_granges)$end_cpg_ind_abs <- match(rownames(beta_locs)[end_cpg_inds], rownames(sorted_locs))
        all_seeds_inds <- as.numeric(unlist(strsplit(as.character(mcols(annotated_dmrs_granges)$seeds_inds), ","), use.names = FALSE))
        all_seeds_abs <- match(rownames(beta_locs)[all_seeds_inds], rownames(sorted_locs))
        seeds_lengths <- lengths(strsplit(as.character(mcols(annotated_dmrs_granges)$seeds_inds), ","))
        seeds_groups <- rep(seq_along(seeds_lengths), seeds_lengths)
        mcols(annotated_dmrs_granges)$seeds_inds_abs <- vapply(split(all_seeds_abs, seeds_groups), function(x) paste(x, collapse = ","), character(1), USE.NAMES = FALSE)

        .log_success("DMR CpG indices conversion completed.", level = 2)
    } else {
        mcols(annotated_dmrs_granges)$start_cpg_ind_abs <- mcols(annotated_dmrs_granges)$start_cpg_ind
        mcols(annotated_dmrs_granges)$end_cpg_ind_abs <- mcols(annotated_dmrs_granges)$end_cpg_ind
        mcols(annotated_dmrs_granges)$seeds_inds_abs <- mcols(annotated_dmrs_granges)$seeds_inds
    }
    final_dmrs_granges <- annotated_dmrs_granges
    if (rank_dmrs) {
        .log_step("Ranking DMRs...", level = 1)
        ranked_dmrs_granges <- rankDMRs(
            final_dmrs_granges,
            beta = beta,
            pheno = pheno_all,
            genome = genome,
            array = array,
            sorted_locs = sorted_locs,
            sample_group_col = sample_group_col,
            covariates = covariates
        )
        .log_success("DMR ranking completed.", level = 1)
        final_dmrs_granges <- ranked_dmrs_granges
    }


    final_dmrs <- as.data.frame(final_dmrs_granges)
    colnames(final_dmrs)[colnames(final_dmrs) == "seqnames"] <- "chr"
    .log_info("Final number of DMRs: ", nrow(final_dmrs), level = 1)
    .log_info("Summary of final DMRs:\n\t", paste(capture.output(summary(final_dmrs)), collapse = "\n\t"), level = 3)
    if (!is.null(output_prefix)) {
        dmrs_file <- paste0(output_prefix, "dmrs.tsv.gz")
        .log_step("Saving DMRs to ", dmrs_file, "..", level = 2)
        gz <- gzfile(dmrs_file, "w")
        write.table(
            final_dmrs,
            gz,
            sep = "\t",
            quote = FALSE,
            qmethod = "double",
            col.names = TRUE,
            row.names = FALSE
        )
        close(gz)
        .log_success("DMRs saved.", level = 2)
    }

    gc()
    invisible(final_dmrs_granges)
}
