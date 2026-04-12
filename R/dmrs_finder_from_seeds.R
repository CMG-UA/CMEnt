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
#' @param beta Character. Path to the beta value file, or a tabix file, or a beta matrix, or a BetaHandler object, or a bed file, or a BSseq object. If a bed file is provided, it must at least contain bed_chrom_col and bed_chrom_start, followed by samples names in the given pheno, with corresponging beta values, and it will be converted to a tabix-indexed beta file internall, with the locations separately saved and queried as a DelayedDataFrame. If a BSseq object is provided, genomic locations and methylation values will be extracted using bsseq methods.
#' @param seeds Character. Path to the seeds TSV file or the seeds dataframe, in a format like the one produced by dmpFinder.
#' @param pheno Data frame. Phenotype data.
#' @param seeds_id_col Character. Column name or index for Seed identifiers in the seeds TSV file. Default is 1.
#' @param sample_group_col Character. Column name for sample group information in the phenotype data. Default is NULL.
#' @param casecontrol_col Boolean Column in pheno for case (TRUE/1) / control (FALSE/0) status . If NULL, controls will be assumed to be the first level of sample_group_col. Default is NULL.
#' @param covariates Character vector of column names in pheno to adjust for (e.g. "age", "sex"). When provided, correlations are computed on residuals after regressing M-values on these covariates within each group
#' @param min_cpg_delta_beta Numeric. Minimum delta beta value for CpGs. Default is 0.1.
#' @param expansion_step Numeric. Index-specific step size for expanding DMRs. Increasing it means higher memory usage and faster computation. Default is 500.
#' @param array Character. Type of array used (e.g., "450K", "EPIC", "EPICv2", "27K"). Ignored if using mm10 genome.
#' @param genome Character. Genome version (e.g., "hg38", "hg19", "hs1", "mm10"). Default is NULL and inferred as "hg19" for 450K, 27K, and EPIC arrays, otherwise "hg38".
#' @param max_pval Numeric. Maximum p-value to assume seeds correlation is significant. Default is 0.05.
#' @param pval_mode Character. "auto" (default) selects between t-based correlation p-values and empirical p-values per sample group using data diagnostics. You can also force "parametric" for t-based correlation p-values or "empirical" for permutation-based p-values.
#' @param empirical_strategy Character. When pval_mode = "empirical": "auto" (default) uses Monte Carlo for groups with <6 samples and permutations for groups with >=6 samples; "montecarlo" always uses Monte Carlo; "permutations" always uses permutations.
#' @param ntries Integer. Number of permutations when pval_mode = "empirical". Default is 0 (disabled).
#' @param max_lookup_dist Numeric. Maximum distance to look up for adjacent seeds belonging to the same DMR during Stage 1. Default is 10000 (10 kb).
#' @param expansion_window Numeric. Stage 2 connectivity is computed only in windows centered on seed-derived Stage 1 DMR neighborhoods, with this total window width in bp. Set <=0 for genome-wide connectivity. Default is -1 for microarrays and 10000 (10 kb) for NGS datasets.
#' @param max_bridge_seeds_gaps Integer. Maximum number of consecutive failed seed-to-seed edges to bridge during Stage 1 when both flanking edges are connected and failures are p-value driven. Set to 0 to disable. Default is 1.
#' @param max_bridge_extension_gaps Integer. Maximum gap size to consider during Stage 2 extension. Default is 1 (i.e., at most 1 consecutive failing CpG to bridge).
#' @param min_seeds Numeric. Minimum number of connected seeds in a DMR. Minimum is 2. Default is 2.
#' @param min_adj_seeds Numeric. Minimum number of seeds, adjusted by array CpG density, in a DMR after extension. Minimum is 2. Default is 2.
#' @param min_cpgs Numeric. Minimum number of CpGs in a DMR after extension, including the seeds. Minimum is 2. Default is 3.
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
#' @importFrom bedr tabix
#' @importFrom BSgenome getSeq
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rtracklayer import.chain
#' @importFrom GenomeInfoDb Seqinfo seqnames
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
.buildConnectivityWindowsFromDMRs <- function(dmrs, expansion_window) {
    if (is.null(dmrs) || nrow(dmrs) == 0L || !is.finite(expansion_window) || expansion_window <= 0) {
        return(data.frame(chr = character(0), start = numeric(0), end = numeric(0)))
    }
    flank <- as.numeric(expansion_window) / 2
    windows <- data.frame(
        chr = as.character(dmrs$chr),
        start = pmax(1, as.numeric(dmrs$start_seed_pos) - flank),
        end = as.numeric(dmrs$end_seed_pos) + flank
    )
    .mergeGenomicWindows(windows)
}


#' @keywords internal
#' @noRd
.normalizeFindDMRsArray <- function(array) {
    if (is.null(array)) {
        return(NULL)
    }
    array <- as.character(array)
    if (length(array) == 0L) {
        return(NULL)
    }
    if (length(array) > 1L) {
        array <- array[1]
    }
    if (is.na(array) || !nzchar(array) || tolower(array) == "null") {
        return(NULL)
    }
    supported_arrays <- c("450K", "27K", "EPIC", "EPICv2")
    matched_idx <- match(tolower(array), tolower(supported_arrays))
    if (is.na(matched_idx)) {
        stop(
            "Unsupported array: ", array,
            ". Supported arrays are: ", paste(c(supported_arrays, "NULL"), collapse = ", ")
        )
    }
    supported_arrays[matched_idx]
}


#' @keywords internal
#' @noRd
.normalizeFindDMRsGenome <- function(genome) {
    if (is.null(genome)) {
        return(NULL)
    }
    genome <- as.character(genome)
    if (length(genome) == 0L) {
        return(NULL)
    }
    if (length(genome) > 1L) {
        genome <- genome[1]
    }
    if (is.na(genome) || !nzchar(genome) || tolower(genome) == "null") {
        return(NULL)
    }
    supported_genomes <- c("hg38", "hg19", "hs1", "mm10", "mm39")
    matched_idx <- match(tolower(genome), tolower(supported_genomes))
    if (is.na(matched_idx)) {
        stop(
            "Unsupported genome: ", genome,
            ". Supported genomes are: ", paste(supported_genomes, collapse = ", ")
        )
    }
    supported_genomes[matched_idx]
}


#' @keywords internal
#' @noRd
.resolveFindDMRsGenome <- function(beta, array = NULL, genome = NULL, bed_provided = FALSE) {
    genome <- .normalizeFindDMRsGenome(genome)
    if (!is.null(genome)) {
        return(genome)
    }

    array <- .normalizeFindDMRsArray(array)

    if (inherits(beta, "BetaHandler")) {
        beta_genome <- tryCatch(.normalizeFindDMRsGenome(beta$genome), error = function(e) NULL)
        if (!is.null(beta_genome)) {
            return(beta_genome)
        }
        beta_array <- tryCatch(.normalizeFindDMRsArray(beta$array), error = function(e) NULL)
        if (!is.null(beta_array) && beta_array %in% c("450K", "27K", "EPIC")) {
            return("hg19")
        }
        return("hg38")
    }

    if (is_bsseq(beta)) {
        gr <- GenomicRanges::granges(beta)
        bsseq_genome <- unique(as.character(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(gr))))
        bsseq_genome <- bsseq_genome[!is.na(bsseq_genome) & nzchar(bsseq_genome)]
        bsseq_genome <- bsseq_genome[tolower(bsseq_genome) %in% c("hg38", "hg19", "hs1", "mm10", "mm39")]
        if (length(bsseq_genome) > 0L) {
            return(.normalizeFindDMRsGenome(bsseq_genome[1]))
        }
        return("hg38")
    }

    is_bed_input <- isTRUE(bed_provided) ||
        (is.character(beta) && length(beta) == 1L && file.exists(beta) && identical(tolower(tools::file_ext(beta)), "bed"))
    if (is_bed_input) {
        return("hg38")
    }

    if (!is.null(array) && array %in% c("450K", "27K", "EPIC")) {
        return("hg19")
    }

    "hg38"
}


#' @keywords internal
#' @noRd
.subsetStage2BetaToWindows <- function(beta_handler, beta_locs, col_names, expansion_windows, njobs = 1L) {
    wins <- .mergeGenomicWindows(expansion_windows)
    if (nrow(wins) == 0L) {
        return(NULL)
    }
    n_sites <- nrow(beta_locs)
    if (n_sites == 0L) {
        return(NULL)
    }
    beta_chr <- as.character(beta_locs[, "chr"])
    beta_start <- as.integer(beta_locs[, "start"])
    keep <- rep(FALSE, n_sites)
    for (chr in intersect(unique(wins$chr), unique(beta_chr))) {
        chr_inds <- which(beta_chr == chr)
        if (length(chr_inds) == 0L) {
            next
        }
        chr_wins <- wins[wins$chr == chr, , drop = FALSE]
        if (nrow(chr_wins) == 0L) {
            next
        }
        win_ir <- IRanges::IRanges(
            start = as.integer(chr_wins$start),
            end = as.integer(chr_wins$end)
        )
        site_ir <- IRanges::IRanges(
            start = beta_start[chr_inds],
            width = 1L
        )
        ov <- IRanges::findOverlaps(site_ir, win_ir)
        if (length(ov) > 0L) {
            keep[chr_inds[unique(S4Vectors::queryHits(ov))]] <- TRUE
        }
    }
    subset_to_full_idx <- which(keep)
    if (length(subset_to_full_idx) == 0L) {
        return(NULL)
    }
    full_to_subset_idx <- rep(NA_integer_, n_sites)
    full_to_subset_idx[subset_to_full_idx] <- seq_along(subset_to_full_idx)

    subset_locs <- as.data.frame(beta_locs[subset_to_full_idx, , drop = FALSE])
    subset_ids <- rownames(beta_locs)[subset_to_full_idx]
    rownames(subset_locs) <- subset_ids

    subset_beta <- beta_handler$getBeta(
        row_names = subset_ids,
        col_names = col_names
    )
    if (!inherits(subset_beta, "DelayedDataFrame")) {
        subset_beta <- DelayedDataFrame::DelayedDataFrame(subset_beta)
    }
    subset_beta_ddf <- subset_beta
    rownames(subset_beta_ddf) <- subset_ids
    if (!is.null(col_names)) {
        colnames(subset_beta_ddf) <- col_names
    }

    subset_beta_handler <- getBetaHandler(
        beta = subset_beta_ddf,
        sorted_locs = subset_locs,
        njobs = njobs
    )
    list(
        beta_handler = subset_beta_handler,
        beta_locs = subset_locs,
        expansion_windows = wins,
        subset_to_full_idx = subset_to_full_idx,
        full_to_subset_idx = full_to_subset_idx
    )
}


.chooseTestingOptions <- function(group, mat, mask, pval_mode, empirical_strategy) {
    n_sites <- nrow(mat)
    x_mat_full <- mat[1:(n_sites - 1), , drop = FALSE] # Sites i
    y_mat_full <- mat[2:n_sites, , drop = FALSE] # Sites i+1
    x_mat <- x_mat_full[mask, , drop = FALSE]
    y_mat <- y_mat_full[mask, , drop = FALSE]
    valid_pairs <- !is.na(x_mat) & !is.na(y_mat)
    n_valid <- rowSums(valid_pairs)

    if (pval_mode == "auto") {
        auto_diag <- .summarizeCorrelationAssumptions(
            x_mat = x_mat,
            y_mat = y_mat,
            n_valid = n_valid
        )
        pval_mode <- if (auto_diag$use_parametric) "parametric" else "empirical"
        .log_info(
            "Auto-selected pval_mode='", pval_mode, "' for group '", group,
            "' (q10_n_valid=", signif(auto_diag$q10_n_valid, 3),
            ", median|pearson-spearman|=", signif(auto_diag$median_abs_delta_spearman, 3),
            ", median|pearson-winsorized|=", signif(auto_diag$median_abs_delta_winsorized, 3),
            ", median|skew|=", signif(auto_diag$median_abs_skew, 3),
            ", median_excess_kurtosis=", signif(auto_diag$median_excess_kurtosis, 3),
            ", pilot_pairs=", auto_diag$n_pairs_used,
            ", reason=", auto_diag$reason, ").",
            level = 2
        )
    }
    if (pval_mode == "empirical" && empirical_strategy == "auto") {
        if (ncol(x_mat) >= 6) {
            .log_info(
                "Group '", group, "' has ", ncol(x_mat), " samples. Using 'montecarlo' empirical strategy for faster computation with sufficient sample size.",
                level = 2
            )
            empirical_strategy <- "montecarlo"
        } else {
            .log_info(
                "Group '", group, "' has only ", ncol(x_mat), " samples. Using 'permutations' empirical strategy for more accurate p-value estimation with small sample size.",
                level = 2
            )
            empirical_strategy <- "permutations"
        }
    }
    list(pval_mode = pval_mode, empirical_strategy = empirical_strategy)
}


#' @keywords internal
#' @noRd
.connectivityChunkWorker <- function(
    split,
    beta_handler,
    beta_chr_vec,
    beta_start_vec,
    group_inds,
    pheno,
    pval_mode_per_group,
    empirical_strategy_per_group,
    col_names = NULL,
    max_pval = 0.05,
    min_delta_beta = 0,
    covariates = NULL,
    max_lookup_dist = 1000,
    entanglement = "strong",
    aggfun = median,
    ntries = 500,
    mid_p = TRUE,
    checked_pairs = NULL,
    use_numeric_row_index = TRUE,
    beta_row_ids = NULL,
    beta_row_ids_offset = 0L
) {
    pair_start <- as.integer(split[[1]])
    pair_end <- as.integer(split[[2]])
    site_start <- pair_start
    site_end <- pair_end + 1L
    inds <- site_start:site_end

    if (!is.null(checked_pairs)) {
        pair_mask <- checked_pairs$before >= site_start & checked_pairs$after <= site_end
        if (!any(pair_mask)) {
            return(list(pair_start = pair_start, pair_end = pair_end, result = data.frame()))
        }
        chunk_pairs <- checked_pairs[pair_mask, c("before", "after"), drop = FALSE]
        inds <- as.vector(t(as.matrix(chunk_pairs)))
    }

    locs <- data.frame(
        chr = beta_chr_vec[inds],
        start = beta_start_vec[inds],
        stringsAsFactors = FALSE
    )
    if (use_numeric_row_index) {
        row_ids <- inds
    } else {
        local_inds <- as.integer(inds - as.integer(beta_row_ids_offset))
        if (any(local_inds < 1L) || any(local_inds > length(beta_row_ids))) {
            stop(
                "Row index offset mismatch while resolving batch row IDs: [",
                min(local_inds), ", ", max(local_inds), "] outside [1, ",
                length(beta_row_ids), "]"
            )
        }
        row_ids <- beta_row_ids[local_inds]
    }
    chunk_beta <- beta_handler$getBeta(
        row_names = row_ids,
        col_names = col_names
    )
    x <- .testConnectivityBatch(
        sites_beta = as.matrix(chunk_beta),
        group_inds = group_inds,
        pheno = pheno,
        covariates = covariates,
        max_pval = max_pval,
        min_delta_beta = min_delta_beta,
        max_lookup_dist = max_lookup_dist,
        sites_locs = locs,
        entanglement = entanglement,
        aggfun = aggfun,
        pval_mode_per_group = pval_mode_per_group,
        empirical_strategy_per_group = empirical_strategy_per_group,
        ntries = ntries,
        mid_p = mid_p,
        check_non_overlapping = !is.null(checked_pairs)
    )
    rm(chunk_beta)
    if (!is.null(checked_pairs) && length(inds) > 0L) {
        # attach to result so outer loop can map back exactly
        # keep only mod 2 = 1
        inds <- inds[seq_along(inds) %% 2 == 1]
        attr(x, "recomputed_pairs") <- inds
    }
    list(pair_start = pair_start, pair_end = pair_end, result = x)
}


#' @keywords internal
#' @noRd
.buildConnectivityArraySinglePass <- function(
    beta_handler,
    beta_locs = NULL,
    pheno,
    group_inds,
    pval_mode_per_group,
    empirical_strategy_per_group,
    col_names = NULL,
    max_pval = 0.05,
    min_delta_beta = 0,
    covariates = NULL,
    max_lookup_dist = 1000,
    chunk_size = getOption("DMRsegal.chunk_size", 1000),
    entanglement = "strong",
    aggfun = median,
    ntries = 500,
    mid_p = TRUE,
    njobs = 1,
    expansion_windows = NULL,
    connectivity_array = NULL,
    ugap = 0L,
    dgap = 0L,
    recheck = NULL,
    splits = NULL
) {
    if (is.null(beta_locs)) {
        beta_locs <- beta_handler$getBetaLocs()
    }
    if (ugap > 0L || dgap > 0L) {
        if (is.null(connectivity_array) || is.null(splits)) {
            stop("ugap and dgap parameters require providing the current connectivity_array and splits data frames.")
        }
        if (ugap > 0L && dgap > 0L) {
            stop("The bridging is either upstream or downstream, but not both at the same time.")
        }
    }
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
        return(list(connectivity_array = ret, splits = NULL, pval_mode_per_group = pval_mode_per_group, empirical_strategy_per_group = empirical_strategy_per_group))
    }
    beta_chr_factor <- as.factor(beta_locs[, "chr"])
    beta_chr_levels <- levels(beta_chr_factor)
    beta_chr_ids <- as.integer(beta_chr_factor)
    chr_runs <- rle(beta_chr_ids)
    chr_run_ends <- cumsum(chr_runs$lengths)
    chr_run_starts <- chr_run_ends - chr_runs$lengths + 1L
    chr_ends <- as.integer(chr_run_ends)
    window_mode <- !is.null(expansion_windows) && nrow(expansion_windows) > 0L
    default_reason <- if (window_mode) "outside_connectivity_window" else ""

    # Estimate a safe upper bound for pair chunks to avoid worker OOM when
    # users request very large chunk_size values on sequencing-scale inputs.
    n_cols_for_chunk <- if (!is.null(col_names)) {
        length(col_names)
    } else {
        length(beta_handler$getBetaColNames())
    }
    chunk_mem_mb <- getOption("DMRsegal.max_chunk_memory_mb", 256)
    if (!is.numeric(chunk_mem_mb) || length(chunk_mem_mb) != 1L || is.na(chunk_mem_mb) || chunk_mem_mb <= 0) {
        chunk_mem_mb <- Inf
    }
    chunk_mem_multiplier <- getOption("DMRsegal.chunk_memory_multiplier", 12)
    if (!is.numeric(chunk_mem_multiplier) || length(chunk_mem_multiplier) != 1L || is.na(chunk_mem_multiplier) || chunk_mem_multiplier <= 0) {
        chunk_mem_multiplier <- 12
    }
    bytes_per_pair_est <- as.numeric(8 * max(1L, n_cols_for_chunk) * chunk_mem_multiplier)
    max_pairs_by_memory <- if (is.finite(chunk_mem_mb)) {
        as.integer(floor((as.numeric(chunk_mem_mb) * 1024^2) / max(1, bytes_per_pair_est)))
    } else {
        .Machine$integer.max
    }
    if (!is.finite(max_pairs_by_memory) || is.na(max_pairs_by_memory)) {
        max_pairs_by_memory <- .Machine$integer.max
    }
    max_pairs_by_memory <- max(100L, max_pairs_by_memory)

    requested_chunk_size <- suppressWarnings(as.integer(chunk_size))
    if (is.finite(requested_chunk_size) && requested_chunk_size > max_pairs_by_memory) {
        .log_info(
            "Capping chunk_size from ", requested_chunk_size,
            " to ", max_pairs_by_memory,
            " to fit the connectivity memory budget (",
            as.numeric(chunk_mem_mb), " MB per worker, ", n_cols_for_chunk,
            " sample columns, multiplier=", chunk_mem_multiplier, ").",
            level = 2
        )
    }

    .applyChunkSizeSafetyCap <- function(chunk_size_eff) {
        chunk_size_eff <- as.integer(chunk_size_eff)
        if (!is.finite(chunk_size_eff) || is.na(chunk_size_eff) || chunk_size_eff < 1L) {
            chunk_size_eff <- 1L
        }
        min(chunk_size_eff, max_pairs_by_memory)
    }

    .makeOutputTemplate <- function(nrows, reason_default) {
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

    .chunkPairRanges <- function(pair_ranges_df) {
        if (is.null(pair_ranges_df) || nrow(pair_ranges_df) == 0L) {
            return(matrix(numeric(0), ncol = 2))
        }
        chunk_size_eff <- as.integer(chunk_size)
        total_pairs <- sum(pair_ranges_df$end_pair - pair_ranges_df$start_pair + 1L)
        if (njobs > 1L && is.finite(total_pairs) && total_pairs > 0L) {
            # Keep at least one chunk per worker even when chunk_size is set very large.
            max_chunk_size_for_parallel <- as.integer(ceiling(total_pairs / njobs))
            chunk_size_eff <- min(chunk_size_eff, max_chunk_size_for_parallel)
        }
        chunk_size_eff <- .applyChunkSizeSafetyCap(chunk_size_eff)
        out <- vector("list", nrow(pair_ranges_df) * 2L)
        out_n <- 0L
        for (i in seq_len(nrow(pair_ranges_df))) {
            ps <- as.integer(pair_ranges_df$start_pair[i])
            pe <- as.integer(pair_ranges_df$end_pair[i])
            if (pe < ps) {
                next
            }
            for (chunk_ps in seq(ps, pe, by = chunk_size_eff)) {
                chunk_pe <- min(chunk_ps + chunk_size_eff - 1L, pe)
                out_n <- out_n + 1L
                out[[out_n]] <- c(chunk_ps, chunk_pe)
            }
        }
        do.call(rbind, out[seq_len(out_n)])
    }

    .poolSplits <- function(splits_mat, split_weights = NULL) {
        if (is.null(splits_mat) || nrow(splits_mat) <= 1L) {
            return(splits_mat)
        }
        if (is.null(split_weights)) {
            split_weights <- as.integer(splits_mat[, 2] - splits_mat[, 1] + 1L)
        } else {
            split_weights <- as.integer(split_weights)
        }
        split_weights[is.na(split_weights) | split_weights < 1L] <- 1L

        chunk_size_eff <- as.integer(chunk_size)
        if (!is.finite(chunk_size_eff) || chunk_size_eff < 1L) {
            chunk_size_eff <- 1L
        }
        if (njobs > 1L) {
            max_chunk_size_for_parallel <- as.integer(ceiling(sum(split_weights) / njobs))
            chunk_size_eff <- min(chunk_size_eff, max(1L, max_chunk_size_for_parallel))
        }
        chunk_size_eff <- .applyChunkSizeSafetyCap(chunk_size_eff)
        if (chunk_size_eff <= 1L) {
            return(splits_mat)
        }

        split_groups <- (cumsum(split_weights) - 1L) %/% chunk_size_eff
        split_idx_groups <- split(seq_len(nrow(splits_mat)), split_groups)
        pooled <- lapply(
            split_idx_groups,
            function(ix) {
                c(
                    min(splits_mat[ix, 1]),
                    max(splits_mat[ix, 2])
                )
            }
        )
        pooled <- do.call(rbind, pooled)
        if (is.null(dim(pooled))) {
            pooled <- matrix(pooled, ncol = 2L, byrow = TRUE)
        }
        storage.mode(pooled) <- "integer"
        pooled
    }

    .buildAllPairRanges <- function() {
        keep_runs <- chr_runs$lengths >= 2L
        if (!any(keep_runs)) {
            return(data.frame(start_pair = integer(0), end_pair = integer(0)))
        }
        data.frame(
            start_pair = as.integer(chr_run_starts[keep_runs]),
            end_pair = as.integer(chr_run_ends[keep_runs] - 1L)
        )
    }

    .build_window_pair_ranges <- function() {
        wins <- .mergeGenomicWindows(expansion_windows)
        if (nrow(wins) == 0L) {
            return(data.frame(start_pair = integer(0), end_pair = integer(0)))
        }
        pair_ranges <- vector("list", length(unique(beta_chr_ids)))
        pair_n <- 0L
        beta_start <- as.integer(beta_locs[, "start"])
        win_chr_ids <- match(unique(wins$chr), beta_chr_levels)
        win_chr_ids <- unique(win_chr_ids[!is.na(win_chr_ids)])
        for (chr_id in win_chr_ids) {
            chr_inds <- which(beta_chr_ids == chr_id)
            if (length(chr_inds) < 2L) {
                next
            }
            chr_wins <- wins[wins$chr == beta_chr_levels[[chr_id]], , drop = FALSE]
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

    if (is.null(connectivity_array)) {
        connectivity_array <- .makeOutputTemplate(n_sites, default_reason)
        connectivity_array[chr_ends, "connected"] <- FALSE
        connectivity_array[chr_ends, "reason"] <- "end-of-input"
        checked_pairs <- NULL
        pair_ranges <- if (window_mode) .build_window_pair_ranges() else .buildAllPairRanges()
        if (nrow(pair_ranges) == 0L) {
            return (
                list(
                    connectivity_array = connectivity_array,
                    splits = NULL,
                    pval_mode_per_group = pval_mode_per_group,
                    empirical_strategy_per_group = empirical_strategy_per_group
                )
            )
        }
        connectivity_array[pair_ranges$start_pair[pair_ranges$start_pair > 1] - 1, "reason"] <- "end-of-input"
        connectivity_array[pair_ranges$end_pair, "reason"] <- "end-of-input"
        splits <- .chunkPairRanges(pair_ranges)
        if (window_mode && nrow(splits) > 1L) {
            old_n <- nrow(splits)
            splits <- .poolSplits(splits)
            .log_info(
                "Pooled initial window chunks from ", old_n, " to ", nrow(splits),
                " (target chunk size: ", as.integer(chunk_size), ").",
                level = 3
            )
        }
    } else {
        # If connectivity array is provided, we assume it has already been filled for all sites up to windows/chromosomes ends.
        # Instead, we re-assess the disconnected sites on the edges of the connected regions,
        # comparing i with i + 2 instead of i with i + 1, to see if we can connect them by bridging the gap of one site.
        # The connectivity array is then updated with the bridged connections.
        connected <- connectivity_array[, "connected"]
        runs <- rle(connected)
        run_ends <- cumsum(runs$lengths)
        run_starts <- run_ends - runs$lengths + 1L
        values <- runs$values
        run_mask <- values == 1
        if (!is.null(recheck)) {
            recheck_inds <- which(recheck)
            run_mask <- run_mask & run_starts %in% recheck_inds
        }
        run_ends <- run_ends[run_mask]
        run_starts <- run_starts[run_mask]
        run_mask <- rep(TRUE, length(run_starts))

        if (ugap > 0L) {
            boundary_mask <- run_starts - ugap > 0L
            run_starts <- run_starts[boundary_mask]
            run_ends <- run_ends[boundary_mask]
            run_mask <- rep(TRUE, length(run_starts))
            for (i in seq_len(ugap)) {
                run_mask <- run_mask & connectivity_array[run_starts - i, "reason"] != "end-of-input"
            }
        }
        if (dgap > 0L) {
            boundary_mask <- run_ends + dgap <= n_sites
            run_starts <- run_starts[boundary_mask]
            run_ends <- run_ends[boundary_mask]
            run_mask <- rep(TRUE, length(run_starts))
            for (i in seq_len(dgap)) {
                run_mask <- run_mask & connectivity_array[run_ends + i, "reason"] != "end-of-input"
            }
        }
        if (!any(run_mask)) {
            # No runs to bridge, return the existing connectivity array and splits
            return(list(connectivity_array = connectivity_array, splits = splits, pval_mode_per_group = pval_mode_per_group, empirical_strategy_per_group = empirical_strategy_per_group))
        }
        run_ends <- run_ends[run_mask]
        run_starts <- run_starts[run_mask]
        if (ugap > 0L) {
            # The following indices will be re-checked
            checked_pairs <- data.frame(before = run_starts - ugap, after = run_starts)
            .log_info(
                "Re-assessing connectivity for ", nrow(checked_pairs), " site pairs at the upstream edges of existing connected regions to see if we can bridge small gaps.",
                level = 3
            )
        } else if (dgap > 0L) {
            checked_pairs <- data.frame(before = run_ends + 1, after = run_ends + dgap + 1)
            .log_info(
                "Re-assessing connectivity for ", nrow(checked_pairs), " site pairs at the downstream edges of existing connected regions to see if we can bridge small gaps.",
                level = 3
            )
        }

        # Updating the splits to include only chunks with complete re-check pairs.
        ninds_in_splits <- apply(splits, 1, function(r) {
            site_start <- as.integer(r[1])
            site_end <- as.integer(r[2]) + 1L
            sum(checked_pairs$before >= site_start & checked_pairs$after <= site_end)
        })
        splits <- splits[ninds_in_splits > 0L, , drop = FALSE]
        ninds_in_splits <- ninds_in_splits[ninds_in_splits > 0L]
        # accumulate bridge re-check chunks to keep a considerable chunk size
        if (nrow(splits) > 1L) {
            old_n <- nrow(splits)
            splits <- .poolSplits(splits, split_weights = ninds_in_splits)
            .log_info(
                "Pooled bridge re-check chunks from ", old_n, " to ", nrow(splits),
                " (target chunk size: ", as.integer(chunk_size), ").",
                level = 3
            )
        }
    }
    .log_info(
        "Connectivity computation mode: ",
        if (window_mode) "windowed" else "genome-wide",
        "; chunks to evaluate: ", nrow(splits), ".",
        level = 3
    )
    verbose <- getOption("DMRsegal.verbose", 1)
    p_ext <- NULL
    if (verbose > 0) {
        p_ext <- progressr::progressor(steps = nrow(splits), message = "Computing connectivity array...")
    }
    beta_row_ids_full <- rownames(beta_locs)
    # Numeric indices in this function are relative to beta_locs, which may be a subset
    # of the full beta matrix. Prefer stable row IDs whenever they are available so that
    # all backends query the same CpGs during connectivity estimation.
    use_numeric_row_index <- is.null(beta_row_ids_full) ||
        length(beta_row_ids_full) != n_sites ||
        anyNA(beta_row_ids_full) ||
        any(!nzchar(beta_row_ids_full))
    beta_row_ids <- if (use_numeric_row_index) NULL else beta_row_ids_full
    if (any(pval_mode_per_group == "auto") || any(empirical_strategy_per_group[pval_mode_per_group == "empirical"] == "auto")) {
        .log_info(
            "Selecting p-value computation mode for each group using the first chunk as a pilot.",
            level = 2
        )
        first_chunk_rows <- splits[1, 1]:(splits[1, 2] + 1L)
        if (!use_numeric_row_index) {
            first_chunk_rows <- beta_row_ids_full[first_chunk_rows]
        }
        # select testing settings using the first chunk as a pilot
        first_chunk <- beta_handler$getBeta(
            row_names = first_chunk_rows,
            col_names = col_names
        )
        sites_locs <- as.data.frame(beta_locs[splits[1, 1]:(splits[1, 2] + 1L), , drop = FALSE])
        s <- nrow(sites_locs)
        if (!is.null(max_lookup_dist) && !is.null(sites_locs)) {
            dists <- as.numeric(sites_locs[2:s, "start"]) - as.numeric(sites_locs[1:(s - 1), "start"])
            exceeded_dist <- dists > max_lookup_dist | sites_locs[2:s, "chr"] != sites_locs[1:(s - 1), "chr"]
        } else {
            exceeded_dist <- rep(FALSE, n_sites - 1L)
        }
        nexdist_mask <- !exceeded_dist
        groups_options <- lapply(
            names(group_inds),
            function(x) {
                idx <- group_inds[[x]]
                chunk_m <- .transformBeta(first_chunk[, idx, drop = FALSE], pheno = pheno[idx, , drop = FALSE], covariates = covariates)
                .chooseTestingOptions(
                    group = x,
                    mat = chunk_m,
                    mask = nexdist_mask,
                    pval_mode = pval_mode_per_group[[x]],
                    empirical_strategy = empirical_strategy_per_group[[x]]
                )
            }
        )
        empirical_strategy_per_group <- sapply(groups_options, function(opt) opt$empirical_strategy)
        names(empirical_strategy_per_group) <- names(group_inds)
        pval_mode_per_group <- sapply(groups_options, function(opt) opt$pval_mode)
        names(pval_mode_per_group) <- names(group_inds)
        rm(first_chunk, sites_locs, exceeded_dist, nexdist_mask, groups_options)
    }

    gc()
    beta_chr_vec <- beta_chr_ids
    beta_start_vec <- as.integer(beta_locs[, "start"])
    .runConnectivityChunk <- function(
        split,
        checked_pairs_local = checked_pairs,
        worker_beta_handler = beta_handler,
        worker_use_numeric_row_index = use_numeric_row_index,
        worker_beta_row_ids = beta_row_ids,
        worker_beta_row_ids_offset = 0L
    ) {
        .connectivityChunkWorker(
            split = split,
            beta_handler = worker_beta_handler,
            beta_chr_vec = beta_chr_vec,
            beta_start_vec = beta_start_vec,
            group_inds = group_inds,
            pheno = pheno,
            pval_mode_per_group = pval_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            col_names = col_names,
            max_pval = max_pval,
            min_delta_beta = min_delta_beta,
            covariates = covariates,
            max_lookup_dist = max_lookup_dist,
            entanglement = entanglement,
            aggfun = aggfun,
            ntries = ntries,
            mid_p = mid_p,
            checked_pairs = checked_pairs_local,
            use_numeric_row_index = worker_use_numeric_row_index,
            beta_row_ids = worker_beta_row_ids,
            beta_row_ids_offset = worker_beta_row_ids_offset
        )
    }
    # Work with local vectors to avoid repeated data.frame copy-on-modify in the hot loop.
    connected_vec <- connectivity_array$connected
    pval_vec <- connectivity_array$pval
    reason_vec <- connectivity_array$reason
    fail_col <- if ("first_failing_group" %in% names(connectivity_array)) "first_failing_group" else if ("failing_groups" %in% names(connectivity_array)) "failing_groups" else NULL
    fail_vec <- if (!is.null(fail_col)) connectivity_array[[fail_col]] else NULL
    delta_vec <- if ("delta_beta" %in% names(connectivity_array)) connectivity_array$delta_beta else NULL

    bridge_mask <- rep(FALSE, n_sites)
    recheck <- rep(FALSE, n_sites)

    .applyChunkResult <- function(item) {
        if (!is.null(p_ext)) {
            p_ext()
        }
        x <- item$result
        if (nrow(x) == 0L) {
            return(invisible(NULL))
        }
        if (is.null(checked_pairs)) {
            # First pass: full overwrite
            idx <- item$pair_start:item$pair_end
            connected_vec[idx] <<- x$connected
            pval_vec[idx] <<- x$pval
            reason_vec[idx] <<- x$reason
            if (!is.null(fail_col) && fail_col %in% names(x)) {
                fail_vec[idx] <<- x[[fail_col]]
            }
            if (!is.null(delta_vec) && "delta_beta" %in% names(x)) {
                delta_vec[idx] <<- x$delta_beta
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
                recheck[masked_idx[update_m]] <- TRUE
                update_idx <- masked_idx[update_m]
                connected_vec[update_idx] <<- x$connected[update_m]
                pval_vec[update_idx] <<- x$pval[update_m]
                reason_vec[update_idx] <<- x$reason[update_m]
                if (!is.null(fail_col) && fail_col %in% names(x)) {
                    fail_vec[update_idx] <<- x[[fail_col]][update_m]
                }
                if (!is.null(delta_vec) && "delta_beta" %in% names(x)) {
                    delta_vec[update_idx] <<- x$delta_beta[update_m]
                }
                gap <- if (ugap > 0L) ugap else dgap
                for (i in 0:(gap - 1L)) {
                    bridge_mask[masked_idx[update_m] + i] <- TRUE
                }
            }
        }
        invisible(NULL)
    }

    .futureBatchConnectivity <- function(
        batch_splits,
        batch_checked_pairs = checked_pairs,
        batch_beta_handler = beta_handler,
        batch_use_numeric_row_index = use_numeric_row_index,
        batch_beta_row_ids = beta_row_ids,
        batch_beta_row_ids_offset = 0L
    ) {
        future.apply::future_lapply(
            X = batch_splits,
            future.seed = TRUE,
            future.stdout = NA,
            future.globals = c(
                ".connectivityChunkWorker",
                ".testConnectivityBatch"
            ),
            FUN = .connectivityChunkWorker,
            beta_handler = batch_beta_handler,
            beta_chr_vec = beta_chr_vec,
            beta_start_vec = beta_start_vec,
            group_inds = group_inds,
            pheno = pheno,
            pval_mode_per_group = pval_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            col_names = col_names,
            max_pval = max_pval,
            min_delta_beta = min_delta_beta,
            covariates = covariates,
            max_lookup_dist = max_lookup_dist,
            entanglement = entanglement,
            aggfun = aggfun,
            ntries = ntries,
            mid_p = mid_p,
            checked_pairs = batch_checked_pairs,
            use_numeric_row_index = batch_use_numeric_row_index,
            beta_row_ids = batch_beta_row_ids,
            beta_row_ids_offset = batch_beta_row_ids_offset
        )
    }

    .subsetCheckedPairsForBatch <- function(batch_inds) {
        if (is.null(checked_pairs) || length(batch_inds) == 0L) {
            return(checked_pairs)
        }
        batch_pair_start <- min(splits[batch_inds, 1])
        batch_pair_end <- max(splits[batch_inds, 2]) + 1L
        checked_pairs[
            checked_pairs$before >= batch_pair_start & checked_pairs$after <= batch_pair_end,
            c("before", "after"),
            drop = FALSE
        ]
    }

    batch_subset_enabled <- isTRUE(getOption("DMRsegal.parallel_batch_beta_subset", TRUE))
    batch_subset_min_rows <- as.integer(getOption("DMRsegal.parallel_batch_beta_subset_min_rows", 5000L))
    if (!is.finite(batch_subset_min_rows) || is.na(batch_subset_min_rows) || batch_subset_min_rows < 1L) {
        batch_subset_min_rows <- 1L
    }

    .resolveBatchWorkerPayload <- function(batch_inds) {
        payload <- list(
            beta_handler = beta_handler,
            use_numeric_row_index = use_numeric_row_index,
            beta_row_ids = beta_row_ids,
            beta_row_ids_offset = 0L
        )
        if (!batch_subset_enabled || length(batch_inds) == 0L) {
            return(payload)
        }
        if (use_numeric_row_index || is.null(beta_row_ids_full)) {
            return(payload)
        }

        batch_pair_start <- min(splits[batch_inds, 1])
        batch_pair_end <- max(splits[batch_inds, 2]) + 1L
        batch_nrows <- batch_pair_end - batch_pair_start + 1L
        if (!is.finite(batch_nrows) || is.na(batch_nrows) || batch_nrows < batch_subset_min_rows) {
            return(payload)
        }

        batch_row_ids <- beta_row_ids_full[batch_pair_start:batch_pair_end]
        subset_handler <- tryCatch(
            beta_handler$subset(row_names = batch_row_ids, col_names = col_names),
            error = function(e) e
        )
        if (inherits(subset_handler, "error")) {
            .log_warn(
                "Failed to create per-batch BetaHandler subset for chunk range ",
                batch_pair_start, "-", batch_pair_end,
                ". Proceeding with the full handler. Error: ",
                conditionMessage(subset_handler)
            )
            return(payload)
        }

        payload$beta_handler <- subset_handler
        payload$use_numeric_row_index <- FALSE
        payload$beta_row_ids <- batch_row_ids
        payload$beta_row_ids_offset <- batch_pair_start - 1L
        payload
    }

    .isRetryableFutureBatchError <- function(e) {
        if (inherits(e, "FutureInterruptError")) {
            return(TRUE)
        }
        msg <- conditionMessage(e)
        grepl(
            "did not deliver a result|was interrupted|multicorefuture|worker.*(died|terminated)|failed to retrieve|connection.*closed|sigkill",
            msg,
            ignore.case = TRUE
        )
    }

    if (njobs == 1L) {
        for (split_ind in seq_len(nrow(splits))) {
            .applyChunkResult(.runConnectivityChunk(splits[split_ind, ], checked_pairs_local = checked_pairs))
        }
    } else {
        .setupParallel()
        on.exit(.finalizeParallel(), add = TRUE)
        all_split_inds <- seq_len(nrow(splits))
        batch_size <- as.integer(getOption("DMRsegal.parallel_result_batch_size", max(njobs * 4L, 16L)))
        batch_size <- max(1L, batch_size)
        max_parallel_retries <- as.integer(getOption("DMRsegal.parallel_batch_retries", 1L))
        if (!is.finite(max_parallel_retries) || is.na(max_parallel_retries) || max_parallel_retries < 0L) {
            max_parallel_retries <- 1L
        }
        max_attempts <- max_parallel_retries + 1L
        fallback_to_sequential <- getOption("DMRsegal.parallel_fallback_to_sequential", TRUE)
        if (is.null(fallback_to_sequential)) {
            fallback_to_sequential <- TRUE
        }
        fallback_to_sequential <- isTRUE(fallback_to_sequential)
        use_parallel <- TRUE
        for (batch_start in seq.int(1L, length(all_split_inds), by = batch_size)) {
            batch_end <- min(batch_start + batch_size - 1L, length(all_split_inds))
            batch_inds <- all_split_inds[batch_start:batch_end]
            batch_splits <- lapply(batch_inds, function(i) as.integer(splits[i, ]))
            batch_checked_pairs <- .subsetCheckedPairsForBatch(batch_inds)
            batch_worker_payload <- .resolveBatchWorkerPayload(batch_inds)
            if (!use_parallel) {
                for (split in batch_splits) {
                    .applyChunkResult(
                        .runConnectivityChunk(
                            split,
                            checked_pairs_local = batch_checked_pairs,
                            worker_beta_handler = batch_worker_payload$beta_handler,
                            worker_use_numeric_row_index = batch_worker_payload$use_numeric_row_index,
                            worker_beta_row_ids = batch_worker_payload$beta_row_ids,
                            worker_beta_row_ids_offset = batch_worker_payload$beta_row_ids_offset
                        )
                    )
                }
                next
            }

            attempt <- 1L
            ret <- NULL
            repeat {
                ret <- tryCatch(
                    .futureBatchConnectivity(
                        batch_splits,
                        batch_checked_pairs = batch_checked_pairs,
                        batch_beta_handler = batch_worker_payload$beta_handler,
                        batch_use_numeric_row_index = batch_worker_payload$use_numeric_row_index,
                        batch_beta_row_ids = batch_worker_payload$beta_row_ids,
                        batch_beta_row_ids_offset = batch_worker_payload$beta_row_ids_offset
                    ),
                    error = function(e) e
                )
                if (!inherits(ret, "error")) {
                    break
                }

                retryable_error <- .isRetryableFutureBatchError(ret)
                if (!retryable_error) {
                    stop(ret)
                }

                if (attempt < max_attempts) {
                    .log_warn(
                        "Parallel connectivity batch ", batch_start, "-", batch_end,
                        " failed with recoverable future error (attempt ", attempt,
                        "/", max_attempts,
                        "). Reinitializing parallel backend and retrying. Error: ",
                        conditionMessage(ret)
                    )
                    .cleanupParallelState()
                    .setupParallel()
                    attempt <- attempt + 1L
                    next
                }

                if (fallback_to_sequential) {
                    .log_warn(
                        "Parallel connectivity batch ", batch_start, "-", batch_end,
                        " failed after ", max_attempts,
                        " attempt(s). Falling back to sequential processing for remaining chunks. Last error: ",
                        conditionMessage(ret)
                    )
                    use_parallel <- FALSE
                    .cleanupParallelState()
                    ret <- NULL
                    break
                }

                stop(ret)
            }

            if (!use_parallel && is.null(ret)) {
                for (split in batch_splits) {
                    .applyChunkResult(
                        .runConnectivityChunk(
                            split,
                            checked_pairs_local = batch_checked_pairs,
                            worker_beta_handler = batch_worker_payload$beta_handler,
                            worker_use_numeric_row_index = batch_worker_payload$use_numeric_row_index,
                            worker_beta_row_ids = batch_worker_payload$beta_row_ids,
                            worker_beta_row_ids_offset = batch_worker_payload$beta_row_ids_offset
                        )
                    )
                }
                next
            }

            for (item in ret) {
                .applyChunkResult(item)
            }
            rm(ret)
        }
    }

    connected_vec[bridge_mask] <- TRUE
    reason_vec[bridge_mask] <- "bridged"

    # Preserve hard window boundaries even when chunk pooling evaluates ranges spanning multiple windows.
    if (is.null(checked_pairs) && window_mode && exists("pair_ranges", inherits = FALSE) && nrow(pair_ranges) > 0L) {
        boundary_left <- pair_ranges$start_pair[pair_ranges$start_pair > 1L] - 1L
        boundary_right <- pair_ranges$end_pair
        boundary_idx <- unique(c(boundary_left, boundary_right))
        if (length(boundary_idx) > 0L) {
            connected_vec[boundary_idx] <- FALSE
            reason_vec[boundary_idx] <- "end-of-input"
        }
    }

    connectivity_array$connected <- connected_vec
    connectivity_array$pval <- pval_vec
    connectivity_array$reason <- reason_vec
    if (!is.null(fail_col)) {
        connectivity_array[[fail_col]] <- fail_vec
    }
    if (!is.null(delta_vec)) {
        connectivity_array$delta_beta <- delta_vec
    }

    list(
        connectivity_array = connectivity_array,
        splits = splits,
        pval_mode_per_group = pval_mode_per_group,
        empirical_strategy_per_group = empirical_strategy_per_group,
        recheck = recheck
    )
}

.buildConnectivityArray <- function(
    beta_handler,
    beta_locs = NULL,
    pheno,
    group_inds,
    pval_mode_per_group,
    empirical_strategy_per_group,
    col_names = NULL,
    max_pval = 0.05,
    min_delta_beta = 0,
    covariates = NULL,
    max_lookup_dist = 1000,
    chunk_size = getOption("DMRsegal.chunk_size", 1000),
    entanglement = "strong",
    aggfun = median,
    ntries = 500,
    mid_p = TRUE,
    njobs = 1,
    expansion_windows = NULL,
    max_bridge_gaps = 0
) {
    connectivity_array <- NULL
    splits <- NULL
    for (gap in seq(0L, max_bridge_gaps)) {
        if (gap == 0L) {
            .log_info("Building initial connectivity array with no gap bridging.", level = 2)
        } else {
            .log_info("Building bridged connectivity array allowing up to ", gap, " gap(s) between connected seeds.", level = 2)
        }
        build_args <- list(
            beta_handler = beta_handler,
            beta_locs = beta_locs,
            group_inds = group_inds,
            col_names = col_names,
            pheno = pheno,
            covariates = covariates,
            max_lookup_dist = max_lookup_dist,
            max_pval = max_pval,
            min_delta_beta = min_delta_beta,
            entanglement = entanglement,
            aggfun = aggfun,
            pval_mode_per_group = pval_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            chunk_size = chunk_size,
            ntries = ntries,
            mid_p = mid_p,
            njobs = njobs,
            expansion_windows = expansion_windows,
            connectivity_array = connectivity_array,
            splits = splits
        )
        .buildConnectivityArraySinglePassWithGaps <- function(build_args, gap) {
            build_args$ugap <- gap
            build_args$dgap <- 0
            build_ret <- do.call(.buildConnectivityArraySinglePass, build_args)
            build_args$connectivity_array <- build_ret$connectivity_array
            urecheck <- build_ret$recheck
            build_args$ugap <- 0
            build_args$dgap <- gap
            build_ret <- do.call(.buildConnectivityArraySinglePass, build_args)
            drecheck <- build_ret$recheck
            list(
                recheck = urecheck | drecheck,
                connectivity_array = build_ret$connectivity_array,
                splits = build_ret$splits,
                pval_mode_per_group = build_ret$pval_mode_per_group,
                empirical_strategy_per_group = build_ret$empirical_strategy_per_group
            )
        }
        .buildConnectivityArraySinglePassWithGapsRecursive <- function(build_args, gap) {
            build_ret <- .buildConnectivityArraySinglePassWithGaps(build_args, gap)
            if (any(build_ret$recheck)) {
                build_args$recheck <- build_ret$recheck
                build_args$connectivity_array <- build_ret$connectivity_array
                if (gap > 1) {
                    for (g in 1:gap) {
                        g_recheck <- rep(FALSE, nrow(build_args$connectivity_array))
                        while (TRUE) {
                            build_ret <- .buildConnectivityArraySinglePassWithGapsRecursive(build_args, g)
                            if (!any(build_ret$recheck)) {
                                break
                            }
                            build_args$recheck <- build_ret$recheck
                            g_recheck <- g_recheck | build_ret$recheck
                            build_args$connectivity_array <- build_ret$connectivity_array
                        }
                        build_args$recheck <- g_recheck
                    }
                } else {
                    build_ret <- .buildConnectivityArraySinglePassWithGapsRecursive(build_args, gap)
                }
            }
            build_ret
        }
        if (gap > 0L) {
            build_ret <- .buildConnectivityArraySinglePassWithGapsRecursive(build_args, gap)
            connectivity_array <- build_ret$connectivity_array
        } else {
            build_ret <- do.call(.buildConnectivityArraySinglePass, build_args)
            connectivity_array <- build_ret$connectivity_array
        }

        splits <- build_ret$splits
        pval_mode_per_group <- build_ret$pval_mode_per_group
        empirical_strategy_per_group <- build_ret$empirical_strategy_per_group
        .log_info("Connectivity array built with gap allowance of ", gap, " (", sum(connectivity_array$connected), " connected CpGs).", level = 2)
    }
    list(connectivity_array = connectivity_array, splits = splits, pval_mode_per_group = pval_mode_per_group, empirical_strategy_per_group = empirical_strategy_per_group)
}


#' @keywords internal
#' @noRd
.expandDMR <- function(dmr,
                       chr_array,
                       chr_locs,
                       min_cpg_delta_beta = 0,
                       min_cpgs = 3,
                       expansion_step = 500,
                       chr_start_base = 0,
                       subset_to_full_idx = NULL,
                       full_locs = NULL,
                       chr_locs_idx_map = NULL) {
    .log_step("Expanding DMR..", level = 4)
    if (is.data.frame(dmr)) {
        if (nrow(dmr) != 1L) {
            stop("dmr must contain exactly one row.")
        }
        dmr <- as.matrix(dmr)[1, ]
    }
    dmr_start <- dmr["start_seed"]
    dmr_end <- dmr["end_seed"]

    # Use pre-built index map for O(1) lookup instead of O(n) match
    if (!is.null(chr_locs_idx_map)) {
        dmr_start_ind <- chr_locs_idx_map[[dmr_start]]
        dmr_end_ind <- chr_locs_idx_map[[dmr_end]]
        if (is.null(dmr_start_ind)) {
            stop("Could not find the start CpG ", dmr_start, " in the beta file row names.")
        }
        if (is.null(dmr_end_ind)) {
            stop("Could not find the end CpG ", dmr_end, " in the beta file row names.")
        }
        chr_locs_rownames <- names(chr_locs_idx_map)
    } else {
        # Fallback to match() if no index map provided
        chr_locs_rownames <- rownames(chr_locs)
        dmr_start_ind <- match(dmr_start, chr_locs_rownames)
        dmr_end_ind <- match(dmr_end, chr_locs_rownames)
        if (is.na(dmr_start_ind)) {
            stop("Could not find the start CpG ", dmr_start, " in the beta file row names.")
        }
        if (is.na(dmr_end_ind)) {
            stop("Could not find the end CpG ", dmr_end, " in the beta file row names.")
        }
    }
    projected_cpg_ids <- chr_locs_rownames
    projected_positions <- as.integer(chr_locs[, "start"])
    if (!is.null(subset_to_full_idx)) {
        if (length(subset_to_full_idx) != nrow(chr_locs)) {
            stop("Length of subset_to_full_idx must match nrow(chr_locs).")
        }
        if (is.null(full_locs)) {
            stop("full_locs must be provided when subset_to_full_idx is used.")
        }
        if (anyNA(subset_to_full_idx)) {
            stop("subset_to_full_idx contains NA values.")
        }
        projected_cpg_ids <- rownames(full_locs)[subset_to_full_idx]
        projected_positions <- as.integer(full_locs[subset_to_full_idx, "start"])
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
    dmr["start_cpg"] <- projected_cpg_ids[ustream_exp]
    dmr["end_cpg"] <- projected_cpg_ids[dstream_exp]
    dmr["start"] <- projected_positions[ustream_exp]
    dmr["end"] <- projected_positions[dstream_exp]

    to_cpg_ids <- function(local_inds) {
        if (length(local_inds) == 0) {
            return(character(0))
        }
        projected_cpg_ids[local_inds]
    }

    dmr["upstream_cpg_expansion_stop_reason"] <- ustream_stop_reason
    upstream_candidate <- if (ustream_exp <= (dmr_start_ind - 1L)) {
        seq.int(ustream_exp, dmr_start_ind - 1L)
    } else {
        integer(0)
    }
    if (length(upstream_candidate) > 0) {
        bridged_upstream_m <- chr_array[upstream_candidate, "reason"] == "bridged"
        upstream_kept <- upstream_candidate[!bridged_upstream_m]
    } else {
        upstream_kept <- integer(0)
    }
    dmr["upstream_cpgs"] <- paste(to_cpg_ids(upstream_kept), collapse = ",")
    dmr["upstream_cpg_expansion"] <- length(upstream_kept)

    dmr["downstream_cpg_expansion_stop_reason"] <- dstream_stop_reason
    downstream_candidate <- if ((dmr_end_ind + 1L) <= dstream_exp) {
        seq.int(dmr_end_ind + 1L, dstream_exp)
    } else {
        integer(0)
    }
    if (length(downstream_candidate) > 0) {
        bridged_downstream_m <- chr_array[downstream_candidate, "reason"] == "bridged"
        downstream_kept <- downstream_candidate[!bridged_downstream_m]
    } else {
        downstream_kept <- integer(0)
    }
    dmr["downstream_cpgs"] <- paste(to_cpg_ids(downstream_kept), collapse = ",")
    dmr["downstream_cpg_expansion"] <- length(downstream_kept)

    .log_success("Expanded DMR finalized: (start_cpg: ", dmr["start_cpg"], ", end_cpg: ", dmr["end_cpg"], ").", level = 4)
    dmr
}


#' @keywords internal
#' @noRd
.expandDMRChunk <- function(dmr_inds,
                            chr_dmrs,
                            chr_array,
                            chr_locs,
                            min_cpg_delta_beta = 0,
                            min_cpgs = 3,
                            expansion_step = 500,
                            chr_start_base = 0,
                            subset_to_full_idx = NULL,
                            full_locs = NULL,
                            chr_locs_idx_map = NULL) {
    if (length(dmr_inds) == 0L) {
        return(list())
    }
    old_warn <- getOption("warn")
    on.exit(options(warn = old_warn), add = TRUE)
    options(warn = 2)

    ret <- vector("list", length(dmr_inds))
    for (i in seq_along(dmr_inds)) {
        ret[[i]] <- .expandDMR(
            dmr = chr_dmrs[dmr_inds[[i]], , drop = FALSE],
            chr_array = chr_array,
            expansion_step = expansion_step,
            min_cpgs = min_cpgs,
            min_cpg_delta_beta = min_cpg_delta_beta,
            chr_locs = chr_locs,
            chr_start_base = chr_start_base,
            subset_to_full_idx = subset_to_full_idx,
            full_locs = full_locs,
            chr_locs_idx_map = chr_locs_idx_map
        )
    }
    ret
}


#' Vectorized connectivity testing for consecutive site pairs, given their beta values
#'
#' @return Data frame with columns: connected, pval, delta_beta, reason, first_failing_group, stop_reason
#' @keywords internal
#' @noRd
.testConnectivityBatch <- function(sites_beta, group_inds, pheno,
                                   pval_mode_per_group,
                                   empirical_strategy_per_group,
                                   max_pval, covariates = NULL,
                                   min_delta_beta = 0,
                                   max_lookup_dist = NULL, sites_locs = NULL,
                                   entanglement = "strong",
                                   aggfun = mean,
                                   ntries = 0, mid_p = FALSE,
                                   check_non_overlapping = FALSE) {
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
    if (check_non_overlapping) {
        n_pairs <- n_sites / 2
        start_pair_inds <- seq(1, n_sites - 1, by = 2)
        end_pair_inds <- seq(2, n_sites, by = 2)
    } else {
        n_pairs <- n_sites - 1
        start_pair_inds <- seq_len(n_sites - 1)
        end_pair_inds <- seq_len(n_sites - 1) + 1L
    }
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
        dists <- as.numeric(sites_locs[end_pair_inds, "start"]) - as.numeric(sites_locs[start_pair_inds, "start"])
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
        group_m <- .transformBeta(group_beta, pheno = pheno[idx, ], covariates = covariates)

        # Extract pairs matrices
        x_mat_full <- group_m[start_pair_inds, , drop = FALSE] # Sites i
        y_mat_full <- group_m[end_pair_inds, , drop = FALSE] # Sites i+1

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
        effective_pval_mode <- pval_mode_per_group[g]

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
                do_permutations <- empirical_strategy_per_group[g] == "permutations"
                skip_empirical <- FALSE
                if (do_permutations) {
                    min_possible_pval <- 1 / (1 + factorial(m))
                    if (min_possible_pval > max_pval_corrected) {
                        .log_warn(
                            "Cannot compute sufficient small empirical p-values for group '", g,
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
        site2_beta_mat <- sites_beta[end_pair_inds, , drop = FALSE]

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

.checkResult <- function(dmrs, stage, start_col = "start", end_col = "end") {
    end_less_than_start <- dmrs[[end_col]] - dmrs[[start_col]] < 0

    if (any(end_less_than_start)) {
        .log_error(
            paste0( "Error in stage ", stage, ": ",
                sum(end_less_than_start),
                " DMRs have been assigned an end larger than start ! (CODE BUG TO BE REPORTED)",
                " Those are: \n\t",
                paste0(capture.output(print(dmrs[end_less_than_start, ])), collapse = "\n\t")
            )
        )
    }
}


#' Find Differentially Methylated Regions (DMRs) from Differentially Methylated Positions (seeds)
#'
#' This function identifies DMRs from a given set of seeds and a beta value file.
#'
#' @param beta Character. Path to the beta value file, or a tabix file, or a beta matrix, or a BetaHandler object, or a bed file. If a bed file is provided, it must at least contain bed_chrom_col and bed_chrom_start, followed by samples names in the given pheno, with corresponging beta values, and it will be converted to a tabix-indexed beta file internall, with the locations separately saved and queried as a DelayedDataFrame. object.
#' @param seeds Character. Path to the seeds (seeds, etc.) TSV file or the seeds dataframe, in a format like the one produced by dmpFinder.
#' @param pheno Data frame. Phenotype data.
#' @param seeds_id_col Character. Column name or index for Seed identifiers in the seeds TSV file. Default is NULL, which corresponds to the rows names if existing, or the first column if not.
#' @param sample_group_col Character. Column name for sample group information in the phenotype data. Default is NULL.
#' @param casecontrol_col Boolean Column in pheno for case (TRUE/1) / control (FALSE/0) status . If NULL, controls will be assumed to be the first level of sample_group_col. Default is NULL.
#' @param covariates Character vector of column names in pheno to adjust for (e.g. "age", "sex"). When provided, correlations are computed on residuals after regressing M-values on these covariates within each group
#' @param min_cpg_delta_beta Numeric. Minimum delta beta value for CpGs. Default is 0.1.
#' @param expansion_step Numeric. Index-specific step size for expanding DMRs. Increasing it means higher memory usage and faster computation. Default is 500.
#' @param array Character. Type of array used (e.g., "450K", "EPIC", "EPICv2", "27K"). Ignored if using a mouse genome. Also ignored if the beta file is provided as a beta values BED file. Default is "450K".
#' @param genome Character. Genome version. Default is NULL and inferred as "hg19" for 450K, 27K, and EPIC arrays, otherwise "hg38".
#' @param max_pval Numeric. Maximum p-value to assume seeds correlation is significant. Default is 0.05.
#' @param entanglement Character. "strong" (default) requires all groups to show significant correlation for connectivity; "weak" requires at least one group to show significant correlation.
#' @param pval_mode Character. "auto" (default) selects between t-based correlation p-values and empirical p-values per sample group using data diagnostics. You can also force "parametric" for t-based correlation p-values or "empirical" for permutation-based p-values.
#' @param empirical_strategy Character. When pval_mode = "empirical": "auto" (default) uses Monte Carlo for groups with <6 samples and permutations for groups with >=6 samples; "montecarlo" always uses Monte Carlo; "permutations" always uses permutations.
#' @param ntries Integer. Number of permutations when pval_mode = "empirical". Default is 0 (disabled).
#' @param max_lookup_dist Numeric. Maximum distance to look up for adjacent seeds belonging to the same DMR during Stage 1. Default is 10000 (10 kb).
#' @param expansion_window Numeric. Stage 2 connectivity is computed only in windows centered on seed-derived Stage 1 DMR neighborhoods, with this total window width in bp. This value sets a maximum effective size of a DMR after stage 2. Set <=0 for genome-wide connectivity. Default is -1 for microarrays and 10000 (10 kb) for NGS datasets.
#' @param max_bridge_seeds_gaps Integer. Maximum number of consecutive failed seed-to-seed edges to bridge during Stage 1 when both flanking edges are connected and failures are p-value driven. Set to 0 to disable. Default is 1.
#' @param max_bridge_extension_gaps Integer. Maximum gap size to consider during Stage 2 extension. Default is 1 (i.e., at most 1 consecutive failing CpG to bridge).
#' @param min_seeds Numeric. Minimum number of connected seeds in a DMR. Minimum is 2. Default is 2.
#' @param min_adj_seeds Numeric. Minimum number of seeds, adjusted by array CpG density, in a DMR after extension. Minimum is 2. Default is 2.
#' @param min_cpgs Numeric. Minimum number of CpGs in a DMR after extension, including the seeds. Minimum is 2. Default is 3.
#' @param aggfun Function or character. Aggregation function to use when calculating delta beta values and p-values of DMRs. Can be "median", "mean", or a function (e.g., median, mean). Default is "median".
#' @param ignored_sample_groups Character vector. Sample groups to ignore during connection and expansion, separated by commas. Can also be "case" or "control". Default is NULL.
#' @param output_prefix Character. Identifier for the output files. If not provided, no output will be saved. Default is NULL.
#' @param njobs Numeric. Number of parallel jobs to use. Default is the number of available cores.
#' @param beta_row_names_file Character. Path to a file containing row names for the beta values. If not provided, row names will be read from the beta file. Default is NULL.
#' @param annotate_with_genes Logical. Whether to annotate DMRs with overlapping genes. Default is TRUE.
#' @param .score_dmrs Logical. Whether to score DMRs based on significance and effect size. Default is TRUE.
#' @param extract_motifs Logical. Whether to compute DMRs seeds motifs. Default is TRUE.
#' @param bed_provided Logical. Whether the beta file is provided as a BED file. Default is FALSE. In case the input has a .bed extension, this will be set to TRUE automatically.
#' @param bed_chrom_col Character. Column name for chromosome in the BED file. Default is "chrom".
#' @param bed_start_col Character. Column name for start position in the BED file. Default is "start".
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
    expansion_step = 500,
    array = c("450K", "27K", "EPIC", "EPICv2", "NULL"),
    genome = NULL,
    max_pval = 0.05,
    entanglement = c("strong", "weak"),
    pval_mode = c("auto", "parametric", "empirical"),
    empirical_strategy = c("auto", "montecarlo", "permutations"),
    ntries = 200L,
    mid_p = FALSE,
    max_lookup_dist = 10000,
    expansion_window = "auto",
    max_bridge_seeds_gaps = 1L,
    max_bridge_extension_gaps = 1L,
    min_seeds = 2,
    min_adj_seeds = 2,
    min_cpgs = 3,
    aggfun = c("median", "mean"),
    ignored_sample_groups = NULL,
    output_prefix = NULL,
    njobs = getOption("DMRsegal.njobs", min(8, future::availableCores() - 1)),
    chunk_size = getOption("DMRsegal.chunk_size", 10000),
    beta_row_names_file = NULL,
    annotate_with_genes = TRUE,
    .score_dmrs = TRUE,
    extract_motifs = TRUE,
    bed_provided = FALSE,
    bed_chrom_col = "chrom",
    bed_start_col = "start",
    verbose = getOption("DMRsegal.verbose", 1),
    .load_debug = FALSE
) {
    pval_mode <- strex::match_arg(pval_mode, ignore_case = TRUE)
    empirical_strategy <- strex::match_arg(empirical_strategy, ignore_case = TRUE)

    pval_mode_per_group <- rep(pval_mode, length.out = length(unique(pheno[[sample_group_col]])))
    names(pval_mode_per_group) <- unique(pheno[[sample_group_col]])
    empirical_strategy_per_group <- rep(empirical_strategy, length.out = length(unique(pheno[[sample_group_col]])))
    names(empirical_strategy_per_group) <- unique(pheno[[sample_group_col]])

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
    .log_info("Resetting parallel state from previous runs...", level = 2)
    .cleanupParallelState()
    withr::defer(.cleanupParallelState(), envir = environment())
    .log_info(
        "Parallel config: requested njobs=", njobs,
        ", available cores=", future::availableCores(),
        ", multicore-capable cores=", future::availableCores("multicore"),
        level = 2
    )


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
    array <- .normalizeFindDMRsArray(array)
    requested_genome <- .normalizeFindDMRsGenome(genome)
    genome <- .resolveFindDMRsGenome(
        beta = beta,
        array = array,
        genome = requested_genome,
        bed_provided = bed_provided
    )
    if (is.null(requested_genome)) {
        .log_info("No genome provided. Using inferred genome: ", genome, ".", level = 2)
    }
    all_cpgs <- NULL
    beta_locs_rownames <- NULL
    if (inherits(beta, "BetaHandler")) {
        beta_handler <- beta
        all_cpgs <- rownames(beta_handler$getGenomicLocs())
        beta_locs_rownames <- beta_handler$getBetaRowNames()
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
                    start_col = bed_start_col
                )
                beta <- ret$tabix_file
                beta_locs <- ret$locations
                beta_locs_rownames <- rownames(beta_locs)
                all_cpgs <- beta_locs_rownames
                .log_success("Conversion to tabix-indexed beta file completed.")
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

    array_based <- beta_handler$isArrayBased()

    if (is.null(all_cpgs)) {
        all_cpgs <- rownames(beta_handler$getGenomicLocs())
        beta_locs_rownames <- beta_handler$getBetaRowNames()
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
    if (expansion_window == "auto") {
        expansion_window <- if (array_based) -1 else 10000
        if (array_based) {
            .log_info("Setting expansion_window to -1 for array-based dataset, meaning genome-wide connectivity will be computed during expansion.", level = 2)
        } else {
            .log_info("Setting expansion_window to 10000 for sequencing-based dataset, meaning connectivity during expansion will only be computed within 10 kb windows around seed-derived DMRs.", level = 2)
        }
    }
    if (!is.numeric(expansion_window) || length(expansion_window) != 1 || is.na(expansion_window)) {
        stop("expansion_window must be a numeric scalar.")
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
        pheno[, "__casecontrol__"] <- as.numeric(pheno[, casecontrol_col])
    }
    if (is.null(ignored_sample_groups)) {
        ignored_sample_groups <- c()
    } else {
        ignored_sample_groups <- trimws(unlist(base::strsplit(ignored_sample_groups, ",")))
        ignored_sample_groups <- ignored_sample_groups[nzchar(ignored_sample_groups)]
    }
    if (!is.null(output_prefix)) {
        output_prefix_base <- output_prefix
        output_dir <- dirname(output_prefix)
        dir.create(output_dir, showWarnings = FALSE)
        output_prefix <- paste0(output_prefix, ".")
    } else {
        output_prefix_base <- NULL
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
    beta_chr <- as.character(beta_locs[, "chr"])
    beta_start <- suppressWarnings(as.numeric(beta_locs[, "start"]))
    if (anyNA(beta_chr) || any(!nzchar(beta_chr))) {
        stop(
            "Beta locations contain missing chromosome labels. Ensure the beta input includes valid chromosome values.",
            call. = FALSE
        )
    }
    if (anyNA(beta_start)) {
        stop(
            "Beta locations contain missing or non-numeric start positions. Ensure the beta input includes numeric genomic start coordinates.",
            call. = FALSE
        )
    }
    if (length(beta_chr) > 1L) {
        chr_runs <- rle(beta_chr)
        if (anyDuplicated(chr_runs$values)) {
            dup_chr <- chr_runs$values[duplicated(chr_runs$values)][1]
            stop(
                "Beta locations are not grouped by chromosome: ", dup_chr,
                " appears in multiple blocks. Ensure the beta input is ordered by chromosome and genomic start position.",
                call. = FALSE
            )
        }

        same_chr_adj <- beta_chr[-1L] == beta_chr[-length(beta_chr)]
        unsorted_adj <- same_chr_adj & (beta_start[-1L] < beta_start[-length(beta_start)])
        if (any(unsorted_adj, na.rm = TRUE)) {
            bad_idx <- which(unsorted_adj)[1L] + 1L
            stop(
                "Beta locations are not sorted within chromosome ", beta_chr[bad_idx],
                " around positions ", beta_start[bad_idx - 1L], " -> ", beta_start[bad_idx],
                ". Ensure the beta input is ordered by genomic start position.",
                call. = FALSE
            )
        }
    }

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
    viewer_beta_output_file <- if (!is.null(output_prefix)) paste0(output_prefix, "seeds_beta.tsv.gz") else NULL
    saveViewerBetaFile <- function(beta_values = NULL) {
        if (is.null(viewer_beta_output_file)) {
            return(invisible(NULL))
        }
        if (is.null(beta_values)) {
            beta_values <- matrix(
                numeric(0),
                nrow = 0,
                ncol = length(beta_col_names_detection),
                dimnames = list(character(0), beta_col_names_detection)
            )
        }
        .log_info("Saving viewer beta values to file: ", viewer_beta_output_file, " ...", level = 2)
        gz <- gzfile(viewer_beta_output_file, "w")
        on.exit(close(gz), add = TRUE)
        write.table(beta_values, gz, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
        invisible(beta_values)
    }
    if (!is.null(output_prefix_base)) {
        viewer_meta_file <- paste0(output_prefix_base, ".meta.rds")
        viewer_meta <- list(
            pheno = pheno_detection,
            genome = genome,
            array = array,
            sample_group_col = sample_group_col
        )
        saveRDS(viewer_meta, file = viewer_meta_file)
        .log_info("Saved viewer metadata to ", viewer_meta_file, ".", level = 2)
    }
    .log_info("Samples to process during DMR detection: ", length(beta_col_names_detection), level = 1)

    sample_groups <- factor(pheno_detection[, sample_group_col])
    group_inds <- split(seq_along(sample_groups), sample_groups)

    .log_step("Reordering seeds by genomic location...", level = 3)


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
    .log_success("Seeds reordered by genomic location.", level = 3)

    # Filter seeds not present in array annotation first (prevents NA logical indices later)

    seeds <- unique(seeds_tsv[, seeds_id_col])
    missing_in_annotation <- setdiff(seeds, all_cpgs)
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
        stop("No seeds remain after filtering against array annotation and beta file.")
    }
    seeds <- seeds[orderByLoc(seeds, genome = genome, genomic_locs = beta_locs)]

    .log_step("Subsetting beta matrix for seeds...", level = 3)
    seeds_locs <- as.data.frame(beta_locs[seeds, , drop = FALSE])
    seeds_beta <- beta_handler$getBeta(row_names = seeds, col_names = beta_col_names_detection)
    rownames(seeds_beta) <- seeds

    if (nrow(seeds_locs) != nrow(seeds_beta)) {
        stop(
            "Number of rows in the queried seeds beta file does not match the number of seeds. Number of rows in beta file: ",
            nrow(seeds_beta),
            " Number of rows in seeds: ",
            nrow(seeds_locs)
        )
    }
    .log_info("Checking for seeds with all NA beta values...", level = 3)
    all.na.rows <- matrixStats::rowAlls(is.na(as.matrix(seeds_beta)))

    if (any(all.na.rows)) {
        stop(
            "Beta extraction failure: the following Seed rows have all NA beta values: ",
            paste(rownames(seeds_beta)[all.na.rows], collapse = ","),
            ". This indicates a mismatch between requested CpG IDs and beta file columns or a parsing issue."
        )
    }

    .log_success("Subset size: ", paste(dim(seeds_beta), collapse = ","), level = 3)
    .log_info("Number of provided seeds: ", length(seeds), level = 2)
    resolved_min_cpg_delta_beta <- as.numeric(min_cpg_delta_beta)

    .log_info(
        "Using min_cpg_delta_beta threshold: ", signif(resolved_min_cpg_delta_beta, 4),
        level = 2
    )

    .log_success("Input preparation complete.", level = 1)
    .log_step("Stage 1: Connecting seeds to form initial DMRs..", level = 1)


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
        .log_step("Building seed connectivity array...", level = 2)
        ret <- .buildConnectivityArray(
            beta_handler = beta_handler,
            beta_locs = seeds_locs,
            pheno = pheno_detection,
            group_inds = group_inds,
            pval_mode_per_group = pval_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            col_names = beta_col_names_detection,
            max_pval = max_pval,
            min_delta_beta = 0, # delta-beta filtering is applied later during the extension
            covariates = covariates,
            max_lookup_dist = max_lookup_dist,
            chunk_size = chunk_size,
            entanglement = entanglement,
            aggfun = aggfun,
            ntries = ntries,
            mid_p = mid_p,
            njobs = njobs,
            expansion_windows = NULL,
            max_bridge_gaps = max_bridge_seeds_gaps
        )
        seeds_connectivity_array <- ret$connectivity_array
        pval_mode_per_group <- ret$pval_mode_per_group
        empirical_strategy_per_group <- ret$empirical_strategy_per_group
        .log_success("Seed connectivity array built.", level = 2)
        # connected_seeds[i] encodes edge i -> i+1
        connected_seeds <- seeds_connectivity_array$connected

        # vector already includes chromosome-end sentinels as FALSE:
        breakpoints <- which(!connected_seeds)

        connected_seeds_segments_starts <- c(1L, head(breakpoints, -1L) + 1L)
        connected_seeds_segments_ends <- breakpoints
        connected_seeds_segments_lengths <- connected_seeds_segments_ends - connected_seeds_segments_starts + 1L
        connected_seeds_segments_chrs <- as.character(seeds_locs[connected_seeds_segments_starts, "chr"])
        connected_seeds_segments_starts_locs <- as.integer(seeds_locs[connected_seeds_segments_starts, "start"])
        connected_seeds_segments_ends_locs <- as.integer(seeds_locs[connected_seeds_segments_ends, "start"])
        stop_reasons <- seeds_connectivity_array$reason[connected_seeds_segments_ends]
        mask <- rep(FALSE, length(seeds))
        mask[connected_seeds_segments_starts] <- TRUE
        seeds_connectivity_array$id <- cumsum(mask)
        seeds_connectivity_array$cid <- seeds_connectivity_array$id
        ids <- unique(seeds_connectivity_array$id)
        seeds_connectivity_array$id[connected_seeds_segments_ends] <- NA
        seeds_connectivity_array$seeds <- seeds
        .log_info("Number of segments (potential DMRs before filtering): ", length(connected_seeds_segments_chrs), level = 2)
        valid_id_mask <- !is.na(seeds_connectivity_array$id)
        if (any(valid_id_mask)) {
            agg_data <- seeds_connectivity_array[valid_id_mask, , drop = FALSE]
            connected_seeds_connection_corr_pval <- aggregate(
                pval ~ id,
                data = agg_data, function(x) if (all(is.na(x))) NA else mean(x, na.rm = TRUE)
            )
        } else {
            connected_seeds_connection_corr_pval <- data.frame(
                id = ids,
                pval = rep(NA_real_, length(ids))
            )
        }
        dmrs_seeds <- data.frame(
            id = ids,
            seeds = mapply(function(start_idx, end_idx) {
                paste(seeds[start_idx:end_idx], collapse = ",")
            }, connected_seeds_segments_starts, connected_seeds_segments_ends, USE.NAMES = FALSE),
            stringsAsFactors = FALSE
        )

        dmrs <- data.frame(
            chr = connected_seeds_segments_chrs,
            start_seed = seeds[connected_seeds_segments_starts],
            end_seed = seeds[connected_seeds_segments_ends],
            start_seed_pos = connected_seeds_segments_starts_locs,
            end_seed_pos = connected_seeds_segments_ends_locs,
            seeds_num = connected_seeds_segments_lengths,
            stop_connection_reason = stop_reasons,
            id = ids,
            stringsAsFactors = FALSE
        )
        dmrs <- dmrs[dmrs$seeds_num >= min_seeds, , drop = FALSE]
        dmrs <- merge(dmrs, connected_seeds_connection_corr_pval, by = "id", all.x = TRUE)
        colnames(dmrs)[colnames(dmrs) == "pval"] <- "connection_corr_pval"
        dmrs <- merge(dmrs, dmrs_seeds, by = "id", all.x = TRUE)


        if (nrow(dmrs) == 0) {
            .log_warn("No DMRs remain after filtering based on min_seeds.")
            saveViewerBetaFile()
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
    .checkResult(dmrs, "1", start_col = "start_seed_pos", end_col = "end_seed_pos")

    .log_success("Initial DMRs formed: ", nrow(dmrs), level = 1)
    .log_step("Stage 2: Expanding DMRs on neighborhood CpGs..", level = 1)

    # Set up progress tracking for DMR expansion
    n_dmrs <- nrow(dmrs)
    stage2_beta_handler <- beta_handler
    stage2_beta_locs <- beta_locs
    stage2_subset_to_full_idx <- NULL
    if (verbose > 1 && .load_debug && file.exists(file.path("debug", "connectivity_array.rds"))) {
        .log_info("Loading debug connectivity array from file...", level = 2)
        connectivity_array <- readRDS(file.path("debug", "connectivity_array.rds"))
    } else {
        expansion_windows <- NULL
        if (is.finite(expansion_window) && expansion_window > 0) {
            expansion_windows <- .buildConnectivityWindowsFromDMRs(
                dmrs = dmrs,
                expansion_window = expansion_window
            )
            if (nrow(expansion_windows) > 0L) {
                .log_info(
                    "Stage 2 connectivity restricted to ", nrow(expansion_windows),
                    " seed-derived windows (total span: ",
                    format(sum(expansion_windows$end - expansion_windows$start + 1), big.mark = ","),
                    " bp).",
                    level = 2
                )
            } else {
                .log_info("No connectivity windows were generated; Stage 2 will return disconnected CpGs outside chromosome termini.", level = 2)
            }
        } else {
            .log_info("Stage 2 connectivity computed genome-wide (expansion_window <= 0).", level = 2)
        }
        if (!is.null(expansion_windows) && nrow(expansion_windows) > 0L) {
            subset_ret <- .subsetStage2BetaToWindows(
                beta_handler = beta_handler,
                beta_locs = beta_locs,
                col_names = beta_col_names_detection,
                expansion_windows = expansion_windows,
                njobs = njobs
            )
            if (is.null(subset_ret)) {
                stop("Stage 2 window subsetting produced an empty beta subset. This indicates inconsistent expansion windows.")
            }
            stage2_beta_handler <- subset_ret$beta_handler
            stage2_beta_locs <- subset_ret$beta_locs
            stage2_subset_to_full_idx <- subset_ret$subset_to_full_idx
            expansion_windows <- subset_ret$expansion_windows
            .log_info(
                "Stage 2 beta subset contains ",
                format(nrow(stage2_beta_locs), big.mark = ","),
                " CpGs across ",
                length(unique(stage2_beta_locs$chr)),
                " chromosome(s).",
                level = 2
            )
        }
        .log_step("Building expansion connectivity array..", level = 2)
        ret <- .buildConnectivityArray(
            beta_handler = stage2_beta_handler,
            beta_locs = stage2_beta_locs,
            pheno = pheno_detection,
            group_inds = group_inds,
            pval_mode_per_group = pval_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            col_names = beta_col_names_detection,
            max_pval = max_pval,
            min_delta_beta = resolved_min_cpg_delta_beta,
            covariates = covariates,
            max_lookup_dist = max_lookup_dist,
            chunk_size = chunk_size,
            entanglement = entanglement,
            aggfun = aggfun,
            ntries = ntries,
            mid_p = mid_p,
            njobs = njobs,
            expansion_windows = expansion_windows,
            max_bridge_gaps = max_bridge_extension_gaps
        )
        connectivity_array <- ret$connectivity_array
    }
    .log_success("Connectivity array built.", level = 2)
    .log_info("Number of underlying connected CpGs found: ", sum(connectivity_array$connected), level = 1)
    if (getOption("DMRsegal.make_debug_dir", FALSE)) {
        .log_info("Saving connectivity array to debug/connectivity_array.rds", level = 1)
        dir.create("debug", showWarnings = FALSE)
        saveRDS(connectivity_array, file = file.path("debug", "connectivity_array.rds"))
    }
    if (nrow(connectivity_array) != nrow(stage2_beta_locs)) {
        stop(
            "Stage 2 connectivity_array row count (", nrow(connectivity_array),
            ") does not match Stage 2 beta_locs row count (", nrow(stage2_beta_locs),
            "). If using .load_debug, rebuild debug artifacts with the current code."
        )
    }
    .log_info("Expanding ", n_dmrs, " DMRs using ", njobs, " jobs...", level = 2)
    if (verbose > 0) {
        p_ext <- NULL
        # check if version of progressr is equal or higher than >= 0.17.0-9002, otherwise p_ext will not be used
        if (utils::packageVersion("progressr") >= "0.17.0-9002") {
            p_ext <- progressr::progressor(steps = n_dmrs, message = "Expanding DMRs to proximal CpGs..")
        }
    }
    ret <- list()
    for (chr in unique(dmrs$chr)) {
        .log_info("Processing ", chr, level = 2)
        chr_dmrs <- dmrs[dmrs$chr == chr, ]
        chr_mask <- stage2_beta_locs[, "chr"] == chr
        if (!any(chr_mask)) {
            stop("No Stage 2 beta rows available for chromosome ", chr, ".")
        }
        first_row <- which(chr_mask)[1]
        chr_start_base <- first_row - 1
        chr_locs <- as.data.frame(stage2_beta_locs[chr_mask, , drop = FALSE])
        chr_array <- connectivity_array[chr_mask, , drop = FALSE]
        chr_subset_to_full_idx <- NULL
        chr_full_locs <- NULL
        if (!is.null(stage2_subset_to_full_idx)) {
            chr_subset_to_full_idx <- stage2_subset_to_full_idx[which(chr_mask)]
            chr_full_locs <- as.data.frame(beta_locs[chr_subset_to_full_idx, , drop = FALSE])
            rownames(chr_full_locs) <- rownames(beta_locs)[chr_subset_to_full_idx]
            chr_subset_to_full_idx <- seq_len(nrow(chr_locs))
        }
        # Create index map once per chromosome for O(1) lookups
        chr_locs_rownames <- rownames(chr_locs)
        chr_locs_idx_map <- setNames(seq_along(chr_locs_rownames), chr_locs_rownames)

        chr_dmr_inds <- seq_len(nrow(chr_dmrs))
        default_dmr_chunk_size <- max(1L, ceiling(length(chr_dmr_inds) / max(njobs * 4L, 1L)))
        dmr_chunk_size <- as.integer(getOption("DMRsegal.parallel_dmr_chunk_size", default_dmr_chunk_size))[1]
        if (!is.finite(dmr_chunk_size) || dmr_chunk_size < 1L) {
            dmr_chunk_size <- default_dmr_chunk_size
        }
        dmr_chunk_size <- min(dmr_chunk_size, length(chr_dmr_inds))
        chr_chunks <- split(chr_dmr_inds, ceiling(chr_dmr_inds / dmr_chunk_size))

        if (njobs == 1L || length(chr_chunks) == 1L) {
            chr_ret <- lapply(
                chr_chunks,
                .expandDMRChunk,
                chr_dmrs = chr_dmrs,
                chr_array = chr_array,
                expansion_step = expansion_step,
                min_cpgs = min_cpgs,
                min_cpg_delta_beta = min_cpg_delta_beta,
                chr_locs = chr_locs,
                chr_start_base = chr_start_base,
                subset_to_full_idx = chr_subset_to_full_idx,
                full_locs = chr_full_locs,
                chr_locs_idx_map = chr_locs_idx_map
            )
        } else {
            .setupParallel()
            chr_ret <- future.apply::future_lapply(
                X = chr_chunks,
                FUN = .expandDMRChunk,
                chr_dmrs = chr_dmrs,
                chr_array = chr_array,
                expansion_step = expansion_step,
                min_cpgs = min_cpgs,
                min_cpg_delta_beta = min_cpg_delta_beta,
                chr_locs = chr_locs,
                chr_start_base = chr_start_base,
                subset_to_full_idx = chr_subset_to_full_idx,
                full_locs = chr_full_locs,
                chr_locs_idx_map = chr_locs_idx_map,
                future.seed = TRUE,
                future.stdout = NA,
                future.globals = c(
                    ".expandDMRChunk",
                    ".expandDMR",
                    "chr_dmrs",
                    "chr_array",
                    "expansion_step",
                    "min_cpgs",
                    "min_cpg_delta_beta",
                    "chr_locs",
                    "chr_start_base",
                    "chr_subset_to_full_idx",
                    "chr_full_locs",
                    "chr_locs_idx_map"
                )
            )
            .finalizeParallel()
        }
        chr_ret <- unlist(chr_ret, recursive = FALSE, use.names = FALSE)
        if (verbose > 0 && !is.null(p_ext) && length(chr_ret) > 0L) {
            p_ext(amount = length(chr_ret))
        }
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


    .checkResult(extended_dmrs, "2", start_col = "start", end_col = "end")
    .log_success("Post-processing complete.", level = 2)
    .log_success("DMR expansion complete.", level = 1)
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

    extended_dmrs_ranges <- GenomicRanges::makeGRangesFromDataFrame(
        extended_dmrs,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqnames.field = "chr",
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
        agg_df[i, "start_seed"] <- cols_vals$start_seed[[1]]
        agg_df[i, "end_seed"] <- cols_vals$end_seed[[length(inds)]]
        agg_df[i, "start_seed_pos"] <- cols_vals$start_seed_pos[[1]]
        agg_df[i, "end_seed_pos"] <- cols_vals$end_seed_pos[[length(inds)]]
        agg_seeds <- unique(unlist(lapply(cols_vals$seeds, .splitCsvValues), use.names = FALSE))
        agg_df[i, "seeds"] <- paste(agg_seeds, collapse = ",")
        agg_df[i, "seeds_num"] <- length(agg_seeds)
        agg_df[i, "connection_corr_pval"] <- aggfun(as.double(cols_vals$connection_corr_pval), na.rm = TRUE)
        agg_df[i, "stop_connection_reason"] <- paste(cols_vals$stop_connection_reason, collapse = ",")
        agg_df[i, "start_cpg"] <- cols_vals$start_cpg[[1]]
        agg_df[i, "end_cpg"] <- cols_vals$end_cpg[[length(inds)]]

        agg_upstream_cpgs <- unique(unlist(lapply(cols_vals$upstream_cpgs, .splitCsvValues), use.names = FALSE))
        agg_df[i, "upstream_cpg_expansion"] <- paste(cols_vals$upstream_cpg_expansion, collapse = ",")
        agg_df[i, "upstream_cpgs"] <- paste(agg_upstream_cpgs, collapse = ",")
        agg_df[i, "upstream_cpg_expansion_stop_reason"] <- paste(cols_vals$upstream_cpg_expansion_stop_reason, collapse = ",")
        agg_downsteam_cpgs <- unique(unlist(lapply(cols_vals$downstream_cpgs, .splitCsvValues), use.names = FALSE))
        agg_df[i, "downstream_cpg_expansion"] <- paste(cols_vals$downstream_cpg_expansion, collapse = ",")
        agg_df[i, "downstream_cpgs"] <- paste(agg_downsteam_cpgs, collapse = ",")
        agg_df[i, "downstream_cpg_expansion_stop_reason"] <- paste(cols_vals$downstream_cpg_expansion_stop_reason, collapse = ",")
        agg_df[i, "merged_dmrs_num"] <- length(inds)
    }
    agg_df[, "cpgs"] <- apply(agg_df[, c("upstream_cpgs", "seeds", "downstream_cpgs")], 1, function(x) {
        vals <- unique(unlist(lapply(x, .splitCsvValues), use.names = FALSE))
        paste(vals, collapse = ",")
    })
    agg_df[, "supporting_cpgs_num"] <- vapply(agg_df$cpgs, function(x) {
        length(.splitCsvValues(x))
    }, integer(1))
    agg_df[, "cpgs_num"] <- match(agg_df$end_cpg, all_cpgs) - match(agg_df$start_cpg, all_cpgs) + 1
    agg_df[, "id"] <- paste0(seqnames(merged_dmrs_ranges), ":", agg_df$start_cpg, "-", agg_df$end_cpg)

    GenomicRanges::mcols(merged_dmrs_ranges) <- agg_df
    .log_success("Overlapping extended DMRs merged: ", length(merged_dmrs_ranges), " resulting DMRs.", level = 1)
    merged_dmrs <- convertToDataFrame(merged_dmrs_ranges)

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
    filtered_dmrs <- convertToDataFrame(filtered_dmrs_ranges)

    if (nrow(filtered_dmrs) == 0) {
        .log_warn("No DMRs passed the filtering step.")
        saveViewerBetaFile()
        if (!is.null(output_prefix)) {
            for (f in c(".methylation.tsv.gz", ".tsv.gz")) {
                gzfile <- gzfile(file.path(output_dir, paste0(output_prefix, f)), "w", compression = 2)
                close(gzfile)
            }
        }
        return(NULL)
    }

    array_based <- !is.null(array)

    if (array_based && min_adj_seeds > min_seeds) {
        .log_step("Calculating CpG content and adjusted seeds number..", level = 2)
        cpgs_num_bg <- getCpGBackgroundCounts(filtered_dmrs_ranges, genome)
        cpgs_num_bg[!is.finite(cpgs_num_bg) | is.na(cpgs_num_bg) | cpgs_num_bg <= 0] <- 1L
        filtered_dmrs$cpgs_num_bg <- cpgs_num_bg

        filtered_dmrs$seeds_num_adj <- ceiling(filtered_dmrs$cpgs_num / filtered_dmrs$cpgs_num_bg * filtered_dmrs$seeds_num)

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
        saveViewerBetaFile()
        if (!is.null(output_prefix)) {
            for (f in c(".methylation.tsv.gz", ".tsv.gz")) {
                gzfile <- gzfile(file.path(output_dir, paste0(output_prefix, f)), "w", compression = 2)
                close(gzfile)
            }
        }
        return(NULL)
    }

    annotated_dmrs <- filtered_dmrs

    all_selected_cpgs <- unique(unlist(base::strsplit(annotated_dmrs$cpgs, ","), use.names = FALSE))
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
    dmrs_seeds <- base::strsplit(annotated_dmrs$seeds, ",")
    .log_step("Adding DMR delta-beta information..", level = 2)

    dmr_seeds_list <- lapply(dmrs_seeds, as.character)
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

    dmr_cpgs_list <- base::strsplit(annotated_dmrs$cpgs, ",")
    dmr_cpgs <- unlist(dmr_cpgs_list, use.names = FALSE)
    dmr_cpgs_groups <- rep(seq_along(dmr_cpgs_list), lengths(dmr_cpgs_list))

    beta_stats_cpgs <- beta_stats[dmr_cpgs, , drop = FALSE]
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

    if (annotate_with_genes) {
        .log_step("Annotating DMRs with gene information...", level = 1)
        annotated_dmrs <- annotateDMRsWithGenes(annotated_dmrs, genome = genome, njobs = njobs)
        .log_success("DMR annotation completed.", level = 1)
    }

    if (.score_dmrs) {
        .log_step("Scoring DMRs...", level = 1)
        annotated_dmrs <- scoreDMRs(
            annotated_dmrs,
            beta = beta_handler,
            pheno = pheno_all,
            genome = genome,
            array = array,
            sorted_locs = beta_handler$getGenomicLocs(),
            sample_group_col = sample_group_col,
            covariates = covariates,
            njobs = njobs
        )
        .log_success("DMR scoring completed.", level = 1)
    }


    if (is.data.frame(annotated_dmrs)) {
        annotated_dmrs <- convertToGRanges(annotated_dmrs, genome = genome)
    }

    if (extract_motifs) {
        .log_step("Extracting DMR motifs...", level = 1)
        annotated_dmrs <- extractDMRMotifs(annotated_dmrs, genome = genome, array = array, beta_locs = beta_locs)
        .log_success("DMR motifs computed.", level = 1)
    }

    final_dmrs_granges <- annotated_dmrs

    final_dmrs <- convertToDataFrame(final_dmrs_granges)

    .log_info("Final number of DMRs: ", nrow(final_dmrs), level = 1)
    if (!is.null(output_prefix)) {
        viewer_cpgs <- unique(unlist(lapply(S4Vectors::mcols(final_dmrs_granges)$cpgs, .splitCsvValues), use.names = FALSE))
        viewer_beta <- all_selected_cpgs_beta[viewer_cpgs, beta_col_names_detection, drop = FALSE]
        saveViewerBetaFile(viewer_beta)
        dmrs_file <- paste0(output_prefix, "dmrs.tsv.gz")
        .log_step("Saving DMRs to ", dmrs_file, "..", level = 2)
        encoded_dmrs <- .encodeNonTabularColumns(final_dmrs)
        final_dmrs_export <- encoded_dmrs$data
        if (length(encoded_dmrs$encoded_columns) > 0L) {
            .log_info(
                "Serialized non-tabular columns for TSV output: ",
                paste(encoded_dmrs$encoded_columns, collapse = ","),
                level = 2
            )
        }
        gz <- gzfile(dmrs_file, "w")
        write.table(
            final_dmrs_export,
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
