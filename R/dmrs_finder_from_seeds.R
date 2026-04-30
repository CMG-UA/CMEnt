#' Find Differentially Methylated Regions (DMRs) from Pre-computed seeds
#'
#' @name findDMRsFromSeeds
#' @description This function identifies Differentially Methylated Regions (DMRs) from pre-computed
#' Differentially Methylated Positions (seeds) using a correlation-based approach. It expands
#' seeds into regions, considering both statistical significance and biological
#' relevance of methylation changes.
#'
#' @section Important Note on Input Data:
#' Do not apply heavy filtering to your seeds prior to using this function, particularly based on
#' beta values or effect sizes. The function works by expanding regions around seeds
#' and connecting nearby sites into larger regions. Filtering out seeds with smaller effect sizes
#' may remove important sites that could serve as "bridges" to connect more seeds into
#' larger, biologically meaningful DMRs. For optimal results, include all statistically
#' significant seeds (e.g., adjusted p-value < 0.05) and let the function handle region expansion
#' and letting the function reconnect proximal sites during expansion using the
#' ext_site_delta_beta parameter if needed.
#'
#' @param beta Character. Path to the beta value file, or a tabix file, or a beta matrix, or a BetaHandler object, or a bed file, or a BSseq object. If a bed file is provided, it must at least contain bed_chrom_col and bed_chrom_start, followed by samples names in the given pheno, with corresponging beta values, and it will be converted to a tabix-indexed beta file internall, with the locations separately saved and queried as a DelayedDataFrame. If a BSseq object is provided, genomic locations and methylation values will be extracted using bsseq methods.
#' @param seeds Character. Path to the seeds TSV file or the seeds dataframe, in a format like the one produced by dmpFinder.
#' @param pheno Data frame. Phenotype data.
#' @param seeds_id_col Character. Column name or index for Seed identifiers in the seeds TSV file. Default is 1.
#' @param sample_group_col Character. Column name for sample group information in the phenotype data. Default is NULL.
#' @param casecontrol_col Boolean Column in pheno for case (TRUE/1) / control (FALSE/0) status . If NULL, controls will be assumed to be the first level of sample_group_col. Default is NULL.
#' @param covariates Character vector of column names in pheno to adjust for (e.g. "age", "sex"). When provided, correlations are computed on residuals after regressing M-values on these covariates within each group
#' @param ext_site_delta_beta Numeric. Minimum absolute delta beta value that will
#' force proximal sites to be treated as connected during Stage 2 expansion,
#' regardless of their correlation p-value. Set to `NA`, `NULL`, or `Inf` to
#' disable this shortcut. A value of `0` means any proximal site with a
#' non-missing case-control delta beta can be force-connected. Default is 0.2.
#' @param array Character. Type of array used (e.g., "450K", "EPIC", "EPICv2", "27K"). Ignored if using mm10 genome.
#' @param genome Character. Genome version (e.g., "hg38", "hg19", "hs1", "mm10"). Default is NULL and inferred as "hg19" for 450K, 27K, and EPIC arrays, otherwise "hg38".
#' @param max_pval Numeric. Maximum p-value to assume seeds correlation is significant. Default is 0.05.
#' @param pval_mode Character. "auto" (default) selects between t-based correlation p-values and empirical p-values per sample group using data diagnostics. You can also force "parametric" for t-based correlation p-values or "empirical" for permutation-based p-values.
#' @param empirical_strategy Character. When pval_mode = "empirical": "auto" (default) uses Monte Carlo for groups with <6 samples and permutations for groups with >=6 samples; "montecarlo" always uses Monte Carlo; "permutations" always uses permutations.
#' @param ntries Integer. Number of permutations when pval_mode = "empirical". Default is 0 (disabled).
#' @param max_lookup_dist Numeric. Maximum distance to look up for adjacent seeds belonging to the same DMR during Stage 1. Default is 10000 (10 kb).
#' @param expansion_window Numeric. Stage 2 connectivity is computed only in windows centered on seed-derived Stage 1 DMR neighborhoods, with this total window width in bp. Set <=0 for genome-wide connectivity. Default is -1 for microarrays and 10000 (10 kb) for NGS datasets.
#' @param max_bridge_seeds_gaps Integer. Maximum number of consecutive failed seed-to-seed edges to bridge during Stage 1 when both flanking edges are connected and failures are p-value driven. Set to 0 to disable. Default is 1.
#' @param max_bridge_extension_gaps Integer. Maximum gap size to consider during Stage 2 extension. Default is 1 (i.e., at most 1 consecutive failing site to bridge).
#' @param min_seeds Numeric. Minimum number of connected seeds in a DMR. Minimum is 2. Default is 2.
#' @param min_adj_seeds Numeric. Minimum number of seeds, adjusted by array site density, in a DMR after extension. Minimum is 2. Default is 2.
#' @param min_sites Numeric. Minimum number of sites in a DMR after extension, including the seeds. Minimum is 2. Default is 3.
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
#'   \item sites_num: Number of sites in the region
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
#'   \item start_site: ID of the starting site after expansion
#'   \item end_site: ID of the ending site after expansion
#'   \item upstream_expansion_length: Number of sites expanded upstream
#'   \item upstream_expansion_stop_reason: Reason for stopping upstream expansion
#'   \item downstream_expansion_length: Number of sites expanded downstream
#'   \item downstream_expansion_stop_reason: Reason for stopping downstream expansion
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
#' @importFrom stringr str_count str_order
#' @importFrom readr read_tsv
#' @importFrom data.table fread fwrite
#' @importFrom bedr tabix
#' @importFrom BSgenome getSeq
#' @importFrom rtracklayer import.chain
#' @importFrom GenomeInfoDb Seqinfo seqnames
#' @importFrom utils write.table read.table
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom GenomicFeatures genes promoters
#' @importFrom AnnotationDbi mapIds
#' @importFrom S4Vectors queryHits subjectHits


.CASE_CONTROL_COL <- "__casecontrol__"

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
.normalizeForceConnectDeltaBeta <- function(force_connect_delta_beta, arg_name = "force_connect_delta_beta") {
    if (is.null(force_connect_delta_beta)) {
        return(NA_real_)
    }
    if (!is.numeric(force_connect_delta_beta) || length(force_connect_delta_beta) != 1L) {
        stop("'", arg_name, "' must be NULL, NA, Inf, or a numeric scalar in [0, 1].")
    }
    force_connect_delta_beta <- as.numeric(force_connect_delta_beta)[1]
    if (is.na(force_connect_delta_beta) || is.infinite(force_connect_delta_beta)) {
        return(NA_real_)
    }
    if (force_connect_delta_beta < 0 || force_connect_delta_beta > 1) {
        stop("'", arg_name, "' must be NULL, NA, Inf, or a numeric scalar in [0, 1].")
    }
    force_connect_delta_beta
}


#' @keywords internal
#' @noRd
.forceConnectDeltaBetaEnabled <- function(force_connect_delta_beta) {
    !is.na(force_connect_delta_beta) && is.finite(force_connect_delta_beta)
}


#' @keywords internal
#' @noRd
.matchSequencingIdsToBeta <- function(ids, beta_chr, beta_start) {
    ids <- as.character(ids)
    parsed <- regexec("^([^:]+):([0-9]+)$", ids)
    pieces <- regmatches(ids, parsed)
    ok <- lengths(pieces) == 3L
    idx <- rep(NA_integer_, length(ids))
    if (!any(ok)) {
        return(idx)
    }
    id_chr <- vapply(pieces[ok], `[[`, character(1), 2L)
    id_start <- suppressWarnings(as.integer(vapply(pieces[ok], `[[`, character(1), 3L)))
    ok_pos <- which(ok)
    for (chr in unique(id_chr[!is.na(id_start)])) {
        query_pos <- which(id_chr == chr & !is.na(id_start))
        beta_pos <- which(beta_chr == chr)
        if (length(beta_pos) == 0L && startsWith(chr, "chr")) {
            beta_pos <- which(beta_chr == sub("^chr", "", chr))
        } else if (length(beta_pos) == 0L) {
            beta_pos <- which(beta_chr == paste0("chr", chr))
        }
        if (length(beta_pos) == 0L) {
            next
        }
        hit <- match(id_start[query_pos], beta_start[beta_pos])
        idx[ok_pos[query_pos]] <- beta_pos[hit]
    }
    idx
}


#' @keywords internal
#' @noRd
.explicitRowNames <- function(x) {
    rn_info <- .row_names_info(x, type = 0L)
    if (is.integer(rn_info) && length(rn_info) == 2L && is.na(rn_info[1L]) && rn_info[2L] < 0L) {
        return(NULL)
    }
    rownames(x)
}


#' @keywords internal
#' @noRd
.availableRamBytes <- function(default_gb = 2) {
    if (file.exists("/proc/meminfo")) {
        mem_available <- grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)
        if (length(mem_available) == 1L) {
            kb <- suppressWarnings(as.numeric(strsplit(mem_available, "\\s+")[[1]][2L]))
            if (is.finite(kb) && !is.na(kb) && kb > 0) {
                return(kb * 1024)
            }
        }
    }
    as.numeric(default_gb) * 1024^3
}


#' @keywords internal
#' @noRd
.connectivityChunkSize <- function(n_samples, njobs, n_pairs, available_ram_bytes = .availableRamBytes()) {
    n_pairs <- as.integer(n_pairs)
    if (!is.finite(n_pairs) || is.na(n_pairs) || n_pairs < 1L) {
        return(1L)
    }
    denom <- max(1, as.integer(njobs)) * max(1, as.integer(n_samples)) * 8 * 12
    chunk_size <- floor(0.9 * as.numeric(available_ram_bytes) / denom)
    if (!is.finite(chunk_size) || is.na(chunk_size) || chunk_size < 1) {
        chunk_size <- 1
    }
    as.integer(max(1L, min(n_pairs, chunk_size)))
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
    active_chr <- unique(beta_chr)
    if (length(active_chr) != 1L) {
        stop(".subsetStage2BetaToWindows expects beta_locs from exactly one chromosome.")
    }
    beta_start <- as.integer(beta_locs[, "start"])
    keep <- rep(FALSE, n_sites)
    active_wins <- wins[wins$chr == active_chr, , drop = FALSE]
    if (nrow(active_wins) > 0L) {
        win_ir <- IRanges::IRanges(
            start = as.integer(active_wins$start),
            end = as.integer(active_wins$end)
        )
        site_ir <- IRanges::IRanges(
            start = beta_start,
            width = 1L
        )
        ov <- IRanges::findOverlaps(site_ir, win_ir)
        if (length(ov) > 0L) {
            keep[unique(S4Vectors::queryHits(ov))] <- TRUE
        }
    }
    subset_idx <- which(keep)
    if (length(subset_idx) == 0L) {
        return(NULL)
    }
    subset_locs <- as.data.frame(beta_locs[subset_idx, , drop = FALSE])
    beta_locs_rownames_info <- .row_names_info(beta_locs, type = 0L)
    beta_locs_rownames <- if (is.integer(beta_locs_rownames_info) &&
        length(beta_locs_rownames_info) == 2L &&
        is.na(beta_locs_rownames_info[1L]) &&
        beta_locs_rownames_info[2L] < 0L) {
        NULL
    } else {
        rownames(beta_locs)
    }
    if (is.null(beta_locs_rownames)) {
        subset_ids <- subset_idx
        subset_names <- paste0(as.character(subset_locs[, "chr"]), ":", as.integer(subset_locs[, "start"]))
    } else {
        subset_ids <- beta_locs_rownames[subset_idx]
        subset_names <- subset_ids
    }
    rownames(subset_locs) <- subset_names

    if (is.null(col_names)) {
        col_names <- beta_handler$getBetaColNames()
    }
    subset_beta_handler <- beta_handler$subset(
        row_names = subset_ids,
        col_names = col_names
    )
    list(
        beta_handler = subset_beta_handler,
        beta_locs = subset_locs,
        expansion_windows = wins
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

.permutationIndexMatrix <- local({
    cache <- new.env(parent = emptyenv())
    build <- function(n) {
        if (n <= 1L) {
            return(matrix(1L, nrow = 1L, ncol = 1L))
        }
        prev <- build(n - 1L)
        out <- matrix(NA_integer_, nrow = nrow(prev) * n, ncol = n)
        out_row <- 1L
        for (i in seq_len(nrow(prev))) {
            base_perm <- prev[i, ]
            for (pos in seq_len(n)) {
                before <- if (pos > 1L) base_perm[seq_len(pos - 1L)] else integer(0)
                after <- if (pos <= length(base_perm)) base_perm[pos:length(base_perm)] else integer(0)
                out[out_row, ] <- c(before, n, after)
                out_row <- out_row + 1L
            }
        }
        out
    }
    function(n) {
        n <- as.integer(n)[1]
        if (!is.finite(n) || is.na(n) || n < 1L) {
            stop("n must be a positive integer.")
        }
        key <- as.character(n)
        if (!exists(key, envir = cache, inherits = FALSE)) {
            assign(key, build(n), envir = cache)
        }
        get(key, envir = cache, inherits = FALSE)
    }
})


#' @keywords internal
#' @noRd
.connectivityChunkWorker <- function(
    split,
    beta_handler,
    beta_start_vec,
    group_inds,
    pheno,
    pval_mode_per_group,
    empirical_strategy_per_group,
    col_names = NULL,
    max_pval = 0.05,
    force_connect_delta_beta = NA_real_,
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
    force_connect_delta_beta <- .normalizeForceConnectDeltaBeta(force_connect_delta_beta)
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

    site_starts <- beta_start_vec[inds]
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
        force_connect_delta_beta = force_connect_delta_beta,
        max_lookup_dist = max_lookup_dist,
        site_starts = site_starts,
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
    force_connect_delta_beta = NA_real_,
    covariates = NULL,
    max_lookup_dist = 1000,
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
    splits = NULL,
    verbose = 1
) {
    force_connect_delta_beta <- .normalizeForceConnectDeltaBeta(force_connect_delta_beta)
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
        if (.forceConnectDeltaBetaEnabled(force_connect_delta_beta)) {
            ret$delta_beta <- rep(NA_real_, n_sites)
        }
        return(list(connectivity_array = ret, splits = NULL, pval_mode_per_group = pval_mode_per_group, empirical_strategy_per_group = empirical_strategy_per_group))
    }
    chromosomes <- as.character(beta_locs[, "chr"])
    chromosome_levels <- unique(chromosomes)
    if (length(chromosome_levels) != 1L) {
        stop(".buildConnectivityArraySinglePass expects beta_locs from exactly one chromosome.")
    }
    chr_ends <- as.integer(n_sites)
    window_mode <- !is.null(expansion_windows) && nrow(expansion_windows) > 0L
    default_reason <- if (window_mode) "outside_connectivity_window" else ""

    n_cols_for_chunk <- if (!is.null(col_names)) {
        length(col_names)
    } else {
        length(beta_handler$getBetaColNames())
    }
    chunk_size <- .connectivityChunkSize(
        n_samples = n_cols_for_chunk,
        njobs = njobs,
        n_pairs = max(1L, n_sites - 1L)
    )
    .log_info(
        "Connectivity chunk size: ", chunk_size,
        " pair(s), derived from available RAM, ", n_cols_for_chunk,
        " sample column(s), and ", max(1L, as.integer(njobs)), " job(s).",
        level = 3
    )

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
        if (.forceConnectDeltaBetaEnabled(force_connect_delta_beta)) {
            ret$delta_beta <- rep(NA_real_, nrows)
        }
        ret
    }

    .asRecheckIndices <- function(recheck, nrows) {
        if (is.null(recheck)) {
            return(NULL)
        }
        if (is.logical(recheck)) {
            recheck <- which(recheck)
        } else {
            recheck <- suppressWarnings(as.integer(recheck))
        }
        recheck <- recheck[!is.na(recheck) & recheck >= 1L & recheck <= nrows]
        sort(unique(recheck))
    }

    .runsOverlapIndices <- function(run_starts, run_ends, inds) {
        if (length(run_starts) == 0L || length(inds) == 0L) {
            return(rep(FALSE, length(run_starts)))
        }
        next_ind_pos <- findInterval(run_starts - 1L, inds) + 1L
        next_ind_pos <= length(inds) & inds[next_ind_pos] <= run_ends
    }

    .chunkPairRanges <- function(pair_ranges_df) {
        if (is.null(pair_ranges_df) || nrow(pair_ranges_df) == 0L) {
            return(matrix(numeric(0), ncol = 2))
        }
        chunk_size_eff <- as.integer(chunk_size)
        if (!is.finite(chunk_size_eff) || is.na(chunk_size_eff) || chunk_size_eff < 1L) {
            chunk_size_eff <- 1L
        }
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
        if (chunk_size_eff <= 1L) {
            return(splits_mat)
        }

        pooled <- vector("list", nrow(splits_mat))
        pooled_n <- 0L
        cur_start <- as.integer(splits_mat[1L, 1L])
        cur_end <- as.integer(splits_mat[1L, 2L])
        cur_weight <- as.integer(split_weights[1L])
        for (i in seq.int(2L, nrow(splits_mat))) {
            next_start <- as.integer(splits_mat[i, 1L])
            next_end <- as.integer(splits_mat[i, 2L])
            next_weight <- as.integer(split_weights[i])
            contiguous <- next_start <= cur_end + 1L
            if (!contiguous && next_start == cur_end + 2L) {
                contiguous <- TRUE
            }
            fits_budget <- cur_weight + next_weight <= chunk_size_eff
            if (contiguous && fits_budget) {
                cur_end <- max(cur_end, next_end)
                cur_weight <- cur_weight + next_weight
            } else {
                pooled_n <- pooled_n + 1L
                pooled[[pooled_n]] <- c(cur_start, cur_end)
                cur_start <- next_start
                cur_end <- next_end
                cur_weight <- next_weight
            }
        }
        pooled_n <- pooled_n + 1L
        pooled[[pooled_n]] <- c(cur_start, cur_end)
        pooled <- do.call(rbind, pooled[seq_len(pooled_n)])
        if (is.null(dim(pooled))) {
            pooled <- matrix(pooled, ncol = 2L, byrow = TRUE)
        }
        storage.mode(pooled) <- "integer"
        pooled
    }

    .buildAllPairRanges <- function() {
        if (n_sites < 2L) {
            return(data.frame(start_pair = integer(0), end_pair = integer(0)))
        }
        data.frame(
            start_pair = 1L,
            end_pair = as.integer(n_sites - 1L)
        )
    }

    .build_window_pair_ranges <- function() {
        wins <- .mergeGenomicWindows(expansion_windows)
        if (nrow(wins) == 0L) {
            return(data.frame(start_pair = integer(0), end_pair = integer(0)))
        }
        beta_start <- as.integer(beta_locs[, "start"])
        active_wins <- wins[wins$chr == chromosome_levels[[1L]], , drop = FALSE]
        if (nrow(active_wins) == 0L || n_sites < 2L) {
            return(data.frame(start_pair = integer(0), end_pair = integer(0)))
        }
        win_ir <- IRanges::IRanges(
            start = as.integer(round(as.numeric(active_wins$start))),
            end = as.integer(round(as.numeric(active_wins$end)))
        )
        site_ir <- IRanges::IRanges(
            start = beta_start,
            width = 1L
        )
        ov <- IRanges::findOverlaps(win_ir, site_ir)
        if (length(ov) == 0L) {
            return(data.frame(start_pair = integer(0), end_pair = integer(0)))
        }
        qh <- S4Vectors::queryHits(ov)
        sh <- S4Vectors::subjectHits(ov)
        min_rel <- as.integer(tapply(sh, qh, min))
        max_rel <- as.integer(tapply(sh, qh, max))
        pair_start <- min_rel
        pair_end <- max_rel - 1L
        keep <- !is.na(pair_start) & !is.na(pair_end) & pair_end >= pair_start
        if (!any(keep)) {
            return(data.frame(start_pair = integer(0), end_pair = integer(0)))
        }
        pair_ranges <- data.frame(
            start_pair = pair_start[keep],
            end_pair = pair_end[keep]
        )
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
            recheck_inds <- .asRecheckIndices(recheck, n_sites)
            run_mask <- run_mask & .runsOverlapIndices(run_starts, run_ends, recheck_inds)
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
            return(list(connectivity_array = connectivity_array, splits = splits, pval_mode_per_group = pval_mode_per_group, empirical_strategy_per_group = empirical_strategy_per_group, recheck = integer(0)))
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
    beta_row_ids_full <- rownames(beta_locs)
    # Numeric indices in this function are relative to beta_locs, which may be a subset
    # of the full beta matrix. Prefer stable row IDs whenever they are available so that
    # all backends query the same sites during connectivity estimation.
    handler_row_names_for_numeric <- tryCatch(beta_handler$getBetaRowNames(), error = function(e) NULL)
    numeric_row_index_matches_locs <- is.null(handler_row_names_for_numeric) ||
        (length(handler_row_names_for_numeric) >= n_sites &&
            identical(handler_row_names_for_numeric[seq_len(n_sites)], beta_row_ids_full))
    prefer_numeric_row_index <- numeric_row_index_matches_locs
    use_numeric_row_index <- prefer_numeric_row_index ||
        is.null(beta_row_ids_full) ||
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
        site_starts <- as.numeric(beta_locs[splits[1, 1]:(splits[1, 2] + 1L), "start"])
        s <- length(site_starts)
        if (!is.null(max_lookup_dist) && !is.null(site_starts)) {
            dists <- site_starts[2:s] - site_starts[1:(s - 1)]
            exceeded_dist <- dists > max_lookup_dist
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
        rm(first_chunk, site_starts, exceeded_dist, nexdist_mask, groups_options)
        gc()
    }
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
            beta_start_vec = beta_start_vec,
            group_inds = group_inds,
            pheno = pheno,
            pval_mode_per_group = pval_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            col_names = col_names,
            max_pval = max_pval,
            force_connect_delta_beta = force_connect_delta_beta,
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
    recheck <- integer(0)

    .applyChunkResult <- function(item) {
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
                bridge_idx <- rep(update_idx, each = gap) + rep.int(seq.int(0L, gap - 1L), length(update_idx))
                bridge_idx <- bridge_idx[bridge_idx >= 1L & bridge_idx <= n_sites]
                if (length(bridge_idx) > 0L) {
                    bridge_mask[bridge_idx] <- TRUE
                    recheck <<- c(recheck, bridge_idx)
                }
            }
        }
        invisible(NULL)
    }

    .futureBatchConnectivity <- function(batch_splits, batch_checked_pairs = checked_pairs) {
        future.apply::future_lapply(
            X = batch_splits,
            future.seed = TRUE,
            future.stdout = NA,
            future.globals = c(
                ".connectivityChunkWorker",
                ".testConnectivityBatch",
                ".permutationIndexMatrix"
            ),
            FUN = .connectivityChunkWorker,
            beta_handler = beta_handler,
            beta_start_vec = beta_start_vec,
            group_inds = group_inds,
            pheno = pheno,
            pval_mode_per_group = pval_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            col_names = col_names,
            max_pval = max_pval,
            force_connect_delta_beta = force_connect_delta_beta,
            covariates = covariates,
            max_lookup_dist = max_lookup_dist,
            entanglement = entanglement,
            aggfun = aggfun,
            ntries = ntries,
            mid_p = mid_p,
            checked_pairs = batch_checked_pairs,
            use_numeric_row_index = use_numeric_row_index,
            beta_row_ids = beta_row_ids,
            beta_row_ids_offset = 0L
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

    if (njobs == 1L || nrow(splits) == 1L) {
        for (split_ind in seq_len(nrow(splits))) {
            .applyChunkResult(.runConnectivityChunk(splits[split_ind, ], checked_pairs_local = checked_pairs))
        }
    } else {
        .setupParallel()
        on.exit(.finalizeParallel(), add = TRUE)
        all_split_inds <- seq_len(nrow(splits))
        batch_size <- max(1L, as.integer(njobs))
        for (batch_start in seq.int(1L, length(all_split_inds), by = batch_size)) {
            batch_end <- min(batch_start + batch_size - 1L, length(all_split_inds))
            batch_inds <- all_split_inds[batch_start:batch_end]
            batch_splits <- lapply(batch_inds, function(i) as.integer(splits[i, ]))
            batch_checked_pairs <- .subsetCheckedPairsForBatch(batch_inds)
            ret <- .futureBatchConnectivity(batch_splits, batch_checked_pairs = batch_checked_pairs)
            for (item in ret) {
                .applyChunkResult(item)
            }
            rm(ret, batch_checked_pairs, batch_splits)
            gc(verbose = FALSE)
        }
    }

    connected_vec[bridge_mask] <- TRUE
    reason_vec[bridge_mask] <- "bridged"
    recheck <- sort(unique(recheck))

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
    force_connect_delta_beta = NA_real_,
    covariates = NULL,
    max_lookup_dist = 1000,
    entanglement = "strong",
    aggfun = median,
    ntries = 500,
    mid_p = TRUE,
    njobs = 1,
    expansion_windows = NULL,
    max_bridge_gaps = 0,
    verbose = getOption("CMEnt.verbose", 1L)
) {
    force_connect_delta_beta <- .normalizeForceConnectDeltaBeta(force_connect_delta_beta)
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
            force_connect_delta_beta = force_connect_delta_beta,
            entanglement = entanglement,
            aggfun = aggfun,
            pval_mode_per_group = pval_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            ntries = ntries,
            mid_p = mid_p,
            njobs = njobs,
            expansion_windows = expansion_windows,
            connectivity_array = connectivity_array,
            splits = splits,
            verbose = verbose
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
                recheck = sort(unique(c(urecheck, drecheck))),
                connectivity_array = build_ret$connectivity_array,
                splits = build_ret$splits,
                pval_mode_per_group = build_ret$pval_mode_per_group,
                empirical_strategy_per_group = build_ret$empirical_strategy_per_group
            )
        }
        .buildConnectivityArrayBridgeFixedPoint <- function(build_args, max_gap) {
            active_recheck <- NULL
            base_splits <- build_args$splits
            pass <- 0L
            repeat {
                cycle_recheck <- integer(0)
                for (bridge_gap in seq_len(max_gap)) {
                    pass <- pass + 1L
                    build_args$ugap <- NULL
                    build_args$dgap <- NULL
                    build_args$recheck <- active_recheck
                    build_args$splits <- base_splits
                    .log_info(
                        "Bridge fixed-point pass ", pass,
                        " for gap ", bridge_gap, "/", max_gap,
                        if (is.null(active_recheck)) " over all connected runs." else paste0(" over ", length(active_recheck), " touched edge(s)."),
                        level = 3
                    )
                    build_ret <- .buildConnectivityArraySinglePassWithGaps(build_args, bridge_gap)
                    build_args$connectivity_array <- build_ret$connectivity_array
                    build_args$pval_mode_per_group <- build_ret$pval_mode_per_group
                    build_args$empirical_strategy_per_group <- build_ret$empirical_strategy_per_group
                    cycle_recheck <- sort(unique(c(cycle_recheck, build_ret$recheck)))
                }
                if (length(cycle_recheck) == 0L) {
                    build_ret$splits <- base_splits
                    return(build_ret)
                }
                active_recheck <- cycle_recheck
            }
        }
        if (gap > 0L) {
            build_ret <- .buildConnectivityArrayBridgeFixedPoint(build_args, gap)
            connectivity_array <- build_ret$connectivity_array
        } else {
            build_ret <- do.call(.buildConnectivityArraySinglePass, build_args)
            connectivity_array <- build_ret$connectivity_array
        }

        splits <- build_ret$splits
        pval_mode_per_group <- build_ret$pval_mode_per_group
        empirical_strategy_per_group <- build_ret$empirical_strategy_per_group
        .log_info("Connectivity array built with gap allowance of ", gap, " (", sum(connectivity_array$connected), " connected sites).", level = 2)
    }
    list(connectivity_array = connectivity_array, splits = splits, pval_mode_per_group = pval_mode_per_group, empirical_strategy_per_group = empirical_strategy_per_group)
}


#' @keywords internal
#' @noRd
.buildExpansionBoundaryLookup <- function(connectivity_array) {
    n_edges <- nrow(connectivity_array)
    if (n_edges == 0L) {
        return(list(
            previous_failed_edge = integer(0),
            next_failed_edge = integer(0),
            reason = character(0)
        ))
    }

    failed_edges <- which(!is.na(connectivity_array$connected) & !connectivity_array$connected)

    previous_failed_edge <- integer(n_edges)
    if (length(failed_edges) > 0L) {
        previous_failed_edge[failed_edges] <- failed_edges
        previous_failed_edge <- cummax(previous_failed_edge)
    }

    next_failed_edge <- rep.int(n_edges + 1L, n_edges)
    if (length(failed_edges) > 0L) {
        next_failed_edge[failed_edges] <- failed_edges
        next_failed_edge <- rev(cummin(rev(next_failed_edge)))
    }

    list(
        previous_failed_edge = previous_failed_edge,
        next_failed_edge = next_failed_edge,
        reason = as.character(connectivity_array$reason)
    )
}


#' @keywords internal
#' @noRd
.expandDMR <- function(dmr,
                       connectivity_array,
                       locs,
                       min_sites = 3,
                       locs_idx_map,
                       expansion_boundaries) {
    .log_step("Expanding DMR..", level = 4)
    if (is.data.frame(dmr)) {
        if (nrow(dmr) != 1L) {
            stop("dmr must contain exactly one row.")
        }
        dmr <- dmr[1L, , drop = FALSE]
    } else {
        dmr <- as.data.frame(as.list(dmr), stringsAsFactors = FALSE)
    }
    dmr_start <- as.character(dmr[["start_seed"]][[1]])
    dmr_end <- as.character(dmr[["end_seed"]][[1]])

    dmr_start_ind <- locs_idx_map[[dmr_start]]
    dmr_end_ind <- locs_idx_map[[dmr_end]]
    if (is.null(dmr_start_ind)) {
        stop("Could not find the start site ", dmr_start, " in the beta file row names.")
    }
    if (is.null(dmr_end_ind)) {
        stop("Could not find the end site ", dmr_end, " in the beta file row names.")
    }
    dmr_start_ind <- as.integer(dmr_start_ind)
    dmr_end_ind <- as.integer(dmr_end_ind)
    locs_rownames <- names(locs_idx_map)
    projected_site_ids <- locs_rownames
    projected_positions <- as.integer(locs[, "start"])

    .check_upstream <- function(ustream_exp, exp_step) {
        ustream_end_lookup_site_ind <- as.integer(ustream_exp - 1L)
        if (ustream_end_lookup_site_ind < 1L) {
            return(list(
                ustream_stop_reason = "end-of-input",
                ustream_exp = ustream_exp
            ))
        }
        fail_edge <- expansion_boundaries$previous_failed_edge[ustream_end_lookup_site_ind]
        if (is.na(fail_edge) || fail_edge < 1L) {
            return(list(
                ustream_stop_reason = "end-of-input",
                ustream_exp = 1L
            ))
        }
        list(
            ustream_stop_reason = expansion_boundaries$reason[fail_edge],
            ustream_exp = as.integer(fail_edge + 1L)
        )
    }

    .check_downstream <- function(dstream_exp, exp_step) {
        dstream_start_lookup_site_ind <- as.integer(dstream_exp)
        if (dstream_start_lookup_site_ind >= nrow(locs)) {
            return(list(
                dstream_stop_reason = "end-of-input",
                dstream_exp = dstream_exp
            ))
        }
        fail_edge <- expansion_boundaries$next_failed_edge[dstream_start_lookup_site_ind]
        if (is.na(fail_edge) || fail_edge > nrow(locs)) {
            return(list(
                dstream_stop_reason = "end-of-input",
                dstream_exp = nrow(locs)
            ))
        }
        list(
            dstream_stop_reason = expansion_boundaries$reason[fail_edge],
            dstream_exp = as.integer(fail_edge)
        )
    }

    ustream_exp <- dmr_start_ind
    ustream_stop_reason <- NULL
    dstream_exp <- dmr_end_ind
    dstream_stop_reason <- NULL

    t <- 0
    while (TRUE) {
        exp_step <- nrow(locs)
        if (t == 0) { # first iteration, use min_sites and remove the DMRs that are not long enough
            csites <- dstream_exp - ustream_exp + 1
            if (csites < (min_sites)) {
                .log_info("DMR  too short (", csites, " sites). Expanding to reach min_sites=", min_sites, ".", level = 4)
                exp_step <- min_sites - csites
            }
            .log_info("Number of sites in DMR: ", csites, level = 5)
        }
        .log_info("Expansion step size: ", exp_step, " bp.", level = 5)
        .log_step("Checking upstream expansion...", level = 5)
        if (is.null(ustream_stop_reason)) {
            res <- .check_upstream(ustream_exp, exp_step)
            ustream_stop_reason <- res$ustream_stop_reason
            if (res$ustream_exp > ustream_exp) {
                .log_info("Upstream expanded by ", ustream_exp - res$ustream_exp, " sites.", level = 5)
            }
            ustream_exp <- res$ustream_exp
        }
        .log_success("Upstream expansion checked.", level = 5)
        .log_step("Checking downstream expansion...", level = 5)
        if (is.null(dstream_stop_reason)) {
            res <- .check_downstream(dstream_exp, exp_step)
            dstream_stop_reason <- res$dstream_stop_reason
            if (res$dstream_exp > dstream_exp) {
                .log_info("Downstream expanded by ", res$dstream_exp - dstream_exp, " sites.", level = 5)
            }
            dstream_exp <- res$dstream_exp
        }
        .log_success("Downstream expansion checked.", level = 5)
        if (t == 0) {
            new_csites <- dstream_exp - ustream_exp + 1
            .log_info("Number of sites in expanded DMR after first iteration: ", new_csites, " from ", csites, level = 4)
            if (new_csites < min_sites) {
                ustream_stop_reason <- "min-sites-not-reached"
                dstream_stop_reason <- "min-sites-not-reached"
                .log_info("DMR could not reach min_sites=", min_sites, " after expansion (", new_csites, "). Stopping expansion.", level = 4)
            }
            t <- 1
        }
        if (!is.null(ustream_stop_reason) && !is.null(dstream_stop_reason)) {
            break
        }
    }
    .log_step("Finalizing expanded DMR.", level = 4)
    dmr[["start_site"]] <- projected_site_ids[ustream_exp]
    dmr[["end_site"]] <- projected_site_ids[dstream_exp]
    dmr[["start"]] <- projected_positions[ustream_exp]
    dmr[["end"]] <- projected_positions[dstream_exp]

    to_site_ids <- function(local_inds) {
        if (length(local_inds) == 0) {
            return(character(0))
        }
        projected_site_ids[local_inds]
    }

    dmr[["upstream_expansion_stop_reason"]] <- ustream_stop_reason
    upstream_candidate <- if (ustream_exp <= (dmr_start_ind - 1L)) {
        seq.int(ustream_exp, dmr_start_ind - 1L)
    } else {
        integer(0)
    }
    if (length(upstream_candidate) > 0) {
        bridged_upstream_m <- connectivity_array[upstream_candidate, "reason"] == "bridged"
        upstream_kept <- upstream_candidate[!bridged_upstream_m]
    } else {
        upstream_kept <- integer(0)
    }
    dmr[["upstream_sites"]] <- paste(to_site_ids(upstream_kept), collapse = ",")
    dmr[["upstream_expansion_length"]] <- length(upstream_kept)

    dmr[["downstream_expansion_stop_reason"]] <- dstream_stop_reason
    downstream_candidate <- if ((dmr_end_ind + 1L) <= dstream_exp) {
        seq.int(dmr_end_ind + 1L, dstream_exp)
    } else {
        integer(0)
    }
    if (length(downstream_candidate) > 0) {
        bridged_downstream_m <- connectivity_array[downstream_candidate, "reason"] == "bridged"
        downstream_kept <- downstream_candidate[!bridged_downstream_m]
    } else {
        downstream_kept <- integer(0)
    }
    dmr[["downstream_sites"]] <- paste(to_site_ids(downstream_kept), collapse = ",")
    dmr[["downstream_expansion_length"]] <- length(downstream_kept)

    .log_success("Expanded DMR finalized: (start_site: ", dmr[["start_site"]], ", end_site: ", dmr[["end_site"]], ").", level = 4)
    dmr
}


#' @keywords internal
#' @noRd
.expandDMRChunk <- function(dmr_inds,
                            dmrs,
                            connectivity_array,
                            locs,
                            min_sites = 3,
                            locs_idx_map,
                            expansion_boundaries) {
    if (length(dmr_inds) == 0L) {
        return(list())
    }
    old_warn <- getOption("warn")
    on.exit(options(warn = old_warn), add = TRUE)
    options(warn = 2)

    ret <- vector("list", length(dmr_inds))
    for (i in seq_along(dmr_inds)) {
        ret[[i]] <- .expandDMR(
            dmr = dmrs[dmr_inds[[i]], , drop = FALSE],
            connectivity_array = connectivity_array,
            min_sites = min_sites,
            locs = locs,
            locs_idx_map = locs_idx_map,
            expansion_boundaries = expansion_boundaries
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
                                   force_connect_delta_beta = NA_real_,
                                   max_lookup_dist = NULL, site_starts = NULL,
                                   entanglement = "strong",
                                   aggfun = mean,
                                   ntries = 0, mid_p = FALSE,
                                   check_non_overlapping = FALSE) {
    force_connect_delta_beta <- .normalizeForceConnectDeltaBeta(force_connect_delta_beta)
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
    if (!is.null(max_lookup_dist) && !is.null(site_starts)) {
        dists <- as.numeric(site_starts[end_pair_inds]) - as.numeric(site_starts[start_pair_inds])
        exceeded_dist <- dists > max_lookup_dist
        connected[exceeded_dist] <- FALSE
        reasons[exceeded_dist] <- "exceeded max distance"
    } else {
        exceeded_dist <- rep(FALSE, n_pairs)
    }
    nexdist_mask <- !exceeded_dist
    .log_info(sum(exceeded_dist), " out of ", n_pairs, " site pairs exceeded the maximum lookup distance and will be marked as not connected.", level = 4)

    high_delta <- rep(FALSE, n_pairs)
    if (.forceConnectDeltaBetaEnabled(force_connect_delta_beta) && length(unique(pheno[, .CASE_CONTROL_COL])) > 1) {
        # Compute this inexpensive case-control effect-size screen up front so
        # proximal high-delta pairs can bypass correlation testing entirely.
        site2_beta_mat <- sites_beta[end_pair_inds, , drop = FALSE]
        case_betas <- apply(site2_beta_mat[, pheno[, .CASE_CONTROL_COL] == 1, drop = FALSE], 1, aggfun, na.rm = TRUE)
        control_betas <- apply(site2_beta_mat[, pheno[, .CASE_CONTROL_COL] == 0, drop = FALSE], 1, aggfun, na.rm = TRUE)
        delta_betas <- case_betas - control_betas
        high_delta <- !is.na(delta_betas) & abs(delta_betas) >= force_connect_delta_beta & nexdist_mask
        connected[high_delta] <- TRUE
        reasons[high_delta] <- "abs(delta_beta)>=force_connect_delta_beta"
    }
    corr_mask <- nexdist_mask & !high_delta

    if (!strict_mode && n_groups > 0) {
        per_group_reasons <- matrix("", nrow = n_groups, ncol = n_pairs)
        per_group_reasons[, exceeded_dist] <- "exceeded max distance"
        per_group_reasons[, high_delta] <- "abs(delta_beta)>=force_connect_delta_beta"
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

        # Extract only correlation-eligible pair rows; this avoids retaining
        # full pair matrices for distance-filtered or delta-beta-rescued chunks.
        x_mat <- group_m[start_pair_inds[corr_mask], , drop = FALSE]
        y_mat <- group_m[end_pair_inds[corr_mask], , drop = FALSE]

        sn_pairs <- nrow(x_mat)
        if (sn_pairs == 0L) {
            next
        }
        g_reasons <- rep("", sn_pairs)
        g_mask <- rep(TRUE, sn_pairs)
        if (strict_mode) {
            g_mask <- connected[corr_mask]
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
        denom <- sqrt(sum_x2 * sum_y2)

        # Compute correlations (fully vectorized)
        cors <- sum_xy / denom

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
                set.seed(getOption("CMEnt.random_seed", 42))
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
                    perm_matrix <- NULL
                    if (do_permutations) {
                        n_unique_perms <- factorial(m)
                        if (is.finite(n_unique_perms) && n_unique_perms <= ntries) {
                            perm_matrix <- .permutationIndexMatrix(m)
                            ntries <- nrow(perm_matrix)
                            .log_info(
                                "Using all ", ntries,
                                " unique sample-label permutations for group '", g,
                                "' instead of redundant random draws.",
                                level = 4
                            )
                        }
                    }
                    abs_cors <- abs(cors)
                    compare_tol <- suppressWarnings(as.numeric(
                        getOption("CMEnt.permutation_compare_tolerance", sqrt(.Machine$double.eps))
                    )[1])
                    if (!is.finite(compare_tol) || is.na(compare_tol) || compare_tol < 0) {
                        compare_tol <- sqrt(.Machine$double.eps)
                    }
                    maxval <- if (do_permutations) NA_real_ else max(y_mat, na.rm = TRUE)
                    minval <- if (do_permutations) NA_real_ else min(y_mat, na.rm = TRUE)
                    for (b in seq_len(ntries)) {
                        # Permute sample labels (columns) only for y; x remains fixed
                        if (do_permutations) {
                            perm <- if (is.null(perm_matrix)) {
                                sample.int(m, size = m, replace = FALSE)
                            } else {
                                perm_matrix[b, ]
                            }
                            yc <- y_centered[, perm, drop = FALSE]
                            sxy <- rowSums(x_centered * yc, na.rm = TRUE)
                            rperm <- sxy / denom
                        } else {
                            yp <- matrix(stats::runif(n = nrow(y_mat) * m, min = minval, max = maxval), nrow = nrow(y_mat), ncol = m)
                            yc <- yp - rowMeans(yp, na.rm = TRUE)
                            sxy <- rowSums(x_centered * yc, na.rm = TRUE)
                            sy2 <- rowSums(yc^2, na.rm = TRUE)
                            rperm <- sxy / sqrt(sum_x2 * sy2)
                        }
                        comp_mask <- is.finite(rperm)
                        if (any(comp_mask)) {
                            ap <- abs(rperm[comp_mask])
                            ao <- abs_cors[comp_mask]
                            counts_ge[comp_mask] <- counts_ge[comp_mask] + (ap > ao + compare_tol)
                            counts_eq[comp_mask] <- counts_eq[comp_mask] + (abs(ap - ao) <= compare_tol)
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
            sconnected <- connected[corr_mask]
            broad_mask <- corr_mask & connected
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
            per_group_p[g_index, corr_mask] <- ps
            per_group_reasons[g_index, corr_mask] <- g_reasons
        }
        .log_success("Finished processing chunk for group '", g, "'.", level = 4)
    }
    if (!strict_mode) {
        not_failed <- per_group_reasons == ""
        connected[] <- FALSE
        connected[high_delta] <- TRUE
        connected[corr_mask] <- colSums(per_group_reasons[, corr_mask, drop = FALSE] == "") > 0
        pvals[corr_mask & connected] <- as.vector(apply(
            per_group_p[, corr_mask & connected, drop = FALSE], 2, function(v) {
                if (all(is.na(v))) {
                    return(NA_real_)
                }
                max(v, na.rm = TRUE)
            }
        ))
        reasons[corr_mask & !connected] <- apply(
            per_group_reasons[, corr_mask & !connected, drop = FALSE], 2, function(v) paste(v, collapse = ";")
        )
        failing_groups[corr_mask & !connected] <- apply(
            !not_failed[, corr_mask & !connected, drop = FALSE], 2, function(v) paste(names(group_inds)[v], collapse = ";")
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
    if (exists("delta_betas", inherits = FALSE)) {
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

#' @keywords internal
#' @noRd
.aggregateMergedDMRRow <- function(hit_idx, qh, sh, orig_mcols, aggfun) {
    inds <- sh[qh == hit_idx]
    cols_vals <- orig_mcols[inds, , drop = FALSE]
    agg_seeds <- unique(unlist(lapply(cols_vals$seeds, .splitCsvValues), use.names = FALSE))
    agg_upstream_sites <- unique(unlist(lapply(cols_vals$upstream_sites, .splitCsvValues), use.names = FALSE))
    agg_downsteam_sites <- unique(unlist(lapply(cols_vals$downstream_sites, .splitCsvValues), use.names = FALSE))
    list(
        idx = as.integer(hit_idx),
        start_seed = as.character(cols_vals$start_seed[[1]]),
        end_seed = as.character(cols_vals$end_seed[[length(inds)]]),
        start_seed_pos = as.numeric(cols_vals$start_seed_pos[[1]]),
        end_seed_pos = as.numeric(cols_vals$end_seed_pos[[length(inds)]]),
        seeds = paste(agg_seeds, collapse = ","),
        seeds_num = as.integer(length(agg_seeds)),
        connection_corr_pval = aggfun(as.double(cols_vals$connection_corr_pval), na.rm = TRUE),
        stop_connection_reason = paste(cols_vals$stop_connection_reason, collapse = ","),
        start_site = as.character(cols_vals$start_site[[1]]),
        end_site = as.character(cols_vals$end_site[[length(inds)]]),
        upstream_expansion_length = paste(cols_vals$upstream_expansion_length, collapse = ","),
        upstream_sites = paste(agg_upstream_sites, collapse = ","),
        upstream_expansion_stop_reason = paste(cols_vals$upstream_expansion_stop_reason, collapse = ","),
        downstream_expansion_length = paste(cols_vals$downstream_expansion_length, collapse = ","),
        downstream_sites = paste(agg_downsteam_sites, collapse = ","),
        downstream_expansion_stop_reason = paste(cols_vals$downstream_expansion_stop_reason, collapse = ","),
        merged_dmrs_num = as.integer(length(inds))
    )
}

#' @keywords internal
#' @noRd
.aggregateMergedDMRChunk <- function(hit_indices, qh, sh, orig_mcols, aggfun) {
    if (length(hit_indices) == 0L) {
        return(list())
    }
    lapply(
        hit_indices,
        .aggregateMergedDMRRow,
        qh = qh,
        sh = sh,
        orig_mcols = orig_mcols,
        aggfun = aggfun
    )
}

#' @keywords internal
#' @noRd
.collapseMergedDMRsitesChunk <- function(components_df) {
    if (is.null(components_df) || nrow(components_df) == 0L) {
        return(character(0))
    }
    vapply(
        seq_len(nrow(components_df)),
        function(i) {
            vals <- unique(unlist(lapply(components_df[i, ], .splitCsvValues), use.names = FALSE))
            paste(vals, collapse = ",")
        },
        character(1)
    )
}

#' @keywords internal
#' @noRd
.aggregateDMRBetaStatsChunk <- function(beta_stats_chunk, aggfun) {
    if (is.null(beta_stats_chunk) || nrow(beta_stats_chunk) == 0L) {
        return(data.frame(
            dmr_id = integer(0),
            cases_beta = numeric(0),
            controls_beta = numeric(0),
            cases_beta_sd = numeric(0),
            controls_beta_sd = numeric(0),
            cases_beta_min = numeric(0),
            cases_beta_max = numeric(0),
            controls_beta_min = numeric(0),
            controls_beta_max = numeric(0)
        ))
    }
    as.data.frame(data.table::as.data.table(beta_stats_chunk)[, .(
        cases_beta = aggfun(abs(cases_beta)) * sign(sum(sign(cases_beta))),
        controls_beta = aggfun(abs(controls_beta)) * sign(sum(sign(controls_beta))),
        cases_beta_sd = aggfun(cases_beta_sd, na.rm = TRUE),
        controls_beta_sd = aggfun(controls_beta_sd, na.rm = TRUE),
        cases_beta_min = min(cases_beta, na.rm = TRUE),
        cases_beta_max = max(cases_beta, na.rm = TRUE),
        controls_beta_min = min(controls_beta, na.rm = TRUE),
        controls_beta_max = max(controls_beta, na.rm = TRUE)
    ), by = dmr_id])
}

#' @keywords internal
#' @importFrom stats aggregate
#' @noRd
.aggregateDMRBetaStats <- function(beta_stats_df,
                                   aggfun,
                                   njobs = 1L,
                                   parallel_enabled = FALSE,
                                   min_groups_for_parallel = 1000L) {
    ret <- .aggregateDMRBetaStatsChunk(data.table::as.data.table(beta_stats_df), aggfun)
    ret <- ret[order(ret$dmr_id), , drop = FALSE]
    rownames(ret) <- NULL
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


#' @keywords internal
#' @noRd
.findDMRsFromSeedsChr <- function(
    beta_handler,
    seed_ids,
    seed_beta_index,
    pheno_detection,
    group_inds,
    pval_mode_per_group,
    empirical_strategy_per_group,
    beta_col_names_detection,
    ext_site_delta_beta,
    array,
    genome,
    max_pval,
    entanglement,
    ntries,
    mid_p,
    max_lookup_dist,
    expansion_window,
    max_bridge_seeds_gaps,
    max_bridge_extension_gaps,
    min_seeds,
    min_adj_seeds,
    min_sites,
    aggfun,
    njobs,
    verbose,
    .load_debug,
    pheno_all,
    beta_col_names,
    sample_group_col,
    covariates,
    annotate_with_genes,
    .score_dmrs,
    extract_motifs
) {
    if (!inherits(beta_handler, "BetaHandler")) {
        stop(".findDMRsFromSeedsChr expects a chromosome-scoped BetaHandler.")
    }
    array_based <- beta_handler$isArrayBased()
    beta_locs <- beta_handler$getBetaLocs()
    all_sites <- .explicitRowNames(beta_locs)
    chromosome_levels <- unique(as.character(beta_locs[, "chr"]))
    chromosome_progress <- NULL
    chromosome_progress_step <- 0L
    chromosome_progress_total <- 6L +
        as.integer(isTRUE(annotate_with_genes)) +
        as.integer(isTRUE(.score_dmrs)) +
        as.integer(isTRUE(extract_motifs))
    if (verbose >= 1L && length(chromosome_levels) == 1L && chromosome_progress_total > 0L) {
        chromosome_progress <- utils::txtProgressBar(
            min = 0,
            max = chromosome_progress_total,
            style = 3,
            file = stderr()
        )
        on.exit(
            {
                if (!is.null(chromosome_progress)) {
                    close(chromosome_progress)
                }
            },
            add = TRUE
        )
    }
    .advanceChromosomeProgress <- function() {
        if (is.null(chromosome_progress)) {
            return(invisible(NULL))
        }
        chromosome_progress_step <<- chromosome_progress_step + 1L
        utils::setTxtProgressBar(
            chromosome_progress,
            min(chromosome_progress_step, chromosome_progress_total)
        )
        invisible(NULL)
    }

    .advanceChromosomeProgress()
    .log_step("Preparing input..", level = 2)


    seed_ids <- unique(as.character(seed_ids))
    if (length(seed_ids) == 0L) {
        .log_warn("No seeds provided for chromosome-specific DMR detection.")
        return(NULL)
    }
    if (length(seed_beta_index) == 0L) {
        .log_warn("No chromosome-scoped seed indices available for DMR detection.")
        return(NULL)
    }
    if (length(seed_ids) != length(seed_beta_index)) {
        stop("seed_ids and seed_beta_index must have the same length for chromosome-specific DMR detection.")
    }

    .log_step("Subsetting beta matrix for seeds...", level = 3)
    seeds_locs <- as.data.frame(beta_locs[seed_beta_index, , drop = FALSE])
    rownames(seeds_locs) <- seed_ids
    seeds_beta <- beta_handler$getBeta(row_names = seed_beta_index, col_names = beta_col_names_detection)
    rownames(seeds_beta) <- seed_ids

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
            ". This indicates a mismatch between requested site IDs and beta file columns or a parsing issue."
        )
    }

    .log_success("Subset size: ", paste(dim(seeds_beta), collapse = ","), level = 3)
    seeds_beta_handler <- getBetaHandler(
        beta = seeds_beta,
        array = array,
        genome = genome,
        sorted_locs = seeds_locs,
        njobs = njobs
    )
    .log_info("Number of provided chromosome-scoped seeds: ", length(seed_ids), level = 2)
    rm(seeds_beta)
    gc(verbose = FALSE)

    .log_success("Input preparation complete.", level = 2)
    .advanceChromosomeProgress()
    .log_step("Stage 1: Connecting seeds to form initial DMRs..", level = 2)


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
        .log_step("Building seed connectivity array...", level = 3)
        ret <- .buildConnectivityArray(
            beta_handler = seeds_beta_handler,
            beta_locs = seeds_locs,
            pheno = pheno_detection,
            group_inds = group_inds,
            pval_mode_per_group = pval_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            col_names = beta_col_names_detection,
            max_pval = max_pval,
            force_connect_delta_beta = NA_real_, # delta-beta based rescue is applied later during the extension
            covariates = covariates,
            max_lookup_dist = max_lookup_dist,
            entanglement = entanglement,
            aggfun = aggfun,
            ntries = ntries,
            mid_p = mid_p,
            njobs = njobs,
            expansion_windows = NULL,
            max_bridge_gaps = max_bridge_seeds_gaps,
            verbose = verbose
        )
        rm(seeds_beta_handler)
        gc(verbose = FALSE)
        seeds_connectivity_array <- ret$connectivity_array
        pval_mode_per_group <- ret$pval_mode_per_group
        empirical_strategy_per_group <- ret$empirical_strategy_per_group
        .log_success("Seed connectivity array built.", level = 3)
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
        mask <- rep(FALSE, length(seed_ids))
        mask[connected_seeds_segments_starts] <- TRUE
        seeds_connectivity_array$id <- cumsum(mask)
        seeds_connectivity_array$cid <- seeds_connectivity_array$id
        ids <- unique(seeds_connectivity_array$id)
        seeds_connectivity_array$id[connected_seeds_segments_ends] <- NA
        seeds_connectivity_array$seeds <- seed_ids
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
                paste(seed_ids[start_idx:end_idx], collapse = ",")
            }, connected_seeds_segments_starts, connected_seeds_segments_ends, USE.NAMES = FALSE),
            stringsAsFactors = FALSE
        )

        dmrs <- data.frame(
            chr = connected_seeds_segments_chrs,
            start_seed = seed_ids[connected_seeds_segments_starts],
            end_seed = seed_ids[connected_seeds_segments_ends],
            start_seed_pos = connected_seeds_segments_starts_locs,
            end_seed_pos = connected_seeds_segments_ends_locs,
            seeds_num = connected_seeds_segments_lengths,
            stop_connection_reason = stop_reasons,
            id = ids,
            stringsAsFactors = FALSE
        )
        dmrs <- dmrs[dmrs$seeds_num >= min_seeds, , drop = FALSE]
        if (min_seeds > 0) {
            .log_info("Number of DMRs after filtering by min_seeds: ", nrow(dmrs), level = 2)
        }
        dmrs <- merge(dmrs, connected_seeds_connection_corr_pval, by = "id", all.x = TRUE)
        colnames(dmrs)[colnames(dmrs) == "pval"] <- "connection_corr_pval"
        dmrs <- merge(dmrs, dmrs_seeds, by = "id", all.x = TRUE)


        if (nrow(dmrs) == 0) {
            .log_warn("No DMRs remain after filtering based on min_seeds.")
            return(NULL)
        }
        if (getOption("CMEnt.make_debug_dir", FALSE)) {
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

    .log_success("Initial DMRs formed: ", nrow(dmrs), level = 2)
    .advanceChromosomeProgress()
    .log_step("Stage 2: Expanding DMRs on neighborhood sites..", level = 2)

    # Set up progress tracking for DMR expansion
    n_dmrs <- nrow(dmrs)
    stage2_beta_handler <- beta_handler
    stage2_beta_locs <- beta_locs
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
                    level = 3
                )
            } else {
                .log_info("No connectivity windows were generated; Stage 2 will return disconnected sites outside chromosome termini.", level = 2)
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
            expansion_windows <- subset_ret$expansion_windows
            if (!isTRUE(extract_motifs)) {
                beta_locs <- stage2_beta_locs
            }
            .log_info(
                "Stage 2 beta subset contains ",
                format(nrow(stage2_beta_locs), big.mark = ","),
                " sites on chromosome ", unique(stage2_beta_locs$chr),
                level = 2
            )
        }
        .log_step("Building expansion connectivity array..", level = 3)
        ret <- .buildConnectivityArray(
            beta_handler = stage2_beta_handler,
            beta_locs = stage2_beta_locs,
            pheno = pheno_detection,
            group_inds = group_inds,
            pval_mode_per_group = pval_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            col_names = beta_col_names_detection,
            max_pval = max_pval,
            force_connect_delta_beta = ext_site_delta_beta,
            covariates = covariates,
            max_lookup_dist = max_lookup_dist,
            entanglement = entanglement,
            aggfun = aggfun,
            ntries = ntries,
            mid_p = mid_p,
            njobs = njobs,
            expansion_windows = expansion_windows,
            max_bridge_gaps = max_bridge_extension_gaps,
            verbose = verbose
        )
        connectivity_array <- ret$connectivity_array
    }
    .log_success("Connectivity array built.", level = 3)
    .log_info("Number of underlying correlated sites found: ", sum(connectivity_array$connected), level = 2)
    if (getOption("CMEnt.make_debug_dir", FALSE)) {
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
    .log_step("Expanding ", n_dmrs, " DMRs using ", njobs, " jobs...", level = 3)
    dmr_chromosomes <- unique(as.character(dmrs$chr))
    if (length(dmr_chromosomes) != 1L) {
        stop(".findDMRsFromSeedsChr expects Stage 1 DMRs from exactly one chromosome.")
    }
    stage2_chromosomes <- unique(as.character(stage2_beta_locs[, "chr"]))
    if (length(stage2_chromosomes) != 1L || !identical(stage2_chromosomes, dmr_chromosomes)) {
        stop(".findDMRsFromSeedsChr expects Stage 2 beta rows from the same single chromosome as the DMRs.")
    }
    chromosome <- dmr_chromosomes[[1L]]
    .log_info("Processing ", chromosome, level = 2)
    dmrs_to_expand <- dmrs
    locs <- as.data.frame(stage2_beta_locs)
    connectivity <- connectivity_array
    locs_rownames <- rownames(locs)
    locs_idx_map <- setNames(seq_along(locs_rownames), locs_rownames)
    expansion_boundaries <- .buildExpansionBoundaryLookup(connectivity)

    dmr_inds <- seq_len(nrow(dmrs_to_expand))
    default_dmr_chunk_size <- max(1L, ceiling(length(dmr_inds) / max(njobs * 4L, 1L)))
    dmr_chunk_size <- min(default_dmr_chunk_size, length(dmr_inds))
    dmr_chunks <- split(dmr_inds, ceiling(dmr_inds / dmr_chunk_size))

    if (njobs == 1L || length(dmr_chunks) == 1L) {
        ret <- lapply(
            dmr_chunks,
            .expandDMRChunk,
            dmrs = dmrs_to_expand,
            connectivity_array = connectivity,
            min_sites = min_sites,
            locs = locs,
            locs_idx_map = locs_idx_map,
            expansion_boundaries = expansion_boundaries
        )
    } else {
        .setupParallel()
        ret <- future.apply::future_lapply(
            X = dmr_chunks,
            FUN = .expandDMRChunk,
            dmrs = dmrs_to_expand,
            connectivity_array = connectivity,
            min_sites = min_sites,
            locs = locs,
            locs_idx_map = locs_idx_map,
            expansion_boundaries = expansion_boundaries,
            future.seed = TRUE,
            future.stdout = NA,
            future.globals = c(
                ".expandDMRChunk",
                ".expandDMR",
                ".buildExpansionBoundaryLookup",
                "dmrs_to_expand",
                "connectivity",
                "min_sites",
                "locs",
                "locs_idx_map",
                "expansion_boundaries"
            )
        )
        .finalizeParallel()
    }
    ret <- unlist(ret, recursive = FALSE, use.names = FALSE)
    .log_info("Chromosome ", chromosome, ": Number of DMRs processed: ", length(ret), level = 2)
    if (inherits(ret, "try-error")) {
        stop(ret)
    }
    upstream_expansion_length_table <- table(sapply(ret, function(x) x[["upstream_expansion_length"]]))
    downstream_expansion_length_table <- table(sapply(ret, function(x) x[["downstream_expansion_length"]]))
    # sort tables by names (expansion sizes)
    upstream_expansion_length_table <- upstream_expansion_length_table[order(as.integer(names(upstream_expansion_length_table)))]
    downstream_expansion_length_table <- downstream_expansion_length_table[order(as.integer(names(downstream_expansion_length_table)))]
    .log_info("Table of upstream_expansion_length:\n\t", paste(capture.output(upstream_expansion_length_table), collapse = "\n\t"), level = 2)
    .log_info("Table of downstream_expansion_length:\n\t", paste(capture.output(downstream_expansion_length_table), collapse = "\n\t"), level = 2)

    .advanceChromosomeProgress()
    .log_step("Post-processing extended DMRs..", level = 2)

    extended_dmrs <- as.data.frame(do.call(rbind, ret))
    extended_dmrs$end <- as.numeric(extended_dmrs$end)
    extended_dmrs$start <- as.numeric(extended_dmrs$start)
    extended_dmrs$start_seed_pos <- as.numeric(extended_dmrs$start_seed_pos)
    extended_dmrs$end_seed_pos <- as.numeric(extended_dmrs$end_seed_pos)
    extended_dmrs$seeds_num <- as.numeric(extended_dmrs$seeds_num)
    extended_dmrs$connection_corr_pval <- as.numeric(extended_dmrs$connection_corr_pval)


    .checkResult(extended_dmrs, "2", start_col = "start", end_col = "end")
    .log_success("Post-processing complete.", level = 3)
    .log_success("DMR expansion complete.", level = 2)
    if (getOption("CMEnt.make_debug_dir", FALSE)) {
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
    .advanceChromosomeProgress()
    .log_step("Stage 3: Merging overlapping extended DMRs..", level = 2)

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
    .log_info("Frequency of N-DMRs overlap:\n", paste(capture.output(print(table(tqh))), collapse = "\n"), level = 3)
    single_hits <- which(tqh == 1)
    if (length(single_hits) > 0) {
        # get the corresponding indices in the original extended_dmrs_ranges
        .log_info("Copying over ", length(single_hits), " non-overlapping extended DMRs...", level = 3)
        # do it vectorizedly
        agg_df[single_hits, ] <- orig_mcols[sh[qh %in% single_hits], ]
        agg_df[single_hits, "merged_dmrs_num"] <- 1
    }

    multiple_hits <- which(tqh > 1)
    .log_info("Merging ", length(multiple_hits), " overlapping extended DMRs...", level = 3)
    if (length(multiple_hits) > 0L) {
        merge_rows <- .aggregateMergedDMRChunk(
            hit_indices = multiple_hits,
            qh = qh,
            sh = sh,
            orig_mcols = orig_mcols,
            aggfun = aggfun
        )

        merge_idx <- vapply(merge_rows, function(x) x$idx, integer(1))
        if (!identical(merge_idx, as.integer(multiple_hits))) {
            ord <- match(multiple_hits, merge_idx)
            if (anyNA(ord)) {
                stop("Merged overlap aggregation lost rows during chunk combination.")
            }
            merge_rows <- merge_rows[ord]
        }

        char_cols <- c(
            "start_seed", "end_seed", "seeds", "stop_connection_reason",
            "start_site", "end_site", "upstream_expansion_length", "upstream_sites",
            "upstream_expansion_stop_reason", "downstream_expansion_length",
            "downstream_sites", "downstream_expansion_stop_reason"
        )
        num_cols <- c(
            "start_seed_pos", "end_seed_pos", "seeds_num",
            "connection_corr_pval", "merged_dmrs_num"
        )
        for (col in char_cols) {
            agg_df[multiple_hits, col] <- vapply(merge_rows, function(x) as.character(x[[col]]), character(1))
        }
        for (col in num_cols) {
            agg_df[multiple_hits, col] <- vapply(merge_rows, function(x) as.numeric(x[[col]]), numeric(1))
        }
    }

    agg_site_components <- agg_df[, c("upstream_sites", "seeds", "downstream_sites"), drop = FALSE]
    agg_df[, "sites"] <- .collapseMergedDMRsitesChunk(agg_site_components)
    agg_df[, "supporting_sites_num"] <- vapply(agg_df$sites, function(x) {
        length(.splitCsvValues(x))
    }, integer(1))
    if (is.null(all_sites)) {
        agg_df[, "sites_num"] <- vapply(agg_df$sites, function(x) {
            length(.splitCsvValues(x))
        }, integer(1))
    } else {
        agg_df[, "sites_num"] <- match(agg_df$end_site, all_sites) - match(agg_df$start_site, all_sites) + 1
    }
    agg_df[, "id"] <- paste0(seqnames(merged_dmrs_ranges), ":", agg_df$start_site, "-", agg_df$end_site)

    GenomicRanges::mcols(merged_dmrs_ranges) <- agg_df
    .log_success("Overlapping extended DMRs merged: ", length(merged_dmrs_ranges), " resulting DMRs.", level = 2)
    merged_dmrs <- convertToDataFrame(merged_dmrs_ranges)

    if (getOption("CMEnt.make_debug_dir", FALSE)) {
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

    .advanceChromosomeProgress()
    .log_step("Stage 4: Filtering and annotating resulting DMRs..", level = 2)

    if (min_sites > 1 || min_seeds > 1) {
        filtered_dmrs_ranges <- merged_dmrs_ranges[
            GenomicRanges::mcols(merged_dmrs_ranges)$seeds_num >= min_seeds &
                GenomicRanges::mcols(merged_dmrs_ranges)$sites_num >= min_sites
        ]
        .log_info(
            "Keeping ",
            length(filtered_dmrs_ranges),
            " out of ",
            length(merged_dmrs_ranges),
            " with at least ",
            min_seeds,
            " supporting seeds and at least ",
            min_sites,
            " sites in the DMR interval.",
            level = 2
        )
    } else {
        filtered_dmrs_ranges <- merged_dmrs_ranges
    }
    filtered_dmrs <- convertToDataFrame(filtered_dmrs_ranges)

    if (nrow(filtered_dmrs) == 0) {
        .log_warn("No DMRs passed the filtering step.")
        return(NULL)
    }

    if (array_based && min_adj_seeds > min_seeds) {
        .log_step("Calculating site content and adjusted seeds number..", level = 3)
        sites_num_bg <- getSiteBackgroundCounts(filtered_dmrs_ranges, genome)
        sites_num_bg[!is.finite(sites_num_bg) | is.na(sites_num_bg) | sites_num_bg <= 0] <- 1L
        filtered_dmrs$sites_num_bg <- sites_num_bg

        filtered_dmrs$seeds_num_adj <- ceiling(filtered_dmrs$sites_num / filtered_dmrs$sites_num_bg * filtered_dmrs$seeds_num)

        .log_success("site content calculated.", level = 3)
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
            level = 2
        )
        filtered_dmrs <- adj_filtered_dmrs
    } else {
        .log_info("Skipping adjusted seeds number calculation as min_adj_seeds <= min_seeds.", level = 2)
        filtered_dmrs$sites_num_bg <- NA
        filtered_dmrs$seeds_num_adj <- NA
    }
    if (nrow(filtered_dmrs) == 0) {
        .log_warn("No DMRs passed the filtering step.")
        return(NULL)
    }

    annotated_dmrs <- filtered_dmrs

    all_selected_sites <- unique(unlist(base::strsplit(annotated_dmrs$sites, ","), use.names = FALSE))
    all_selected_sites_beta <- beta_handler$getBeta(row_names = all_selected_sites, col_names = beta_col_names)
    .log_step("Calculating per-site beta statistics..", level = 3)
    beta_stats <- .calculateBetaStats(
        beta_values = all_selected_sites_beta,
        pheno = pheno_all,
        aggfun = aggfun
    )
    rm(all_selected_sites_beta)
    gc(verbose = FALSE)
    .log_success("Per-site beta statistics calculated.", level = 3)

    beta_stats <- as.data.frame(beta_stats)
    rownames(beta_stats) <- all_selected_sites
    dmrs_seeds <- base::strsplit(annotated_dmrs$seeds, ",")
    .log_step("Adding DMR delta-beta information..", level = 3)

    dmr_seeds_list <- lapply(dmrs_seeds, as.character)
    dmr_seeds_indices <- unlist(dmr_seeds_list, use.names = FALSE)
    dmr_seeds_groups <- rep(seq_along(dmr_seeds_list), lengths(dmr_seeds_list))

    beta_stats_seeds <- beta_stats[dmr_seeds_indices, , drop = FALSE]
    beta_stats_seeds$dmr_id <- dmr_seeds_groups

    seeds_agg <- .aggregateDMRBetaStats(
        beta_stats_df = beta_stats_seeds,
        aggfun = aggfun,
        njobs = njobs,
        parallel_enabled = FALSE
    )
    seeds_match <- match(seq_len(nrow(annotated_dmrs)), seeds_agg$dmr_id)
    if (anyNA(seeds_match)) {
        stop("Internal error while aggregating seed beta statistics: missing dmr_id rows.")
    }
    seeds_agg <- seeds_agg[seeds_match, , drop = FALSE]

    annotated_dmrs$cases_beta <- seeds_agg$cases_beta
    annotated_dmrs$controls_beta <- seeds_agg$controls_beta
    annotated_dmrs$delta_beta <- annotated_dmrs$cases_beta - annotated_dmrs$controls_beta
    annotated_dmrs$cases_beta_sd <- seeds_agg$cases_beta_sd
    annotated_dmrs$controls_beta_sd <- seeds_agg$controls_beta_sd
    annotated_dmrs$cases_beta_min <- seeds_agg$cases_beta_min
    annotated_dmrs$cases_beta_max <- seeds_agg$cases_beta_max
    annotated_dmrs$controls_beta_min <- seeds_agg$controls_beta_min
    annotated_dmrs$controls_beta_max <- seeds_agg$controls_beta_max

    dmr_sites_list <- base::strsplit(annotated_dmrs$sites, ",")
    dmr_sites <- unlist(dmr_sites_list, use.names = FALSE)
    dmr_sites_groups <- rep(seq_along(dmr_sites_list), lengths(dmr_sites_list))

    beta_stats_sites <- beta_stats[dmr_sites, , drop = FALSE]
    beta_stats_sites$dmr_id <- dmr_sites_groups

    sites_agg <- .aggregateDMRBetaStats(
        beta_stats_df = beta_stats_sites,
        aggfun = aggfun,
        njobs = njobs,
        parallel_enabled = FALSE
    )
    sites_match <- match(seq_len(nrow(annotated_dmrs)), sites_agg$dmr_id)
    if (anyNA(sites_match)) {
        stop("Internal error while aggregating site beta statistics: missing dmr_id rows.")
    }
    sites_agg <- sites_agg[sites_match, , drop = FALSE]

    annotated_dmrs$sites_cases_beta <- sites_agg$cases_beta
    annotated_dmrs$sites_controls_beta <- sites_agg$controls_beta
    annotated_dmrs$sites_delta_beta <- annotated_dmrs$sites_cases_beta - annotated_dmrs$sites_controls_beta
    annotated_dmrs$sites_cases_beta_sd <- sites_agg$cases_beta_sd
    annotated_dmrs$sites_controls_beta_sd <- sites_agg$controls_beta_sd
    annotated_dmrs$sites_cases_beta_min <- sites_agg$cases_beta_min
    annotated_dmrs$sites_cases_beta_max <- sites_agg$cases_beta_max
    annotated_dmrs$sites_controls_beta_min <- sites_agg$controls_beta_min
    annotated_dmrs$sites_controls_beta_max <- sites_agg$controls_beta_max

    .log_success("DMR delta-beta information added.", level = 3)

    if (annotate_with_genes) {
        .advanceChromosomeProgress()
        .log_step("Annotating DMRs with gene information...", level = 2)
        annotated_dmrs <- annotateDMRsWithGenes(annotated_dmrs, genome = genome, njobs = njobs)
        .log_success("DMR annotation completed.", level = 2)
    }

    if (.score_dmrs) {
        .advanceChromosomeProgress()
        .log_step("Scoring DMRs...", level = 2)
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
        .log_success("DMR scoring completed.", level = 2)
    }


    if (is.data.frame(annotated_dmrs)) {
        annotated_dmrs <- convertToGRanges(annotated_dmrs, genome = genome)
    }

    if (extract_motifs) {
        .advanceChromosomeProgress()
        .log_step("Extracting DMR motifs...", level = 2)
        annotated_dmrs <- extractDMRMotifs(annotated_dmrs, genome = genome, array = array, beta_locs = beta_locs)
        .log_success("DMR motifs computed.", level = 2)
    }

    final_dmrs_granges <- annotated_dmrs

    final_dmrs <- convertToDataFrame(final_dmrs_granges)

    .log_info("Final number of chromosome-scoped DMRs: ", nrow(final_dmrs), level = 2)

    gc()
    invisible(final_dmrs_granges)
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
#' @param ext_site_delta_beta Numeric. Minimum absolute delta beta value that will
#' force proximal sites to be treated as connected during Stage 2 expansion,
#' regardless of their correlation p-value. Set to `NA`, `NULL`, or `Inf` to
#' disable this shortcut. A value of `0` means any proximal site with a
#' non-missing case-control delta beta can be force-connected. Default is 0.2.
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
#' @param max_bridge_extension_gaps Integer. Maximum gap size to consider during Stage 2 extension. Default is 1 (i.e., at most 1 consecutive failing site to bridge).
#' @param min_seeds Numeric. Minimum number of connected seeds in a DMR. Minimum is 2. Default is 2.
#' @param min_adj_seeds Numeric. Minimum number of seeds, adjusted by array site density, in a DMR after extension. Minimum is 2. Default is 2.
#' @param min_sites Numeric. Minimum number of sites in a DMR after extension, including the seeds. Minimum is 2. Default is 3.
#' @param aggfun Function or character. Aggregation function to use when calculating delta beta values and p-values of DMRs. Can be "median", "mean", or a function (e.g., median, mean). Default is "median".
#' @param ignored_sample_groups Character vector. Sample groups to ignore during connection and expansion, separated by commas. Can also be "case" or "control". Default is NULL.
#' @param output_prefix Character. Identifier for the output files. If not provided, no output will be saved. Default is NULL.
#' @param njobs Numeric. Number of parallel jobs to use. Default is the number of available cores.
#' @param beta_row_names_file Character. Path to a file containing row names for the beta values. If not provided, row names will be read from the beta file. Default is NULL.
#' @param annotate_with_genes Logical. Whether to annotate DMRs with overlapping genes. Default is TRUE.
#' @param .score_dmrs Logical. Whether to score DMRs based on cross-validated SVM predictions. Default is TRUE.
#' @param extract_motifs Logical. Whether to compute DMRs seeds motifs. Default is TRUE.
#' @param bed_provided Logical. Whether the beta file is provided as a BED file. Default is FALSE. In case the input has a .bed extension, this will be set to TRUE automatically.
#' @param bed_chrom_col Character. Column name for chromosome in the BED file. Default is "chrom".
#' @param bed_start_col Character. Column name for start position in the BED file. Default is "start".
#' @param verbose Numeric. Level of verbosity for logging messages, from 0 (not verbose) to 5 (very very verbose). Default is retrieved from option "CMEnt.verbose".
#' @param .load_debug Logical. If TRUE, enables debug mode for loading beta files. Default is FALSE.
#'
#' @return Data frame of identified DMRs.
#' 
#' @examples
#' \dontrun{
#' beta <- loadExampleInputDataChr5And11("beta")
#' dmps <- loadExampleInputDataChr5And11("dmps")
#' pheno <- loadExampleInputDataChr5And11("pheno")
#' array_type <- loadExampleInputDataChr5And11("array_type")
#' dmrs <- findDMRsFromSeeds(
#'   beta = beta,
#'   seeds = seeds,
#'   pheno = pheno,
#'   array = array_type,
#'   sample_group_col = "Sample_Group"
#' )
#' }
#' @export
findDMRsFromSeeds <- function(
    beta,
    seeds,
    pheno,
    seeds_id_col = NULL,
    sample_group_col = "Sample_Group",
    casecontrol_col = NULL,
    covariates = NULL,
    ext_site_delta_beta = 0.2,
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
    min_sites = 3,
    aggfun = c("median", "mean"),
    ignored_sample_groups = NULL,
    output_prefix = NULL,
    njobs = getOption("CMEnt.njobs", min(8, future::availableCores() - 1)),
    beta_row_names_file = NULL,
    annotate_with_genes = TRUE,
    .score_dmrs = TRUE,
    extract_motifs = TRUE,
    bed_provided = FALSE,
    bed_chrom_col = "chrom",
    bed_start_col = "start",
    verbose = getOption("CMEnt.verbose", 1),
    .load_debug = FALSE
) {
    .emptyOutputs <- function(prefix_base) {
        if (is.null(prefix_base)) {
            return(invisible(NULL))
        }
        for (suffix in c(".dmrs.tsv.gz", ".seeds_beta.tsv.gz")) {
            con <- gzfile(paste0(prefix_base, suffix), "w", compression = 2)
            close(con)
        }
        invisible(NULL)
    }
    .readSeeds <- function(seeds, seeds_id_col) {
        if (is.character(seeds) && length(seeds) == 1) {
            seeds_tsv <- try(as.data.frame(read.table(
                seeds, header = TRUE, sep = "\t", check.names = FALSE,
                quote = "", comment.char = "", row.names = NULL
            )))
        } else if (is.data.frame(seeds)) {
            seeds_tsv <- as.data.frame(seeds)
        } else {
            stop("seeds must be either a file path or a data frame")
        }
        if (inherits(seeds_tsv, "try-error")) {
            return(list(data = NULL, id_col = seeds_id_col))
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
        list(data = seeds_tsv, id_col = seeds_id_col)
    }
    pval_mode <- strex::match_arg(pval_mode, ignore_case = TRUE)
    empirical_strategy <- strex::match_arg(empirical_strategy, ignore_case = TRUE)
    entanglement <- strex::match_arg(entanglement, ignore_case = TRUE)
    options(CMEnt.verbose = verbose, future.globals.maxSize = Inf, "CMEnt.njobs" = njobs)
    if (Sys.info()[["sysname"]] != "Windows") {
        includes <- "#include <sys/wait.h>"
        code <- "int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};"
        wait <- inline::cfunction(body = code, includes = includes, convention = ".C")
        withr::defer(wait())
    }
    .cleanupParallelState()
    withr::defer(.cleanupParallelState(), envir = environment())

    .log_step("Preparing chromosome-sequential DMR input...")
    seeds_ret <- .readSeeds(seeds, seeds_id_col)
    seeds_df <- seeds_ret$data
    seeds_id_col <- seeds_ret$id_col
    if (is.null(seeds_df) || nrow(seeds_df) == 0L) {
        .log_warn("Provided seeds file has no data rows. Not proceeding.")
        .emptyOutputs(output_prefix)
        return(NULL)
    }
    if (!is.null(covariates)) {
        missing_covars <- covariates[!covariates %in% colnames(pheno)]
        if (length(missing_covars) > 0) {
            stop("The following covariates are not present in pheno: ", paste(missing_covars, collapse = ", "))
        }
    }

    array <- .normalizeFindDMRsArray(array)
    requested_genome <- .normalizeFindDMRsGenome(genome)
    genome <- .resolveFindDMRsGenome(beta, array = array, genome = requested_genome, bed_provided = bed_provided)
    .assertDependencyRequirements(
        requirements = .findDMRsDependencyRequirements(
            beta = beta,
            array = array,
            genome = genome,
            annotate_with_genes = annotate_with_genes,
            extract_motifs = extract_motifs,
            bed_provided = bed_provided
        ),
        context = "findDMRsFromSeeds()"
    )
    if (is.null(requested_genome)) {
        .log_info("No genome provided. Using inferred genome: ", genome, ".", level = 2)
    }
    beta_locs <- NULL
    if (!inherits(beta, "BetaHandler") && is.character(beta) && length(beta) == 1 && file.exists(beta)) {
        beta_file_ext <- tools::file_ext(beta)
        if (beta_file_ext == "bed" || bed_provided) {
            bed_provided <- TRUE
            seed_ids <- seeds_df[, seeds_id_col]
            if (!all(grepl("^(chr)?[0-9XYM]+:[0-9]+$", seed_ids))) {
                stop("When providing a bed file as beta input, seed IDs must be in 'chr:pos' format.")
            }
            ret <- readCustomMethylationBedData(
                bed_file = beta, pheno = pheno, genome = genome,
                chrom_col = bed_chrom_col, start_col = bed_start_col
            )
            beta <- ret$tabix_file
            beta_locs <- ret$locations
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
    beta <- NULL
    array_based <- beta_handler$isArrayBased()

    if (!is.function(aggfun)) {
        aggfun_choice <- strex::match_arg(aggfun, ignore_case = TRUE)
        aggfun <- switch(aggfun_choice,
            median = stats::median,
            mean = mean
        )
    }
    stopifnot(sample_group_col %in% colnames(pheno))
    stopifnot(!is.null(max_pval))
    stopifnot(!is.null(min_seeds))
    stopifnot(!is.null(min_sites))
    if (min_seeds < 2 && min_sites < 2) {
        stop("min_seeds or min_sites must be at least 2, to define a DMR")
    }
    stopifnot(!is.null(min_adj_seeds))
    if (min_adj_seeds < 2) {
        stop("min_adj_seeds must be at least 2, to define a DMR")
    }
    stopifnot(!is.null(max_lookup_dist))
    ext_site_delta_beta <- .normalizeForceConnectDeltaBeta(
        ext_site_delta_beta,
        arg_name = "ext_site_delta_beta"
    )
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
    if (is.null(casecontrol_col)) {
        pheno_all[, .CASE_CONTROL_COL] <- ifelse(
            pheno_all[, sample_group_col] == levels(as.factor(pheno_all[, sample_group_col]))[1],
            0, 1
        )
    } else {
        pheno_all[, .CASE_CONTROL_COL] <- as.numeric(pheno_all[, casecontrol_col])
    }
    ignored_sample_groups_chr <- if (is.null(ignored_sample_groups)) character(0) else {
        x <- trimws(unlist(base::strsplit(ignored_sample_groups, ",")))
        x[nzchar(x)]
    }
    samples_selection_mask <- !(pheno_all[, sample_group_col] %in% ignored_sample_groups_chr)
    if ("case" %in% ignored_sample_groups_chr) {
        samples_selection_mask <- samples_selection_mask & (pheno_all[, .CASE_CONTROL_COL] != 1)
    }
    if ("control" %in% ignored_sample_groups_chr) {
        samples_selection_mask <- samples_selection_mask & (pheno_all[, .CASE_CONTROL_COL] != 0)
    }
    beta_col_names_detection <- beta_col_names[samples_selection_mask]
    if (length(beta_col_names_detection) < 2) {
        stop("At least two samples are required after applying ignored_sample_groups.")
    }
    pheno_detection <- pheno_all[beta_col_names_detection, , drop = FALSE]
    sample_groups <- factor(pheno_detection[, sample_group_col])
    group_inds <- split(seq_along(sample_groups), sample_groups)
    pval_mode_per_group <- rep(pval_mode, length.out = length(unique(pheno_detection[[sample_group_col]])))
    names(pval_mode_per_group) <- unique(pheno_detection[[sample_group_col]])
    empirical_strategy_per_group <- rep(empirical_strategy, length.out = length(unique(pheno_detection[[sample_group_col]])))
    names(empirical_strategy_per_group) <- unique(pheno_detection[[sample_group_col]])

    output_prefix_base <- output_prefix
    output_prefix_dot <- NULL
    if (!is.null(output_prefix_base)) {
        dir.create(dirname(output_prefix_base), showWarnings = FALSE, recursive = TRUE)
        output_prefix_dot <- paste0(output_prefix_base, ".")
        saveRDS(
            list(pheno = pheno_detection, genome = genome, array = array, sample_group_col = sample_group_col),
            file = paste0(output_prefix_base, ".meta.rds")
        )
    }

    beta_locs <- beta_handler$getBetaLocs()
    beta_chr <- as.character(beta_locs[, "chr"])
    beta_start <- suppressWarnings(as.numeric(beta_locs[, "start"]))
    if (anyNA(beta_chr) || any(!nzchar(beta_chr))) {
        stop("Beta locations contain missing chromosome labels.", call. = FALSE)
    }
    if (anyNA(beta_start)) {
        stop("Beta locations contain missing or non-numeric start positions.", call. = FALSE)
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
            stop("Beta locations are not sorted within chromosome ", beta_chr[bad_idx], call. = FALSE)
        }
    }

    beta_locs_rownames <- .explicitRowNames(beta_locs)
    use_numeric_sequencing_rows <- !array_based && is.null(beta_locs_rownames)
    beta_row_names <- if (use_numeric_sequencing_rows) NULL else beta_handler$getBetaRowNames()
    if (use_numeric_sequencing_rows) {
        seed_beta_index <- .matchSequencingIdsToBeta(seeds_df[, seeds_id_col], beta_chr, beta_start)
        if (all(is.na(seed_beta_index))) {
            stop("None of the IDs in seeds_id_col match the beta genomic locations.")
        }
        if (anyNA(seed_beta_index)) {
            seeds_df <- seeds_df[!is.na(seed_beta_index), , drop = FALSE]
            seed_beta_index <- seed_beta_index[!is.na(seed_beta_index)]
        }
        seeds_df$.__beta_row_index__ <- seed_beta_index
        seeds_df <- seeds_df[order(seeds_df$.__beta_row_index__), , drop = FALSE]
        keep_unique <- !duplicated(seeds_df[, seeds_id_col])
        seed_ids <- as.character(seeds_df[keep_unique, seeds_id_col])
        seed_beta_index <- seeds_df[keep_unique, ".__beta_row_index__"]
    } else {
        if (!all(seeds_df[, seeds_id_col] %in% beta_row_names)) {
            if (!any(seeds_df[, seeds_id_col] %in% beta_row_names)) {
                seeds_id_col_found <- NULL
                for (col in colnames(seeds_df)) {
                    if (all(seeds_df[, col] %in% beta_row_names)) {
                        seeds_id_col_found <- col
                        break
                    }
                }
                if (is.null(seeds_id_col_found)) {
                    stop("None of the IDs in seeds_id_col match the beta file row names.")
                }
                .log_warn("Switching seeds_id_col from '", seeds_id_col, "' to '", seeds_id_col_found, "'.")
                seeds_id_col <- seeds_id_col_found
            }
            seeds_df <- seeds_df[seeds_df[, seeds_id_col] %in% beta_row_names, , drop = FALSE]
        }
        seeds_df <- seeds_df[orderByLoc(seeds_df[, seeds_id_col], genomic_locs = beta_locs), , drop = FALSE]
        seed_ids <- unique(as.character(seeds_df[, seeds_id_col]))
        seed_ids <- seed_ids[orderByLoc(seed_ids, genome = genome, genomic_locs = beta_locs)]
        seed_beta_index <- seed_ids
    }
    if (length(seed_ids) == 0L) {
        stop("No seeds remain after filtering against beta locations.")
    }
    seeds_locs <- as.data.frame(beta_locs[seed_beta_index, , drop = FALSE])
    rownames(seeds_locs) <- seed_ids
    seed_chr <- as.character(seeds_locs[, "chr"])
    chromosomes <- unique(seed_chr)
    .log_info("Processing ", length(chromosomes), " chromosome(s): ", paste(chromosomes, collapse = ", "), level = 1)

    beta_row_ids_all <- if (is.null(beta_locs_rownames)) {
        seq_len(nrow(beta_locs))
    } else {
        beta_locs_rownames
    }
    chr_results <- vector("list", length(chromosomes))
    names(chr_results) <- chromosomes
    for (chr in chromosomes) {
        .log_step("Processing chromosome ", chr, "...", level = 1)
        chr_beta_idx <- which(beta_chr == chr)
        chr_row_ids <- if (is.null(beta_locs_rownames)) chr_beta_idx else beta_row_ids_all[chr_beta_idx]
        chr_handler <- beta_handler$subset(row_names = chr_row_ids, col_names = beta_col_names)
        chr_seed_mask <- seed_chr == chr
        chr_seed_ids <- seed_ids[chr_seed_mask]
        chr_seed_beta_index <- if (use_numeric_sequencing_rows) {
            as.integer(seed_beta_index[chr_seed_mask] - chr_beta_idx[1L] + 1L)
        } else {
            chr_seed_ids
        }
        chr_ret <- withCallingHandlers(
            .findDMRsFromSeedsChr(
                beta_handler = chr_handler,
                seed_ids = chr_seed_ids,
                seed_beta_index = chr_seed_beta_index,
                pheno_detection = pheno_detection,
                group_inds = group_inds,
                pval_mode_per_group = pval_mode_per_group,
                empirical_strategy_per_group = empirical_strategy_per_group,
                beta_col_names_detection = beta_col_names_detection,
                ext_site_delta_beta = ext_site_delta_beta,
                array = array,
                genome = genome,
                max_pval = max_pval,
                entanglement = entanglement,
                ntries = ntries,
                mid_p = mid_p,
                max_lookup_dist = max_lookup_dist,
                expansion_window = expansion_window,
                max_bridge_seeds_gaps = max_bridge_seeds_gaps,
                max_bridge_extension_gaps = max_bridge_extension_gaps,
                min_seeds = min_seeds,
                min_adj_seeds = min_adj_seeds,
                min_sites = min_sites,
                aggfun = aggfun,
                njobs = njobs,
                verbose = verbose,
                .load_debug = .load_debug,
                pheno_all = pheno_all,
                beta_col_names = beta_col_names,
                sample_group_col = sample_group_col,
                covariates = covariates,
                annotate_with_genes = annotate_with_genes,
                .score_dmrs = .score_dmrs,
                extract_motifs = extract_motifs
            ),
            warning = function(w) {
                if (grepl("No DMRs|No connectivity windows", conditionMessage(w))) {
                    invokeRestart("muffleWarning")
                }
            }
        )
        if (!is.null(chr_ret) && length(chr_ret) > 0L) {
            chr_results[[chr]] <- chr_ret
        }
        rm(chr_handler, chr_ret)
        gc(verbose = FALSE)
    }
    chr_results <- chr_results[lengths(chr_results) > 0L]
    if (length(chr_results) == 0L) {
        .log_warn("No DMRs remain after filtering based on min_seeds.")
        .emptyOutputs(output_prefix_base)
        return(NULL)
    }
    final_dmrs_granges <- if (length(chr_results) == 1L) {
        chr_results[[1L]]
    } else {
        Reduce(c, unname(chr_results))
    }
    final_ord <- order(as.character(GenomicRanges::seqnames(final_dmrs_granges)), GenomicRanges::start(final_dmrs_granges), GenomicRanges::end(final_dmrs_granges))
    final_dmrs_granges <- final_dmrs_granges[final_ord]

    if (!is.null(output_prefix_dot)) {
        viewer_sites <- unique(unlist(lapply(S4Vectors::mcols(final_dmrs_granges)$sites, .splitCsvValues), use.names = FALSE))
        viewer_beta <- beta_handler$getBeta(row_names = viewer_sites, col_names = beta_col_names_detection)
        viewer_file <- paste0(output_prefix_dot, "seeds_beta.tsv.gz")
        gz <- gzfile(viewer_file, "w")
        write.table(viewer_beta, gz, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
        close(gz)

        final_dmrs <- convertToDataFrame(final_dmrs_granges)
        encoded_dmrs <- .encodeNonTabularColumns(final_dmrs)
        dmrs_file <- paste0(output_prefix_dot, "dmrs.tsv.gz")
        gz <- gzfile(dmrs_file, "w")
        write.table(
            encoded_dmrs$data,
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

    invisible(final_dmrs_granges)
}
