.CASE_CONTROL_COL <- "__casecontrol__" # nolint

#' @keywords internal
#' @noRd
.mergeGenomicWindows <- function(windows, return_granges = FALSE) {
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
    if (return_granges) {
        return(gr)
    }
    data.frame(
        chr = as.character(GenomicRanges::seqnames(gr)),
        start = GenomicRanges::start(gr),
        end = GenomicRanges::end(gr)
    )
}

#' @keywords internal
#' @noRd
.buildWindowsFromDMRs <- function(dmrs, expansion_window) {
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
.normalizeBuildDMRsArray <- function(array) {
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
    supported_arrays <- c("450K", "27K", "EPIC", "EPICv2", "Mouse")
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
.normalizeBuildDMRsGenome <- function(genome) {
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
.resolveBuildDMRsGenome <- function(beta, array = NULL, genome = NULL, bed_provided = FALSE) {
    genome <- .normalizeBuildDMRsGenome(genome)
    if (!is.null(genome)) {
        return(genome)
    }

    array <- .normalizeBuildDMRsArray(array)

    if (inherits(beta, "BetaHandler")) {
        beta_genome <- tryCatch(.normalizeBuildDMRsGenome(beta$genome), error = function(e) NULL)
        if (!is.null(beta_genome)) {
            return(beta_genome)
        }
        beta_array <- tryCatch(.normalizeBuildDMRsArray(beta$array), error = function(e) NULL)
        if (!is.null(beta_array) && beta_array %in% c("450K", "27K", "EPIC")) {
            return("hg19")
        }
        if (identical(beta_array, "Mouse")) {
            return("mm10")
        }
        return("hg38")
    }

    if (is_bsseq(beta)) {
        gr <- GenomicRanges::granges(beta)
        bsseq_genome <- unique(as.character(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(gr))))
        bsseq_genome <- bsseq_genome[!is.na(bsseq_genome) & nzchar(bsseq_genome)]
        bsseq_genome <- bsseq_genome[tolower(bsseq_genome) %in% c("hg38", "hg19", "hs1", "mm10", "mm39")]
        if (length(bsseq_genome) > 0L) {
            return(.normalizeBuildDMRsGenome(bsseq_genome[1]))
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
    if (identical(array, "Mouse")) {
        return("mm10")
    }

    "hg38"
}


#' @keywords internal
#' @noRd
.forceConnectDeltaBetaEnabled <- function(ext_site_delta_beta) {
    !is.na(ext_site_delta_beta) && is.finite(ext_site_delta_beta)
}


#' @keywords internal
#' @noRd
.usesFileBackedBSseq <- function(beta_handler) {
    if (!inherits(beta_handler, "BetaHandler")) {
        return(FALSE)
    }
    bsseq_object <- tryCatch(
        beta_handler$.__enclos_env__$private$.bsseq_object,
        error = function(e) NULL
    )
    !is.null(bsseq_object) && !.bsseqIsInMemory(bsseq_object)
}


#' @keywords internal
#' @noRd
.usesBSseqBackend <- function(beta_handler) {
    if (!inherits(beta_handler, "BetaHandler")) {
        return(FALSE)
    }
    bsseq_object <- tryCatch(
        beta_handler$.__enclos_env__$private$.bsseq_object,
        error = function(e) NULL
    )
    !is.null(bsseq_object)
}


#' @keywords internal
#' @noRd
.bsseqBackendGRanges <- function(beta_handler) {
    bsseq_object <- tryCatch(
        beta_handler$.__enclos_env__$private$.bsseq_object,
        error = function(e) NULL
    )
    if (is.null(bsseq_object)) {
        return(NULL)
    }
    GenomicRanges::granges(bsseq_object)
}


#' @keywords internal
#' @noRd
.seedCoordinatesFromSeeds <- function(seeds_df, seeds_id_col) {
    if (all(c("chr", "start") %in% colnames(seeds_df))) {
        return(list(
            chr = as.character(seeds_df$chr),
            start = suppressWarnings(as.integer(as.numeric(seeds_df$start)))
        ))
    }
    ids <- as.character(seeds_df[, seeds_id_col])
    parsed <- regexec("^([^:]+):([0-9]+)$", ids)
    pieces <- regmatches(ids, parsed)
    ok <- lengths(pieces) == 3L
    seed_chr <- rep(NA_character_, length(ids))
    seed_start <- rep(NA_integer_, length(ids))
    if (any(ok)) {
        seed_chr[ok] <- vapply(pieces[ok], `[[`, character(1), 2L)
        seed_start[ok] <- suppressWarnings(as.integer(vapply(pieces[ok], `[[`, character(1), 3L)))
    }
    list(chr = seed_chr, start = seed_start)
}


#' @keywords internal
#' @noRd
.chromosomeRunIndex <- function(beta_chr, require_unique = FALSE) {
    if (inherits(beta_chr, "Rle")) {
        chr_values <- as.character(S4Vectors::runValue(beta_chr))
        chr_lengths <- as.integer(S4Vectors::runLength(beta_chr))
    } else {
        beta_chr <- as.character(beta_chr)
        if (length(beta_chr) == 0L) {
            return(data.frame(chr = character(0), start = integer(0), end = integer(0)))
        }
        chr_runs <- rle(beta_chr)
        chr_values <- as.character(chr_runs$values)
        chr_lengths <- as.integer(chr_runs$lengths)
    }

    if (length(chr_values) == 0L) {
        return(data.frame(chr = character(0), start = integer(0), end = integer(0)))
    }
    if (isTRUE(require_unique) && anyDuplicated(chr_values)) {
        dup_chr <- chr_values[duplicated(chr_values)][1L]
        stop(
            "Beta locations are not grouped by chromosome: ", dup_chr,
            " appears in multiple blocks. Ensure the beta input is ordered by chromosome and genomic start position.",
            call. = FALSE
        )
    }
    chr_end <- cumsum(chr_lengths)
    data.frame(
        chr = chr_values,
        start = as.integer(chr_end - chr_lengths + 1L),
        end = as.integer(chr_end),
        stringsAsFactors = FALSE
    )
}


#' @keywords internal
#' @noRd
.chromosomeRunFor <- function(chr, chromosome_runs) {
    chr <- as.character(chr)[1L]
    if (is.na(chr) || !nzchar(chr)) {
        return(NULL)
    }
    candidates <- chr
    if (startsWith(chr, "chr")) {
        candidates <- c(candidates, sub("^chr", "", chr))
    } else {
        candidates <- c(candidates, paste0("chr", chr))
    }
    run_idx <- match(candidates, chromosome_runs$chr)
    run_idx <- run_idx[!is.na(run_idx)]
    if (length(run_idx) == 0L) {
        return(NULL)
    }
    chromosome_runs[run_idx, , drop = FALSE]
}


#' @keywords internal
#' @noRd
.buildDMRsChromosomeTasks <- function(beta_handler,
                                      beta_chr,
                                      beta_locs_rownames,
                                      chromosomes,
                                      seed_ids,
                                      seed_beta_index,
                                      seed_chr,
                                      beta_col_names,
                                      use_numeric_sequencing_rows,
                                      chromosome_runs = .chromosomeRunIndex(beta_chr, require_unique = TRUE)) {
    tasks <- vector("list", length(chromosomes))
    names(tasks) <- chromosomes

    for (i in seq_along(chromosomes)) {
        chr <- chromosomes[[i]]
        chr_run <- .chromosomeRunFor(chr, chromosome_runs)
        if (is.null(chr_run)) {
            stop("No beta rows found for chromosome ", chr, call. = FALSE)
        }
        if (nrow(chr_run) != 1L) {
            stop("Beta rows for chromosome ", chr, " are not contiguous.", call. = FALSE)
        }
        chr_row_start <- chr_run$start[[1L]]
        chr_row_end <- chr_run$end[[1L]]

        chr_row_ids <- if (isTRUE(use_numeric_sequencing_rows) || is.null(beta_locs_rownames)) {
            NULL
        } else {
            beta_locs_rownames[seq.int(chr_row_start, chr_row_end)]
        }
        chr_seed_mask <- seed_chr == chr
        chr_seed_ids <- seed_ids[chr_seed_mask]
        chr_seed_beta_index <- if (isTRUE(use_numeric_sequencing_rows)) {
            as.integer(seed_beta_index[chr_seed_mask] - chr_row_start + 1L)
        } else {
            chr_seed_ids
        }

        tasks[[i]] <- list(
            chr = chr,
            row_start = chr_row_start,
            row_end = chr_row_end,
            seed_ids = chr_seed_ids,
            seed_beta_index = chr_seed_beta_index
        )
        if (!is.null(chr_row_ids)) {
            tasks[[i]]$row_ids <- chr_row_ids
        }
    }

    tasks
}


#' @keywords internal
#' @noRd
.matchSequencingIdsToBeta <- function(ids, beta_chr, beta_start,
                                      chromosome_runs = .chromosomeRunIndex(beta_chr)) {
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
        chr_run <- .chromosomeRunFor(chr, chromosome_runs)
        if (is.null(chr_run)) {
            next
        }
        beta_pos <- unlist(
            Map(seq.int, chr_run$start, chr_run$end),
            use.names = FALSE
        )
        hit <- match(id_start[query_pos], beta_start[beta_pos])
        idx[ok_pos[query_pos]] <- as.integer(beta_pos[hit])
    }
    idx
}


#' @keywords internal
#' @noRd
.matchSequencingCoordinatesToBeta <- function(seed_chr, seed_start, beta_chr, beta_start,
                                             chromosome_runs = .chromosomeRunIndex(beta_chr)) {
    seed_chr <- as.character(seed_chr)
    seed_start <- suppressWarnings(as.integer(seed_start))
    idx <- rep(NA_integer_, length(seed_chr))
    ok <- !is.na(seed_chr) & nzchar(seed_chr) & !is.na(seed_start)
    if (!any(ok)) {
        return(idx)
    }
    for (chr in unique(seed_chr[ok])) {
        query_pos <- which(ok & seed_chr == chr)
        chr_run <- .chromosomeRunFor(chr, chromosome_runs)
        if (is.null(chr_run)) {
            next
        }
        beta_pos <- unlist(
            Map(seq.int, chr_run$start, chr_run$end),
            use.names = FALSE
        )
        hit <- match(seed_start[query_pos], beta_start[beta_pos])
        idx[query_pos] <- as.integer(beta_pos[hit])
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
.readMemoryBytesFile <- function(path) {
    if (!file.exists(path)) {
        return(NA_real_)
    }
    value <- tryCatch(readLines(path, n = 1L, warn = FALSE), error = function(e) character(0))
    if (length(value) != 1L) {
        return(NA_real_)
    }
    value <- trimws(value)
    if (identical(value, "max")) {
        return(Inf)
    }
    value <- suppressWarnings(as.numeric(value))
    if (!is.na(value) && value > 0) value else NA_real_
}


#' @keywords internal
#' @noRd
.meminfoAvailableRamBytes <- function(path = "/proc/meminfo") {
    if (!file.exists(path)) {
        return(NA_real_)
    }
    meminfo <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
    mem_available <- grep("^MemAvailable:", meminfo, value = TRUE)
    if (length(mem_available) != 1L) {
        return(NA_real_)
    }
    kb <- suppressWarnings(as.numeric(strsplit(mem_available, "\\s+")[[1]][2L]))
    if (is.finite(kb) && !is.na(kb) && kb > 0) kb * 1024 else NA_real_
}


#' @keywords internal
#' @noRd
.cgroupAvailableRamBytes <- function(
    limit_paths = c("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"),
    current_paths = c("/sys/fs/cgroup/memory.current", "/sys/fs/cgroup/memory/memory.usage_in_bytes")
) {
    limits <- vapply(limit_paths, .readMemoryBytesFile, numeric(1))
    currents <- vapply(current_paths, .readMemoryBytesFile, numeric(1))
    limits <- limits[is.finite(limits) & limits > 0 & limits < 2^60]
    currents <- currents[is.finite(currents) & currents >= 0]
    if (length(limits) == 0L || length(currents) == 0L) {
        return(NA_real_)
    }
    available <- min(limits) - max(currents)
    if (is.finite(available) && !is.na(available) && available > 0) available else 1
}


#' @keywords internal
#' @noRd
.availableRamBytes <- function(
    default_gb = 2,
    meminfo_path = "/proc/meminfo",
    cgroup_limit_paths = c("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"),
    cgroup_current_paths = c("/sys/fs/cgroup/memory.current", "/sys/fs/cgroup/memory/memory.usage_in_bytes")
) {
    candidates <- c(
        .cgroupAvailableRamBytes(cgroup_limit_paths, cgroup_current_paths),
        .meminfoAvailableRamBytes(meminfo_path)
    )
    candidates <- candidates[is.finite(candidates) & !is.na(candidates) & candidates > 0]
    if (length(candidates) > 0L) {
        return(min(candidates))
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
.poolConnectivitySplits <- function(splits_mat,
                                    split_weights = NULL,
                                    chunk_size,
                                    min_splits = 1L) {
    .rowsBefore <- function(x, i) {
        if (i <= 1L) {
            return(x[integer(0), , drop = FALSE])
        }
        x[seq_len(i - 1L), , drop = FALSE]
    }
    .rowsAfter <- function(x, i) {
        if (i >= nrow(x)) {
            return(x[integer(0), , drop = FALSE])
        }
        x[seq.int(i + 1L, nrow(x)), , drop = FALSE]
    }
    .valsBefore <- function(x, i) {
        if (i <= 1L) {
            return(x[integer(0)])
        }
        x[seq_len(i - 1L)]
    }
    .valsAfter <- function(x, i) {
        if (i >= length(x)) {
            return(x[integer(0)])
        }
        x[seq.int(i + 1L, length(x))]
    }

    if (is.null(splits_mat)) {
        return(matrix(numeric(0), ncol = 2L))
    }
    if (is.null(dim(splits_mat))) {
        splits_mat <- matrix(splits_mat, ncol = 2L, byrow = TRUE)
    } else {
        splits_mat <- as.matrix(splits_mat)
    }
    if (nrow(splits_mat) == 0L) {
        return(matrix(numeric(0), ncol = 2L))
    }
    if (ncol(splits_mat) != 2L) {
        stop("splits_mat must have exactly two columns.")
    }
    storage.mode(splits_mat) <- "integer"
    n_input_splits <- nrow(splits_mat)
    if (is.null(split_weights)) {
        split_weights <- as.integer(splits_mat[, 2L] - splits_mat[, 1L] + 1L)
    } else {
        split_weights <- as.integer(split_weights)
        if (length(split_weights) != n_input_splits) {
            stop("split_weights must have one value per split.")
        }
    }
    split_weights[is.na(split_weights) | split_weights < 1L] <- 1L

    chunk_size_eff <- as.integer(chunk_size)
    if (!is.finite(chunk_size_eff) || is.na(chunk_size_eff) || chunk_size_eff < 1L) {
        chunk_size_eff <- 1L
    }
    min_splits_eff <- as.integer(min_splits)[1L]
    if (!is.finite(min_splits_eff) || is.na(min_splits_eff) || min_splits_eff < 1L) {
        min_splits_eff <- 1L
    }
    min_splits_eff <- min(min_splits_eff, n_input_splits)

    pooled <- vector("list", n_input_splits)
    pooled_weights <- integer(n_input_splits)
    pooled_n <- 0L
    cur_start <- as.integer(splits_mat[1L, 1L])
    cur_end <- as.integer(splits_mat[1L, 2L])
    cur_weight <- as.integer(split_weights[1L])
    merges_remaining <- max(0L, n_input_splits - min(min_splits_eff, n_input_splits))
    if (n_input_splits > 1L) {
        for (i in seq.int(2L, n_input_splits)) {
            next_start <- as.integer(splits_mat[i, 1L])
            next_end <- as.integer(splits_mat[i, 2L])
            next_weight <- as.integer(split_weights[i])
            contiguous <- next_start <= cur_end + 1L
            if (!contiguous && next_start == cur_end + 2L) {
                contiguous <- TRUE
            }
            fits_budget <- cur_weight + next_weight <= chunk_size_eff
            if (contiguous && fits_budget && merges_remaining > 0L) {
                cur_end <- max(cur_end, next_end)
                cur_weight <- cur_weight + next_weight
                merges_remaining <- merges_remaining - 1L
            } else {
                pooled_n <- pooled_n + 1L
                pooled[[pooled_n]] <- c(cur_start, cur_end)
                pooled_weights[[pooled_n]] <- cur_weight
                cur_start <- next_start
                cur_end <- next_end
                cur_weight <- next_weight
            }
        }
    }
    pooled_n <- pooled_n + 1L
    pooled[[pooled_n]] <- c(cur_start, cur_end)
    pooled_weights[[pooled_n]] <- cur_weight

    pooled <- do.call(rbind, pooled[seq_len(pooled_n)])
    pooled_weights <- pooled_weights[seq_len(pooled_n)]
    if (is.null(dim(pooled))) {
        pooled <- matrix(pooled, ncol = 2L, byrow = TRUE)
    }
    storage.mode(pooled) <- "integer"

    while (nrow(pooled) < min_splits_eff) {
        spans <- as.integer(pooled[, 2L] - pooled[, 1L] + 1L)
        candidates <- which(spans > 1L & pooled_weights > 1L)
        if (length(candidates) == 0L) {
            break
        }
        split_ind <- candidates[which.max(pooled_weights[candidates])]
        left_span <- as.integer(floor(spans[split_ind] / 2L))
        if (left_span < 1L) {
            break
        }
        left_start <- pooled[split_ind, 1L]
        left_end <- left_start + left_span - 1L
        right_start <- left_end + 1L
        right_end <- pooled[split_ind, 2L]
        left_weight <- as.integer(floor(pooled_weights[split_ind] * left_span / spans[split_ind]))
        left_weight <- max(1L, min(pooled_weights[split_ind] - 1L, left_weight))
        right_weight <- pooled_weights[split_ind] - left_weight

        pooled <- rbind(
            .rowsBefore(pooled, split_ind),
            c(left_start, left_end),
            c(right_start, right_end),
            .rowsAfter(pooled, split_ind)
        )
        pooled_weights <- c(
            .valsBefore(pooled_weights, split_ind),
            left_weight,
            right_weight,
            .valsAfter(pooled_weights, split_ind)
        )
        storage.mode(pooled) <- "integer"
    }

    pooled
}


#' @keywords internal
#' @noRd
.subsetStage2BetaToWindows <- function(beta_handler, beta_locs, col_names, expansion_windows, njobs = 1L) {
    if (nrow(expansion_windows) == 0L) {
        return(NULL)
    }
    n_sites <- nrow(beta_locs)
    if (n_sites == 0L) {
        return(NULL)
    }
    beta_start <- as.integer(beta_locs[, "start"])
    keep <- rep(FALSE, n_sites)
    if (nrow(expansion_windows) > 0L) {
        win_ir <- IRanges::IRanges(
            start = as.integer(expansion_windows$start),
            end = as.integer(expansion_windows$end)
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
    beta_locs_rownames <- {
        if (
            is.integer(beta_locs_rownames_info) &&
                length(beta_locs_rownames_info) == 2L &&
                is.na(beta_locs_rownames_info[1L]) &&
                beta_locs_rownames_info[2L] < 0L
        ) {
            NULL
        } else {
            rownames(beta_locs)
        }
    }
    if (is.null(beta_locs_rownames) || .usesBSseqBackend(beta_handler)) {
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
        beta_locs = subset_locs
    )
}


.chooseTestingOptions <- function(group, mat, mask, testing_mode, empirical_strategy) {
    n_sites <- nrow(mat)
    x_mat_full <- mat[seq_len(n_sites - 1L), , drop = FALSE] # Sites i
    y_mat_full <- mat[seq.int(2L, n_sites), , drop = FALSE] # Sites i+1
    x_mat <- x_mat_full[mask, , drop = FALSE]
    y_mat <- y_mat_full[mask, , drop = FALSE]
    valid_pairs <- !is.na(x_mat) & !is.na(y_mat)
    n_valid <- rowSums(valid_pairs)

    if (testing_mode == "auto") {
        auto_diag <- .summarizeCorrAssumptions(
            x_mat = x_mat,
            y_mat = y_mat,
            n_valid = n_valid
        )
        testing_mode <- if (auto_diag$use_parametric) "parametric" else "empirical"
        .log_info(
            "Auto-selected testing_mode='", testing_mode, "' for group '", group,
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
    if (testing_mode == "empirical" && empirical_strategy == "auto") {
        if (ncol(x_mat) >= 6) {
            .log_info(
                "Group '", group, "' has ", ncol(x_mat),
                " samples. Using 'montecarlo' empirical ",
                "strategy for faster computation with sufficient sample size.",
                level = 2
            )
            empirical_strategy <- "montecarlo"
        } else {
            .log_info(
                "Group '", group, "' has only ", ncol(x_mat),
                " samples. Using 'permutations' empirical strategy for more",
                " accurate p-value estimation with small sample size.",
                level = 2
            )
            empirical_strategy <- "permutations"
        }
    }
    list(testing_mode = testing_mode, empirical_strategy = empirical_strategy)
}

#' @keywords internal
#' @noRd
.rowPairCorrelationStats <- function(mat, start_inds, end_inds) {
    n_pairs <- length(start_inds)
    if (n_pairs == 0L) {
        return(list(
            x_means = numeric(0),
            y_means = numeric(0),
            sum_xy = numeric(0),
            sum_x2 = numeric(0),
            sum_y2 = numeric(0),
            n_valid = integer(0)
        ))
    }

    x_sum <- numeric(n_pairs)
    y_sum <- numeric(n_pairs)
    x_count <- integer(n_pairs)
    y_count <- integer(n_pairs)
    n_valid <- integer(n_pairs)

    for (j in seq_len(ncol(mat))) {
        x <- mat[start_inds, j]
        y <- mat[end_inds, j]
        x_has_na <- anyNA(x)
        y_has_na <- anyNA(y)

        if (x_has_na) {
            x_ok <- !is.na(x)
            x_sum[x_ok] <- x_sum[x_ok] + x[x_ok]
            x_count <- x_count + as.integer(x_ok)
        } else {
            x_sum <- x_sum + x
            x_count <- x_count + 1L
        }

        if (y_has_na) {
            y_ok <- !is.na(y)
            y_sum[y_ok] <- y_sum[y_ok] + y[y_ok]
            y_count <- y_count + as.integer(y_ok)
        } else {
            y_sum <- y_sum + y
            y_count <- y_count + 1L
        }

        if (x_has_na || y_has_na) {
            if (!x_has_na) {
                x_ok <- rep(TRUE, n_pairs)
            }
            if (!y_has_na) {
                y_ok <- rep(TRUE, n_pairs)
            }
            n_valid <- n_valid + as.integer(x_ok & y_ok)
        } else {
            n_valid <- n_valid + 1L
        }
    }

    x_means <- x_sum / x_count
    y_means <- y_sum / y_count
    sum_xy <- numeric(n_pairs)
    sum_x2 <- numeric(n_pairs)
    sum_y2 <- numeric(n_pairs)

    for (j in seq_len(ncol(mat))) {
        x <- mat[start_inds, j]
        y <- mat[end_inds, j]
        x_centered <- x - x_means
        y_centered <- y - y_means

        if (anyNA(x)) {
            x_centered[is.na(x)] <- 0
        }
        if (anyNA(y)) {
            y_centered[is.na(y)] <- 0
        }

        sum_xy <- sum_xy + x_centered * y_centered
        sum_x2 <- sum_x2 + x_centered * x_centered
        sum_y2 <- sum_y2 + y_centered * y_centered
    }

    list(
        x_means = x_means,
        y_means = y_means,
        sum_xy = sum_xy,
        sum_x2 = sum_x2,
        sum_y2 = sum_y2,
        n_valid = n_valid
    )
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
    testing_mode_per_group,
    empirical_strategy_per_group,
    col_names = NULL,
    max_pval = 0.05,
    ext_site_delta_beta = NA_real_,
    covariates = NULL,
    covariate_models = NULL,
    max_lookup_dist = 1000,
    entanglement = "strong",
    aggfun = stats::median,
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
    .log_step(
        "Testing connectivity for pairs ", pair_start, "-", pair_end,
        if (!is.null(checked_pairs)) paste0(" (", sum(pair_mask), " checked pair(s))") else "",
        " within beta rows [", min(inds), ", ", max(inds), "].",
        level = 3
    )
    x <- .testConnectivityBatch(
        sites_beta = as.matrix(chunk_beta),
        group_inds = group_inds,
        pheno = pheno,
        covariates = covariates,
        covariate_models = covariate_models,
        max_pval = max_pval,
        ext_site_delta_beta = ext_site_delta_beta,
        max_lookup_dist = max_lookup_dist,
        site_starts = site_starts,
        n_sites = length(inds),
        entanglement = entanglement,
        aggfun = aggfun,
        testing_mode_per_group = testing_mode_per_group,
        empirical_strategy_per_group = empirical_strategy_per_group,
        ntries = ntries,
        mid_p = mid_p,
        check_non_overlapping = !is.null(checked_pairs)
    )
    rm(chunk_beta)
    .log_success("Connectivity tested.", level = 3)
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
.prepareGroupCovariateModels <- function(pheno, group_inds, covariates = NULL) {
    if (is.null(covariates) || length(covariates) == 0L) {
        return(NULL)
    }
    ret <- lapply(group_inds, function(idx) {
        .prepareCovariateModel(pheno = pheno[idx, , drop = FALSE], covariates = covariates)
    })
    names(ret) <- names(group_inds)
    ret
}


#' @keywords internal
#' @noRd
.buildCASinglePass <- function(
    beta_handler,
    beta_locs = NULL,
    pheno,
    group_inds,
    testing_mode_per_group,
    empirical_strategy_per_group,
    col_names = NULL,
    max_pval = 0.05,
    ext_site_delta_beta = NA_real_,
    covariates = NULL,
    covariate_models = NULL,
    max_lookup_dist = 1000,
    entanglement = "strong",
    aggfun = stats::median,
    ntries = 500,
    mid_p = TRUE,
    njobs = 1,
    expansion_windows = NULL,
    connectivity_array = NULL,
    ugap = 0L,
    dgap = 0L,
    recheck = NULL,
    splits = NULL,
    memory_njobs = njobs,
    verbose = 1
) {
    if (is.null(covariate_models)) {
        covariate_models <- .prepareGroupCovariateModels(pheno, group_inds, covariates)
    }
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
        if (.forceConnectDeltaBetaEnabled(ext_site_delta_beta)) {
            ret$delta_beta <- rep(NA_real_, n_sites)
        }
        return(
            list(
                connectivity_array = ret, splits = NULL,
                testing_mode_per_group = testing_mode_per_group,
                empirical_strategy_per_group = empirical_strategy_per_group
            )
        )
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
        njobs = memory_njobs,
        n_pairs = max(1L, n_sites - 1L)
    )
    .log_info(
        "Connectivity chunk size: ", chunk_size,
        " pair(s), derived from available RAM, ", n_cols_for_chunk,
        " sample column(s), and ", max(1L, as.integer(memory_njobs)), " budgeted job(s).",
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
        if (.forceConnectDeltaBetaEnabled(ext_site_delta_beta)) {
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
        .poolConnectivitySplits(
            splits_mat = splits_mat,
            split_weights = split_weights,
            chunk_size = chunk_size,
            min_splits = 1
        )
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
        if (nrow(expansion_windows) == 0L) {
            return(data.frame(start_pair = integer(0), end_pair = integer(0)))
        }
        beta_start <- as.integer(beta_locs[, "start"])
        if (nrow(expansion_windows) == 0L || n_sites < 2L) {
            return(data.frame(start_pair = integer(0), end_pair = integer(0)))
        }
        win_ir <- IRanges::IRanges(
            start = as.integer(round(as.numeric(expansion_windows$start))),
            end = as.integer(round(as.numeric(expansion_windows$end)))
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
            return(
                list(
                    connectivity_array = connectivity_array,
                    splits = NULL,
                    testing_mode_per_group = testing_mode_per_group,
                    empirical_strategy_per_group = empirical_strategy_per_group
                )
            )
        }
        connectivity_array[pair_ranges$start_pair[pair_ranges$start_pair > 1] - 1, "reason"] <- "end-of-input"
        connectivity_array[pair_ranges$end_pair, "reason"] <- "end-of-input"
        splits <- .chunkPairRanges(pair_ranges)
        if (nrow(splits) > 0L) {
            old_n <- nrow(splits)
            splits <- .poolSplits(splits)
            if (old_n != nrow(splits)) {
                .log_info(
                    "Adjusted initial connectivity chunks from ", old_n, " to ", nrow(splits),
                    " (target chunk size: ", as.integer(chunk_size),
                    "; available jobs: ", max(1L, as.integer(njobs)), ").",
                    level = 3
                )
            }
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
            return(
                list(
                    connectivity_array = connectivity_array,
                    splits = splits,
                    testing_mode_per_group = testing_mode_per_group,
                    empirical_strategy_per_group = empirical_strategy_per_group,
                    recheck = integer(0)
                )
            )
        }
        run_ends <- run_ends[run_mask]
        run_starts <- run_starts[run_mask]
        if (ugap > 0L) {
            # The following indices will be re-checked
            checked_pairs <- data.frame(before = run_starts - ugap, after = run_starts)
            .log_info(
                "Re-assessing connectivity for ",
                nrow(checked_pairs),
                " site pairs at the upstream edges of existing connected",
                " regions to see if we can bridge small gaps.",
                level = 3
            )
        } else if (dgap > 0L) {
            checked_pairs <- data.frame(before = run_ends + 1, after = run_ends + dgap + 1)
            if (!is.null(max_lookup_dist) && is.finite(max_lookup_dist) && nrow(checked_pairs) > 0L) {
                site_starts_for_shift <- as.integer(beta_locs[, "start"])
                candidate_dist <- site_starts_for_shift[checked_pairs$after] - site_starts_for_shift[checked_pairs$before]
                shift_mask <- !is.na(candidate_dist) & candidate_dist > max_lookup_dist
                if (any(shift_mask)) {
                    shifted_before <- run_ends + 1L - dgap
                    shifted_after <- run_ends + 1L
                    valid_shift <- shift_mask & shifted_before >= 1L & shifted_after <= n_sites
                    if (any(valid_shift)) {
                        shifted_dist <- rep(NA_real_, length(valid_shift))
                        valid_shift_inds <- which(valid_shift)
                        shifted_dist[valid_shift_inds] <- site_starts_for_shift[shifted_after[valid_shift_inds]] -
                            site_starts_for_shift[shifted_before[valid_shift_inds]]
                        valid_shift <- valid_shift & !is.na(shifted_dist) & shifted_dist <= max_lookup_dist
                        shifted_edge_end <- shifted_after - 1L
                        shifted_has_eoi <- rep(FALSE, length(valid_shift))
                        shift_width <- unique(shifted_edge_end[valid_shift] - shifted_before[valid_shift] + 1L)
                        if (length(shift_width) > 0L) {
                            for (offset in seq_len(max(shift_width)) - 1L) {
                                edge_idx <- shifted_before + offset
                                in_range <- valid_shift & edge_idx <= shifted_edge_end
                                shifted_has_eoi[in_range] <- shifted_has_eoi[in_range] |
                                    connectivity_array[edge_idx[in_range], "reason"] == "end-of-input"
                            }
                        }
                        valid_shift <- valid_shift & !shifted_has_eoi
                        if (any(valid_shift)) {
                            checked_pairs$before[valid_shift] <- shifted_before[valid_shift]
                            checked_pairs$after[valid_shift] <- shifted_after[valid_shift]
                            .log_info(
                                "Back-shifted ", sum(valid_shift),
                                " downstream bridge assessment(s) after the forward bridged",
                                " candidate exceeded max_lookup_dist.",
                                level = 3
                            )
                        }
                    }
                }
            }
            .log_info(
                "Re-assessing connectivity for ", nrow(checked_pairs),
                " site pairs at the downstream edges of existing connected regions",
                " to see if we can bridge small gaps.",
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
        # Accumulate sparse bridge re-check chunks without dropping below the worker count.
        if (nrow(splits) > 0L) {
            old_n <- nrow(splits)
            splits <- .poolSplits(splits, split_weights = ninds_in_splits)
            if (old_n != nrow(splits)) {
                .log_info(
                    "Adjusted bridge re-check chunks from ", old_n, " to ", nrow(splits),
                    " (target chunk size: ", as.integer(chunk_size),
                    "; available jobs: ", max(1L, as.integer(njobs)), ").",
                    level = 3
                )
            }
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
    use_bsseq_numeric_index <- .usesBSseqBackend(beta_handler)
    handler_row_names_for_numeric <- if (use_bsseq_numeric_index) {
        NULL
    } else {
        tryCatch(beta_handler$getBetaRowNames(), error = function(e) NULL)
    }
    numeric_row_index_matches_locs <- use_bsseq_numeric_index ||
        is.null(handler_row_names_for_numeric) ||
        (
            length(handler_row_names_for_numeric) >= n_sites &&
                identical(handler_row_names_for_numeric[seq_len(n_sites)], beta_row_ids_full)
        )
    prefer_numeric_row_index <- numeric_row_index_matches_locs
    use_numeric_row_index <- prefer_numeric_row_index ||
        is.null(beta_row_ids_full) ||
        length(beta_row_ids_full) != n_sites ||
        anyNA(beta_row_ids_full) ||
        any(!nzchar(beta_row_ids_full))
    beta_row_ids <- if (use_numeric_row_index) NULL else beta_row_ids_full
    if (any(testing_mode_per_group == "auto") || any(empirical_strategy_per_group[testing_mode_per_group == "empirical"] == "auto")) {
        .log_info(
            "Selecting p-value computation mode for each group using the first",
            " chunk as a pilot.",
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
            dists <- site_starts[seq.int(2L, s)] - site_starts[seq_len(s - 1L)]
            exceeded_dist <- dists > max_lookup_dist
        } else {
            exceeded_dist <- rep(FALSE, n_sites - 1L)
        }
        nexdist_mask <- !exceeded_dist
        groups_options <- lapply(
            names(group_inds),
            function(x) {
                idx <- group_inds[[x]]
                chunk_m <- .transformBeta(
                    first_chunk,
                    pheno = pheno,
                    covariates = covariates,
                    covariate_model = covariate_models[[x]],
                    cols = idx
                )
                .chooseTestingOptions(
                    group = x,
                    mat = chunk_m,
                    mask = nexdist_mask,
                    testing_mode = testing_mode_per_group[[x]],
                    empirical_strategy = empirical_strategy_per_group[[x]]
                )
            }
        )
        empirical_strategy_per_group <- vapply(groups_options, function(opt) opt$empirical_strategy, character(1))
        names(empirical_strategy_per_group) <- names(group_inds)
        testing_mode_per_group <- vapply(groups_options, function(opt) opt$testing_mode, character(1))
        names(testing_mode_per_group) <- names(group_inds)
        rm(first_chunk, site_starts, exceeded_dist, nexdist_mask, groups_options)
    }
    beta_start_vec <- as.integer(beta_locs[, "start"])
    .runConnectivityChunk <- function(split,
                                      checked_pairs_local = checked_pairs,
                                      worker_beta_handler = beta_handler,
                                      worker_use_numeric_row_index = use_numeric_row_index,
                                      worker_beta_row_ids = beta_row_ids,
                                      worker_beta_row_ids_offset = 0L) {
        .connectivityChunkWorker(
            split = split,
            beta_handler = worker_beta_handler,
            beta_start_vec = beta_start_vec,
            group_inds = group_inds,
            pheno = pheno,
            testing_mode_per_group = testing_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            col_names = col_names,
            max_pval = max_pval,
            ext_site_delta_beta = ext_site_delta_beta,
            covariates = covariates,
            covariate_models = covariate_models,
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
    connectivity_state <- new.env(parent = emptyenv())
    connectivity_state$connected_vec <- connectivity_array$connected
    connectivity_state$pval_vec <- connectivity_array$pval
    connectivity_state$reason_vec <- connectivity_array$reason
    fail_col <- if ("first_failing_group" %in% names(connectivity_array))
        "first_failing_group" else if ("failing_groups" %in% names(connectivity_array))
        "failing_groups" else NULL
    connectivity_state$fail_vec <- if (!is.null(fail_col)) connectivity_array[[fail_col]] else NULL
    connectivity_state$delta_vec <- if ("delta_beta" %in% names(connectivity_array)) connectivity_array$delta_beta else NULL

    connectivity_state$bridge_mask <- rep(FALSE, n_sites)
    connectivity_state$recheck <- integer(0)

    .updatePvalMin <- function(idx, values) {
        if (length(idx) == 0L || length(values) == 0L) {
            return(invisible(NULL))
        }
        pval_vec <- connectivity_state$pval_vec
        keep <- !is.na(idx) & idx >= 1L & idx <= length(pval_vec) & !is.na(values)
        idx <- as.integer(idx[keep])
        values <- as.numeric(values[keep])
        if (length(idx) == 0L) {
            return(invisible(NULL))
        }
        for (i in unique(idx)) {
            new_pval <- min(values[idx == i], na.rm = TRUE)
            old_pval <- pval_vec[i]
            pval_vec[i] <- if (is.na(old_pval)) new_pval else min(old_pval, new_pval)
        }
        connectivity_state$pval_vec <- pval_vec
        invisible(NULL)
    }

    .applyChunkResult <- function(item) {
        x <- item$result
        if (nrow(x) == 0L) {
            return(invisible(NULL))
        }
        if (is.null(checked_pairs)) {
            # First pass: full overwrite
            idx <- item$pair_start:item$pair_end
            connectivity_state$connected_vec[idx] <- x$connected
            connectivity_state$pval_vec[idx] <- x$pval
            connectivity_state$reason_vec[idx] <- x$reason
            if (!is.null(fail_col) && fail_col %in% names(x)) {
                fail_vec <- connectivity_state$fail_vec
                fail_vec[idx] <- x[[fail_col]]
                connectivity_state$fail_vec <- fail_vec
            }
            if (!is.null(connectivity_state$delta_vec) && "delta_beta" %in% names(x)) {
                delta_vec <- connectivity_state$delta_vec
                delta_vec[idx] <- x$delta_beta
                connectivity_state$delta_vec <- delta_vec
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
                update_pval <- x$pval[update_m]
                update_was_connected <- connectivity_state$connected_vec[update_idx]
                newly_updated_mask <- !update_was_connected
                if (any(newly_updated_mask)) {
                    newly_updated_idx <- update_idx[newly_updated_mask]
                    connectivity_state$connected_vec[newly_updated_idx] <- x$connected[update_m][newly_updated_mask]
                    connectivity_state$reason_vec[newly_updated_idx] <- x$reason[update_m][newly_updated_mask]
                    if (!is.null(fail_col) && fail_col %in% names(x)) {
                        fail_vec <- connectivity_state$fail_vec
                        fail_vec[newly_updated_idx] <- x[[fail_col]][update_m][newly_updated_mask]
                        connectivity_state$fail_vec <- fail_vec
                    }
                    if (!is.null(connectivity_state$delta_vec) && "delta_beta" %in% names(x)) {
                        delta_vec <- connectivity_state$delta_vec
                        delta_vec[newly_updated_idx] <- x$delta_beta[update_m][newly_updated_mask]
                        connectivity_state$delta_vec <- delta_vec
                    }
                }
                gap <- if (ugap > 0L) ugap else dgap
                bridge_idx <- rep(update_idx, each = gap) + rep.int(seq.int(0L, gap - 1L), length(update_idx))
                bridge_pval <- rep(update_pval, each = gap)
                bridge_keep <- bridge_idx >= 1L & bridge_idx <= n_sites
                bridge_idx <- bridge_idx[bridge_keep]
                bridge_pval <- bridge_pval[bridge_keep]
                if (length(bridge_idx) > 0L) {
                    bridge_was_connected <- connectivity_state$connected_vec[bridge_idx]
                    update_bridge_match <- match(bridge_idx, update_idx)
                    matched_update_bridge <- !is.na(update_bridge_match)
                    bridge_was_connected[matched_update_bridge] <- update_was_connected[update_bridge_match[matched_update_bridge]]
                    newly_connected_mask <- !bridge_was_connected
                    newly_connected_idx <- bridge_idx[newly_connected_mask]
                    if (length(newly_connected_idx) > 0L) {
                        .updatePvalMin(newly_connected_idx, bridge_pval[newly_connected_mask])
                        connectivity_state$bridge_mask <- replace(connectivity_state$bridge_mask, newly_connected_idx, TRUE)
                        connectivity_state$recheck <- c(connectivity_state$recheck, newly_connected_idx)
                    }
                }
            }
        }
        invisible(NULL)
    }

    .bpBatchConnectivity <- function(batch_splits, batch_checked_pairs = checked_pairs) {
        BiocParallel::bplapply(
            X = batch_splits,
            FUN = .connectivityChunkWorker,
            beta_handler = beta_handler,
            beta_start_vec = beta_start_vec,
            group_inds = group_inds,
            pheno = pheno,
            testing_mode_per_group = testing_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            col_names = col_names,
            max_pval = max_pval,
            ext_site_delta_beta = ext_site_delta_beta,
            covariates = covariates,
            covariate_models = covariate_models,
            max_lookup_dist = max_lookup_dist,
            entanglement = entanglement,
            aggfun = aggfun,
            ntries = ntries,
            mid_p = mid_p,
            checked_pairs = batch_checked_pairs,
            use_numeric_row_index = use_numeric_row_index,
            beta_row_ids = beta_row_ids,
            beta_row_ids_offset = 0L,
            BPPARAM = .makeBiocParallelParam(
                njobs,
                n_tasks = length(batch_splits)
            )
        )
    }

    pairs_to_eval <- if (is.null(checked_pairs)) {
        sum(as.integer(splits[, 2]) - as.integer(splits[, 1]) + 1L)
    } else {
        nrow(checked_pairs)
    }
    min_pairs_for_parallel <- suppressWarnings(as.integer(getOption("CMEnt.min_pairs_for_parallel", 5000L)))
    if (is.na(min_pairs_for_parallel) || min_pairs_for_parallel < 1L) {
        min_pairs_for_parallel <- 1L
    }
    use_parallel <- njobs > 1L && nrow(splits) > 1L && pairs_to_eval >= min_pairs_for_parallel

    if (!use_parallel) {
        for (split_ind in seq_len(nrow(splits))) {
            .applyChunkResult(.runConnectivityChunk(splits[split_ind, ], checked_pairs_local = checked_pairs))
        }
    } else {
        batch_splits <- lapply(seq_len(nrow(splits)), function(i) as.integer(splits[i, ]))
        ret <- .bpBatchConnectivity(batch_splits, batch_checked_pairs = checked_pairs)
        for (item in ret) {
            .applyChunkResult(item)
        }
        rm(ret, batch_splits)
    }

    connected_vec <- connectivity_state$connected_vec
    pval_vec <- connectivity_state$pval_vec
    reason_vec <- connectivity_state$reason_vec
    fail_vec <- connectivity_state$fail_vec
    delta_vec <- connectivity_state$delta_vec
    bridge_mask <- connectivity_state$bridge_mask
    recheck <- sort(unique(connectivity_state$recheck))

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
        testing_mode_per_group = testing_mode_per_group,
        empirical_strategy_per_group = empirical_strategy_per_group,
        recheck = recheck
    )
}

.buildConnectivityArray <- function(
    beta_handler,
    beta_locs = NULL,
    pheno,
    group_inds,
    testing_mode_per_group,
    empirical_strategy_per_group,
    col_names = NULL,
    max_pval = 0.05,
    ext_site_delta_beta = NA_real_,
    covariates = NULL,
    covariate_models = NULL,
    max_lookup_dist = 1000,
    entanglement = "strong",
    aggfun = stats::median,
    ntries = 500,
    mid_p = TRUE,
    njobs = 1,
    expansion_windows = NULL,
    max_bridge_gaps = 0,
    memory_njobs = njobs,
    verbose = getOption("CMEnt.verbose", 1L)
) {
    if (is.null(covariate_models)) {
        covariate_models <- .prepareGroupCovariateModels(pheno, group_inds, covariates)
    }
    connectivity_array <- NULL
    splits <- NULL
    for (gap in seq(0L, max_bridge_gaps)) {
        if (gap == 0L) {
            .log_info("Building connectivity array with no gap bridging.", level = 2)
        } else {
            .log_info(
                "Building bridged connectivity array allowing up to ", gap,
                " gap(s) between connected seeds.",
                level = 2
            )
        }
        build_args <- list(
            beta_handler = beta_handler,
            beta_locs = beta_locs,
            group_inds = group_inds,
            col_names = col_names,
            pheno = pheno,
            covariates = covariates,
            covariate_models = covariate_models,
            max_lookup_dist = max_lookup_dist,
            max_pval = max_pval,
            ext_site_delta_beta = ext_site_delta_beta,
            entanglement = entanglement,
            aggfun = aggfun,
            testing_mode_per_group = testing_mode_per_group,
            empirical_strategy_per_group = empirical_strategy_per_group,
            ntries = ntries,
            mid_p = mid_p,
            njobs = njobs,
            expansion_windows = expansion_windows,
            connectivity_array = connectivity_array,
            splits = splits,
            memory_njobs = memory_njobs,
            verbose = verbose
        )
        .buildCASinglePassWithGaps <- function(build_args, gap) {
            build_args$ugap <- gap
            build_args$dgap <- 0
            build_ret <- do.call(.buildCASinglePass, build_args)
            build_args$connectivity_array <- build_ret$connectivity_array
            urecheck <- build_ret$recheck
            build_args$ugap <- 0
            build_args$dgap <- gap
            build_ret <- do.call(.buildCASinglePass, build_args)
            drecheck <- build_ret$recheck
            list(
                recheck = sort(unique(c(urecheck, drecheck))),
                connectivity_array = build_ret$connectivity_array,
                splits = build_ret$splits,
                testing_mode_per_group = build_ret$testing_mode_per_group,
                empirical_strategy_per_group = build_ret$empirical_strategy_per_group
            )
        }
        .buildCABridgeFixedPoint <- function(build_args, max_gap) {
            active_recheck <- NULL
            base_splits <- build_args$splits
            pass <- 0L
            max_passes <- max(1L, nrow(build_args$connectivity_array) * max(1L, as.integer(max_gap)) * 2L)
            repeat {
                cycle_recheck <- integer(0)
                connected_before <- sum(build_args$connectivity_array$connected, na.rm = TRUE)
                for (bridge_gap in seq_len(max_gap)) {
                    pass <- pass + 1L
                    if (pass > max_passes) {
                        stop(
                            "Bridge fixed-point connectivity did not converge after ",
                            max_passes, " pass(es). This indicates a cycling recheck set."
                        )
                    }
                    build_args$ugap <- NULL
                    build_args$dgap <- NULL
                    build_args$recheck <- active_recheck
                    build_args$splits <- base_splits
                    .log_info(
                        "Bridge fixed-point pass ", pass,
                        " for gap ", bridge_gap, "/", max_gap,
                        if (is.null(active_recheck)) {
                            " over all connected runs."
                        } else {
                            paste0(
                                " over ", length(active_recheck), " touched edge(s)."
                            )
                        },
                        level = 3
                    )
                    build_ret <- .buildCASinglePassWithGaps(build_args, bridge_gap)
                    build_args$connectivity_array <- build_ret$connectivity_array
                    build_args$testing_mode_per_group <- build_ret$testing_mode_per_group
                    build_args$empirical_strategy_per_group <- build_ret$empirical_strategy_per_group
                    cycle_recheck <- sort(unique(c(cycle_recheck, build_ret$recheck)))
                }
                if (length(cycle_recheck) == 0L) {
                    build_ret$splits <- base_splits
                    return(build_ret)
                }
                connected_after <- sum(build_args$connectivity_array$connected, na.rm = TRUE)
                if (connected_after <= connected_before) {
                    .log_warn(
                        "Bridge fixed-point connectivity stopped because a recheck pass made no progress."
                    )
                    build_ret$recheck <- integer(0)
                    build_ret$splits <- base_splits
                    return(build_ret)
                }
                active_recheck <- cycle_recheck
            }
        }
        if (gap > 0L) {
            build_ret <- .buildCABridgeFixedPoint(build_args, gap)
            connectivity_array <- build_ret$connectivity_array
        } else {
            build_ret <- do.call(.buildCASinglePass, build_args)
            connectivity_array <- build_ret$connectivity_array
        }

        splits <- build_ret$splits
        testing_mode_per_group <- build_ret$testing_mode_per_group
        empirical_strategy_per_group <- build_ret$empirical_strategy_per_group
        .log_info("Connectivity array built with gap allowance of ", gap, " (", sum(connectivity_array$connected), " connected sites).", level = 2)
    }
    list(
        connectivity_array = connectivity_array, splits = splits,
        testing_mode_per_group = testing_mode_per_group,
        empirical_strategy_per_group = empirical_strategy_per_group
    )
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

    dmrs_chunk <- dmrs[dmr_inds, , drop = FALSE]
    n_dmrs <- nrow(dmrs_chunk)
    n_locs <- nrow(locs)
    locs_rownames <- names(locs_idx_map)
    projected_positions <- as.integer(locs[, "start"])
    connectivity_reason <- as.character(connectivity_array[["reason"]])

    start_seed <- as.character(dmrs_chunk[["start_seed"]])
    end_seed <- as.character(dmrs_chunk[["end_seed"]])
    start_idx <- unname(as.integer(locs_idx_map[start_seed]))
    end_idx <- unname(as.integer(locs_idx_map[end_seed]))
    if (anyNA(start_idx)) {
        stop("Could not find the start site ", start_seed[which(is.na(start_idx))[1L]], " in the beta file row names.")
    }
    if (anyNA(end_idx)) {
        stop("Could not find the end site ", end_seed[which(is.na(end_idx))[1L]], " in the beta file row names.")
    }

    upstream_exp <- start_idx
    upstream_stop_reason <- rep(NA_character_, n_dmrs)
    upstream_lookup <- start_idx - 1L
    has_upstream_lookup <- upstream_lookup >= 1L
    upstream_stop_reason[!has_upstream_lookup] <- "end-of-input"
    upstream_fail_edge <- rep(NA_integer_, n_dmrs)
    upstream_fail_edge[has_upstream_lookup] <- expansion_boundaries$previous_failed_edge[upstream_lookup[has_upstream_lookup]]
    upstream_no_failure <- has_upstream_lookup & (is.na(upstream_fail_edge) | upstream_fail_edge < 1L)
    upstream_exp[upstream_no_failure] <- 1L
    upstream_stop_reason[upstream_no_failure] <- "end-of-input"
    upstream_has_failure <- has_upstream_lookup & !upstream_no_failure
    upstream_exp[upstream_has_failure] <- upstream_fail_edge[upstream_has_failure] + 1L
    upstream_stop_reason[upstream_has_failure] <- expansion_boundaries$reason[upstream_fail_edge[upstream_has_failure]]

    downstream_exp <- end_idx
    downstream_stop_reason <- rep(NA_character_, n_dmrs)
    downstream_lookup <- end_idx
    has_downstream_lookup <- downstream_lookup < n_locs
    downstream_stop_reason[!has_downstream_lookup] <- "end-of-input"
    downstream_fail_edge <- rep(NA_integer_, n_dmrs)
    downstream_fail_edge[has_downstream_lookup] <- expansion_boundaries$next_failed_edge[downstream_lookup[has_downstream_lookup]]
    downstream_no_failure <- has_downstream_lookup & (is.na(downstream_fail_edge) | downstream_fail_edge > n_locs)
    downstream_exp[downstream_no_failure] <- n_locs
    downstream_stop_reason[downstream_no_failure] <- "end-of-input"
    downstream_has_failure <- has_downstream_lookup & !downstream_no_failure
    downstream_exp[downstream_has_failure] <- downstream_fail_edge[downstream_has_failure]
    downstream_stop_reason[downstream_has_failure] <- expansion_boundaries$reason[downstream_fail_edge[downstream_has_failure]]

    too_short <- (downstream_exp - upstream_exp + 1L) < min_sites
    upstream_stop_reason[too_short] <- "min-sites-not-reached"
    downstream_stop_reason[too_short] <- "min-sites-not-reached"

    expansion_kept_indices <- function(local_inds) {
        if (length(local_inds) == 0L) {
            return(integer(0))
        }
        bridged <- connectivity_reason[local_inds] == "bridged"
        local_inds[!bridged]
    }

    expansion_site_ids <- function(local_inds) {
        kept <- expansion_kept_indices(local_inds)
        if (length(kept) == 0L) {
            return("")
        }
        paste(locs_rownames[kept], collapse = ",")
    }

    upstream_sites <- character(n_dmrs)
    downstream_sites <- character(n_dmrs)
    upstream_expansion_length <- integer(n_dmrs)
    downstream_expansion_length <- integer(n_dmrs)
    for (i in seq_len(n_dmrs)) {
        upstream_candidate <- if (upstream_exp[[i]] <= (start_idx[[i]] - 1L)) {
            seq.int(upstream_exp[[i]], start_idx[[i]] - 1L)
        } else {
            integer(0)
        }
        downstream_candidate <- if ((end_idx[[i]] + 1L) <= downstream_exp[[i]]) {
            seq.int(end_idx[[i]] + 1L, downstream_exp[[i]])
        } else {
            integer(0)
        }
        upstream_sites[[i]] <- expansion_site_ids(upstream_candidate)
        downstream_sites[[i]] <- expansion_site_ids(downstream_candidate)
        upstream_expansion_length[[i]] <- length(expansion_kept_indices(upstream_candidate))
        downstream_expansion_length[[i]] <- length(expansion_kept_indices(downstream_candidate))
    }

    dmrs_chunk[["start_site"]] <- locs_rownames[upstream_exp]
    dmrs_chunk[["end_site"]] <- locs_rownames[downstream_exp]
    dmrs_chunk[["start"]] <- projected_positions[upstream_exp]
    dmrs_chunk[["end"]] <- projected_positions[downstream_exp]
    dmrs_chunk[["upstream_expansion_stop_reason"]] <- upstream_stop_reason
    dmrs_chunk[["upstream_sites"]] <- upstream_sites
    dmrs_chunk[["upstream_expansion_length"]] <- upstream_expansion_length
    dmrs_chunk[["downstream_expansion_stop_reason"]] <- downstream_stop_reason
    dmrs_chunk[["downstream_sites"]] <- downstream_sites
    dmrs_chunk[["downstream_expansion_length"]] <- downstream_expansion_length

    dmrs_chunk
}


#' Vectorized connectivity testing for consecutive site pairs, given their beta values
#'
#' @return Data frame with columns: connected, pval, delta_beta, reason, first_failing_group, stop_reason
#' @keywords internal
#' @noRd
.testConnectivityBatch <- function(sites_beta, group_inds, pheno,
                                   testing_mode_per_group,
                                   empirical_strategy_per_group,
                                   max_pval, covariates = NULL,
                                   covariate_models = NULL,
                                   ext_site_delta_beta = NA_real_,
                                   max_lookup_dist = NULL,
                                   site_starts = NULL,
                                   n_sites = NULL,
                                   entanglement = "weak",
                                   aggfun = mean,
                                   ntries = 0, mid_p = FALSE,
                                   check_non_overlapping = FALSE) {
    if (is.null(n_sites)) {
        n_sites <- nrow(sites_beta)
    } else {
        n_sites <- suppressWarnings(as.integer(n_sites)[1L])
        if (!is.finite(n_sites) || is.na(n_sites) || n_sites < 0L) {
            stop("n_sites must be a non-negative integer.", call. = FALSE)
        }
    }
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
        .log_info(
            "Using max p-value of ",
            max_pval_corrected,
            "(group multi-testing corrected) for connectivity testing.",
            level = 4
        )
    } else {
        # null hypothesis: at least one group must be significant -> independent testing
        max_pval_corrected <- max_pval
        .log_info(
            "Using max p-value of ",
            max_pval_corrected,
            " for connectivity testing.",
            level = 4
        )
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
    .log_info(
        sum(exceeded_dist),
        " out of ",
        n_pairs,
        " site pairs exceeded the maximum lookup distance and will be marked as not connected.",
        level = 4
    )

    high_delta <- rep(FALSE, n_pairs)
    if (.forceConnectDeltaBetaEnabled(ext_site_delta_beta) && length(unique(pheno[, .CASE_CONTROL_COL])) > 1) {
        # Compute this inexpensive case-control effect-size screen up front so
        # proximal high-delta pairs can bypass correlation testing entirely.
        site2_beta_mat <- sites_beta[end_pair_inds, , drop = FALSE]
        case_betas <- apply(site2_beta_mat[, pheno[, .CASE_CONTROL_COL] == 1, drop = FALSE], 1, aggfun, na.rm = TRUE)
        control_betas <- apply(site2_beta_mat[, pheno[, .CASE_CONTROL_COL] == 0, drop = FALSE], 1, aggfun, na.rm = TRUE)
        delta_betas <- case_betas - control_betas
        high_delta <- is.finite(delta_betas) & abs(delta_betas) >= ext_site_delta_beta & nexdist_mask
        connected[high_delta] <- TRUE
        reasons[high_delta] <- "abs(delta_beta)>=ext_site_delta_beta"
    }
    corr_mask <- nexdist_mask & !high_delta

    if (!strict_mode && n_groups > 0) {
        per_group_reasons <- matrix("", nrow = n_groups, ncol = n_pairs)
        per_group_reasons[, exceeded_dist] <- "exceeded max distance"
        per_group_reasons[, high_delta] <- "abs(delta_beta)>=ext_site_delta_beta"
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

        group_m <- .transformBeta(
            sites_beta,
            pheno = pheno,
            covariates = covariates,
            covariate_model = covariate_models[[g]],
            cols = idx
        )

        # Compute only correlation-eligible pair rows; this avoids retaining
        # full pair matrices for distance-filtered or delta-beta-rescued chunks.
        corr_start_inds <- start_pair_inds[corr_mask]
        corr_end_inds <- end_pair_inds[corr_mask]
        sn_pairs <- length(corr_start_inds)
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

        corr_stats <- .rowPairCorrelationStats(group_m, corr_start_inds, corr_end_inds)
        x_means <- corr_stats$x_means
        y_means <- corr_stats$y_means
        sum_xy <- corr_stats$sum_xy
        sum_x2 <- corr_stats$sum_x2
        sum_y2 <- corr_stats$sum_y2
        n_valid <- corr_stats$n_valid
        denom <- sqrt(sum_x2 * sum_y2)

        # Compute correlations (fully vectorized)
        cors <- sum_xy / denom
        zero_variance_cor <- g_mask & is.finite(sum_x2) & is.finite(sum_y2) &
            (sum_x2 <= .Machine$double.eps | sum_y2 <= .Machine$double.eps)

        dfs <- n_valid - 2L

        low_df <- (dfs < 1) & g_mask
        g_reasons[low_df] <- ifelse(g_reasons[low_df] == "", "df<1", g_reasons[low_df])
        g_mask[low_df] <- FALSE

        cors[zero_variance_cor] <- 1

        na_r <- is.na(cors) & g_mask
        g_reasons[na_r] <- ifelse(g_reasons[na_r] == "", "na r", g_reasons[na_r])
        g_mask[na_r] <- FALSE

        ps <- rep(NA_real_, sn_pairs)
        ps[zero_variance_cor] <- 0
        effective_testing_mode <- testing_mode_per_group[g]

        # Precompute parametric p-values as fallback when empirical is not feasible/resolved
        if (effective_testing_mode == "parametric") {
            tstats <- cors * sqrt(dfs / pmax(1e-12, 1 - cors * cors))
            # Compute t-statistics (vectorized)
            na_tstat <- is.na(tstats) & g_mask
            g_reasons[na_tstat] <- "na tstat"
            g_mask[na_tstat] <- FALSE
            pval_mask <- g_mask & !zero_variance_cor
            ps[pval_mask] <- -2 * expm1(pt(abs(tstats[pval_mask]), df = dfs[pval_mask], log.p = TRUE))
        } else {
            # Empirical p-values via permutations of sample labels within group
            # Only compute for rows that are still connected and have finite cors
            mask <- is.finite(cors) & g_mask & !zero_variance_cor
            if (any(mask)) {
                x_mat <- group_m[corr_start_inds, , drop = FALSE]
                y_mat <- group_m[corr_end_inds, , drop = FALSE]
                x_centered <- x_mat - x_means
                y_centered <- y_mat - y_means
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
                        ps[g_mask & !zero_variance_cor] <- 1
                        skip_empirical <- TRUE
                    }
                }
                if (!skip_empirical) {
                    .log_info(
                        "Computing empirical p-values for group '", g, "' using ",
                        if (do_permutations) "permutations" else "Monte Carlo",
                        " with ", ntries, " tries.",
                        level = 4
                    )
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
                        ps[mask] <- (counts_ge[mask] + 0.5 * counts_eq[mask] + 1) / (ntries + 1)
                    } else {
                        ps[mask] <- (counts_ge[mask] + counts_eq[mask] + 1) / (ntries + 1)
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
                min(v, na.rm = TRUE)
            }
        ))
        disconnected_pval_mask <- corr_mask & !connected
        if (any(disconnected_pval_mask)) {
            disconnected_pvals <- per_group_p[, disconnected_pval_mask, drop = FALSE]
            disconnected_reasons <- per_group_reasons[, disconnected_pval_mask, drop = FALSE]
            pvals[disconnected_pval_mask] <- vapply(seq_len(ncol(disconnected_pvals)), function(i) {
                failed_by_pval <- disconnected_reasons[, i] == "pval>max_pval" & !is.na(disconnected_pvals[, i])
                if (!any(failed_by_pval)) {
                    return(NA_real_)
                }
                min(disconnected_pvals[failed_by_pval, i], na.rm = TRUE)
            }, numeric(1))
        }
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
.safeFiniteAgg <- function(x, aggfun) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(NA_real_)
    }
    as.numeric(aggfun(x, na.rm = TRUE))
}

#' @keywords internal
#' @noRd
.safeFiniteMin <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(NA_real_)
    }
    min(x)
}

#' @keywords internal
#' @noRd
.safeFiniteMax <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(NA_real_)
    }
    max(x)
}

#' @keywords internal
#' @noRd
.signedFiniteAgg <- function(x, aggfun) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
        return(NA_real_)
    }
    .safeFiniteAgg(abs(x), aggfun) * sign(sum(sign(x)))
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
        cases_beta = .signedFiniteAgg(cases_beta, aggfun),
        controls_beta = .signedFiniteAgg(controls_beta, aggfun),
        cases_beta_sd = .safeFiniteAgg(cases_beta_sd, aggfun),
        controls_beta_sd = .safeFiniteAgg(controls_beta_sd, aggfun),
        cases_beta_min = .safeFiniteMin(cases_beta),
        cases_beta_max = .safeFiniteMax(cases_beta),
        controls_beta_min = .safeFiniteMin(controls_beta),
        controls_beta_max = .safeFiniteMax(controls_beta)
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

.stoufferPvalue <- function(pvals) {
    pvals <- suppressWarnings(as.numeric(pvals))
    pvals <- pvals[is.finite(pvals) & !is.na(pvals)]
    if (length(pvals) == 0L) {
        return(NA_real_)
    }
    if (length(pvals) == 1L) {
        return(pvals[[1L]])
    }
    pvals <- pmin(pmax(pvals, .Machine$double.xmin), 1 - .Machine$double.eps)
    stats::pnorm(sum(stats::qnorm(pvals, lower.tail = FALSE)) / sqrt(length(pvals)), lower.tail = FALSE)
}

.combineDMRSeedPvalues <- function(seed_lists, seed_pvals) {
    if (is.null(seed_pvals)) {
        return(rep(NA_real_, length(seed_lists)))
    }
    vapply(seed_lists, function(seeds) {
        .stoufferPvalue(seed_pvals[.splitCsvValues(seeds)])
    }, numeric(1))
}

#' @keywords internal
#' @noRd
.seedPvalColumn <- function(seeds_df) {
    candidates <- c("pval", "P.Value", "p.value", "p_value")
    candidates[candidates %in% colnames(seeds_df)][1L]
}

#' @keywords internal
#' @noRd
.seedPvalSource <- function(seeds_df, seeds_id_col) {
    pval_col <- .seedPvalColumn(seeds_df)
    if (is.na(pval_col)) {
        return(NULL)
    }
    list(
        ids = as.character(seeds_df[[seeds_id_col]]),
        pvals = seeds_df[[pval_col]],
        pval_col = pval_col
    )
}

#' @keywords internal
#' @noRd
.seedPvaluesForSelectedSeeds <- function(seed_pval_source, selected_seed_ids) {
    if (is.null(seed_pval_source) || length(selected_seed_ids) == 0L) {
        return(NULL)
    }
    selected_seed_ids <- unique(as.character(selected_seed_ids))
    ids <- seed_pval_source$ids
    keep <- !is.na(ids) & nzchar(ids) & ids %in% selected_seed_ids
    if (!any(keep)) {
        return(NULL)
    }
    raw_pvals <- seed_pval_source$pvals[keep]
    if (is.list(raw_pvals)) {
        raw_pvals <- vapply(raw_pvals, function(x) {
            if (length(x) == 0L) NA_character_ else as.character(x[[1L]])
        }, character(1))
    }
    pvals <- suppressWarnings(as.numeric(as.character(raw_pvals)))
    invalid <- !is.na(pvals) & (pvals < 0 | pvals > 1)
    if (any(invalid)) {
        stop("Seed p-values in column '", seed_pval_source$pval_col, "' must be between 0 and 1.")
    }
    ids <- ids[keep]
    pvals_by_id <- tapply(pvals, ids, function(x) {
        if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
    })
    stats::setNames(as.numeric(pvals_by_id), names(pvals_by_id))
}

#' @keywords internal
#' @noRd
.stage1SeedChunkRanges <- function(n_seeds, max_chunk_seeds = getOption("CMEnt.max_stage1_seeds_per_chunk", 100000)) {
    n_seeds <- as.integer(n_seeds)
    max_chunk_seeds <- suppressWarnings(as.integer(max_chunk_seeds))
    if (length(max_chunk_seeds) == 0L || is.na(max_chunk_seeds) || max_chunk_seeds < 2L) {
        max_chunk_seeds <- 100000
    }
    if (n_seeds <= 0L) {
        out <- matrix(integer(0), ncol = 2L)
        colnames(out) <- c("start", "end")
        return(out)
    }
    if (n_seeds <= max_chunk_seeds) {
        out <- matrix(c(1L, n_seeds), ncol = 2L)
        colnames(out) <- c("start", "end")
        return(out)
    }
    starts <- seq.int(1L, n_seeds - 1L, by = max_chunk_seeds - 1L)
    ends <- pmin(starts + max_chunk_seeds - 1L, n_seeds)
    cbind(start = starts, end = ends)
}

#' @keywords internal
#' @noRd
.dmrsFromSeedConnectivity <- function(seeds_connectivity_array, seed_ids, seeds_locs, min_seeds) {
    connected_seeds <- seeds_connectivity_array$connected
    breakpoints <- which(!connected_seeds)
    if (length(breakpoints) == 0L || breakpoints[[length(breakpoints)]] != length(seed_ids)) {
        breakpoints <- c(breakpoints, length(seed_ids))
    }

    segment_starts <- c(1L, head(breakpoints, -1L) + 1L)
    segment_ends <- breakpoints
    segment_lengths <- segment_ends - segment_starts + 1L
    ids <- seq_along(segment_starts)
    edge_ids <- rep(ids, segment_lengths)
    edge_ids[segment_ends] <- NA_integer_

    valid_edge <- !is.na(edge_ids)
    if (any(valid_edge)) {
        pvals <- tapply(
            seeds_connectivity_array$pval[valid_edge],
            edge_ids[valid_edge],
            function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
        )
    } else {
        pvals <- numeric(0)
    }

    dmrs <- data.frame(
        chr = as.character(seeds_locs[segment_starts, "chr"]),
        start_seed = seed_ids[segment_starts],
        end_seed = seed_ids[segment_ends],
        start_seed_pos = as.integer(seeds_locs[segment_starts, "start"]),
        end_seed_pos = as.integer(seeds_locs[segment_ends, "start"]),
        seeds_num = segment_lengths,
        stop_connection_reason = seeds_connectivity_array$reason[segment_ends],
        id = ids,
        stringsAsFactors = FALSE
    )
    dmrs <- dmrs[dmrs$seeds_num >= min_seeds, , drop = FALSE]
    if (nrow(dmrs) == 0L) {
        dmrs$connection_corr_pval <- numeric(0)
        dmrs$seeds <- character(0)
        return(dmrs)
    }

    dmrs$connection_corr_pval <- as.numeric(pvals[as.character(dmrs$id)])
    dmrs$seeds <- mapply(function(start_idx, end_idx) {
        paste(seed_ids[start_idx:end_idx], collapse = ",")
    }, segment_starts[dmrs$id], segment_ends[dmrs$id], USE.NAMES = FALSE)
    rownames(dmrs) <- NULL
    dmrs
}

.checkResult <- function(dmrs, stage, start_col = "start", end_col = "end") {
    end_less_than_start <- dmrs[[end_col]] - dmrs[[start_col]] < 0

    if (any(end_less_than_start)) {
        .log_error(
            paste0(
                "Error in stage ", stage, ": ",
                sum(end_less_than_start),
                " DMRs have been assigned an end larger than start ! (CODE BUG TO BE REPORTED)",
                " Those are: \n\t",
                paste0(capture.output(dmrs[end_less_than_start, ]), collapse = "\n\t")
            )
        )
    }
}


#' @keywords internal
#' @noRd
.buildDMRsChr <- function(
    beta_handler,
    seed_ids,
    seed_beta_index,
    pheno_detection,
    group_inds,
    testing_mode_per_group,
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
    memory_njobs,
    verbose,
    .load_debug,
    pheno,
    beta_col_names,
    sample_group_col,
    covariates,
    covariate_models,
    annotate_with_genes,
    score_dmrs,
    extract_motifs
) {
    if (!inherits(beta_handler, "BetaHandler")) {
        stop(".buildDMRsChr expects a chromosome-scoped BetaHandler.")
    }
    array_based <- beta_handler$isArrayBased()
    bsseq_gr <- if (.usesBSseqBackend(beta_handler)) .bsseqBackendGRanges(beta_handler) else NULL
    if (is.null(bsseq_gr)) {
        beta_locs <- beta_handler$getBetaLocs()
        all_sites <- .explicitRowNames(beta_locs)
        chromosome <- unique(as.character(beta_locs[, "chr"]))
    } else {
        beta_locs <- NULL
        all_sites <- NULL
        chromosome <- unique(as.character(GenomeInfoDb::seqnames(bsseq_gr)))
    }
    .log_step("Building DMRs for ", chromosome, "..", level = 1)

    chromosome_progress <- NULL
    chromosome_progress_state <- new.env(parent = emptyenv())
    chromosome_progress_state$step <- 0L
    chromosome_progress_total <- 6L +
        as.integer(isTRUE(annotate_with_genes)) +
        as.integer(isTRUE(score_dmrs)) +
        as.integer(isTRUE(extract_motifs))
    if (verbose >= 1L) {
        chromosome_progress <- utils::txtProgressBar(
            min = 0,
            max = chromosome_progress_total,
            style = 3,
            title = paste0("Processing ", chromosome),
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
        chromosome_progress_state$step <- chromosome_progress_state$step + 1L
        utils::setTxtProgressBar(
            chromosome_progress,
            min(chromosome_progress_state$step, chromosome_progress_total)
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

    .log_info("Number of provided chromosome-scoped seeds: ", length(seed_ids), level = 2)

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
        seed_chunk_ranges <- .stage1SeedChunkRanges(length(seed_ids))
        .log_step(
            "Building seed connectivity array in ", nrow(seed_chunk_ranges),
            " chunk(s) with up to ",
            getOption("CMEnt.max_stage1_seeds_per_chunk", 100000),
            " seeds per chunk...",
            level = 3
        )
        dmrs_chunks <- vector("list", nrow(seed_chunk_ranges))
        for (chunk_i in seq_len(nrow(seed_chunk_ranges))) {
            chunk_idx <- seq.int(seed_chunk_ranges[chunk_i, "start"], seed_chunk_ranges[chunk_i, "end"])
            chunk_seed_ids <- seed_ids[chunk_idx]
            .log_info(
                "Processing seed chunk ", chunk_i, "/", nrow(seed_chunk_ranges),
                " (", length(chunk_idx), " seeds).",
                level = 3
            )
            seeds_locs <- if (is.null(bsseq_gr)) {
                as.data.frame(beta_locs[seed_beta_index[chunk_idx], , drop = FALSE])
            } else {
                .bsseqLocsFromGRanges(bsseq_gr, seed_beta_index[chunk_idx])
            }
            rownames(seeds_locs) <- chunk_seed_ids
            seeds_beta <- beta_handler$getBeta(row_names = seed_beta_index[chunk_idx], col_names = beta_col_names_detection)
            rownames(seeds_beta) <- chunk_seed_ids

            if (nrow(seeds_locs) != nrow(seeds_beta)) {
                stop(
                    "Number of rows in the queried seeds beta file does not match the number of seeds. Number of rows in beta file: ",
                    nrow(seeds_beta),
                    " Number of rows in seeds: ",
                    nrow(seeds_locs)
                )
            }
            all.na.rows <- matrixStats::rowAlls(is.na(as.matrix(seeds_beta)))
            if (any(all.na.rows)) {
                stop(
                    "Beta extraction failure: the following Seed rows have all NA beta values: ",
                    paste(rownames(seeds_beta)[all.na.rows], collapse = ","),
                    ". This indicates a mismatch between requested site IDs and beta file columns or a parsing issue."
                )
            }

            seeds_beta_handler <- getBetaHandler(
                beta = seeds_beta,
                array = array,
                genome = genome,
                sorted_locs = seeds_locs,
                njobs = njobs
            )
            rm(seeds_beta)

            ret <- .buildConnectivityArray(
                beta_handler = seeds_beta_handler,
                beta_locs = seeds_locs,
                pheno = pheno_detection,
                group_inds = group_inds,
                testing_mode_per_group = testing_mode_per_group,
                empirical_strategy_per_group = empirical_strategy_per_group,
                col_names = beta_col_names_detection,
                max_pval = max_pval,
                ext_site_delta_beta = NA_real_, # delta-beta based rescue is applied later during the extension
                covariates = covariates,
                covariate_models = covariate_models,
                max_lookup_dist = max_lookup_dist,
                entanglement = entanglement,
                aggfun = aggfun,
                ntries = ntries,
                mid_p = mid_p,
                njobs = njobs,
                memory_njobs = memory_njobs,
                expansion_windows = NULL,
                max_bridge_gaps = max_bridge_seeds_gaps,
                verbose = verbose
            )
            testing_mode_per_group <- ret$testing_mode_per_group
            empirical_strategy_per_group <- ret$empirical_strategy_per_group
            dmrs_chunks[[chunk_i]] <- .dmrsFromSeedConnectivity(
                seeds_connectivity_array = ret$connectivity_array,
                seed_ids = chunk_seed_ids,
                seeds_locs = seeds_locs,
                min_seeds = min_seeds
            )
            rm(ret, seeds_beta_handler, seeds_locs)
            gc(FALSE)
        }

        dmrs_chunks <- Filter(function(x) !is.null(x) && nrow(x) > 0L, dmrs_chunks)
        dmrs <- if (length(dmrs_chunks) > 0L) {
            do.call(rbind, dmrs_chunks)
        } else {
            data.frame()
        }
        if (nrow(dmrs) > 0L) {
            rownames(dmrs) <- NULL
            dmrs$id <- seq_len(nrow(dmrs))
        }
        if (min_seeds > 0) {
            .log_info("Number of DMRs after filtering by min_seeds: ", nrow(dmrs), level = 2)
        }

        if (nrow(dmrs) == 0) {
            .log_warn("No DMRs remain after filtering based on min_seeds.")
            return(NULL)
        }
        if (getOption("CMEnt.make_debug_dir", FALSE)) {
            debug_path <- file.path("debug", paste0("01_dmrs_from_connected_seeds_", chromosome, ".tsv"))
            .log_info("Saving initial DMRs from connected seeds to ", debug_path, level = 1)
            dir.create("debug", showWarnings = FALSE)
            write.table(dmrs,
                file = debug_path,
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
    if (is.null(beta_locs)) {
        beta_locs <- beta_handler$getBetaLocs()
        all_sites <- .explicitRowNames(beta_locs)
    }
    .expandStage2DMRs <- function(dmrs_to_expand, stage2_beta_locs, connectivity_array) {
        if (nrow(connectivity_array) != nrow(stage2_beta_locs)) {
            stop(
                "Stage 2 connectivity_array row count (", nrow(connectivity_array),
                ") does not match Stage 2 beta_locs row count (", nrow(stage2_beta_locs),
                "). If using .load_debug, rebuild debug artifacts with the current code."
            )
        }
        .log_step("Expanding ", nrow(dmrs_to_expand), " DMRs using ", njobs, " jobs...", level = 3)
        locs <- as.data.frame(stage2_beta_locs)
        connectivity <- connectivity_array
        locs_rownames <- rownames(locs)
        locs_idx_map <- setNames(seq_along(locs_rownames), locs_rownames)
        expansion_boundaries <- .buildExpansionBoundaryLookup(connectivity)

        dmr_inds <- seq_len(nrow(dmrs_to_expand))
        default_dmr_chunk_size <- max(1L, ceiling(length(dmr_inds) / max(njobs * 4L, 1L)))
        dmr_chunk_size <- min(default_dmr_chunk_size, length(dmr_inds))
        dmr_chunks <- split(dmr_inds, ceiling(dmr_inds / dmr_chunk_size))
        min_dmrs_for_parallel <- suppressWarnings(as.integer(getOption("CMEnt.min_dmrs_for_parallel", 1000L)))
        if (is.na(min_dmrs_for_parallel) || min_dmrs_for_parallel < 1L) {
            min_dmrs_for_parallel <- 1L
        }

        if (njobs == 1L || length(dmr_chunks) == 1L || length(dmr_inds) < min_dmrs_for_parallel) {
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
            ret <- .safeBiocParallelApply(
                X = dmr_chunks,
                FUN = .expandDMRChunk,
                dmrs = dmrs_to_expand,
                connectivity_array = connectivity,
                min_sites = min_sites,
                locs = locs,
                locs_idx_map = locs_idx_map,
                expansion_boundaries = expansion_boundaries,
                BPPARAM = .makeBiocParallelParam(njobs, n_tasks = length(dmr_chunks))
            )
        }
        .log_success("", level = 3) 
        if (inherits(ret, "try-error")) {
            stop(ret)
        }
        as.data.frame(do.call(rbind, ret))
    }

    connectivity_array <- NULL
    connectivity_connected_total <- 0L
    if (verbose > 1 && .load_debug && file.exists(file.path("debug", "connectivity_array.rds"))) {
        .log_info("Loading debug connectivity array from file...", level = 2)
        connectivity_array <- readRDS(file.path("debug", "connectivity_array.rds"))
        extended_dmrs <- .expandStage2DMRs(dmrs, beta_locs, connectivity_array)
        connectivity_connected_total <- sum(connectivity_array$connected)
    } else {
        expansion_windows <- NULL
        if (is.finite(expansion_window) && expansion_window > 0) {
            expansion_windows <- .buildWindowsFromDMRs(
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
            stage2_site_gr <- GenomicRanges::GRanges(
                seqnames = as.character(beta_locs[, "chr"]),
                ranges = IRanges::IRanges(start = as.integer(beta_locs[, "start"]), width = 1L)
            )
            expansion_window_gr <- GenomicRanges::GRanges(
                seqnames = as.character(expansion_windows$chr),
                ranges = IRanges::IRanges(
                    start = as.integer(expansion_windows$start),
                    end = as.integer(expansion_windows$end)
                )
            )
            window_hits <- GenomicRanges::findOverlaps(expansion_window_gr, stage2_site_gr, ignore.strand = TRUE)
            window_site_counts <- tabulate(S4Vectors::queryHits(window_hits), nbins = nrow(expansion_windows))
            window_pair_counts <- pmax(1L, window_site_counts - 1L)
            target_pairs_per_chunk <- .connectivityChunkSize(
                n_samples = length(beta_col_names_detection),
                njobs = memory_njobs,
                n_pairs = sum(window_pair_counts)
            )
            window_chunk_ranges <- vector("list", nrow(expansion_windows))
            window_chunk_n <- 0L
            window_chunk_start <- 1L
            window_chunk_pairs <- 0L
            for (window_i in seq_len(nrow(expansion_windows))) {
                next_pairs <- as.integer(window_pair_counts[[window_i]])
                if (window_i > window_chunk_start && window_chunk_pairs + next_pairs > target_pairs_per_chunk) {
                    window_chunk_n <- window_chunk_n + 1L
                    window_chunk_ranges[[window_chunk_n]] <- c(window_chunk_start, window_i - 1L)
                    window_chunk_start <- window_i
                    window_chunk_pairs <- 0L
                }
                window_chunk_pairs <- window_chunk_pairs + next_pairs
            }
            window_chunk_n <- window_chunk_n + 1L
            window_chunk_ranges[[window_chunk_n]] <- c(window_chunk_start, nrow(expansion_windows))
            window_chunk_ranges <- window_chunk_ranges[seq_len(window_chunk_n)]
            .log_info(
                "Stage 2 expansion windows split into ", length(window_chunk_ranges),
                " chunk(s) before connectivity array construction.",
                level = 2
            )

            dmr_seed_gr <- GenomicRanges::GRanges(
                seqnames = as.character(dmrs$chr),
                ranges = IRanges::IRanges(
                    start = as.integer(dmrs$start_seed_pos),
                    end = as.integer(dmrs$end_seed_pos)
                )
            )
            extended_chunks <- vector("list", length(window_chunk_ranges))
            expanded_dmr_mask <- rep(FALSE, nrow(dmrs))
            for (window_chunk_i in seq_along(window_chunk_ranges)) {
                window_chunk_range <- window_chunk_ranges[[window_chunk_i]]
                window <- expansion_windows[seq.int(window_chunk_range[[1L]], window_chunk_range[[2L]]), , drop = FALSE]
                window_gr <- GenomicRanges::GRanges(
                    seqnames = as.character(window$chr),
                    ranges = IRanges::IRanges(
                        start = as.integer(window$start),
                        end = as.integer(window$end)
                    )
                )
                dmr_hits <- GenomicRanges::findOverlaps(dmr_seed_gr, window_gr, type = "within", ignore.strand = TRUE)
                dmr_chunk_inds <- sort(unique(S4Vectors::queryHits(dmr_hits)))
                if (length(dmr_chunk_inds) == 0L) {
                    next
                }
                if (any(expanded_dmr_mask[dmr_chunk_inds])) {
                    stop("Stage 2 expansion window chunking assigned a DMR to more than one chunk.")
                }
                expanded_dmr_mask[dmr_chunk_inds] <- TRUE
                subset_ret <- .subsetStage2BetaToWindows(
                    beta_handler = beta_handler,
                    beta_locs = beta_locs,
                    col_names = beta_col_names_detection,
                    expansion_windows = window,
                    njobs = njobs
                )
                if (is.null(subset_ret)) {
                    stop("Stage 2 window subsetting produced an empty beta subset. This indicates inconsistent expansion windows.")
                }
                .log_info(
                    "Stage 2 beta subset contains ",
                    format(nrow(subset_ret$beta_locs), big.mark = ","),
                    " sites for expansion window chunk ", window_chunk_i, "/", length(window_chunk_ranges),
                    " (", nrow(window), " window(s)).",
                    level = 3
                )
                .log_step(
                    "Building expansion connectivity array for window chunk ",
                    window_chunk_i, "/", length(window_chunk_ranges), "..",
                    level = 3
                )
                ret <- .buildConnectivityArray(
                    beta_handler = subset_ret$beta_handler,
                    beta_locs = subset_ret$beta_locs,
                    pheno = pheno_detection,
                    group_inds = group_inds,
                    testing_mode_per_group = testing_mode_per_group,
                    empirical_strategy_per_group = empirical_strategy_per_group,
                    col_names = beta_col_names_detection,
                    max_pval = max_pval,
                    ext_site_delta_beta = ext_site_delta_beta,
                    covariates = covariates,
                    covariate_models = covariate_models,
                    max_lookup_dist = max_lookup_dist,
                    entanglement = entanglement,
                    aggfun = aggfun,
                    ntries = ntries,
                    mid_p = mid_p,
                    njobs = njobs,
                    memory_njobs = memory_njobs,
                    expansion_windows = window,
                    max_bridge_gaps = max_bridge_extension_gaps,
                    verbose = verbose
                )
                testing_mode_per_group <- ret$testing_mode_per_group
                empirical_strategy_per_group <- ret$empirical_strategy_per_group
                connectivity_connected_total <- connectivity_connected_total + sum(ret$connectivity_array$connected)
                if (getOption("CMEnt.make_debug_dir", FALSE)) {
                    debug_path <- file.path("debug", paste0("connectivity_array_", chromosome, "_window_chunk", window_chunk_i, ".rds"))
                    .log_info("Saving connectivity array to ", debug_path, level = 1)
                    dir.create("debug", showWarnings = FALSE)
                    saveRDS(ret$connectivity_array, file = debug_path)
                }
                extended_chunks[[window_chunk_i]] <- .expandStage2DMRs(
                    dmrs_to_expand = dmrs[dmr_chunk_inds, , drop = FALSE],
                    stage2_beta_locs = subset_ret$beta_locs,
                    connectivity_array = ret$connectivity_array
                )
                rm(ret, subset_ret)
                gc(FALSE)
                .log_success("Stage 2 expansion window chunk ", window_chunk_i, "/", length(window_chunk_ranges), " complete.", level = 3)
            }
            if (!all(expanded_dmr_mask)) {
                stop("Stage 2 expansion window assignment missed one or more DMRs.")
            }
            extended_chunks <- Filter(function(x) !is.null(x) && nrow(x) > 0L, extended_chunks)
            extended_dmrs <- as.data.frame(do.call(rbind, extended_chunks))
        } else {
            .log_step("Building expansion connectivity array..", level = 3)
            ret <- .buildConnectivityArray(
                beta_handler = beta_handler,
                beta_locs = beta_locs,
                pheno = pheno_detection,
                group_inds = group_inds,
                testing_mode_per_group = testing_mode_per_group,
                empirical_strategy_per_group = empirical_strategy_per_group,
                col_names = beta_col_names_detection,
                max_pval = max_pval,
                ext_site_delta_beta = ext_site_delta_beta,
                covariates = covariates,
                covariate_models = covariate_models,
                max_lookup_dist = max_lookup_dist,
                entanglement = entanglement,
                aggfun = aggfun,
                ntries = ntries,
                mid_p = mid_p,
                njobs = njobs,
                memory_njobs = memory_njobs,
                expansion_windows = NULL,
                max_bridge_gaps = max_bridge_extension_gaps,
                verbose = verbose
            )
            testing_mode_per_group <- ret$testing_mode_per_group
            empirical_strategy_per_group <- ret$empirical_strategy_per_group
            connectivity_array <- ret$connectivity_array
            connectivity_connected_total <- sum(connectivity_array$connected)
            extended_dmrs <- .expandStage2DMRs(dmrs, beta_locs, connectivity_array)
            .log_success("Connectivity array built.", level = 3)
        }
    }
    .log_info("Number of underlying correlated sites found: ", connectivity_connected_total, level = 3)
    if (!is.null(connectivity_array) && getOption("CMEnt.make_debug_dir", FALSE)) {
        debug_path <- file.path("debug", paste0("connectivity_array_", chromosome, ".rds"))
        .log_info("Saving connectivity array to ", debug_path, level = 1)
        dir.create("debug", showWarnings = FALSE)
        saveRDS(connectivity_array, file = debug_path)
    }
    rownames(extended_dmrs) <- NULL
    u_exp_len_table <- table(
        extended_dmrs$upstream_expansion_length
    )
    d_exp_len_table <- table(
        extended_dmrs$downstream_expansion_length
    )
    # sort tables by names (expansion sizes)
    u_exp_len_table <- u_exp_len_table[order(as.integer(names(u_exp_len_table)))]
    d_exp_len_table <- d_exp_len_table[order(as.integer(names(d_exp_len_table)))]
    .log_info(
        "Table of upstream_expansion_length:\n\t",
        paste(
            capture.output(u_exp_len_table),
            collapse = "\n\t"
        ),
        level = 3
    )
    .log_info(
        "Table of downstream_expansion_length:\n\t",
        paste(
            capture.output(d_exp_len_table),
            collapse = "\n\t"
        ),
        level = 3
    )

    .log_success("DMR expansion complete.", level = 2)
    .advanceChromosomeProgress()
    .log_step("Post-processing extended DMRs..", level = 2)

    extended_dmrs$end <- as.numeric(extended_dmrs$end)
    extended_dmrs$start <- as.numeric(extended_dmrs$start)
    extended_dmrs$start_seed_pos <- as.numeric(extended_dmrs$start_seed_pos)
    extended_dmrs$end_seed_pos <- as.numeric(extended_dmrs$end_seed_pos)
    extended_dmrs$seeds_num <- as.numeric(extended_dmrs$seeds_num)
    extended_dmrs$connection_corr_pval <- as.numeric(extended_dmrs$connection_corr_pval)


    .checkResult(extended_dmrs, "2", start_col = "start", end_col = "end")
    .log_success("Post-processing complete.", level = 2)
    if (getOption("CMEnt.make_debug_dir", FALSE)) {
        debug_path <- file.path("debug", paste0("02_extended_dmrs_", chromosome, ".tsv"))
        .log_info("Saving extended DMRs prior to filtering to ", debug_path, level = 1)
        dir.create("debug", showWarnings = FALSE)
        write.table(extended_dmrs,
            file = debug_path,
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
    .log_info("Frequency of N-DMRs overlap:\n", paste(capture.output(table(tqh)), collapse = "\n"), level = 3)
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
    agg_df[, "pval"] <- NA_real_
    if (is.null(all_sites)) {
        agg_df[, "sites_num"] <- vapply(agg_df$sites, function(x) {
            length(.splitCsvValues(x))
        }, integer(1))
    } else {
        agg_df[, "sites_num"] <- match(agg_df$end_site, all_sites) - match(agg_df$start_site, all_sites) + 1
    }
    start_site_pos <- sub("^[^:]+:", "", agg_df$start_site)
    end_site_pos <- sub("^[^:]+:", "", agg_df$end_site)
    agg_df[, "id"] <- paste0(
        as.character(GenomeInfoDb::seqnames(merged_dmrs_ranges)), ":",
        start_site_pos,
        "-",
        end_site_pos
    )

    GenomicRanges::mcols(merged_dmrs_ranges) <- agg_df
    .log_success("Overlapping extended DMRs merged: ", length(merged_dmrs_ranges), " resulting DMRs.", level = 2)
    merged_dmrs <- .convertToDataFrame(merged_dmrs_ranges)

    if (getOption("CMEnt.make_debug_dir", FALSE)) {
        debug_path <- file.path("debug", paste0("03_merged_dmrs_", chromosome, ".tsv"))
        .log_info("Saving merged DMRs prior to filtering to ", debug_path, level = 1)
        dir.create("debug", showWarnings = FALSE)
        write.table(merged_dmrs,
            file = debug_path,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
    }

    .advanceChromosomeProgress()
    .log_step("Stage 4: Filtering resulting DMRs..", level = 2)

    if (min_sites > 1) {
        filtered_dmrs_ranges <- merged_dmrs_ranges[
            GenomicRanges::mcols(merged_dmrs_ranges)$sites_num >= min_sites
        ]
        .log_info(
            "Keeping ",
            length(filtered_dmrs_ranges),
            " out of ",
            length(merged_dmrs_ranges),
            " with at least ",
            min_sites,
            " sites in the DMR interval.",
            level = 2
        )
    } else {
        filtered_dmrs_ranges <- merged_dmrs_ranges
    }
    filtered_dmrs <- .convertToDataFrame(filtered_dmrs_ranges)

    if (nrow(filtered_dmrs) == 0) {
        .log_warn("No DMRs passed the filtering step.")
        return(NULL)
    }

    if (array_based && !is.null(min_adj_seeds) && min_adj_seeds > min_seeds) {
        .log_step("Calculating site content and adjusted seeds number..", level = 3)
        sites_num_bg <- .getSiteBackgroundCounts(filtered_dmrs_ranges, genome)
        sites_num_bg[!is.finite(sites_num_bg) | is.na(sites_num_bg) | sites_num_bg <= 0] <- 1L
        filtered_dmrs$sites_num_bg <- sites_num_bg

        filtered_dmrs$seeds_num_adj <- ceiling(filtered_dmrs$sites_num_bg / filtered_dmrs$sites_num * filtered_dmrs$seeds_num)

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
        pheno = pheno,
        aggfun = aggfun
    )
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
        stop("Internal failure while aggregating seed beta statistics: missing dmr_id rows.")
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
        stop("Internal failure while aggregating site beta statistics: missing dmr_id rows.")
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
    .log_success("DMR filtering complete.", level = 2)
    if (annotate_with_genes) {
        .log_step("Annotating DMRs with gene information...", level = 2)
        .advanceChromosomeProgress()
        site_delta_beta <- beta_stats$cases_beta - beta_stats$controls_beta
        names(site_delta_beta) <- rownames(beta_stats)
        annotated_dmrs <- annotateDMRsWithGenes(
            annotated_dmrs,
            genome = genome,
            njobs = njobs,
            site_locs = beta_handler$getBetaLocs()[all_selected_sites, , drop = FALSE],
            site_delta_beta = site_delta_beta,
            aggfun = aggfun
        )
        .log_success("DMR annotation complete.", level = 2)
    }

    if (score_dmrs) {
        .advanceChromosomeProgress()
        .log_step("Scoring DMRs...", level = 2)
        annotated_dmrs <- scoreDMRs(
            annotated_dmrs,
            beta = beta_handler,
            pheno = pheno,
            genome = genome,
            array = array,
            sorted_locs = beta_handler$getGenomicLocs(),
            sample_group_col = sample_group_col,
            covariates = covariates,
            njobs = njobs,
            show_progress = verbose >=2,
            .dmr_beta = all_selected_sites_beta,
            .memory_njobs = memory_njobs
        )
        .log_success("DMR scoring complete.", level = 2)
    }
    rm(all_selected_sites_beta)


    if (is.data.frame(annotated_dmrs)) {
        annotated_dmrs <- .convertToGRanges(annotated_dmrs, genome = genome)
    }

    if (extract_motifs) {
        .advanceChromosomeProgress()
        .log_step("Extracting DMR motifs...", level = 2)
        annotated_dmrs <- extractDMRMotifs(annotated_dmrs, genome = genome, array = array, beta_locs = beta_locs)
        .log_success("DMR motifs computed.", level = 2)
    }

    final_dmrs_granges <- annotated_dmrs

    final_dmrs <- .convertToDataFrame(final_dmrs_granges)

    .log_info("Final number of chromosome-scoped DMRs: ", nrow(final_dmrs), level = 2)
    .log_success("DMR building complete for ", chromosome, ".", level = 1)
    invisible(final_dmrs_granges)
}


#' @keywords internal
#' @noRd
.buildDMRsChromosomeTask <- function(
    task,
    beta_handler,
    pheno_detection,
    group_inds,
    testing_mode_per_group,
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
    memory_njobs,
    verbose,
    .load_debug,
    pheno,
    beta_col_names,
    sample_group_col,
    covariates,
    covariate_models,
    annotate_with_genes,
    score_dmrs,
    extract_motifs
) {
    .log_step("Processing chromosome ", task$chr, "...", level = 1)
    task_row_ids <- if (!is.null(task$row_ids)) {
        task$row_ids
    } else {
        seq.int(task$row_start, task$row_end)
    }
    chr_handler <- beta_handler$subset(row_names = task_row_ids, col_names = beta_col_names)
    on.exit(
        {
            rm(chr_handler, task_row_ids)
            gc(FALSE)
        },
        add = TRUE
    )

    withCallingHandlers(
        .buildDMRsChr(
            beta_handler = chr_handler,
            seed_ids = task$seed_ids,
            seed_beta_index = task$seed_beta_index,
            pheno_detection = pheno_detection,
            group_inds = group_inds,
            testing_mode_per_group = testing_mode_per_group,
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
            memory_njobs = memory_njobs,
            verbose = verbose,
            .load_debug = .load_debug,
            pheno = pheno,
            beta_col_names = beta_col_names,
            sample_group_col = sample_group_col,
            covariates = covariates,
            covariate_models = covariate_models,
            annotate_with_genes = annotate_with_genes,
            score_dmrs = score_dmrs,
            extract_motifs = extract_motifs
        ),
        warning = function(w) {
            if (grepl("No DMRs|No connectivity windows", conditionMessage(w))) {
                invokeRestart("muffleWarning")
            }
        }
    )
}

#' Build Differentially Methylated Regions (DMRs) from Differentially Methylated Positions (seeds)
#'
#' This function assembles DMRs from a given set of seeds and a beta value file. It operates in three main stages:
#' 1. **Seed Connectivity**: It builds a connectivity array based on the correlation of beta values between seeds and their proximal sites,
#'  connecting seeds into preliminary DMRs based on significant correlations.
#' 2. **DMR Expansion**: It expands the preliminary DMRs by including nearby sites that show significant correlation with the seeds,
#'  allowing for a specified delta-beta threshold to connect sites that may not meet the correlation p-value cutoff but have a strong effect size.
#' 3. **DMR Merging and Filtering**: It merges overlapping extended DMRs into final DMRs and applies filtering based on the number of seeds and sites,
#'  as well as optional adjustments for array-based analyses.
#'
#' @section Note on Input Data:
#' Do not apply heavy filtering to your seeds prior to using this function, particularly based on
#' beta values or effect sizes. The function works by expanding regions around seeds
#' and connecting nearby sites into larger regions. Filtering out seeds with smaller effect sizes
#' may remove important sites that could serve as "bridges" to connect more seeds into
#' larger, biologically meaningful DMRs. For optimal results, include all statistically
#' significant seeds (e.g., adjusted p-value < 0.05) and let the function handle region expansion
#' and letting the function reconnect proximal sites during expansion using the
#' ext_site_delta_beta parameter if needed. The p-value adjustment can be done using [combinePvalues()] or other methods,
#' but avoid filtering based on beta value thresholds or effect size cutoffs before running this function. For BSSeq data,
#' [findDMPsBSSeq()] performs seed finding using DSS.
#'
#' @param beta Character. Path to the beta value file, or a tabix file, or a beta matrix, or a BetaHandler object, or a bed file.
#'  If a bed file is provided, it must at least contain bed_chrom_col and bed_chrom_start,
#'  followed by samples names in the given pheno, with corresponging beta values,
#'  and it will be converted to a tabix-indexed beta file internall,
#'  with the locations separately saved and queried as a DelayedDataFrame. object.
#' @param seeds Character. Path to the seeds (seeds, etc.) TSV file or the seeds dataframe,
#'  in a format like the one produced by dmpFinder. If a `pval`, `P.Value`, `p.value`,
#'  or `p_value` column is present, DMR-level `pval` is computed from final
#'  supporting seed p-values using Stouffer's method and `qval` is FDR-corrected globally.
#' @param pheno Character. Path to the phenotype TSV file or the phenotype dataframe,
#'  containing sample information including group labels and optionally covariates.
#' @param seeds_id_col Character. Column name or index for Seed identifiers in the seeds TSV file.
#'  Default is NULL, which corresponds to the rows names if existing, or the first column if not.
#' @param sample_group_col Character. Required column name for sample group information in the phenotype data.
#' @param casecontrol_col Boolean Column in pheno for case (TRUE/1) / control (FALSE/0) status .
#'  If NULL, controls will be assumed to be the first level of sample_group_col. Default is NULL.
#' @param covariates Character vector of column names in pheno to adjust for (e.g. "age", "sex").
#'  When provided, correlations are computed on residuals after regressing M-values on these covariates within each group
#' @param ext_site_delta_beta Numeric. Minimum absolute delta beta value that will
#'  force proximal sites to be treated as connected during Stage 2 expansion,
#'  regardless of their correlation p-value. Set to `NA`, `NULL`, or `Inf` to
#'  disable this shortcut. A value of `0` means any proximal site with a
#'  non-missing case-control delta beta can be force-connected. Default is 0.2.
#' @param array Character. Type of array used (e.g., "450K", "EPIC",
#'   "EPICv2", "27K", "Mouse"). Required when beta row names are array probe
#'   IDs.
#'  Ignored if the beta file is provided as a beta values BED file.
#' @param genome Optional character genome version. If `NULL`, inferred from
#'   `beta` when possible, otherwise from `array`; BED inputs default to
#'   `"hg38"`.
#' @param max_pval Numeric. Maximum p-value to assume seeds correlation is significant. Default is 0.05.
#' @param entanglement Character. "weak" (default) requires at least one group to show significant correlation;
#'  "strong" requires all groups to show significant correlation for connectivity.
#' @param testing_mode Character. "auto" (default) selects between t-based correlation p-values and empirical p-values
#'  per sample group using data diagnostics. You can also force "parametric" for t-based correlation p-values
#'  or "empirical" for permutation-based p-values.
#' @param empirical_strategy Character. When testing_mode = "empirical": "auto" (default) uses Monte Carlo for
#'  groups with <6 samples and permutations for groups with >=6 samples; "montecarlo" always uses Monte Carlo;
#'  "permutations" always uses permutations.
#' @param ntries Integer. Number of permutations when testing_mode = "empirical". Default is 0 (disabled).
#' @param mid_p Logical. Whether to use mid-p values for empirical correlation tests. Default is FALSE.
#' @param max_lookup_dist Numeric. Maximum distance to look up for adjacent seeds belonging to the same DMR during Stage 1. Default is 10000 (10 kb).
#' @param expansion_window Numeric. Stage 2 connectivity is computed only in windows centered on seed-derived Stage 1 DMR neighborhoods,
#'  with this total window width in bp. This value sets a maximum effective size of a DMR after stage 2.
#'   Set <=0 for genome-wide connectivity. Default is -1 for microarrays and 10000 (10 kb) for NGS datasets.
#' @param max_bridge_seeds_gaps Integer. Maximum number of consecutive failed seed-to-seed edges to bridge during Stage 1
#'  when both flanking edges are connected and failures are p-value driven. Set to 0 to disable. Default is 1.
#' @param max_bridge_extension_gaps Integer. Maximum gap size to consider during Stage 2 extension.
#'  Default is 1 (i.e., at most 1 consecutive failing site to bridge).
#' @param min_seeds Numeric. Minimum number of connected seeds in a DMR. Minimum is 1. Default is 2.
#' @param min_adj_seeds Numeric. Minimum number of seeds, adjusted by array site density, in a DMR after extension.
#'  It serves as a less stringent cutoff for arrays with variable site density, allowing regions in sparse areas to be retained
#'  if they have enough seeds relative to the local site density. Default is NULL (disabled).
#' @param min_sites Numeric. Minimum number of sites in a DMR after extension, including the seeds. Minimum is 2. Default is 3.
#' @param aggfun Function or character. Aggregation function to use when calculating delta beta values and p-values of DMRs.
#'  Can be "median", "mean", or a function (e.g., median, mean). Default is "median".
#' @param ignored_sample_groups Character vector. Sample groups to ignore during connection and expansion, separated by commas.
#'  Can also be "case" or "control". Default is NULL.
#' @param output_prefix Character. Identifier for the output files. If not provided, no output will be saved. Default is NULL.
#' @param njobs Numeric. Number of parallel jobs to use. Default is the number of available cores.
#' @param beta_row_names_file Character. Path to a file containing row names for the beta values.
#'  If not provided, row names will be read from the beta file. Default is NULL.
#' @param annotate_with_genes Logical. Whether to annotate DMRs with overlapping genes. Default is TRUE.
#' @param score_dmrs Logical. Whether to add complementary cross-validated
#'  SVM discrimination scores to DMRs. Default is TRUE.
#' @param extract_motifs Logical. Whether to compute DMRs seeds motifs. Default is TRUE.
#' @param bed_provided Logical. Whether the beta file is provided as a BED file. Default is FALSE.
#'  In case the input has a .bed extension, this will be set to TRUE automatically.
#' @param bed_chrom_col Character. Column name for chromosome in the BED file. Default is "chrom".
#' @param bed_start_col Character. Column name for start position in the BED file. Default is "start".
#' @param verbose Numeric. Level of verbosity for logging messages, from 0 (not verbose) to 5 (very very verbose).
#'  Default is retrieved from option "CMEnt.verbose".
#' @param .load_debug Logical. If TRUE, enables debug mode for loading beta files. Default is FALSE.
#'
#' @return GRanges object of identified DMRs with metadata including DMR-level
#'  `pval` and FDR-adjusted `qval` when seed p-value aggregation is enabled.
#'
#' @examples
#' loadExampleInputDataChr21And22("beta", "dmps", "pheno", "array_type")
#' \donttest{
#' dmrs <- buildDMRs(
#'     beta = beta,
#'     seeds = dmps,
#'     pheno = pheno,
#'     array = array_type,
#'     sample_group_col = "Sample_Group"
#' )
#' }
#' @export
buildDMRs <- function(
    beta,
    seeds,
    pheno,
    seeds_id_col = NULL,
    sample_group_col = NULL,
    casecontrol_col = NULL,
    covariates = NULL,
    ext_site_delta_beta = 0.2,
    array = NULL,
    genome = NULL,
    max_pval = 0.05,
    entanglement = c("weak", "strong"),
    testing_mode = c("auto", "parametric", "empirical"),
    empirical_strategy = c("auto", "montecarlo", "permutations"),
    ntries = 200L,
    mid_p = FALSE,
    max_lookup_dist = 10000,
    expansion_window = "auto",
    max_bridge_seeds_gaps = 1L,
    max_bridge_extension_gaps = 1L,
    min_seeds = 2,
    min_adj_seeds = NULL,
    min_sites = 3,
    aggfun = c("median", "mean"),
    ignored_sample_groups = NULL,
    output_prefix = NULL,
    njobs = getOption("CMEnt.njobs", .defaultNJobs()),
    beta_row_names_file = NULL,
    annotate_with_genes = TRUE,
    score_dmrs = TRUE,
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
            seeds_tsv <- try(as.data.frame(data.table::fread(
                seeds,
                header = TRUE, check.names = FALSE,
                quote = "", comment.char = "", data.table = FALSE
            )))
        } else if (is.data.frame(seeds)) {
            seeds_tsv <- as.data.frame(seeds)
        } else if (is.vector(seeds) && length(seeds) > 1) {
            seeds_tsv <- as.data.frame(seeds)
            colnames(seeds_tsv) <- "seeds"
            seeds_id_col <- "seeds"
        } else {
            stop("seeds must be either a file path, a vector, or a data frame")
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
    .readPheno <- function(pheno) {
        if (is.character(pheno) && length(pheno) == 1) {
            pheno_df <- try(as.data.frame(data.table::fread(
                pheno,
                header = TRUE, check.names = FALSE,
                quote = "", comment.char = ""
            )))
        } else if (is.data.frame(pheno)) {
            pheno_df <- as.data.frame(pheno)
        } else {
            stop("pheno must be either a file path or a data frame")
        }
        if (inherits(pheno_df, "try-error")) {
            stop("Could not read pheno file: ", pheno)
        }
        if (all(rownames(pheno_df) == as.character(seq_len(nrow(pheno_df))))) {
            pheno_row_names <- as.character(pheno_df[, 1])
            if (!anyNA(pheno_row_names) && all(nzchar(pheno_row_names)) && !anyDuplicated(pheno_row_names)) {
                rownames(pheno_df) <- pheno_row_names
            }
        }
        if (!is.null(covariates)) {
            missing_covars <- covariates[!covariates %in% colnames(pheno_df)]
            if (length(missing_covars) > 0) {
                stop(
                    "The following covariates are not present in pheno: ",
                    paste(missing_covars, collapse = ", ")
                )
            }
        }
        if (!sample_group_col %in% colnames(pheno_df)) {
            stop("sample_group_col '", sample_group_col, "' is not present in pheno.")
        }
        if (!is.null(casecontrol_col) && !casecontrol_col %in% colnames(pheno_df)) {
            stop("casecontrol_col '", casecontrol_col, "' is not present in pheno.")
        }
        sample_group_values <- .coercePhenoColumn(pheno_df[[sample_group_col]], sample_group_col)
        pheno_df[[sample_group_col]] <- sample_group_values
        if (!is.null(casecontrol_col)) {
            pheno_df[[casecontrol_col]] <- .coercePhenoColumn(pheno_df[[casecontrol_col]], casecontrol_col)
        }
        pheno_df
    }
    .addCaseControlColumn <- function(pheno_df) {
        sample_group_values <- pheno_df[[sample_group_col]]
        if (is.null(casecontrol_col)) {
            control_group <- levels(as.factor(sample_group_values))[1]
            pheno_df[, .CASE_CONTROL_COL] <- ifelse(sample_group_values == control_group, 0, 1)
        } else {
            pheno_df[, .CASE_CONTROL_COL] <- as.numeric(pheno_df[[casecontrol_col]])
        }
        pheno_df
    }
    testing_mode <- strex::match_arg(testing_mode, ignore_case = TRUE)
    empirical_strategy <- strex::match_arg(empirical_strategy, ignore_case = TRUE)
    entanglement <- strex::match_arg(entanglement, ignore_case = TRUE)
    options(CMEnt.verbose = verbose, "CMEnt.njobs" = njobs)
    if (Sys.info()[["sysname"]] != "Windows") {
        includes <- "#include <sys/wait.h>"
        code <- "int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};"
        wait <- inline::cfunction(body = code, includes = includes, convention = ".C")
        withr::defer(wait())
    }

    .log_step("Preparing DMR input...")
    .log_step("Reading seeds and phenotype files...", level = 3)
    seeds_ret <- .readSeeds(seeds, seeds_id_col)
    .log_step("Validating sample group column in phenotype data...", level = 3)
    sample_group_col <- .requireSampleGroupCol(sample_group_col, "buildDMRs()")
    .log_step("Reading phenotype file...", level = 3)
    pheno <- .readPheno(pheno)
    seeds_df <- seeds_ret$data
    seeds_id_col <- seeds_ret$id_col
    if (is.null(seeds_df) || nrow(seeds_df) == 0L) {
        .log_warn("Provided seeds file has no data rows. Not proceeding.")
        .emptyOutputs(output_prefix)
        return(NULL)
    }

    array <- .normalizeBuildDMRsArray(array)
    requested_genome <- .normalizeBuildDMRsGenome(genome)
    genome <- .resolveBuildDMRsGenome(beta, array = array, genome = requested_genome, bed_provided = bed_provided)
    .assertDependencyRequirements(
        requirements = .buildDMRsDependencyReqs(
            beta = beta,
            array = array,
            genome = genome,
            annotate_with_genes = annotate_with_genes,
            extract_motifs = extract_motifs,
            bed_provided = bed_provided
        ),
        context = "buildDMRs()"
    )
    if (is.null(requested_genome)) {
        .log_info("No genome provided. Using inferred genome: ", genome, ".", level = 2)
    }
    output_prefix_base <- output_prefix
    output_prefix_dot <- NULL
    if (!is.null(output_prefix_base)) {
        dir.create(dirname(output_prefix_base), showWarnings = FALSE, recursive = TRUE)
        output_prefix_dot <- paste0(output_prefix_base, ".")
    }
    .log_step("Preparing beta data handler...", level = 3)
    beta_locs <- NULL
    if (!inherits(beta, "BetaHandler") && is.character(beta) && length(beta) == 1 && file.exists(beta)) {
        beta_file_ext <- tools::file_ext(beta)
        if (beta_file_ext == "bed" || bed_provided) {
            bed_provided <- TRUE
            seed_ids <- seeds_df[, seeds_id_col]
            if (!all(grepl("^(chr)?[0-9XYM]+:[0-9]+$", seed_ids))) {
                stop(
                    "When providing a bed file as beta input,",
                    " seed IDs must be in 'chr:pos' format."
                )
            }
            ret <- readCustomMethylationBedData(
                bed_file = beta, pheno = pheno, genome = genome,
                chrom_col = bed_chrom_col, start_col = bed_start_col,
                output_prefix = output_prefix_base,
                njobs = njobs
            )
            beta <- ret$beta_file
            beta_locs <- ret$locations
        }
    }

    beta_handler <- getBetaHandler(
        beta = beta,
        beta_row_names_file = beta_row_names_file,
        njobs = njobs,
        sorted_locs = beta_locs,
        array = array,
        genome = genome,
        output_prefix = output_prefix_base
    )
    beta <- NULL
    array_based <- beta_handler$isArrayBased()
    .log_step("Validating parameters...", level = 3)

    if (!is.function(aggfun)) {
        aggfun_choice <- strex::match_arg(aggfun, ignore_case = TRUE)
        aggfun <- switch(aggfun_choice,
            median = stats::median,
            mean = mean
        )
    }
    stopifnot(!is.null(max_pval))
    stopifnot(!is.null(min_seeds))
    stopifnot(!is.null(min_sites))
    if (min_seeds < 2 && min_sites < 2) {
        stop("min_seeds or min_sites must be at least 2, to define a DMR")
    }
    stopifnot(!is.null(max_lookup_dist))
    ext_site_delta_beta <- {
        if (is.null(ext_site_delta_beta) || is.na(ext_site_delta_beta)) {
            NA_real_
        } else {
            try(ext_site_delta_beta <- as.numeric(ext_site_delta_beta), silent = TRUE)
            if (!is.numeric(ext_site_delta_beta) || length(ext_site_delta_beta) != 1L) {
                stop(
                    "ext_site_delta_beta must be NULL, NA, Inf,",
                    " or a numeric scalar in [0, 1]."
                )
            }
            ext_site_delta_beta <- ext_site_delta_beta[1]
            if (is.infinite(ext_site_delta_beta)) {
                NA_real_
            } else {
                if (ext_site_delta_beta < 0 || ext_site_delta_beta > 1) {
                    stop(
                        "ext_site_delta_beta must be NULL, NA, Inf,",
                        " or a numeric scalar in [0, 1]."
                    )
                }
                ext_site_delta_beta
            }
        }
    }

    if (expansion_window == "auto") {
        expansion_window <- if (array_based) -1 else 10000
        if (array_based) {
            .log_info(
                "Setting expansion_window to -1 for array-based dataset,",
                " meaning genome-wide connectivity will be computed during",
                " expansion.",
                level = 2
            )
        } else {
            .log_info(
                "Setting expansion_window to 10000 for sequencing-based dataset,",
                " meaning connectivity during expansion will only be computed",
                " within 10 kb windows around seed-derived DMRs.",
                level = 2
            )
        }
    }
    if (!is.numeric(expansion_window) || length(expansion_window) != 1 || is.na(expansion_window)) {
        stop("expansion_window must be a numeric scalar.")
    }
    if (!is.numeric(max_bridge_seeds_gaps) || length(max_bridge_seeds_gaps) != 1 || is.na(max_bridge_seeds_gaps) || max_bridge_seeds_gaps < 0) {
        stop("max_bridge_seeds_gaps must be a non-negative integer.")
    }
    max_bridge_seeds_gaps <- as.integer(max_bridge_seeds_gaps)
    if (
        !is.numeric(max_bridge_extension_gaps) ||
            length(max_bridge_extension_gaps) != 1 ||
            is.na(max_bridge_extension_gaps) ||
            max_bridge_extension_gaps < 0
    ) {
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
    pheno <- pheno[beta_col_names, , drop = FALSE]
    pheno <- .addCaseControlColumn(pheno)

    ignored_sample_groups_chr <- {
        if (is.null(ignored_sample_groups)) {
            character(0)
        } else {
            x <- trimws(unlist(base::strsplit(ignored_sample_groups, ",")))
            x[nzchar(x)]
        }
    }
    samples_selection_mask <- !(pheno[[sample_group_col]] %in% ignored_sample_groups_chr)
    if ("case" %in% ignored_sample_groups_chr) {
        samples_selection_mask <- samples_selection_mask & (pheno[, .CASE_CONTROL_COL] != 1)
    }
    if ("control" %in% ignored_sample_groups_chr) {
        samples_selection_mask <- samples_selection_mask & (pheno[, .CASE_CONTROL_COL] != 0)
    }
    beta_col_names_detection <- beta_col_names[samples_selection_mask]
    if (length(beta_col_names_detection) < 2) {
        stop("At least two samples are required after applying ignored_sample_groups.")
    }
    pheno_detection <- pheno[beta_col_names_detection, , drop = FALSE]
    sample_groups <- factor(pheno_detection[[sample_group_col]])
    group_inds <- split(seq_along(sample_groups), sample_groups)
    covariate_models <- .prepareGroupCovariateModels(pheno_detection, group_inds, covariates)
    testing_mode_per_group <- rep(testing_mode, length.out = length(unique(pheno_detection[[sample_group_col]])))
    names(testing_mode_per_group) <- unique(pheno_detection[[sample_group_col]])
    empirical_strategy_per_group <- rep(
        empirical_strategy,
        length.out = length(unique(pheno_detection[[sample_group_col]]))
    )
    names(empirical_strategy_per_group) <- unique(pheno_detection[[sample_group_col]])

    if (!is.null(output_prefix_base)) {
        saveRDS(
            list(
                pheno = pheno_detection,
                genome = genome,
                array = array,
                sample_group_col = sample_group_col
            ),
            file = paste0(output_prefix_base, ".meta.rds")
        )
    }

    bsseq_gr <- if (.usesBSseqBackend(beta_handler)) .bsseqBackendGRanges(beta_handler) else NULL
    if (!is.null(bsseq_gr)) {
        beta_locs <- NULL
        beta_chr <- GenomeInfoDb::seqnames(bsseq_gr)
        beta_start <- GenomicRanges::start(bsseq_gr)
        beta_locs_rownames <- NULL
    } else {
        beta_locs <- beta_handler$getBetaLocs()
        beta_chr <- as.character(beta_locs[, "chr"])
        beta_start <- suppressWarnings(as.numeric(beta_locs[, "start"]))
        if (anyNA(beta_chr) || any(!nzchar(beta_chr))) {
            stop("Beta locations contain missing chromosome labels.", call. = FALSE)
        }
        if (anyNA(beta_start)) {
            stop("Beta locations contain missing or non-numeric start positions.", call. = FALSE)
        }
        beta_chr_runs <- .chromosomeRunIndex(beta_chr, require_unique = TRUE)
        if (length(beta_chr) > 1L) {
            same_chr_adj <- beta_chr[-1L] == beta_chr[-length(beta_chr)]
            unsorted_adj <- same_chr_adj & (beta_start[-1L] < beta_start[-length(beta_start)])
            if (any(unsorted_adj, na.rm = TRUE)) {
                bad_idx <- which(unsorted_adj)[1L] + 1L
                stop("Beta locations are not sorted within chromosome ", beta_chr[bad_idx], call. = FALSE)
            }
        }
        beta_locs_rownames <- .explicitRowNames(beta_locs)
    }
    if (!exists("beta_chr_runs", inherits = FALSE)) {
        beta_chr_runs <- .chromosomeRunIndex(beta_chr, require_unique = TRUE)
    }
    use_numeric_sequencing_rows <- !array_based &&
        (is.null(beta_locs_rownames) || .usesBSseqBackend(beta_handler))
    beta_row_names <- if (use_numeric_sequencing_rows) NULL else beta_handler$getBetaRowNames()
    if (use_numeric_sequencing_rows) {
        .log_step("Matching seed IDs to beta genomic locations...", level = 3)
        seed_coords <- .seedCoordinatesFromSeeds(seeds_df, seeds_id_col)
        seed_beta_index <- .matchSequencingCoordinatesToBeta(
            seed_coords$chr,
            seed_coords$start,
            beta_chr,
            beta_start,
            chromosome_runs = beta_chr_runs
        )
        if (all(is.na(seed_beta_index))) {
            stop("None of the IDs in seeds_id_col match the beta genomic locations.")
        }
        if (anyNA(seed_beta_index)) {
            keep_matched <- !is.na(seed_beta_index)
            seeds_df <- seeds_df[keep_matched, , drop = FALSE]
            seed_coords$chr <- seed_coords$chr[keep_matched]
            seed_coords$start <- seed_coords$start[keep_matched]
            seed_beta_index <- seed_beta_index[keep_matched]
        }
        seeds_df$.__beta_row_index__ <- seed_beta_index
        seed_order <- order(seeds_df$.__beta_row_index__)
        seeds_df <- seeds_df[seed_order, , drop = FALSE]
        seed_coords$chr <- seed_coords$chr[seed_order]
        seed_coords$start <- seed_coords$start[seed_order]
        keep_unique <- !duplicated(seeds_df[, seeds_id_col])
        seed_ids <- as.character(seeds_df[keep_unique, seeds_id_col])
        seed_beta_index <- seeds_df[keep_unique, ".__beta_row_index__"]
        seed_chr <- as.character(beta_chr[seed_beta_index])
    } else {
        .log_info("Matching seed IDs to beta row names...", level = 2)
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
    .log_success("", level = 3)
    if (length(seed_ids) == 0L) {
        stop("No seeds remain after filtering against beta locations.")
    }
    if (!use_numeric_sequencing_rows) {
        seeds_locs <- as.data.frame(beta_locs[seed_beta_index, , drop = FALSE])
        rownames(seeds_locs) <- seed_ids
        seed_chr <- as.character(seeds_locs[, "chr"])
    }
    chromosomes <- unique(seed_chr)
    seed_pval_source <- .seedPvalSource(seeds_df, seeds_id_col)
    .log_info("Processing ", length(chromosomes), " chromosome(s): ", paste(chromosomes, collapse = ", "), level = 1)

    chromosome_tasks <- .buildDMRsChromosomeTasks(
        beta_handler = beta_handler,
        beta_chr = beta_chr,
        beta_locs_rownames = beta_locs_rownames,
        chromosomes = chromosomes,
        seed_ids = seed_ids,
        seed_beta_index = seed_beta_index,
        seed_chr = seed_chr,
        beta_col_names = beta_col_names,
        use_numeric_sequencing_rows = use_numeric_sequencing_rows,
        chromosome_runs = beta_chr_runs
    )
    task_chromosomes <- names(chromosome_tasks)
    rm(
        beta_locs,
        beta_chr,
        beta_start,
        beta_locs_rownames,
        beta_chr_runs,
        seed_chr,
        seed_ids,
        seed_beta_index,
        use_numeric_sequencing_rows,
        seeds_df,
        chromosomes
    )
    if (exists("seeds_locs", inherits = FALSE)) {
        rm(seeds_locs)
    }
    if (exists("seed_coords", inherits = FALSE)) {
        rm(seed_coords)
    }
    if (exists("bsseq_gr", inherits = FALSE)) {
        rm(bsseq_gr)
    }
    if (exists("beta_row_names", inherits = FALSE)) {
        rm(beta_row_names)
    }
    gc(FALSE)

    chr_parallel_jobs <- suppressWarnings(as.integer(getOption("CMEnt.chr_njobs", 3L)))
    if (is.na(chr_parallel_jobs) || chr_parallel_jobs < 1L) {
        chr_parallel_jobs <- 1L
    }
    chr_parallel_jobs <- min(chr_parallel_jobs, length(chromosome_tasks), max(1L, njobs))
    chr_njobs <- if (chr_parallel_jobs > 1L) max(1L, floor(njobs / chr_parallel_jobs)) else njobs
    memory_njobs <- max(1L, chr_parallel_jobs * chr_njobs)
    if (chr_parallel_jobs > 1L) {
        .log_info(
            "Using ", chr_parallel_jobs, " chromosome job(s) with ", chr_njobs,
            " inner job(s) per chromosome.",
            level = 2
        )
    }
    if (chr_parallel_jobs > 1L) {
        chr_results <- .safeBiocParallelApply(
            chromosome_tasks,
            .buildDMRsChromosomeTask,
            beta_handler = beta_handler,
            pheno_detection = pheno_detection,
            group_inds = group_inds,
            testing_mode_per_group = testing_mode_per_group,
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
            njobs = chr_njobs,
            memory_njobs = memory_njobs,
            verbose = verbose,
            .load_debug = .load_debug,
            pheno = pheno,
            beta_col_names = beta_col_names,
            sample_group_col = sample_group_col,
            covariates = covariates,
            covariate_models = covariate_models,
            annotate_with_genes = annotate_with_genes,
            score_dmrs = score_dmrs,
            extract_motifs = extract_motifs,
            BPPARAM = .makeBiocParallelParam(
                chr_parallel_jobs,
                n_tasks = length(chromosome_tasks)
            )
        )
    } else {
        chr_results <- lapply(
            chromosome_tasks,
            .buildDMRsChromosomeTask,
            beta_handler = beta_handler,
            pheno_detection = pheno_detection,
            group_inds = group_inds,
            testing_mode_per_group = testing_mode_per_group,
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
            njobs = chr_njobs,
            memory_njobs = memory_njobs,
            verbose = verbose,
            .load_debug = .load_debug,
            pheno = pheno,
            beta_col_names = beta_col_names,
            sample_group_col = sample_group_col,
            covariates = covariates,
            covariate_models = covariate_models,
            annotate_with_genes = annotate_with_genes,
            score_dmrs = score_dmrs,
            extract_motifs = extract_motifs
        )
    }
    names(chr_results) <- task_chromosomes
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
    .log_info("Total DMRs assembled: ", length(final_dmrs_granges), level = 2)
    final_ord <- order(
        as.character(GenomicRanges::seqnames(final_dmrs_granges)),
        GenomicRanges::start(final_dmrs_granges),
        GenomicRanges::end(final_dmrs_granges)
    )
    final_dmrs_granges <- final_dmrs_granges[final_ord]
    final_seed_ids <- unique(unlist(
        lapply(S4Vectors::mcols(final_dmrs_granges)$seeds, .splitCsvValues),
        use.names = FALSE
    ))
    seed_pvals <- .seedPvaluesForSelectedSeeds(seed_pval_source, final_seed_ids)
    S4Vectors::mcols(final_dmrs_granges)$pval <- .combineDMRSeedPvalues(
        S4Vectors::mcols(final_dmrs_granges)$seeds,
        seed_pvals
    )
    S4Vectors::mcols(final_dmrs_granges)$qval <- stats::p.adjust(
        S4Vectors::mcols(final_dmrs_granges)$pval,
        method = "BH"
    )
    final_ids <- as.character(S4Vectors::mcols(final_dmrs_granges)$id)
    if (length(final_ids) == length(final_dmrs_granges)) {
        names(final_dmrs_granges) <- make.unique(ifelse(is.na(final_ids) | !nzchar(final_ids), "dmr", final_ids))
    }

    if (!is.null(output_prefix_dot)) {
        viewer_sites <- unique(unlist(lapply(S4Vectors::mcols(final_dmrs_granges)$sites, .splitCsvValues), use.names = FALSE))
        viewer_beta <- beta_handler$getBeta(row_names = viewer_sites, col_names = beta_col_names_detection)
        viewer_file <- paste0(output_prefix_dot, "seeds_beta.tsv.gz")
        gz <- gzfile(viewer_file, "w")
        write.table(viewer_beta, gz, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
        close(gz)

        final_dmrs <- .convertToDataFrame(final_dmrs_granges)
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
