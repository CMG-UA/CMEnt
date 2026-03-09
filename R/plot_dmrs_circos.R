.filterDMRsForCircos <- function(dmrs, max_dmrs_per_chr) {
    dmrs_split <- split(seq_along(dmrs), as.character(GenomicRanges::seqnames(dmrs)))
    filtered_dmrs_list <- lapply(dmrs_split, function(inds) {
        chr_dmrs <- dmrs[inds]
        if (length(chr_dmrs) <= max_dmrs_per_chr) {
            return(chr_dmrs)
        }
        abs_delta_beta <- abs(chr_dmrs$delta_beta)
        score <- minmaxscale(abs_delta_beta)
        ord <- order(score, decreasing = TRUE)
        selected_indices <- ord[1:max_dmrs_per_chr]
        chr_dmrs[selected_indices]
    })
    filtered_dmrs <- unlist(GenomicRanges::GRangesList(filtered_dmrs_list))
    filtered_dmrs
}

.hasNonEmptyString <- function(x) {
    x <- as.character(x)
    !is.na(x) & nzchar(trimws(x))
}

.parseCoordinateNum <- function(x) {
    x <- trimws(as.character(x))
    is_kb <- grepl("k", x, ignore.case = TRUE)
    is_mb <- grepl("m", x, ignore.case = TRUE)
    is_gb <- grepl("g", x, ignore.case = TRUE)
    x <- sub("[k|m|g]*b?p?$", "", x, ignore.case = TRUE)
    e <- try(x <- as.numeric(x), silent = TRUE)
    if (inherits(e, "try-error") || !is.finite(x) || x < 1) {
        stop("Coordinate must be a finite numeric value >= 1. Invalid coordinate: '", x, "'")
    }
    if (is_kb) {
        x <- x * 1e3
    } else if (is_mb) {
        x <- x * 1e6
    } else if (is_gb) {
        x <- x * 1e9
    }
    if (round(x) != x) {
        warning("Coordinate value '", x, "' is not an integer. Rounding to nearest integer.")
    }
    as.integer(round(x))
}

.parseRegionString <- function(region_str) {
    clean <- trimws(unlist(strsplit(region_str, ",")))
    matched <- list()
    for (i in seq_along(clean)) {
        chr <- strsplit(clean[i], ":", fixed = TRUE)[[1]]
        if (length(chr) != 2) {
            stop("Each region must be in the format 'chr:start-end' (e.g., 'chr7:100000-200000'). Invalid region: '", clean[i], "'")
        }
        start_end <- strsplit(chr[2], "-", fixed = TRUE)[[1]]
        chr <- trimws(chr[[1]])
        if (length(start_end) != 2) {
            stop("Each region must be in the format 'chr:start-end' (e.g., 'chr7:100000-200000'). Invalid region: '", clean[i], "'")
        }
        e <- try(start <- .parseCoordinateNum(start_end[1]))
        if (inherits(e, "try-error") || !is.finite(start) || start < 1) {
            stop("Start position must be a valid genomic number in the format 'chr:start-end'. Invalid start: '", start_end[1], "' in region '", clean[i], "'")
        }
        e <- try(end <- .parseCoordinateNum(start_end[2]))
        if (inherits(e, "try-error") || !is.finite(end) || end < 1) {
            stop("End position must be a valid genomic number in the format 'chr:start-end'. Invalid end: '", start_end[2], "' in region '", clean[i], "'")
        }
        matched[[i]] <- c(chr, start, end)
    }
    matched <- do.call(rbind, matched)
    data.frame(
        chr = matched[, 1],
        start = as.numeric(matched[, 2]),
        end = as.numeric(matched[, 3]),
        stringsAsFactors = FALSE
    )
}

.normalizeCircosRegion <- function(region, cytoband = NULL) {
    if (is.null(region)) {
        return(NULL)
    }

    region_df <- NULL
    if (inherits(region, "GRanges")) {
        if (length(region) == 0) {
            return(NULL)
        }
        region_df <- data.frame(
            chr = as.character(GenomeInfoDb::seqnames(region)),
            start = GenomicRanges::start(region),
            end = GenomicRanges::end(region),
            stringsAsFactors = FALSE
        )
    } else if (is.character(region)) {
        region_df <- .parseRegionString(region)
    } else if (is.data.frame(region)) {
        req_cols <- c("chr", "start", "end")
        if (!all(req_cols %in% colnames(region))) {
            stop("region data frame must contain columns: chr, start, end.")
        }
        region_df <- data.frame(
            chr = as.character(region$chr),
            start = as.numeric(region$start),
            end = as.numeric(region$end),
            stringsAsFactors = FALSE
        )
    } else if (is.list(region) && all(c("chr", "start", "end") %in% names(region))) {
        region_df <- data.frame(
            chr = as.character(region$chr),
            start = as.numeric(region$start),
            end = as.numeric(region$end),
            stringsAsFactors = FALSE
        )
    } else {
        stop("region must be NULL, GRanges, a 'chr:start-end' string, or a data.frame/list with chr/start/end.")
    }

    region_df <- region_df[.hasNonEmptyString(region_df$chr), , drop = FALSE]
    if (nrow(region_df) == 0) {
        return(NULL)
    }
    if (any(!is.finite(region_df$start)) || any(!is.finite(region_df$end))) {
        stop("region start/end must be finite numeric values.")
    }
    region_df$start <- floor(region_df$start)
    region_df$end <- ceiling(region_df$end)
    if (any(region_df$start < 1L) || any(region_df$end < 1L)) {
        stop("region start/end must be >= 1.")
    }
    if (any(region_df$start > region_df$end)) {
        stop("region start must be <= end.")
    }
    if (!is.null(cytoband) && nrow(cytoband) > 0) {
        cytoband_gr <- GenomicRanges::GRanges(
            seqnames = cytoband$V1,
            ranges = IRanges::IRanges(start = cytoband$V2, end = cytoband$V3)
        )
        region_gr <- GenomicRanges::GRanges(
            seqnames = region_df$chr,
            ranges = IRanges::IRanges(start = region_df$start, end = region_df$end)
        )
        region_gr <- intersect(region_gr, cytoband_gr, ignore.strand = TRUE)
        if (length(region_gr) == 0) {
            return(NULL)
        }
        region_df <- as.data.frame(region_gr)[, c("seqnames", "start", "end")]
        colnames(region_df) <- c("chr", "start", "end")
    }
    unique(region_df)
}

.filterDMRsByScopeForCircos <- function(dmrs, chromosomes = NULL, region_df = NULL) {
    ret <- dmrs
    if (!is.null(chromosomes)) {
        keep_chr <- as.character(GenomicRanges::seqnames(ret)) %in% chromosomes
        ret <- ret[keep_chr]
    }
    if (!is.null(region_df) && nrow(region_df) > 0) {
        ret <- lapply(seq_len(nrow(region_df)), function(i) {
            chr <- region_df$chr[i]
            reg_start <- region_df$start[i]
            reg_end <- region_df$end[i]
            chr_dmrs <- ret[as.character(GenomicRanges::seqnames(ret)) == chr]
            if (length(chr_dmrs) == 0) {
                return(NULL)
            }
            overlaps <- GenomicRanges::findOverlaps(
                chr_dmrs,
                GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = reg_start, end = reg_end)),
                ignore.strand = TRUE
            )
            if (length(overlaps) == 0) {
                return(NULL)
            }
            overlap <- chr_dmrs[S4Vectors::queryHits(overlaps)]
            GenomicRanges::start(overlap) <- pmax(GenomicRanges::start(overlap), reg_start)
            GenomicRanges::end(overlap) <- pmin(GenomicRanges::end(overlap), reg_end)
            overlap <- overlap[GenomicRanges::start(overlap) <= GenomicRanges::end(overlap)]
            if (length(overlap) == 0) {
                return(NULL)
            }
            new_seqlevel <- sprintf("%s:%d-%d", region_df$chr[i], region_df$start[i], region_df$end[i])
            names(new_seqlevel) <- chr
            GenomeInfoDb::renameSeqlevels(overlap, new_seqlevel)
        })
        ret <- suppressWarnings(do.call(c, ret))
    }
    ret
}

.orderChromosomesNaturally <- function(chromosomes) {
    chromosomes <- unique(as.character(chromosomes))
    if (length(chromosomes) == 0) {
        return(chromosomes)
    }
    chr_clean <- gsub("^chr", "", chromosomes, ignore.case = TRUE)
    chr_num <- suppressWarnings(as.numeric(chr_clean))
    chr_special <- match(toupper(chr_clean), c("X", "Y", "M", "MT"))
    chr_special[is.na(chr_special)] <- Inf
    ord <- order(is.na(chr_num), chr_num, chr_special, chr_clean)
    chromosomes[ord]
}

.subsetCytobandForCircos <- function(cytoband, unique_chrs = NULL, region_df = NULL) {
    if (is.null(cytoband) || nrow(cytoband) == 0) {
        return(NULL)
    }
    cb <- cytoband
    if (!is.null(unique_chrs)) {
        unique_chrs <- as.character(unique_chrs)
        cb <- cb[cb$V1 %in% lapply(strsplit(unique_chrs, ":"), function(x) x[1]), , drop = FALSE]
    }
    if (nrow(cb) == 0) {
        return(NULL)
    }
    cb <- cb[order(match(cb$V1, unique_chrs), cb$V2, cb$V3), , drop = FALSE]

    if (!is.null(region_df) && nrow(region_df) > 0) {
        trimmed <- do.call(rbind, lapply(seq_len(nrow(region_df)), function(i) {
            chr <- region_df$chr[i]
            chr_cb <- cb[cb$V1 == chr, , drop = FALSE]
            if (nrow(chr_cb) == 0) {
                return(NULL)
            }
            reg_start <- region_df$start[i]
            reg_end <- region_df$end[i]
            overlap <- chr_cb[chr_cb$V3 >= reg_start & chr_cb$V2 <= reg_end, , drop = FALSE]
            if (nrow(overlap) == 0) {
                return(NULL)
            }
            overlap$V1 <- sprintf("%s:%d-%d", region_df$chr[i], reg_start, reg_end)
            overlap$V2 <- pmax(overlap$V2, reg_start)
            overlap$V3 <- pmin(overlap$V3, reg_end)
            overlap[overlap$V2 < overlap$V3, , drop = FALSE]
        }))
        if (is.null(trimmed) || nrow(trimmed) == 0) {
            return(NULL)
        }
        cb <- unique(trimmed)
    }

    if (nrow(cb) == 0) {
        return(NULL)
    }
    cb
}

.orderCircosInteractionRows <- function(link_data) {
    if (is.null(link_data) || nrow(link_data) == 0) {
        return(link_data)
    }
    if ("component_best_rank" %in% colnames(link_data) && any(is.finite(link_data$component_best_rank))) {
        rank_vec <- link_data$component_best_rank
        rank_vec[!is.finite(rank_vec)] <- Inf
        ord <- order(rank_vec, -link_data$sim, link_data$component_id)
    } else {
        ord <- order(-link_data$sim, link_data$component_id)
    }
    link_data[ord, , drop = FALSE]
}

.selectCircosInteractions <- function(link_data, max_components) {
    link_data <- .orderCircosInteractionRows(link_data)
    if (is.null(link_data) || nrow(link_data) == 0) {
        return(link_data)
    }
    if (is.null(max_components) || !is.finite(max_components)) {
        return(link_data)
    }
    max_components <- as.integer(max_components)
    if (max_components < 1L) {
        return(link_data[0, , drop = FALSE])
    }
    if (nrow(link_data) <= max_components) {
        return(link_data)
    }

    matched_mask <- !is.na(link_data$has_jaspar_match) & link_data$has_jaspar_match
    matched <- .orderCircosInteractionRows(link_data[matched_mask, , drop = FALSE])
    unmatched <- .orderCircosInteractionRows(link_data[!matched_mask, , drop = FALSE])
    unique_matched <- unique(matched$component_id)
    n_matched_take <- min(length(unique_matched), max_components)
    matched_comps <- unique_matched[seq_len(n_matched_take)]

    ret <- if (n_matched_take > 0) {
        matched[matched$component_id %in% matched_comps, , drop = FALSE]
    } else {
        matched[0, , drop = FALSE]
    }
    if (length(matched_comps) < max_components && nrow(unmatched) > 0) {
        unique_unmatched <- unique(unmatched$component_id)
        n_unmatched_take <- min(length(unique_unmatched), max_components - length(matched_comps))
        unmatched_comps <- unique_unmatched[seq_len(n_unmatched_take)]
        ret <- rbind(ret, unmatched[unmatched$component_id %in% unmatched_comps, , drop = FALSE])
    }
    ret
}

.wrapCircosLegendLabel <- function(label, width = 88L) {
    paste(strwrap(label, width = width), collapse = "\n")
}

.emptyCircosCandidateFrame <- function() {
    data.frame(
        chr = character(0),
        start = numeric(0),
        end = numeric(0),
        priority = numeric(0),
        annotation_priority = numeric(0),
        size = numeric(0),
        mean_abs_delta = numeric(0),
        source = character(0),
        source_id = character(0),
        stringsAsFactors = FALSE
    )
}

.safeFiniteMean <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0) {
        NA_real_
    } else {
        mean(x)
    }
}

.getCircosDMRPriority <- function(dmrs) {
    n <- length(dmrs)
    mcols_df <- S4Vectors::mcols(dmrs)
    ranks <- if ("rank" %in% colnames(mcols_df)) {
        suppressWarnings(as.numeric(mcols_df$rank))
    } else {
        rep(NA_real_, n)
    }
    scores <- if ("score" %in% colnames(mcols_df)) {
        suppressWarnings(as.numeric(mcols_df$score))
    } else {
        rep(NA_real_, n)
    }
    abs_delta <- if ("delta_beta" %in% colnames(mcols_df)) {
        abs(suppressWarnings(as.numeric(mcols_df$delta_beta)))
    } else {
        rep(NA_real_, n)
    }
    ord <- order(
        ifelse(is.finite(ranks), ranks, Inf),
        ifelse(is.finite(scores), -scores, Inf),
        ifelse(is.finite(abs_delta), -abs_delta, Inf),
        seq_len(n)
    )
    priority <- integer(n)
    priority[ord] <- seq_len(n)
    priority
}

.buildCircosBlockCandidates <- function(dmrs, region_flank_bp) {
    if (!"block_id" %in% colnames(S4Vectors::mcols(dmrs))) {
        return(.emptyCircosCandidateFrame())
    }
    block_ids <- as.character(S4Vectors::mcols(dmrs)$block_id)
    keep <- .hasNonEmptyString(block_ids)
    if (!any(keep)) {
        return(.emptyCircosCandidateFrame())
    }
    abs_delta <- if ("delta_beta" %in% colnames(S4Vectors::mcols(dmrs))) {
        abs(suppressWarnings(as.numeric(S4Vectors::mcols(dmrs)$delta_beta)))
    } else {
        rep(NA_real_, length(dmrs))
    }
    block_df <- data.frame(
        chr = as.character(GenomicRanges::seqnames(dmrs))[keep],
        start = GenomicRanges::start(dmrs)[keep],
        end = GenomicRanges::end(dmrs)[keep],
        block_id = block_ids[keep],
        priority = .getCircosDMRPriority(dmrs)[keep],
        abs_delta = abs_delta[keep],
        stringsAsFactors = FALSE
    )
    split_keys <- split(seq_len(nrow(block_df)), paste(block_df$chr, block_df$block_id, sep = "::"))
    candidates <- do.call(rbind, lapply(split_keys, function(idx) {
        d <- block_df[idx, , drop = FALSE]
        data.frame(
            chr = d$chr[1],
            start = pmax(1, floor(min(d$start) - region_flank_bp)),
            end = ceiling(max(d$end) + region_flank_bp),
            priority = min(d$priority),
            annotation_priority = 0,
            size = nrow(d),
            mean_abs_delta = .safeFiniteMean(d$abs_delta),
            source = "block",
            source_id = d$block_id[1],
            stringsAsFactors = FALSE
        )
    }))
    candidates <- candidates[order(
        candidates$priority,
        -candidates$annotation_priority,
        -candidates$size,
        -ifelse(is.finite(candidates$mean_abs_delta), candidates$mean_abs_delta, -Inf),
        candidates$chr,
        candidates$start
    ), , drop = FALSE]
    rownames(candidates) <- NULL
    candidates
}

.buildCircosQuickCandidates <- function(dmrs, region_flank_bp) {
    abs_delta <- if ("delta_beta" %in% colnames(S4Vectors::mcols(dmrs))) {
        abs(suppressWarnings(as.numeric(S4Vectors::mcols(dmrs)$delta_beta)))
    } else {
        rep(NA_real_, length(dmrs))
    }
    candidates <- data.frame(
        chr = as.character(GenomicRanges::seqnames(dmrs)),
        start = pmax(1, floor(GenomicRanges::start(dmrs) - region_flank_bp)),
        end = ceiling(GenomicRanges::end(dmrs) + region_flank_bp),
        priority = .getCircosDMRPriority(dmrs),
        annotation_priority = 0,
        size = 1,
        mean_abs_delta = abs_delta,
        source = "dmr",
        source_id = as.character(seq_len(length(dmrs))),
        stringsAsFactors = FALSE
    )
    candidates <- candidates[order(
        candidates$priority,
        -candidates$annotation_priority,
        -ifelse(is.finite(candidates$mean_abs_delta), candidates$mean_abs_delta, -Inf),
        candidates$chr,
        candidates$start
    ), , drop = FALSE]
    rownames(candidates) <- NULL
    candidates
}

.buildCircosComponentCandidates <- function(dmrs, components, region_flank_bp) {
    if (is.null(components) || nrow(components) == 0 || !"indices" %in% colnames(components)) {
        return(.emptyCircosCandidateFrame())
    }
    abs_delta <- if ("delta_beta" %in% colnames(S4Vectors::mcols(dmrs))) {
        abs(suppressWarnings(as.numeric(S4Vectors::mcols(dmrs)$delta_beta)))
    } else {
        rep(NA_real_, length(dmrs))
    }
    dmr_priority <- .getCircosDMRPriority(dmrs)
    candidates <- lapply(seq_len(nrow(components)), function(i) {
        idxs <- suppressWarnings(as.integer(components$indices[[i]]))
        idxs <- unique(idxs[is.finite(idxs) & idxs >= 1 & idxs <= length(dmrs)])
        if (length(idxs) == 0) {
            return(NULL)
        }
        has_jaspar_match <- if ("jaspar_names" %in% colnames(components)) {
            .hasNonEmptyString(components$jaspar_names[[i]])
        } else {
            FALSE
        }
        component_size <- if ("size" %in% colnames(components)) {
            suppressWarnings(as.numeric(components$size[[i]]))
        } else {
            length(idxs)
        }
        split_idxs <- split(idxs, as.character(GenomicRanges::seqnames(dmrs)[idxs]))
        do.call(rbind, lapply(split_idxs, function(chr_idxs) {
            data.frame(
                chr = as.character(GenomicRanges::seqnames(dmrs)[chr_idxs[1]]),
                start = pmax(1, floor(min(GenomicRanges::start(dmrs)[chr_idxs]) - region_flank_bp)),
                end = ceiling(max(GenomicRanges::end(dmrs)[chr_idxs]) + region_flank_bp),
                priority = min(dmr_priority[chr_idxs]),
                annotation_priority = if (has_jaspar_match) 1 else 0,
                size = component_size,
                mean_abs_delta = .safeFiniteMean(abs_delta[chr_idxs]),
                source = "component",
                source_id = as.character(components$component_id[[i]]),
                stringsAsFactors = FALSE
            )
        }))
    })
    candidates <- do.call(rbind, candidates)
    if (is.null(candidates) || nrow(candidates) == 0) {
        return(.emptyCircosCandidateFrame())
    }
    candidates <- candidates[order(
        candidates$priority,
        -candidates$annotation_priority,
        -candidates$size,
        -ifelse(is.finite(candidates$mean_abs_delta), candidates$mean_abs_delta, -Inf),
        candidates$chr,
        candidates$start
    ), , drop = FALSE]
    rownames(candidates) <- NULL
    candidates
}

.interleaveCircosCandidates <- function(primary_candidates, secondary_candidates, max_pool) {
    primary_candidates <- primary_candidates[seq_len(min(nrow(primary_candidates), max_pool)), , drop = FALSE]
    secondary_candidates <- secondary_candidates[seq_len(min(nrow(secondary_candidates), max_pool)), , drop = FALSE]
    ret <- list()
    i <- 1L
    j <- 1L
    while (i <= nrow(primary_candidates) || j <= nrow(secondary_candidates)) {
        if (i <= nrow(primary_candidates)) {
            ret[[length(ret) + 1L]] <- primary_candidates[i, , drop = FALSE]
            i <- i + 1L
        }
        if (j <= nrow(secondary_candidates)) {
            ret[[length(ret) + 1L]] <- secondary_candidates[j, , drop = FALSE]
            j <- j + 1L
        }
    }
    if (length(ret) == 0) {
        return(.emptyCircosCandidateFrame())
    }
    do.call(rbind, ret)
}

.greedySelectCircosRegions <- function(candidates, n_regions, max_regions_per_chr, min_inter_region_bp) {
    if (is.null(candidates) || nrow(candidates) == 0) {
        return(data.frame(chr = character(0), start = numeric(0), end = numeric(0)))
    }
    selected <- candidates[0, , drop = FALSE]
    selected$selection_order <- numeric(0)
    for (i in seq_len(nrow(candidates))) {
        cand <- candidates[i, , drop = FALSE]
        same_chr <- selected$chr == cand$chr
        near_hits <- which(
            same_chr &
                cand$start <= (selected$end + min_inter_region_bp) &
                cand$end >= (selected$start - min_inter_region_bp)
        )
        if (length(near_hits) > 0) {
            keep_hit <- near_hits[1]
            merged_start <- min(c(selected$start[near_hits], cand$start))
            merged_end <- max(c(selected$end[near_hits], cand$end))
            merged_priority <- min(c(selected$priority[near_hits], cand$priority))
            merged_annotation_priority <- max(c(selected$annotation_priority[near_hits], cand$annotation_priority), na.rm = TRUE)
            merged_size <- sum(c(selected$size[near_hits], cand$size), na.rm = TRUE)
            merged_abs <- .safeFiniteMean(c(selected$mean_abs_delta[near_hits], cand$mean_abs_delta))
            merged_source <- paste(unique(c(selected$source[near_hits], cand$source)), collapse = ",")
            merged_source_id <- paste(unique(c(selected$source_id[near_hits], cand$source_id)), collapse = ",")
            selected$start[keep_hit] <- merged_start
            selected$end[keep_hit] <- merged_end
            selected$priority[keep_hit] <- merged_priority
            selected$annotation_priority[keep_hit] <- merged_annotation_priority
            selected$size[keep_hit] <- merged_size
            selected$mean_abs_delta[keep_hit] <- merged_abs
            selected$source[keep_hit] <- merged_source
            selected$source_id[keep_hit] <- merged_source_id
            if (length(near_hits) > 1) {
                selected <- selected[-near_hits[-1], , drop = FALSE]
            }
            next
        }
        chr_count <- sum(selected$chr == cand$chr)
        if (nrow(selected) >= n_regions) {
            next
        }
        if (!is.null(max_regions_per_chr) && is.finite(max_regions_per_chr) && chr_count >= max_regions_per_chr) {
            next
        }
        cand$selection_order <- nrow(selected) + 1L
        selected <- rbind(selected, cand)
    }
    if (nrow(selected) == 0) {
        return(data.frame(chr = character(0), start = numeric(0), end = numeric(0)))
    }
    chr_levels <- .orderChromosomesNaturally(selected$chr)
    selected <- selected[order(
        selected$priority,
        -selected$annotation_priority,
        -selected$size,
        factor(selected$chr, levels = chr_levels),
        selected$start
    ), , drop = FALSE]
    rownames(selected) <- NULL
    selected[, c("chr", "start", "end"), drop = FALSE]
}

.selectCircosRegions <- function(
    dmrs,
    method = c("blocks", "components", "hybrid", "quick"),
    n_regions = 6L,
    region_flank_bp = 1e6,
    max_regions_per_chr = 2L,
    min_inter_region_bp = 5e6,
    components = NULL
) {
    method <- match.arg(method)
    if (!is.numeric(n_regions) || length(n_regions) != 1 || is.na(n_regions) || n_regions < 1) {
        stop("n_regions must be a single numeric value >= 1.")
    }
    if (!is.numeric(region_flank_bp) || length(region_flank_bp) != 1 || is.na(region_flank_bp) || region_flank_bp < 0) {
        stop("region_flank_bp must be a single numeric value >= 0.")
    }
    if (!is.null(max_regions_per_chr)) {
        if (!is.numeric(max_regions_per_chr) || length(max_regions_per_chr) != 1 || is.na(max_regions_per_chr) || max_regions_per_chr < 1) {
            stop("max_regions_per_chr must be NULL or a single numeric value >= 1.")
        }
        max_regions_per_chr <- as.integer(max_regions_per_chr)
    }
    if (!is.numeric(min_inter_region_bp) || length(min_inter_region_bp) != 1 || is.na(min_inter_region_bp) || min_inter_region_bp < 0) {
        stop("min_inter_region_bp must be a single numeric value >= 0.")
    }
    block_candidates <- .buildCircosBlockCandidates(dmrs, region_flank_bp = region_flank_bp)
    quick_candidates <- .buildCircosQuickCandidates(dmrs, region_flank_bp = region_flank_bp)
    component_candidates <- .buildCircosComponentCandidates(
        dmrs = dmrs,
        components = components,
        region_flank_bp = region_flank_bp
    )

    candidates <- switch(
        method,
        blocks = block_candidates,
        components = component_candidates,
        hybrid = .interleaveCircosCandidates(
            primary_candidates = block_candidates,
            secondary_candidates = component_candidates,
            max_pool = max(2L, as.integer(n_regions) * 2L)
        ),
        quick = quick_candidates
    )
    if (method == "blocks" && nrow(candidates) == 0) {
        .log_warn("No finite block_id values found. Falling back to quick DMR-based Circos region selection.")
        candidates <- quick_candidates
    }
    if (method == "components" && nrow(candidates) == 0) {
        .log_warn("No retained interaction components found. Falling back to quick DMR-based Circos region selection.")
        candidates <- quick_candidates
    }
    if (method == "hybrid" && nrow(candidates) == 0) {
        .log_warn("No block or component candidates found. Falling back to quick DMR-based Circos region selection.")
        candidates <- quick_candidates
    }
    selected <- .greedySelectCircosRegions(
        candidates = candidates,
        n_regions = as.integer(n_regions),
        max_regions_per_chr = max_regions_per_chr,
        min_inter_region_bp = as.numeric(min_inter_region_bp)
    )
    if (nrow(selected) == 0) {
        stop("Automatic Circos region selection did not yield any regions.")
    }
    selected$chr <- as.character(selected$chr)
    selected$start <- as.integer(round(selected$start))
    selected$end <- as.integer(round(selected$end))
    rownames(selected) <- NULL
    selected
}

.formatCircosRegionSelection <- function(region_df) {
    paste(sprintf("%s:%s-%s", region_df$chr, format(region_df$start, big.mark = ","), format(region_df$end, big.mark = ",")), collapse = ", ")
}

.computeCircosInteractionState <- function(dmrs,
                                          genome,
                                          array,
                                          sorted_locs,
                                          min_similarity,
                                          motif_cpg_flank_size,
                                          min_component_size,
                                          query_components_with_jaspar,
                                          components = NULL,
                                          interactions = NULL) {
    if (xor(is.null(components), is.null(interactions))) {
        stop("components and interactions must both be NULL or both be provided.")
    }
    if (!is.null(components) && !is.null(interactions)) {
        return(list(
            dmrs = dmrs,
            components = components,
            interactions = interactions
        ))
    }
    beta_locs <- NULL
    if (!is.null(sorted_locs)) {
        beta_locs <- sorted_locs
    }
    computeDMRsInteraction(
        dmrs = dmrs,
        genome = genome,
        array = array,
        min_similarity = min_similarity,
        beta_locs = beta_locs,
        motif_cpg_flank_size = motif_cpg_flank_size,
        min_component_size = min_component_size,
        query_components_with_jaspar = query_components_with_jaspar
    )
}

.remapCircosInteractionsToDMRs <- function(dmrs, interactions, components = NULL) {
    if (is.null(interactions) || nrow(interactions) == 0) {
        return(list(interactions = interactions, components = components))
    }
    dmr_keys <- paste(
        as.character(GenomicRanges::seqnames(dmrs)),
        GenomicRanges::start(dmrs),
        GenomicRanges::end(dmrs),
        sep = "::"
    )
    key1 <- paste(interactions$chr1, interactions$start1, interactions$end1, sep = "::")
    key2 <- paste(interactions$chr2, interactions$start2, interactions$end2, sep = "::")
    match1 <- match(key1, dmr_keys)
    match2 <- match(key2, dmr_keys)
    keep <- !is.na(match1) & !is.na(match2)
    interactions <- interactions[keep, , drop = FALSE]
    if (nrow(interactions) == 0) {
        if (!is.null(components) && "component_id" %in% colnames(components)) {
            components <- components[0, , drop = FALSE]
        }
        return(list(interactions = interactions, components = components))
    }
    interactions$index1 <- match1[keep]
    interactions$index2 <- match2[keep]
    if (!is.null(components) && "component_id" %in% colnames(components) && "component_id" %in% colnames(interactions)) {
        used_components <- unique(interactions$component_id[!is.na(interactions$component_id)])
        components <- components[components$component_id %in% used_components, , drop = FALSE]
    }
    list(interactions = interactions, components = components)
}

#' Plot Circos Visualization of DMRs
#'
#' @description Creates a circular genome plot using circlize showing DMRs with multiple tracks:
#' ideogram track with chromosome bands, DMR arcs colored by delta beta,
#' beta value heatmaps for each sample, and motif-based interaction links between DMRs.
#'
#' @param dmrs GRanges object or data frame. DMR results from findDMRsFromSeeds.
#' @param beta BetaHandler object, character path to beta file, or beta values matrix.
#' @param pheno Data frame or character path to phenotype file. Sample information with
#'   rownames matching beta column names (required for beta track).
#' @param genome Character. Genome version (e.g., "hg38").
#' @param array Character. Array platform type (default: "450K"). Ignored if sorted_locs is provided.
#' @param sorted_locs Data frame. Genomic locations sorted by position (optional). If NULL, will be fetched based on array and genome.
#' @param components Data frame. Output from motif component detection (optional, will be computed if missing).
#' @param interactions Data frame. Output from motif interaction detection (optional, will be computed if missing).
#' @param sample_group_col Character. Column in pheno for sample grouping (default: "Sample_Group").
#' @param min_similarity Numeric. Minimum motifs PWM similarity threshold for considering DMRs are related (default: 0.8).
#' @param motif_cpg_flank_size Integer. Flanking region size for motif extraction in bp (default: 5).
#' @param max_num_samples_per_group Integer. Maximum number of samples to show per group in heatmap (default: 5).
#' @param max_dmrs_per_chr Integer. Maximum number of DMRs to use per chromosome (default: 10). The DMRs with highest absolute delta beta will be selected.
#' @param max_cpgs_per_dmr Integer. Maximum number of CpGs to show per DMR in scatter/heatmap (default: 5).
#' @param min_component_size Integer. Minimum motif component size to retain (default: 2).
#' @param max_components Integer. Maximum number of interactions to plot (default: 30).
#' @param chromosomes Character vector. Subset of chromosomes to display (default: NULL, show all available).
#' @param region Genomic region to display. Can be NULL, a GRanges, a string in the form `chr:start-end`,
#'   or a data.frame/list with columns `chr`, `start`, `end`.
#' @param query_components_with_jaspar Logical. Whether computed motif components should be queried against JASPAR
#'   before plotting (default: `TRUE`). Set to `FALSE` to keep link computation cheaper.
#' @param unmatched_interaction_color Character. Color used for interaction components without JASPAR matches.
#'   These links are shown but omitted from the interaction legend (default: `"#B3B3B3"`).
#' @param legend_width_ratio Numeric. Fraction of horizontal canvas reserved for legends (default: 0.34).
#' @param degenerate_resolution Integer. Resolution in base pairs for simplifying narrow glyphs:
#'   link ribbons are drawn as lines when both anchors are below this span, and DMR arcs
#'   are drawn as lines instead of rectangles below this span (default: 1e6).
#' @param output_file Character or NULL. Optional PDF path for the plot.
#' @param verbose Numeric. Optional verbosity override.
#' @param ... Additional arguments (currently unused).
#'
#' @return `NULL` invisibly after drawing the plot.
#'
#' @examples
#' \dontrun{
#' dmrs <- readRDS("dmrs.rds")
#' plotDMRsCircos(dmrs, beta = "beta.txt", pheno = pheno_df)
#' plotDMRsCircos(dmrs, beta = "beta.txt", pheno = pheno_df, genome = "hg38")
#' }
#'
#' @importFrom circlize circos.initializeWithIdeogram circos.trackPlotRegion circos.genomicHeatmap
#' @importFrom circlize circos.genomicLink circos.clear colorRamp2 circos.rect circos.link CELL_META
#' @export
plotDMRsCircos <- function(dmrs,
                           beta,
                           pheno,
                           genome = "hg38",
                           array = "450K",
                           sorted_locs = NULL,
                           components = NULL,
                           interactions = NULL,
                           sample_group_col = "Sample_Group",
                           min_similarity = 0.8,
                           motif_cpg_flank_size = 5,
                           max_num_samples_per_group = 5,
                           max_dmrs_per_chr = 10,
                           max_cpgs_per_dmr = 5,
                           min_component_size = 2,
                           max_components = 30,
                           chromosomes = NULL,
                           region = NULL,
                           query_components_with_jaspar = TRUE,
                           unmatched_interaction_color = "#B3B3B3",
                           legend_width_ratio = 0.34,
                           degenerate_resolution = 1e6,
                           output_file = NULL,
                           verbose = NULL,
                           ...) {
    dmrs <- convertToGRanges(dmrs, genome)
    if (!is.null(max_components)) {
        if (!is.numeric(max_components) || length(max_components) != 1 || is.na(max_components)) {
            stop("max_components must be NULL or a single numeric value.")
        }
        if (max_components < 1) {
            stop("max_components must be >= 1 when provided.")
        }
        max_components <- as.integer(max_components)
    }
    if (!is.logical(query_components_with_jaspar) || length(query_components_with_jaspar) != 1 || is.na(query_components_with_jaspar)) {
        stop("query_components_with_jaspar must be TRUE or FALSE.")
    }
    if (!is.character(unmatched_interaction_color) || length(unmatched_interaction_color) != 1 || !nzchar(unmatched_interaction_color)) {
        stop("unmatched_interaction_color must be a single non-empty color string.")
    }
    if (!is.numeric(legend_width_ratio) || length(legend_width_ratio) != 1 || is.na(legend_width_ratio) ||
            legend_width_ratio <= 0.05 || legend_width_ratio >= 0.80) {
        stop("legend_width_ratio must be a single numeric value in (0.05, 0.80).")
    }

    requested_chrs <- NULL
    if (!is.null(chromosomes)) {
        if (!is.character(chromosomes)) {
            stop("chromosomes must be a character vector when provided.")
        }
        requested_chrs <- unique(trimws(chromosomes))
        requested_chrs <- requested_chrs[nzchar(requested_chrs)]
        if (length(requested_chrs) == 0) {
            stop("chromosomes must contain at least one non-empty chromosome name.")
        }
    }
    cytoband <- .getCytobandData(genome)

    if (!is.null(verbose)) {
        options("DMRsegal.verbose" = verbose)
    }
    verbose <- getOption("DMRsegal.verbose", default = 2)
    beta_handler <- getBetaHandler(
        beta = beta,
        array = array,
        genome = genome,
        sorted_locs = sorted_locs
    )

    beta_locs <- beta_handler$getBetaLocs()
    region_df <- .normalizeCircosRegion(region, cytoband)

    if (!is.null(region_df)) {
        beta_locs <- convertToDataFrame(.filterDMRsByScopeForCircos(
            convertToGRanges(beta_locs, genome),
            chromosomes = requested_chrs,
            region_df = region_df
        ))
        beta_locs <- beta_locs[-c(4, 5)]
    }
    if (!"pwm" %in% colnames(mcols(dmrs))) {
        .log_info("DMR motifs not precomputed. Extracting motifs...", level = 2)
        dmrs <- extractDMRMotifs(dmrs, genome, array, beta_locs = beta_locs, motif_cpg_flank_size = motif_cpg_flank_size)
    }
    dmrs <- .filterDMRsByScopeForCircos(dmrs, chromosomes = requested_chrs, region_df = region_df)
    if (length(dmrs) == 0) {
        stop("No DMRs remain after applying chromosome/region filters.")
    }
    present_chrs <- unique(as.character(GenomicRanges::seqnames(dmrs)))
    unique_chrs <- .orderChromosomesNaturally(present_chrs)

    if (!is.null(region_df)) {
        region_df <- region_df[
            sprintf("%s:%d-%d", region_df$chr, region_df$start, region_df$end) %in% unique_chrs, ,
            drop = FALSE
        ]
        if (nrow(region_df) == 0) {
            stop("No region entries overlap the selected chromosomes/DMRs.")
        }
    }

    if (is.character(pheno) && length(pheno) == 1 && file.exists(pheno)) {
        pheno <- read.table(pheno, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
    }
    if (!is.data.frame(pheno)) {
        stop("pheno must be a data frame or a valid file path to a phenotype table")
    }
    if (!(sample_group_col %in% colnames(pheno))) {
        stop(sprintf("sample_group_col '%s' not found in pheno data frame", sample_group_col))
    }

    .log_step("Preparing data for Circos plot...")

    .log_step("Filtering DMRs for plotting by maximum absolute delta beta...", level = 2)
    heatmap_dmrs <- .filterDMRsForCircos(dmrs, max_dmrs_per_chr)
    .log_success("DMRs filtered for Circos plot heatmap", level = 2)
    .log_info("Total DMRs to plot on heatmap: ", length(heatmap_dmrs), level = 2)

    cytoband_subset <- .subsetCytobandForCircos(cytoband, unique_chrs, region_df = region_df)

    .log_step("Preparing DMRs data...", level = 2)
    arc_data <- .prepareCircosArcData(dmrs)
    .log_success("DMR arcs data prepared", level = 2)
    if (getOption("DMRsegal.make_debug_dir", FALSE)) {
        .log_info("Saving DMR arcs data to debug/circos_arc_data.tsv", level = 1)
        dir.create("debug", showWarnings = FALSE)
        write.table(arc_data, file = "debug/circos_arc_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    }

    .log_step("Preparing heatmap data...", level = 2)
    heatmap_data <- .prepareCircosHeatmapData(
        heatmap_dmrs, beta_handler, pheno, sample_group_col,
        max_cpgs_per_dmr, max_num_samples_per_group
    )

    heatmap_df <- heatmap_data$heatmap_df
    if (!is.null(region_df)) {
        heatmap_df <- convertToDataFrame(
            .filterDMRsByScopeForCircos(
                convertToGRanges(heatmap_df, genome),
                chromosomes = requested_chrs,
                region_df = region_df
            )
        )
        heatmap_df <- heatmap_df[-c(4, 5)]
    }
    reduced_pheno <- heatmap_data$reduced_pheno
    if (!is.null(heatmap_df) && nrow(heatmap_df) > 1) {
        ord_heatmap <- order(match(heatmap_df$chr, unique_chrs), heatmap_df$start, heatmap_df$end)
        heatmap_df <- heatmap_df[ord_heatmap, , drop = FALSE]
    }
    .log_success("Heatmap data prepared", level = 2)
    .log_info("Total heatmap entries: ", nrow(heatmap_df), level = 2)
    if (getOption("DMRsegal.make_debug_dir", FALSE)) {
        .log_info("Saving heatmap data to debug/circos_heatmap_data.tsv", level = 1)
        dir.create("debug", showWarnings = FALSE)
        write.table(heatmap_df, file = "debug/circos_heatmap_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    }

    .log_step("Computing DMRs links...", level = 2)
    link_data <- .prepareCircosLinkData(
        dmrs = dmrs,
        genome = genome,
        array = array,
        beta_locs = beta_locs,
        min_similarity = min_similarity,
        motif_cpg_flank_size = motif_cpg_flank_size,
        max_components = max_components,
        min_component_size = min_component_size,
        query_components_with_jaspar = query_components_with_jaspar,
        components = components,
        interactions = interactions
    )
    .log_success("DMR links data prepared", level = 2)
    if (getOption("DMRsegal.make_debug_dir", FALSE) && !is.null(link_data) && nrow(link_data) > 0) {
        .log_info("Saving DMR links data to debug/circos_link_data.tsv", level = 1)
        dir.create("debug", showWarnings = FALSE)
        write.table(link_data, file = "debug/circos_link_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    }
    if (!is.null(output_file)) {
        grDevices::cairo_pdf(output_file, width = 14, height = 10, fallback_resolution = 300)
    }
    if (.Device == "null device") {
        grDevices::cairo_pdf(width = 14, height = 10, fallback_resolution = 300)
    }

    .log_step("Creating Circos plot...")
    on.exit(try(circlize::circos.clear(), silent = TRUE), add = TRUE)
    graphics::plot.new()
    circle_width_ratio <- 1 - legend_width_ratio
    root_layout <- grid::viewport(layout = grid::grid.layout(
        nrow = 1, ncol = 2,
        widths = grid::unit.c(grid::unit(circle_width_ratio, "npc"), grid::unit(legend_width_ratio, "npc"))
    ))
    grid::pushViewport(root_layout)
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::pushViewport(grid::viewport(width = grid::unit(1, "snpc"), height = grid::unit(1, "snpc")))
    graphics::par(omi = gridBase::gridOMI(), new = TRUE)
    circlize::circos.par(gap.after = c(rep(2, length(unique_chrs) - 1), 10))

    heatmap_legend <- NULL
    arc_legend <- NULL
    link_legend <- NULL
    use_manual_init <- FALSE
    if (!is.null(cytoband_subset) && nrow(cytoband_subset) > 0) {
        circlize::circos.initializeWithIdeogram(
            cytoband = cytoband_subset,
            plotType = c("labels", "axis"),
            tickLabelsStartFromZero = FALSE
        )
    } else {
        circlize::circos.initializeWithIdeogram(
            species = genome,
            chromosome.index = unique_chrs,
            plotType = c("labels", "axis")
        )
    }
    if (!is.null(heatmap_df) && nrow(heatmap_df) > 0) {
        .log_step("Adding heatmap track...", level = 2)
        heatmap_numeric <- as.matrix(heatmap_df[, -(1:3), drop = FALSE])
        storage.mode(heatmap_numeric) <- "numeric"
        valid_vals <- as.numeric(heatmap_numeric)
        valid_vals <- valid_vals[is.finite(valid_vals)]
        if (length(valid_vals) == 0) {
            .log_warn("Heatmap beta values are all missing/non-finite. Skipping heatmap track.")
        } else {
            q <- stats::quantile(valid_vals, probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE, type = 8)
            vals_range <- range(valid_vals, na.rm = TRUE)
            if (vals_range[1] <= 0.5 && vals_range[2] >= 0.5) {
                q <- sort(c(q[1], 0.5, q[2]))
                col_fun <- circlize::colorRamp2(q, c("#2b83ba", "#f7f7f7", "#d7191c"))
            } else if (vals_range[2] < 0.5) {
                col_fun <- circlize::colorRamp2(q, c("#2b83ba", "#f7f7f7"))
            } else {
                col_fun <- circlize::colorRamp2(q, c("#f7f7f7", "#d7191c"))
            }
            heatmap_height <- 0.3
            heatmap_df_plot <- data.frame(
                heatmap_df[, 1:3, drop = FALSE],
                heatmap_numeric,
                check.names = FALSE,
                stringsAsFactors = FALSE
            )
            previous_track_index <- circlize::get.current.track.index()
            circlize::circos.genomicHeatmap(
                bed = heatmap_df_plot,
                col = col_fun,
                border = NA,
                side = "outside",
                border_lwd = 0,
                line_col = "#6E6E6E",
                line_lwd = 0.45,
                heatmap_height = heatmap_height
            )
            heatmap_track_index <- previous_track_index + 1
        }

        if (length(valid_vals) > 0) {
            suppressMessages(circlize::circos.track(track.index = heatmap_track_index, panel.fun = function(x, y) {
                if (circlize::CELL_META$sector.numeric.index == length(unique_chrs)) {
                    groups <- as.character(reduced_pheno[[sample_group_col]])
                    group_runs <- rle(groups)
                    unique_groups <- unique(group_runs$values)
                    group_colors <- colorspace::qualitative_hcl(length(unique_groups), palette = "Pastel 1")
                    names(group_colors) <- unique_groups
                    y_limits <- circlize::CELL_META$cell.ylim
                    n_samples <- sum(group_runs$lengths)
                    if (n_samples <= 0) {
                        return(invisible(NULL))
                    }
                    row_height <- (y_limits[2] - y_limits[1]) / n_samples
                    y_cursor <- y_limits[1]
                    for (i in seq_along(group_runs$values)) {
                        group <- group_runs$values[i]
                        group_size <- group_runs$lengths[i]
                        y_next <- y_cursor + group_size * row_height
                        y_bottom <- min(y_cursor, y_next)
                        y_top <- max(y_cursor, y_next)
                        circlize::circos.rect(
                            circlize::CELL_META$cell.xlim[2] + circlize::convert_x(1, "mm"),
                            y_bottom,
                            circlize::CELL_META$cell.xlim[2] + circlize::convert_x(5, "mm"),
                            y_top,
                            col = group_colors[group],
                            border = NA
                        )
                        circlize::circos.text(
                            circlize::CELL_META$cell.xlim[2] + circlize::convert_x(3, "mm"),
                            y_bottom + (y_top - y_bottom) / 2,
                            group,
                            cex = 0.5,
                            facing = "clockwise"
                        )
                        y_cursor <- y_next
                    }
                }
            }, bg.border = NA))

            legend_at <- signif(q, 2)
            heatmap_legend <- ComplexHeatmap::Legend(
                title = "Samples \u03b2-values",
                at = legend_at,
                col_fun = col_fun,
                title_position = "topleft",
                legend_height = grid::unit(4, "cm"),
                labels_gp = grid::gpar(fontsize = 8),
                title_gp = grid::gpar(fontsize = 10, fontface = "bold")
            )
            .log_success("Heatmap track added", level = 2)
        }
    }

    if (!use_manual_init && !is.null(cytoband_subset) && nrow(cytoband_subset) > 0) {
        circlize::circos.genomicIdeogram(cytoband = cytoband_subset)
    } else if (!use_manual_init) {
        circlize::circos.genomicIdeogram(species = genome)
    }

    if (!is.null(arc_data)) {
        .log_step("Adding DMR arc track...", level = 2)
        delta_beta <- as.numeric(arc_data$delta_beta)
        valid_vals <- delta_beta[is.finite(delta_beta)]
        if (length(valid_vals) == 0) {
            .log_warn("DMR delta_beta values are all missing/non-finite. Skipping arc track.")
            arc_data <- NULL
        } else {
            neg_min <- suppressWarnings(min(valid_vals[valid_vals < 0], na.rm = TRUE))
            pos_max <- suppressWarnings(max(valid_vals[valid_vals > 0], na.rm = TRUE))
            if (!is.finite(neg_min)) {
                neg_min <- -max(1e-6, abs(pos_max) * 0.05)
            }
            if (!is.finite(pos_max)) {
                pos_max <- max(1e-6, abs(neg_min) * 0.05)
            }
            if (neg_min >= 0) {
                neg_min <- -max(1e-6, abs(pos_max) * 0.05)
            }
            if (pos_max <= 0) {
                pos_max <- max(1e-6, abs(neg_min) * 0.05)
            }
            q <- c(neg_min, 0, pos_max)
            col_fun <- circlize::colorRamp2(q, c("#055709", "white", "#801414"))
        }

        if (!is.null(arc_data)) {
            arc_colors <- col_fun(arc_data$delta_beta)
            arc_df <- data.frame(
                chr = arc_data$chr,
                start = arc_data$start,
                end = arc_data$end,
                delta_beta = arc_data$delta_beta,
                color = arc_colors,
                stringsAsFactors = FALSE
            )

            circlize::circos.trackPlotRegion(
                sectors = unique_chrs,
                bg.border = NA,
                bg.col = "white",
                track.height = 0.05,
                ylim = c(0, 1),
                panel.fun = function(x, y) {
                    chr <- circlize::CELL_META$sector.index
                    dmr_chr <- arc_df[arc_df$chr == chr, , drop = FALSE]
                    if (nrow(dmr_chr) > 0) {
                        for (i in seq_len(nrow(dmr_chr))) {
                            if (dmr_chr$end[i] - dmr_chr$start[i] < degenerate_resolution) {
                                circlize::circos.lines(
                                    x = c((dmr_chr$start[i] + dmr_chr$end[i]) / 2, (dmr_chr$start[i] + dmr_chr$end[i]) / 2),
                                    y = c(0, 1),
                                    col = dmr_chr$color[i],
                                    border = NA
                                )
                            } else {
                                circlize::circos.rect(
                                    xleft = dmr_chr$start[i],
                                    xright = dmr_chr$end[i],
                                    ybottom = 0,
                                    ytop = 1,
                                    col = dmr_chr$color[i],
                                    border = NA,
                                    density = NULL
                                )
                            }
                        }
                    }
                }
            )
            arc_legend <- ComplexHeatmap::Legend(
                title = "DMR \u0394\u03b2",
                at = signif(q, 2),
                col_fun = col_fun,
                title_position = "topleft",
                legend_height = grid::unit(4, "cm"),
                labels_gp = grid::gpar(fontsize = 8),
                title_gp = grid::gpar(fontsize = 10, fontface = "bold")
            )
            .log_success("Arc track added", level = 2)
        }
    }

    if (!is.null(link_data) && nrow(link_data) > 0) {
        .log_step("Adding link track...", level = 2)

        link_data$has_jaspar_match <- !is.na(link_data$has_jaspar_match) & link_data$has_jaspar_match
        comp_data <- link_data[
            !duplicated(link_data$component_id),
            c("component_id", "size", "consensus_seq", "jaspar_names", "jaspar_corr", "has_jaspar_match", "component_best_rank"),
            drop = FALSE
        ]

        component_colors <- setNames(
            rep(unmatched_interaction_color, nrow(comp_data)),
            as.character(comp_data$component_id)
        )
        matched_components <- comp_data[comp_data$has_jaspar_match, , drop = FALSE]
        if (nrow(matched_components) > 0) {
            rank_vec <- matched_components$component_best_rank
            rank_vec[!is.finite(rank_vec)] <- Inf
            matched_components <- matched_components[
                order(rank_vec, -matched_components$size, matched_components$component_id),
                ,
                drop = FALSE
            ]
            hues <- seq(15, 375, length.out = nrow(matched_components) + 1)
            matched_cols <- grDevices::hcl(h = hues, l = 62, c = 95)[seq_len(nrow(matched_components))]
            component_colors[as.character(matched_components$component_id)] <- matched_cols
        }

        link_data$colors <- vapply(seq_len(nrow(link_data)), function(i) {
            comp_id <- as.character(link_data$component_id[i])
            if (!isTRUE(link_data$has_jaspar_match[i])) {
                return(grDevices::adjustcolor(unmatched_interaction_color, alpha.f = 0.85))
            }
            base_color <- component_colors[comp_id]
            sim_val <- link_data$sim[i]
            if (!is.finite(sim_val)) {
                sim_val <- 0
            }
            sim_val <- pmax(0, pmin(1, sim_val))
            shade_factor <- 0.45 + 0.55 * sim_val
            base_rgb <- grDevices::col2rgb(base_color) / 255
            grDevices::rgb(
                red = base_rgb[1] * shade_factor,
                green = base_rgb[2] * shade_factor,
                blue = base_rgb[3] * shade_factor,
                alpha = 0.85
            )
        }, character(1))

        link_legend_labels <- NULL
        legend_components <- comp_data[comp_data$has_jaspar_match, , drop = FALSE]
        if (nrow(legend_components) > 0) {
            rank_vec <- legend_components$component_best_rank
            rank_vec[!is.finite(rank_vec)] <- Inf
            legend_components <- legend_components[order(rank_vec, legend_components$component_id), , drop = FALSE]
            link_legend_colors <- component_colors[as.character(legend_components$component_id)]
            link_legend_labels <- vapply(seq_len(nrow(legend_components)), function(i) {
                rank_prefix <- if (is.finite(legend_components$component_best_rank[i])) {
                    paste0("[rank=", as.integer(legend_components$component_best_rank[i]), "] ")
                } else {
                    ""
                }
                label <- paste0(rank_prefix, "[n=", legend_components$size[i], "] ", legend_components$consensus_seq[i])
                if (.hasNonEmptyString(legend_components$jaspar_names[i])) {
                    jas_names <- trimws(strsplit(legend_components$jaspar_names[i], ",", fixed = TRUE)[[1]])
                    jas_corr <- if (.hasNonEmptyString(legend_components$jaspar_corr[i])) {
                        trimws(strsplit(legend_components$jaspar_corr[i], ",", fixed = TRUE)[[1]])
                    } else {
                        rep("", length(jas_names))
                    }
                    n_show <- min(3L, length(jas_names))
                    for (j in seq_len(n_show)) {
                        corr_val <- suppressWarnings(as.numeric(jas_corr[j]))
                        corr_txt <- if (is.finite(corr_val)) paste0(" (", signif(corr_val, 3), ")") else ""
                        label <- paste0(label, " | ", jas_names[j], corr_txt)
                    }
                    if (length(jas_names) > n_show) {
                        label <- paste0(label, " ...")
                    }
                }
                .wrapCircosLegendLabel(label)
            }, character(1))
        }

        if (any(!link_data$has_jaspar_match)) {
            if (is.null(link_legend_labels)) {
                link_legend_labels <- c("DMR Interactions without JASPAR Matches")
                link_legend_colors <- unmatched_interaction_color
            } else if (any(!comp_data$has_jaspar_match)) {
                link_legend_labels <- c(link_legend_labels, "DMR Interactions without JASPAR Matches")
                link_legend_colors <- c(link_legend_colors, unmatched_interaction_color)
            }
        }
        if (!is.null(link_legend_labels)) {
            link_legend <- ComplexHeatmap::Legend(
                labels = link_legend_labels,
                legend_gp = grid::gpar(fill = link_legend_colors),
                title = "DMR Interactions",
                title_position = "topleft",
                title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                labels_gp = grid::gpar(fontsize = 8)
            )
        }

        has_directionality <- "rank" %in% colnames(S4Vectors::mcols(dmrs)) &&
            any(is.finite(suppressWarnings(as.numeric(S4Vectors::mcols(dmrs)$rank))))
        directional_flag <- if (has_directionality) 1 else 0
        for (i in seq_len(nrow(link_data))) {
            point1 <- c(link_data$start1[i], link_data$end1[i])
            point2 <- c(link_data$start2[i], link_data$end2[i])
            if ((point1[2] - point1[1]) < degenerate_resolution) {
                point1 <- mean(point1)
            }
            if ((point2[2] - point2[1]) < degenerate_resolution) {
                point2 <- mean(point2)
            }

            if (has_directionality) {
                circlize::circos.link(
                    sector.index1 = link_data$chr1[i],
                    point1 = point1,
                    sector.index2 = link_data$chr2[i],
                    point2 = point2,
                    col = link_data$colors[i],
                    directional = directional_flag,
                    arr.type = "triangle",
                    arr.length = 0.3,
                    arr.width = 0.05,
                    border = NA
                )
            } else {
                circlize::circos.link(
                    sector.index1 = link_data$chr1[i],
                    point1 = point1,
                    sector.index2 = link_data$chr2[i],
                    point2 = point2,
                    col = link_data$colors[i],
                    directional = 0,
                    border = NA
                )
            }
        }
        .log_success("Link track added", level = 2)
    }

    circlize::circos.clear()
    grid::upViewport(2)
    top_legends <- Filter(Negate(is.null), list(heatmap_legend, arc_legend))
    has_any_legend <- length(top_legends) > 0 || !is.null(link_legend)
    if (has_any_legend) {
        grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
        top_packed <- NULL
        if (length(top_legends) > 0) {
            top_args <- c(top_legends, list(
                direction = "horizontal",
                column_gap = grid::unit(6, "mm"),
                row_gap = grid::unit(6, "mm"),
                max_width = grid::unit(0.98, "npc")
            ))
            top_packed <- do.call(ComplexHeatmap::packLegend, top_args)
        }
        lgd_list <- top_packed
        if (!is.null(link_legend)) {
            if (is.null(top_packed)) {
                lgd_list <- link_legend
            } else {
                lgd_list <- ComplexHeatmap::packLegend(
                    top_packed,
                    link_legend,
                    direction = "vertical",
                    row_gap = grid::unit(3, "mm"),
                    max_height = grid::unit(0.98, "npc")
                )
            }
        }
        ComplexHeatmap::draw(
            lgd_list,
            x = grid::unit(0, "npc"),
            y = grid::unit(0.5, "npc"),
            just = c("left", "center")
        )
        grid::upViewport()
    }
    grid::upViewport()
    .log_success("Circos plot created successfully")
    if (!is.null(output_file)) {
        grDevices::dev.off()
    }

    invisible(NULL)
}

#' Plot DMR Circos Views Using Automatically Selected Regions
#'
#' @description Selects a small set of informative genomic windows from the input DMRs
#' and forwards them to [plotDMRsCircos()]. The default `blocks` mode prefers
#' localized ranked DMR blocks when `block_id` is available, while `quick` mode
#' falls back to top-ranked individual DMR windows. Both selectors are greedy and
#' avoid exhaustive optimization to keep region finding inexpensive.
#'
#' @param method Character. Automatic region selection mode: `"blocks"` (default),
#'   `"components"`, `"hybrid"`, or `"quick"`.
#' @param n_regions Integer. Target maximum number of regions to show (default: `6`).
#' @param region_flank_bp Numeric. Flank in base pairs added around each selected
#'   block or DMR before plotting (default: `1e6`).
#' @param max_regions_per_chr Integer or `NULL`. Per-chromosome cap on automatically
#'   selected regions (default: `2`).
#' @param min_inter_region_bp Numeric. Nearby candidate windows within this gap are
#'   merged instead of being shown as separate Circos sectors (default: `5e6`).
#' @inheritParams plotDMRsCircos
#' @param ... Additional arguments passed to [plotDMRsCircos()]. This is where
#'   plot-only settings such as `max_dmrs_per_chr`, `max_cpgs_per_dmr`,
#'   `max_num_samples_per_group`, `max_components`, `unmatched_interaction_color`,
#'   `legend_width_ratio`, `degenerate_resolution`, `output_file`, and `verbose`
#'   can be supplied. Arguments managed by automatic selection, including `region`,
#'   should not be passed through `...`.
#'
#' @return Invisibly returns the selected region data frame with columns `chr`,
#'   `start`, and `end`.
#'
#' @examples
#' \dontrun{
#' plotAutoDMRsCircos(dmrs, beta = "beta.txt", pheno = pheno_df)
#' plotAutoDMRsCircos(dmrs, beta = "beta.txt", pheno = pheno_df, method = "blocks")
#' plotAutoDMRsCircos(dmrs, beta = "beta.txt", pheno = pheno_df, method = "components")
#' plotAutoDMRsCircos(dmrs, beta = "beta.txt", pheno = pheno_df, method = "hybrid")
#' plotAutoDMRsCircos(dmrs, beta = "beta.txt", pheno = pheno_df, method = "quick", n_regions = 4)
#' }
#' @export
plotAutoDMRsCircos <- function(dmrs,
                               beta,
                               pheno,
                               method = c("blocks", "components", "hybrid", "quick"),
                               n_regions = 6,
                               region_flank_bp = 1e6,
                               max_regions_per_chr = 2,
                               min_inter_region_bp = 5e6,
                               genome = "hg38",
                               array = "450K",
                               sorted_locs = NULL,
                               components = NULL,
                               interactions = NULL,
                               min_similarity = 0.8,
                               motif_cpg_flank_size = 5,
                               min_component_size = 2,
                               chromosomes = NULL,
                               query_components_with_jaspar = FALSE,
                               ...) {
    call_arg_names <- names(as.list(sys.call()))
    if (!is.null(call_arg_names) && "region" %in% call_arg_names[nzchar(call_arg_names)]) {
        stop("plotAutoDMRsCircos() manages `region` directly. Pass only plotDMRsCircos()-specific arguments through `...`.")
    }
    dots <- list(...)
    dot_names <- names(dots)
    if (!is.null(dot_names) && "region" %in% dot_names[nzchar(dot_names)]) {
        stop("plotAutoDMRsCircos() manages `region` directly. Pass only plotDMRsCircos()-specific arguments through `...`.")
    }
    method <- match.arg(method)
    dmrs_gr <- convertToGRanges(dmrs, genome = genome)
    if (!is.null(chromosomes)) {
        dmrs_gr <- dmrs_gr[as.character(GenomicRanges::seqnames(dmrs_gr)) %in% chromosomes]
    }
    if (length(dmrs_gr) == 0) {
        stop("No DMRs remain after applying automatic-selection chromosome filters.")
    }
    interaction_state <- NULL
    if (method %in% c("components", "hybrid") || !is.null(components) || !is.null(interactions)) {
        interaction_state <- .computeCircosInteractionState(
            dmrs = dmrs_gr,
            genome = genome,
            array = array,
            sorted_locs = sorted_locs,
            min_similarity = min_similarity,
            motif_cpg_flank_size = motif_cpg_flank_size,
            min_component_size = min_component_size,
            query_components_with_jaspar = query_components_with_jaspar,
            components = components,
            interactions = interactions
        )
    }
    dmrs_for_plot <- if (!is.null(interaction_state$dmrs)) interaction_state$dmrs else dmrs_gr
    region_df <- .selectCircosRegions(
        dmrs = dmrs_for_plot,
        method = method,
        n_regions = n_regions,
        region_flank_bp = region_flank_bp,
        max_regions_per_chr = max_regions_per_chr,
        min_inter_region_bp = min_inter_region_bp,
        components = interaction_state$components
    )
    .log_info("Automatically selected Circos regions (", method, "): ", .formatCircosRegionSelection(region_df), level = 2)

    plotDMRsCircos(
        dmrs = dmrs_for_plot,
        beta = beta,
        pheno = pheno,
        genome = genome,
        array = array,
        sorted_locs = sorted_locs,
        components = interaction_state$components,
        interactions = interaction_state$interactions,
        min_similarity = min_similarity,
        motif_cpg_flank_size = motif_cpg_flank_size,
        min_component_size = min_component_size,
        chromosomes = chromosomes,
        region = region_df,
        query_components_with_jaspar = query_components_with_jaspar,
        ...
    )

    invisible(region_df)
}

.getCytobandData <- function(genome) {
    cache_dir <- tools::R_user_dir("DMRsegal", which = "cache")
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE)
    }

    cache_file <- file.path(cache_dir, paste0("cytoband_", genome, ".rds"))
    if (file.exists(cache_file)) {
        .log_info("Loading cached cytoband data for ", genome, level = 3)
        return(readRDS(cache_file))
    }

    .log_step("Downloading cytoband data from UCSC for ", genome, "...", level = 3)
    tryCatch(
        {
            url <- paste0("https://hgdownload.cse.ucsc.edu/goldenpath/", genome, "/database/cytoBandIdeo.txt.gz")
            temp_file <- tempfile(fileext = ".txt.gz")
            result <- tryCatch(
                {
                    utils::download.file(url, temp_file, quiet = TRUE, mode = "wb")
                    TRUE
                },
                error = function(e) {
                    .log_warn("cytoBandIdeo not found, trying cytoBand table...")
                    FALSE
                }
            )
            if (!result) {
                url <- paste0("https://hgdownload.cse.ucsc.edu/goldenpath/", genome, "/database/cytoBand.txt.gz")
                utils::download.file(url, temp_file, quiet = TRUE, mode = "wb")
            }
            cytoband <- read.table(temp_file, sep = "\t", stringsAsFactors = FALSE, header = FALSE)
            if (ncol(cytoband) >= 5) {
                colnames(cytoband) <- c("V1", "V2", "V3", "V4", "V5")
            } else if (ncol(cytoband) == 4) {
                colnames(cytoband) <- c("V1", "V2", "V3", "V4")
                cytoband$V5 <- ""
            } else {
                stop("Unexpected cytoband format")
            }
            unlink(temp_file)
            saveRDS(cytoband, cache_file)
            .log_success("Cytoband data downloaded and cached", level = 3)
            cytoband
        },
        error = function(e) {
            .log_warn("Failed to download cytoband data: ", e$message, ". Using default ideogram.")
            NULL
        }
    )
}

.closest_rows_indices_to_centroids <- function(beta_pcs, centers, verbose = FALSE) { # nolint
    stopifnot(nrow(beta_pcs) >= nrow(centers))
    selection <- integer()
    assessed_pcs <- beta_pcs
    assessed_centers <- centers
    assessed_indices <- seq_len(nrow(beta_pcs))
    while (length(selection) < nrow(centers)) {
        nn <- FNN::get.knnx(assessed_pcs, assessed_centers, k = 1)
        center_indices <- as.vector(nn$nn.index)
        duplicated_centers <- duplicated(center_indices)
        chosen_indices <- center_indices[!duplicated_centers]
        selection <- c(selection, assessed_indices[chosen_indices])
        if (any(duplicated_centers)) {
            .log_info(sum(duplicated_centers), " duplicated nearest neighbours, retrying on the duplicates...", level = 3)
            assessed_pcs <- assessed_pcs[-chosen_indices, , drop = FALSE]
            assessed_indices <- assessed_indices[-chosen_indices]
            assessed_centers <- assessed_centers[duplicated_centers, , drop = FALSE]
        }
    }
    sort(selection)
}

.prepareCircosHeatmapData <- function(dmrs, beta_handler, pheno, sample_group_col, max_sup_cpgs_per_dmr_side = 2, max_num_samples_per_group = 10) {
    beta_col_names <- beta_handler$getBetaColNames()
    pheno <- pheno[rownames(pheno) %in% beta_col_names, , drop = FALSE]
    if (nrow(pheno) == 0) {
        .log_warn("No samples in pheno match the samples in beta values. Skipping heatmap track.")
        return(list(heatmap_df = NULL, reduced_pheno = NULL))
    }
    pheno <- pheno[order(pheno[[sample_group_col]]), , drop = FALSE]
    dmrs_cpgs_list <- strsplit(mcols(dmrs)$cpgs, split = ",")
    dmrs_cpgs_inds <- unique(unlist(dmrs_cpgs_list, use.names = FALSE))
    if (length(dmrs_cpgs_inds) == 0) {
        .log_warn("No supporting CpGs available for selected DMRs. Skipping heatmap track.")
        return(list(heatmap_df = NULL, reduced_pheno = NULL))
    }

    shown_locs <- beta_handler$getBetaLocs()[dmrs_cpgs_inds, c("chr", "start", "end"), drop = FALSE]
    shown_locs <- as.data.frame(shown_locs)
    shown_locs$chr <- as.character(shown_locs$chr)
    shown_locs$start <- as.numeric(shown_locs$start)
    shown_locs$end <- as.numeric(shown_locs$end)
    if (nrow(shown_locs) > 1) {
        chr_levels <- .orderChromosomesNaturally(shown_locs$chr)
        ord <- order(
            factor(shown_locs$chr, levels = chr_levels),
            shown_locs$start,
            shown_locs$end,
            rownames(shown_locs)
        )
        shown_locs <- shown_locs[ord, , drop = FALSE]
    }
    beta_data <- beta_handler$getBeta(row_names = rownames(shown_locs), col_names = rownames(pheno))
    if (max_num_samples_per_group > 0) {
        selected_samples <- c()
        groups <- pheno[[sample_group_col]]
        unique_groups <- unique(groups)
        reduced_pheno <- data.frame()
        for (group in unique_groups) {
            group_samples <- rownames(pheno)[groups == group]
            if (length(group_samples) > max_num_samples_per_group) {
                group_beta <- beta_data[, group_samples, drop = FALSE]
                non_zero_cols <- apply(group_beta, 2, function(v) stats::var(v, na.rm = TRUE) != 0)
                group_beta <- group_beta[, non_zero_cols, drop = FALSE]
                group_beta <- t(group_beta)
                pcs <- stats::prcomp(group_beta, center = TRUE, scale. = TRUE, rank. = min(10, ncol(group_beta)))$x
                set.seed(getOption("DMRsegal.random_seed", 42))
                kmeans <- stats::kmeans(
                    pcs,
                    centers = max_num_samples_per_group,
                    algorithm = "Lloyd",
                    iter.max = 1000,
                    nstart = 5
                )
                group_samples <- group_samples[.closest_rows_indices_to_centroids(pcs, kmeans$centers)]
            }
            reduced_pheno <- rbind(reduced_pheno, pheno[group_samples, , drop = FALSE])
            selected_samples <- c(selected_samples, group_samples)
        }
        beta_data <- beta_data[, selected_samples, drop = FALSE]
    } else {
        reduced_pheno <- pheno
    }
    beta_data <- as.matrix(beta_data)
    storage.mode(beta_data) <- "numeric"
    if (nrow(beta_data) > 0 && nrow(shown_locs) > 0 && !is.null(rownames(beta_data))) {
        beta_data <- beta_data[match(rownames(shown_locs), rownames(beta_data)), , drop = FALSE]
    }
    heatmap_df <- data.frame(shown_locs, beta_data, check.names = FALSE, stringsAsFactors = FALSE)
    list(heatmap_df = heatmap_df, reduced_pheno = reduced_pheno)
}

.prepareCircosArcData <- function(dmrs) {
    list(
        chr = as.character(GenomicRanges::seqnames(dmrs)),
        start = GenomicRanges::start(dmrs),
        end = GenomicRanges::end(dmrs),
        delta_beta = as.numeric(S4Vectors::mcols(dmrs)$delta_beta)
    )
}

.prepareCircosLinkData <- function(dmrs,
                                   genome,
                                   array,
                                   beta_locs,
                                   min_similarity,
                                   motif_cpg_flank_size,
                                   max_components,
                                   min_component_size,
                                   query_components_with_jaspar = TRUE,
                                   components = NULL,
                                   interactions = NULL) {
    if (is.null(components) || is.null(interactions)) {
        .log_info("Computing motif-based interactions for DMRs with min_similarity = ", min_similarity, "...")
        ret <- tryCatch(
            {
                computeDMRsInteraction(
                    dmrs = dmrs,
                    genome = genome,
                    array = array,
                    min_similarity = min_similarity,
                    beta_locs = beta_locs,
                    motif_cpg_flank_size = motif_cpg_flank_size,
                    min_component_size = min_component_size,
                    query_components_with_jaspar = query_components_with_jaspar
                )
            },
            error = function(e) {
                .log_warn("Failed to compute motif-based interactions: ", e$message)
                NULL
            }
        )
    } else {
        ret <- list(components = components, interactions = interactions)
        remapped <- .remapCircosInteractionsToDMRs(
            dmrs = dmrs,
            interactions = ret$interactions,
            components = ret$components
        )
        ret$interactions <- remapped$interactions
        ret$components <- remapped$components
    }

    if ((is.null(ret) || nrow(ret$interactions) == 0)) {
        .log_info("No significant interactions found at similarity >=", min_similarity, ". Skipping link track.", level = 2)
        return(NULL)
    }

    .log_success("Found ", nrow(ret$interactions), " motif-based interactions")
    if (!is.null(ret$components) && nrow(ret$components) > 0) {
        largest_component <- max(ret$components$size, na.rm = TRUE)
        if (is.finite(largest_component) && largest_component > 0.8 * length(dmrs)) {
            .log_warn(
                "Largest interaction component contains ", largest_component, " / ", length(dmrs),
                " DMRs. This often indicates dense generic CpG-context motif connectivity."
            )
        }
    }
    link_data <- ret$interactions
    if ("component_id" %in% colnames(link_data)) {
        link_data <- link_data[!is.na(link_data$component_id), , drop = FALSE]
    }
    if (nrow(link_data) == 0) {
        .log_warn("No motif interactions remain after component filtering. Skipping link track.")
        return(NULL)
    }

    if (!is.null(ret$components) && nrow(ret$components) > 0) {
        link_data <- merge(
            link_data,
            ret$components,
            by = "component_id",
            all.x = TRUE,
            suffixes = c("", "_comp"),
            sort = FALSE
        )
    }
    required_defaults <- list(
        size = NA_real_,
        consensus_seq = NA_character_,
        jaspar_names = NA_character_,
        jaspar_corr = NA_character_
    )
    for (nm in names(required_defaults)) {
        if (!(nm %in% colnames(link_data))) {
            link_data[[nm]] <- rep(required_defaults[[nm]], nrow(link_data))
        }
    }

    link_data$has_jaspar_match <- .hasNonEmptyString(link_data$jaspar_names)
    if ("rank" %in% colnames(S4Vectors::mcols(dmrs))) {
        ranks <- suppressWarnings(as.numeric(S4Vectors::mcols(dmrs)$rank))
        pair_rank <- pmin(ranks[link_data$index1], ranks[link_data$index2], na.rm = TRUE)
        pair_rank[!is.finite(pair_rank)] <- NA_real_
        link_data$interaction_rank <- pair_rank
        comp_rank_df <- stats::aggregate(
            interaction_rank ~ component_id,
            data = link_data,
            FUN = function(x) {
                x <- x[is.finite(x)]
                if (length(x) == 0) {
                    NA_real_
                } else {
                    min(x)
                }
            }
        )
        colnames(comp_rank_df)[colnames(comp_rank_df) == "interaction_rank"] <- "component_best_rank"
        link_data <- merge(link_data, comp_rank_df, by = "component_id", all.x = TRUE, sort = FALSE)
    } else {
        link_data$component_best_rank <- NA_real_
    }

    original_n <- nrow(link_data)
    link_data <- .selectCircosInteractions(link_data, max_components)
    if (!is.null(max_components) && nrow(link_data) < original_n) {
        .log_info("Selected ", nrow(link_data), " interactions (max_components=", max_components, ").", level = 2)
    }
    if (nrow(link_data) > 0) {
        matched_n <- sum(link_data$has_jaspar_match, na.rm = TRUE)
        .log_info(
            "Selected interactions: ",
            matched_n, " JASPAR-matched + ",
            nrow(link_data) - matched_n, " unmatched.",
            level = 2
        )
    }
    link_data <- .orderCircosInteractionRows(link_data)
    if (nrow(link_data) == 0) {
        .log_warn("No motif interactions remain after filtering. Skipping link track.")
        return(NULL)
    }
    rownames(link_data) <- NULL
    link_data
}
