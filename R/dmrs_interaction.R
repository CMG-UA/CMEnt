.maxPWMCorr <- function(pwm1, pwm2) {
    widthMin <- min(ncol(pwm1), ncol(pwm2))
    n1 <- ncol(pwm1) - widthMin + 1L
    n2 <- ncol(pwm2) - widthMin + 1L
    pwm1_center <- pwm1 - 0.25
    pwm2_center <- pwm2 - 0.25
    pwm1_center_sq <- pwm1_center^2
    pwm2_center_sq <- pwm2_center^2
    best <- -Inf
    for (i in seq_len(n1)) {
        i_idx <- i:(i + widthMin - 1L)
        pwm1Temp <- pwm1_center[, i_idx, drop = FALSE]
        pwm1TempNorm <- colSums(pwm1_center_sq[, i_idx, drop = FALSE])
        for (j in seq_len(n2)) {
            j_idx <- j:(j + widthMin - 1L)
            pwm2Temp <- pwm2_center[, j_idx, drop = FALSE]
            top <- colSums(pwm1Temp * pwm2Temp)
            bottom <- sqrt(pwm1TempNorm * colSums(pwm2_center_sq[, j_idx, drop = FALSE]))
            curr <- mean(top / (bottom + 1e-10))
            if (curr > best) {
                best <- curr
                if (best >= 0.999999) {
                    return(best)
                }
            }
        }
    }
    best
}

.motifCorr <- function(pwmSubjectMatrixList, pwmQuery, pwmQueryRevcomp = NULL) {
    sapply(pwmSubjectMatrixList, function(pwmSubject) {
        best <- .maxPWMCorr(pwmSubject, pwmQuery)
        if (!is.null(pwmQueryRevcomp)) {
            best <- max(best, .maxPWMCorr(pwmSubject, pwmQueryRevcomp))
        }
        best
    })
}

.parseDMRComponentIndices <- function(x) {
    if (is.null(x) || length(x) == 0) {
        return(integer())
    }
    values <- if (is.list(x) && !is.data.frame(x)) {
        unlist(x, recursive = TRUE, use.names = FALSE)
    } else {
        x
    }
    if (is.character(values)) {
        values <- trimws(unlist(base::strsplit(values, ",", fixed = TRUE), use.names = FALSE))
        values <- values[nzchar(values)]
    }
    idxs <- suppressWarnings(as.integer(values))
    idxs[!is.na(idxs)]
}

.serializeDMRInteractionComponentsForStorage <- function(components_df) {
    if (!is.data.frame(components_df) || nrow(components_df) == 0) {
        return(components_df)
    }
    ret <- components_df
    if ("avg_pwm" %in% colnames(ret)) {
        ret$avg_pwm <- NULL
    }
    if ("indices" %in% colnames(ret)) {
        ret$indices <- vapply(ret$indices, function(x) {
            idxs <- .parseDMRComponentIndices(x)
            if (length(idxs) == 0) {
                NA_character_
            } else {
                paste(idxs, collapse = ",")
            }
        }, character(1))
    }
    list_cols <- vapply(ret, is.list, logical(1))
    for (col in names(ret)[list_cols]) {
        ret[[col]] <- vapply(ret[[col]], function(x) {
            if (is.null(x) || length(x) == 0) {
                NA_character_
            } else {
                paste(as.character(unlist(x, recursive = TRUE, use.names = FALSE)), collapse = ",")
            }
        }, character(1))
    }
    ret
}

.writeTabularOutputAtomic <- function(x, path, sep = "\t", quote = FALSE, row.names = FALSE) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    tmp_path <- tempfile(
        pattern = paste0(".", basename(path), "-"),
        tmpdir = dirname(path)
    )
    on.exit({
        if (file.exists(tmp_path)) {
            unlink(tmp_path)
        }
    }, add = TRUE)

    write.table(
        x,
        tmp_path,
        sep = sep,
        quote = quote,
        row.names = row.names
    )

    if (!file.rename(tmp_path, path)) {
        copied <- file.copy(tmp_path, path, overwrite = TRUE, copy.mode = TRUE)
        if (!copied) {
            stop("Failed to move temporary output into place: ", path, call. = FALSE)
        }
        unlink(tmp_path)
    }

    invisible(path)
}

comparePWMToJaspar <- function(pwm_queries) {
    cache <- getOption(
        "DMRsegal.jaspar_cache_dir",
        .getOSCacheDir(file.path("R", "DMRsegal", "jaspar_cache"))
    )
    tax_group <- getOption("DMRsegal.jaspar_tax_group", "vertebrates")
    jaspar_version <- getOption("DMRsegal.jaspar_version", 2024)
    corr_threshold <- getOption("DMRsegal.jaspar_corr_threshold", 0.9)
    dir.create(cache, showWarnings = FALSE, recursive = TRUE)
    pwms_file <- file.path(cache, paste0("jaspar", jaspar_version, "_", tax_group, "_pwms.rds"))
    if (file.exists(pwms_file)) {
        .log_info("Loading JASPAR PWMs from cache...", level = 3)
        jaspar_pwms <- readRDS(pwms_file)
    } else {
        .log_info("Downloading JASPAR PWMs...", level = 2)
        jaspar_pkg <- paste0("JASPAR", jaspar_version)
        if (!requireNamespace(jaspar_pkg, quietly = TRUE)) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager")
            }
            BiocManager::install(jaspar_pkg, ask = FALSE, update = FALSE)
        }
        db <- getExportedValue(jaspar_pkg, jaspar_pkg)()@db
        opts <- list()
        opts[["tax_group"]] <- tax_group
        vertebrate_pfms <- TFBSTools::getMatrixSet(db, opts)
        vertebrate_pwms <- TFBSTools::toPWM(vertebrate_pfms, type = "prob")
        saveRDS(vertebrate_pwms, pwms_file)
        jaspar_pwms <- vertebrate_pwms
    }
    jaspar_pwm_mats <- lapply(jaspar_pwms, function(pwm) pwm@profileMatrix)
    complimentary_bases <- c("A" = "T", "C" = "G", "G" = "C", "T" = "A")
    revcomp_queries <- lapply(pwm_queries, function(pwm) {
        pwm_revcomp <- pwm[complimentary_bases[rownames(pwm)], rev(seq_len(ncol(pwm))), drop = FALSE]
        rownames(pwm_revcomp) <- rownames(pwm)
        pwm_revcomp
    })
    n_subject <- length(jaspar_pwm_mats)
    similarities <- vapply(pwm_queries, function(query) {
        .motifCorr(jaspar_pwm_mats, query)
    }, numeric(n_subject))
    revcomp_similarities <- vapply(revcomp_queries, function(query) {
        .motifCorr(jaspar_pwm_mats, query)
    }, numeric(n_subject))
    similarities <- pmax(similarities, revcomp_similarities)
    jaspar_names <- vapply(jaspar_pwms, function(x) x@name, character(1))

    found <- which(similarities >= corr_threshold, arr.ind = TRUE)
    if (length(found) == 0) {
        .log_info("No similar motifs found in JASPAR database with correlation >=", corr_threshold, level = 2)
        similarities <-  c()
    } else {
        similar_motifs <- data.frame(
            query_index = found[, 2],
            jaspar_id = names(jaspar_pwms)[found[, 1]],
            jaspar_name = jaspar_names[found[, 1]],
            correlation = similarities[found]
        )
        # group by query_index and return comma-separated jaspar_ids and correlations
        similarities <- lapply(split(similar_motifs, similar_motifs$query_index, drop = FALSE), function(df) {
            if (is.data.frame(df)) {
                df <- df[order(df$correlation, decreasing = TRUE), ]
            }
            list(
                jaspar_names = paste(df$jaspar_name, collapse = ","),
                jaspar_ids = paste(df$jaspar_id, collapse = ","),
                jaspar_corr = paste(round(df$correlation, 3), collapse = ",")
            )
        })
    }
    # Fill in results for all queries
    results <- lapply(seq_along(pwm_queries), function(i) {
        if (as.character(i) %in% names(similarities)) {
            similarities[[as.character(i)]]
        } else {
            list(
                jaspar_names = NA,
                jaspar_ids = NA,
                jaspar_corr = NA
            )
        }
    })
    results <- do.call(rbind, lapply(results, as.data.frame))

    results
}


getBackgroundArrayMotif <- function(genome, array, motif_cpg_flank_size = 5, .sorted_locs = NULL) {
    cache_dir <- getOption("DMRsegal.annotation_cache_dir",
        .getOSCacheDir(file.path("R", "DMRsegal", "annotations"))
    )
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    cache_file <- file.path(cache_dir, paste0("bgpwm_", genome, "_", array, "_", motif_cpg_flank_size, ".rds"))
    if (!file.exists(cache_file)) {
        .log_info("Background array motif pwm not existing in cache, computing it..", level = 2)
        if (is.null(.sorted_locs)) {
            sorted_locs <- getSortedGenomicLocs(array = array, genome = genome)
        } else {
            sorted_locs <- .sorted_locs
        }
        sorted_locs <- convertToGRanges(sorted_locs, genome)
        # Array annotations store CpGs as width-2 "CG" ranges, but motif windows in
        # this workflow are anchored on the CpG start base plus one downstream base.
        # Resize to width 1 so the extracted background windows match extractDMRMotifs().
        sorted_locs <- GenomicRanges::resize(sorted_locs, width = 1, fix = "start")
        cpg_seqs <- getDMRSequences(sorted_locs, genome, uflank_size = motif_cpg_flank_size, dflank_size = motif_cpg_flank_size + 1)
        expected_len <- 2 * motif_cpg_flank_size + 2
        valid_cpg_seqs <- !is.na(cpg_seqs) & nchar(cpg_seqs) == expected_len
        if (!all(valid_cpg_seqs)) {
            .log_warn(sum(!valid_cpg_seqs), " background motif windows had unexpected length and were ignored.")
        }
        cpg_seqs <- cpg_seqs[valid_cpg_seqs]
        if (length(cpg_seqs) == 0) {
            stop("Could not compute background motif PWM: no valid motif windows were extracted.")
        }
        cpg_seqs <- matrix(unlist(base::strsplit(cpg_seqs, split = "")), nrow = 2 * motif_cpg_flank_size + 2, byrow = FALSE)
        bg_frequencies <- as.matrix(apply(cpg_seqs, 1, function(x) table(factor(toupper(x), levels = Biostrings::DNA_BASES))))
        bg_pwm <- bg_frequencies / colSums(bg_frequencies) # row: position, column: base
        tryCatch(
            saveRDS(bg_pwm, cache_file),
            error = function(e) {
                .log_info("Could not cache background motif PWM: ", e$message, level = 3)
            }
        )
    } else {
        bg_pwm <- readRDS(cache_file)
    }
    bg_pwm
}


#' Extract DMR Motif Frequencies
#'
#' @description Extracts motif frequencies around CpG sites within DMRs.
#' For each DMR, retrieves sequences around the start and end CpG sites,
#' calculates base frequencies at each position, and stores the results in
#' the DMR metadata.
#' @param dmrs Dataframe or GRanges object containing DMR coordinates and CpG indices
#' @param genome Character. Genome version to use for sequence extraction. Defaults to hg38.
#' @param array Character. Array platform type (e.g., "450K", "EPIC"). Ignored if input is not array-based. (default: "450K")
#' @param beta_locs Data frame. Optional pre-computed genomic locations. If NULL,
#' locations will be retrieved using getSortedGenomicLocs (default: NULL)
#' @param motif_cpg_flank_size Integer. Number of base pairs to include as flanking regions around each CpG site (default: 5)
#' @return The input Dataframe/GRanges object with an additional metadata column:
#' \itemize{
#'   \item pwm: A matrix of base frequencies (rows: positions relative to CpG, columns: bases A, C, G, T)
#'   \item consensus_seq: A character string representing the consensus sequence derived from the PWM
#' }
#' @examples
#' # Extract motif frequencies for DMRs
#' dmrs <- data.frame(
#'     chr = c("chr16", "chr3"),
#'     start = c(53468112, 37459206),
#'     end = c(53468712, 37493431),
#'     start_cpg = c("cg00000029", "cg00000108"),
#'     start_seed = c("cg00000029", "cg00000108"),
#'     end_cpg = c("cg13426503", "cg08730726"),
#'     end_seed = c("cg13426503", "cg08730726"),
#'     seeds = c("cg00000029,cg13426503", "cg00000108,cg08730726")
#' )
#' dmrs_with_motifs <- extractDMRMotifs(dmrs, genome = "hg38", array = "450K")
#' # Access motif frequencies for the first DMR
#' motif_freqs_dmr1 <- dmrs_with_motifs$pwm[[1]]
#' @export
extractDMRMotifs <- function(
    dmrs, genome="hg38", array = "450k", beta_locs = NULL, motif_cpg_flank_size = 5, plot_dir = NULL
) {
    input_is_df <- is.data.frame(dmrs)
    dmrs <- convertToGRanges(dmrs, genome)
    if (!is.null(array)) {
        if (length(array) > 1) {
            array <- array[[1]]
        }
        array <- strex::match_arg(
            array,
            choices = c("450K", "27K", "EPIC", "EPICv2", "Mouse"),
            ignore_case = TRUE
        )
    }
    if (is.null(beta_locs) || (is.character(beta_locs) && length(beta_locs) == 1 && file.exists(beta_locs))) {
        beta_locs <- getSortedGenomicLocs(array = array, genome = genome, locations_file = beta_locs)
    }
    array_based <- !is.null(array)
    if (array_based) {
        bg_pwm <- getBackgroundArrayMotif(genome, array, motif_cpg_flank_size = motif_cpg_flank_size)
    }
    sequences <- getDMRSequences(
        dmrs, genome, uflank_size = motif_cpg_flank_size, dflank_size = motif_cpg_flank_size + 1
    )
    dmrs_seeds <- base::strsplit(as.character(mcols(dmrs)[, "seeds"]), split = ",", fixed = TRUE)
    all_seeds <- unique(unlist(dmrs_seeds, use.names = FALSE))
    all_seeds <- all_seeds[nzchar(all_seeds)]
    beta_locs_start <- as.integer(beta_locs[all_seeds, "start", drop = TRUE])
    names(beta_locs_start) <- all_seeds
    expected_len <- 2 * motif_cpg_flank_size + 2
    pwms <- vector("list", length(dmrs))
    consensus_seq <- rep(NA_character_, length(dmrs))
    for (i in seq_along(dmrs)) {
        dmr_seeds <- dmrs_seeds[[i]]
        dmr_seeds <- dmr_seeds[nzchar(dmr_seeds)]
        if (length(dmr_seeds) == 0) {
            next
        }

        start_locs <- beta_locs_start[dmr_seeds]
        valid_start_locs <- !is.na(start_locs)
        if (!all(valid_start_locs)) {
            .log_warn(sum(!valid_start_locs), " seed(s) were missing genomic locations and were ignored in DMR motif extraction.")
            start_locs <- start_locs[valid_start_locs]
        }
        if (length(start_locs) == 0) {
            next
        }

        sequence <- sequences[[i]]
        start_loc_base <- start_locs[[1]]
        seq_cpg_inds <- start_locs - start_loc_base + 1 + motif_cpg_flank_size
        cpg_seqs <- substring(sequence, seq_cpg_inds - motif_cpg_flank_size, seq_cpg_inds + motif_cpg_flank_size + 1)
        valid_cpg_seqs <- !is.na(cpg_seqs) & nchar(cpg_seqs) == expected_len
        if (!all(valid_cpg_seqs)) {
            .log_warn(sum(!valid_cpg_seqs), " motif windows had unexpected length and were ignored for one DMR.")
            cpg_seqs <- cpg_seqs[valid_cpg_seqs]
        }
        if (length(cpg_seqs) == 0) {
            next
        }
        # Apply transpose to get each sequence as a column, and then calculate base frequencies per row
        cpg_seqs <- matrix(unlist(base::strsplit(cpg_seqs, split = "")), nrow = 2 * motif_cpg_flank_size + 2, byrow = FALSE)
        frequencies <- as.matrix(apply(cpg_seqs, 1, function(x) table(factor(toupper(x), levels = Biostrings::DNA_BASES)))) # nolint
        if (array_based) {
            frequencies <- frequencies * (1 / max(bg_pwm, 1e-7))
        }
        pwms[[i]] <- frequencies / colSums(frequencies) # row: position, column: base
        consensus_seq[[i]] <- paste(Biostrings::DNA_BASES[apply(frequencies, 2, which.max)], collapse = "")
    }
    mcols(dmrs)$pwm <- pwms
    mcols(dmrs)$consensus_seq <- consensus_seq
    if (input_is_df) {
        dmrs <- convertToDataFrame(dmrs)
    }
    invisible(dmrs)
}


.extractMotifsSimilarity <- function(dmrs, motif_cpg_flank_size = 5) {
    if (inherits(dmrs, "GRanges")) {
        pwms <- mcols(dmrs)$pwm
    } else {
        pwms <- dmrs[, "pwm"]
    }
    if (length(pwms) == 0) {
        return(matrix(0, nrow = 0, ncol = 0))
    }
    cpg_cols <- c(motif_cpg_flank_size + 1, motif_cpg_flank_size + 2)
    # Compare motifs with a position-aware centered cosine:
    # 1) remove central CpG positions,
    # 2) center base frequencies around 0.25,
    # 3) normalize each position independently,
    # 4) average similarity over positions.
    pwm_vectors <- lapply(pwms, function(x) {
        keep_cols <- setdiff(seq_len(ncol(x)), cpg_cols)
        x <- x[, keep_cols, drop = FALSE] - 0.25
        col_norms <- sqrt(colSums(x^2))
        col_norms[col_norms < 1e-10] <- 1e-10
        x <- sweep(x, 2, col_norms, "/")
        as.vector(x)
    })
    pwm_matrix <- do.call(cbind, pwm_vectors)
    num_positions <- max(1, ncol(pwms[[1]]) - 2)
    ret <- crossprod(pwm_matrix) / num_positions
    ret[ret > 1] <- 1
    ret[ret < -1] <- -1
    ret
}


#' Compute Motif-Based DMR Interactions
#'
#' @description Computes motif-based interactions between DMRs based on their
#' motif similarity. Identifies pairs of DMRs with significant motif similarity
#' and returns a data frame of interactions. Assigns directionality based on score, if available.
#' @param dmrs Dataframe or GRanges object containing DMR coordinates and motif information
#' @param genome Character. Genome version to use for sequence extraction (e.g., "hg38")
#' @param array Character. Array platform type (e.g., "450K", "EPIC"). Must be NULL if input is not array-based (default: "450K")
#' @param min_similarity Numeric. Minimum motifs PWM similarity threshold for considering DMRs are related (default: 0.8)
#' @param beta_locs Data frame. Optional pre-computed genomic locations. If NULL,
#' locations will be retrieved using getSortedGenomicLocs (default: NULL)
#' @param motif_cpg_flank_size Integer. Number of base pairs to include as flanking regions around each CpG site (default: 5)
#' @param find_components Logical. Whether to identify connected components of interacting DMRs (default: TRUE)
#' @param min_component_size Integer. Minimum size of connected components to consider (default: 2)
#' @param query_components_with_jaspar Logical. Whether to query connected components average PWMs against JASPAR database (default: TRUE)
#' @param output_prefix Character. Prefix for output files to save interactions and components (optional). If NULL, results are not saved to file (default: NULL)
#' @param plot_dir Character. Directory to save diagnostic plots (optional). If NULL, no plots are saved (default: NULL)
#' @return A list with:
#' \itemize{
#'   \item interactions: Data frame of motif-based DMR interactions with columns
#'   \code{chr1}, \code{start1}, \code{end1}, \code{chr2}, \code{start2},
#'   \code{end2}, \code{sim}, and when components are requested,
#'   \code{component_id}
#'   \item components: Data frame of discovered motif components
#'   \item dmrs: The input DMRs in the same class as provided, with an added
#'   \code{component_ids} column containing comma-separated component IDs for
#'   each DMR (or \code{NA} when the DMR is not part of a retained component)
#' }
#' @examples
#' # Compute motif-based interactions for DMRs
#' dmrs <- data.frame(
#'     chr = c("chr16", "chr3"),
#'     start = c(53468112, 37459206),
#'     end = c(53468712, 37493431),
#'     start_cpg = c("cg00000029", "cg00000108"),
#'     end_cpg = c("cg13426503", "cg08730726")
#' )
#' dmrs_with_motifs <- extractDMRMotifs(dmrs, genome = "hg38", array = "450K")
#' interactions <- computeDMRsInteraction(
#'     dmrs_with_motifs,
#'     genome = "hg38",
#'     array = "450K",
#' )
#' @export
computeDMRsInteraction <- function(
    dmrs,
    genome = "hg38",
    array = "450K",
    min_similarity = getOption("DMRsegal.min_motif_similarity", 0.8),
    beta_locs = NULL,
    motif_cpg_flank_size = 5,
    find_components = TRUE,
    min_component_size = 2,
    query_components_with_jaspar = TRUE,
    plot_dir = NULL,
    output_prefix = NULL
) {
    input_is_df <- is.data.frame(dmrs)
    dmrs <- convertToGRanges(dmrs, genome)
    mcols(dmrs)$component_ids <- rep(NA_character_, length(dmrs))
    if (length(dmrs) == 0) {
        .log_info("No DMRs provided for interaction analysis.", level = 2)
        return(list(
            interactions = data.frame(),
            components = data.frame(),
            dmrs = if (input_is_df) convertToDataFrame(dmrs) else dmrs
        ))
    }
    if (!"pwm" %in% colnames(mcols(dmrs))) {
        .log_info("DMR motifs not precomputed. Extracting motifs...", level = 2)
        dmrs <- extractDMRMotifs(dmrs, genome, array, beta_locs = beta_locs, motif_cpg_flank_size = motif_cpg_flank_size)
    }
    if (length(dmrs) == 1) {
        .log_info("Only one DMR provided, skipping interaction analysis.", level = 2)
        return(list(
            interactions = data.frame(),
            components = data.frame(),
            dmrs = if (input_is_df) convertToDataFrame(dmrs) else dmrs
        ))
    }
    similarity_matrix <- .extractMotifsSimilarity(dmrs, motif_cpg_flank_size = motif_cpg_flank_size)
    if (!is.null(plot_dir)) {
        dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
        to_show <- similarity_matrix
        diag(to_show) <- NA
        to_show[lower.tri(to_show)] <- NA
        to_show[to_show < min_similarity] <- NA
        sim_melt <- reshape2::melt(to_show, na.rm = TRUE)
        colnames(sim_melt) <- c("DMR1", "DMR2", "Similarity")
        p <- ggplot2::ggplot(sim_melt, ggplot2::aes_string(x = "DMR1", y = "DMR2", fill = "Similarity")) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_viridis_c() +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
            ggplot2::labs(title = "DMR Motif Similarity Heatmap", x = "DMR Index", y = "DMR Index")
        ggplot2::ggsave(filename = file.path(plot_dir, "dmr_motif_similarity_heatmap.png"), plot = p, width = 8, height = 6)
    }
    mask <- !is.na(similarity_matrix) & (similarity_matrix >= min_similarity)
    # remove diagonal
    diag(mask) <- FALSE
    interactions_df <- data.frame()
    components_df <- data.frame()
    components <- NULL
    membership_to_component <- rep(NA_integer_, length(dmrs))
    if ((!find_components || min_component_size > 1) && !any(mask)) {
        .log_info("No motif similarities found above the threshold.", level = 2)
        return(list(
            interactions = interactions_df,
            components = components_df,
            dmrs = if (input_is_df) convertToDataFrame(dmrs) else dmrs
        ))
    }
    has_score <- inherits(dmrs, "GRanges") && "score" %in% colnames(mcols(dmrs))
    if (any(mask)) {
        if (has_score) {
            rowcol_df <- which(mask, arr.ind = TRUE)
            scores <- mcols(dmrs)$score
            keep <- scores[rowcol_df[, 1]] >= scores[rowcol_df[, 2]]
            rowcol_df <- rowcol_df[keep, , drop = FALSE]
            oriented_mask <- matrix(FALSE, nrow = nrow(mask), ncol = ncol(mask))
            oriented_mask[rowcol_df] <- TRUE
            mask <- oriented_mask
        } else {
            # remove lower triangle to avoid duplicate interactions
            mask[lower.tri(mask)] <- FALSE
            rowcol_df <- which(mask, arr.ind = TRUE)
        }
        start_dmrs <- dmrs[rowcol_df[, 1], ]
        end_dmrs <- dmrs[rowcol_df[, 2], ]
        interactions_df <- data.frame(
            index1 = rowcol_df[, 1],
            chr1 = as.character(GenomeInfoDb::seqnames(start_dmrs)),
            start1 = GenomicRanges::start(start_dmrs),
            end1 = GenomicRanges::end(start_dmrs),
            index2 = rowcol_df[, 2],
            chr2 = as.character(GenomeInfoDb::seqnames(end_dmrs)),
            start2 = GenomicRanges::start(end_dmrs),
            end2 = GenomicRanges::end(end_dmrs),
            sim = similarity_matrix[rowcol_df]
        )
        interactions_df <- interactions_df[is.finite(interactions_df$sim) & interactions_df$sim >= min_similarity, , drop = FALSE]
    }

    if (find_components) {
        # Components should represent undirected interaction connectivity, even when
        # interaction links are oriented by score
        component_mask <- !is.na(similarity_matrix) & (similarity_matrix >= min_similarity)
        diag(component_mask) <- FALSE
        component_mask <- component_mask | t(component_mask)
        g1 <- igraph::graph_from_adjacency_matrix(component_mask, mode = "undirected", diag = FALSE)
        components <- igraph::components(g1)
        keep_component_ids <- which(components$csize >= min_component_size)

        if (length(keep_component_ids) == 0) {
            .log_info("No connected components found with size >=", min_component_size, level = 2)
        } else {
            membership_to_component <- match(components$membership, keep_component_ids)
            components_df <- data.frame(
                component_id = seq_along(keep_component_ids),
                size = as.numeric(components$csize[keep_component_ids])
            )
            components_df$indices <- lapply(keep_component_ids, function(i) {
                which(components$membership == i)
            })
            # Find the average PWM for each component
            components_df$avg_pwm <- lapply(components_df$indices, function(idxs) {
                pwms <- mcols(dmrs)[idxs, "pwm"]
                mat <- Reduce("+", pwms) / length(pwms)
                mat / colSums(mat)
            })
            components_df$consensus_seq <- sapply(components_df$avg_pwm, function(pwm) {
                paste(Biostrings::DNA_BASES[apply(pwm, 2, which.max)], collapse = "")
            })
            # Order by component size
            components_df <- components_df[order(-components_df$size), ]
            old_component_ids <- components_df$component_id
            components_df$component_id <- seq_len(nrow(components_df))
            id_remap <- rep(NA_integer_, max(old_component_ids))
            id_remap[old_component_ids] <- components_df$component_id
            membership_to_component <- id_remap[membership_to_component]
            largest_component <- max(components_df$size)
            if (largest_component >= ceiling(0.8 * length(dmrs))) {
                .log_info(
                    "Largest motif component spans ",
                    largest_component, "/", length(dmrs),
                    " DMRs. This indicates broad motif similarity; consider a stricter min_similarity.",
                    level = 2
                )
            }
            if (query_components_with_jaspar) {
                # Find similarities to JASPAR motifs
                components_df <- cbind(components_df, comparePWMToJaspar(components_df$avg_pwm))
            }
        }
    }
    if (find_components && nrow(components_df) > 0) {
        component_ids_by_dmr <- vector("list", length(dmrs))
        for (i in seq_len(nrow(components_df))) {
            idxs <- components_df$indices[[i]]
            component_id <- as.character(components_df$component_id[[i]])
            for (idx in idxs) {
                component_ids_by_dmr[[idx]] <- c(component_ids_by_dmr[[idx]], component_id)
            }
        }
        mcols(dmrs)$component_ids <- vapply(component_ids_by_dmr, function(ids) {
            ids <- unique(ids)
            if (length(ids) == 0) {
                NA_character_
            } else {
                paste(ids, collapse = ",")
            }
        }, character(1))
    }
    if (find_components && nrow(interactions_df) > 0) {
        interactions_df$component_id <- membership_to_component[interactions_df$index1]
    }
    if (!is.null(output_prefix)) {
        .writeTabularOutputAtomic(
            interactions_df,
            paste0(output_prefix, ".dmr_interactions.tsv")
        )
        .writeTabularOutputAtomic(
            .serializeDMRInteractionComponentsForStorage(components_df),
            paste0(output_prefix, ".dmr_components.tsv")
        )
    }
    list(
        interactions = interactions_df,
        components = components_df,
        dmrs = if (input_is_df) convertToDataFrame(dmrs) else dmrs
    )
}
