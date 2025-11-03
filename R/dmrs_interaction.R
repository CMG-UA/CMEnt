comparePWMToJaspar <- function(pwm_queries, corr_threshold = 0.7) {
    cache <- getOption("DMRsegal.jaspar_cache_dir", file.path(
        path.expand("~"),
        ".cache", "R", "DMRsegal", "jaspar_cache"
    ))
    tax_group <- getOption("DMRsegal.jaspar_tax_group", "vertebrates")
    jaspar_version <- getOption("DMRsegal.jaspar_version", 2024)
    dir.create(cache, showWarnings = FALSE, recursive = TRUE)
    pwms_file <- file.path(cache, paste0("jaspar", jaspar_version, "_", tax_group, "_pwms.rds"))
    if (file.exists(pwms_file)) {
        .log_info("Loading JASPAR PWMs from cache...", level = 3)
        jaspar_pwms <- readRDS(pwms_file)
    } else {
        .log_info("Downloading JASPAR PWMs...", level = 2)
        jaspar_pkg <- paste0("JASPAR", jaspar_version)
        if (!requireNamespace(jaspar_pkg, quietly = TRUE)) {
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
    similarities <- sapply(pwm_queries, function(query) TFBSTools::PWMSimilarity(jaspar_pwms, query, method = "Pearson"))
    found <- which(similarities >= corr_threshold, arr.ind = TRUE)
    if (nrow(found) == 0) {
        .log_warn("No similar motifs found in JASPAR database with correlation >=", corr_threshold)
        similarities <- c()
    } else {
        similar_motifs <- data.frame(
            query_index = found[, 2],
            jaspar_id = names(jaspar_pwms)[found[, 1]],
            jaspar_name = sapply(found[, 1], function(x) jaspar_pwms[[x]]@name),
            correlation = similarities[found]
        )
        # group by query_index and return comma-separated jaspar_ids and correlations
        similarities <- lapply(split(similar_motifs, similar_motifs$query_index, drop=FALSE), function(df) {
            if (is.data.frame(df)) {
                df <- df[order(df$correlation, decreasing = TRUE),]
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

    return(results)
}


#' Extract DMR Motif Frequencies
#'
#' @description Extracts motif frequencies around CpG sites within DMRs.
#' For each DMR, retrieves sequences around the start and end CpG sites,
#' calculates base frequencies at each position, and stores the results in
#' the DMR metadata.
#' @param dmrs Dataframe or GRanges object containing DMR coordinates and CpG indices
#' @param genome Character. Genome version to use for sequence extraction (e.g., "hg19")
#' @param array Character. Array platform type (e.g., "450K", "EPIC") (default: "450K")
#' @param genomic_locs Data frame. Optional pre-computed genomic locations. If NULL,
#' locations will be retrieved using getSortedGenomicLocs (default: NULL)
#' @param flank_size Integer. Number of base pairs to include as flanking regions around each CpG site (default: 5)
#' @return The input Dataframe/GRanges object with an additional metadata column:
#' \itemize{
#'   \item pwm: A matrix of base frequencies (rows: positions relative to CpG, columns: bases A, C, G, T)
#'   \item consensus_sequence: A character string representing the consensus sequence derived from the PWM
#' }
#' @examples
#' # Extract motif frequencies for DMRs
#' dmrs <- data.frame(
#'     chr = c("chr16", "chr3"),
#'     start = c(53468112, 37459206),
#'     end = c(53468712, 37493431),
#'     start_cpg = c("cg00000029", "cg00000108"),
#'     start_dmp = c("cg00000029", "cg00000108"),
#'     end_cpg = c("cg13426503", "cg08730726"),
#'     end_dmp = c("cg13426503", "cg08730726"),
#'     dmps = c("cg00000029,cg13426503", "cg00000108,cg08730726")
#' )
#' dmrs_with_motifs <- extractDMRMotifs(dmrs, genome = "hg19", array = "450K")
#' # Access motif frequencies for the first DMR
#' motif_freqs_dmr1 <- dmrs_with_motifs$pwm[[1]]
#' @export
extractDMRMotifs <- function(dmrs, genome, array = "450k", genomic_locs = NULL, flank_size = 5, plot.dir = NULL) {
    input_is_df <- is.data.frame(dmrs)
    dmrs <- convertToGRanges(dmrs, genome)
    if (is.null(genomic_locs) || (is.character(genomic_locs) && length(genomic_locs) == 1 && file.exists(genomic_locs))) {
        genomic_locs <- getSortedGenomicLocs(array = array, genome = genome, locations_file = genomic_locs)
    }

    sequences <- getDMRSequences(dmrs, genome, uflank_size = flank_size, dflank_size = flank_size + 1)
    dmrs_cpgs_inds <- getSupportingSites(dmrs, rownames(genomic_locs), ret_index=TRUE, separate_by_section = FALSE)
    start_inds <- idToGenomicLocsIndex(mcols(dmrs)$start_cpg, genomic_locs)
    end_inds <- idToGenomicLocsIndex(mcols(dmrs)$end_cpg, genomic_locs)
    for (i in seq_along(dmrs)) {
        if (is.na(start_inds[i]) || is.na(end_inds[i])) {
            next
        }
        absolute_cpg_inds <- dmrs_cpgs_inds[[i]]
        sequence <- sequences[[i]]
        start_loc_base <- genomic_locs[start_inds[[i]], "start"]
        seq_cpg_inds <- genomic_locs[absolute_cpg_inds, "start"] - start_loc_base + 1 + flank_size
        cpg_seqs <- substring(sequence, seq_cpg_inds - flank_size, seq_cpg_inds + flank_size + 1)
        # Apply transpose to get each sequence as a column, and then calculate base frequencies per row
        cpg_seqs <- matrix(unlist(strsplit(cpg_seqs, split = "")), nrow = 2 * flank_size + 2, byrow = FALSE)
        frequencies <- as.matrix(apply(cpg_seqs, 1, function(x) table(factor(toupper(x), levels = DNA_BASES)))) # nolint
        mcols(dmrs)$pwm[[i]] <- frequencies / colSums(frequencies) # row: position, column: base
        mcols(dmrs)$consensus_sequence[[i]] <- paste(DNA_BASES[apply(frequencies, 2, which.max)], collapse = "")
    }
    if (input_is_df) {
        dmrs <- as.data.frame(dmrs)
        colnames(dmrs)[colnames(dmrs) == "seqnames"] <- "chr"
    }
    invisible(dmrs)
}


.extractMotifsSimilarity <- function(dmrs, flank_size = 5) {
    if (inherits(dmrs, "GRanges")) {
        pwms <- mcols(dmrs)$pwm
    } else {
        pwms <- dmrs[, "pwm"]
    }
    # Remove CpG position
    pwms <- do.call(cbind, lapply(pwms, function(x) unlist(as.list(x[, -c(flank_size + 1, flank_size + 2)]))))
    ret <- lsa::cosine(pwms)
    return(ret)
}


#' Compute Motif-Based DMR Interactions
#'
#' @description Computes motif-based interactions between DMRs based on their
#' motif similarity. Identifies pairs of DMRs with significant motif similarity
#' and returns a data frame of interactions.
#' @param dmrs Dataframe or GRanges object containing DMR coordinates and motif information
#' @param genome Character. Genome version to use for sequence extraction (e.g., "hg19")
#' @param array Character. Array platform type (e.g., "450K", "EPIC") (default: "450K")
#' @param min_sim Numeric. Minimum motifs PWM similarity threshold for considering DMRs are related (default: 0.7).
#' @param genomic_locs Data frame. Optional pre-computed genomic locations. If NULL,
#' locations will be retrieved using getSortedGenomicLocs (default: NULL)
#' @param flank_size Integer. Number of base pairs to include as flanking regions around each CpG site (default: 5)
#' @param find_components Logical. Whether to identify connected components of interacting DMRs (default: TRUE)
#' @param query_components_with_jaspar Logical. Whether to query connected components average PWMs against JASPAR database (default: TRUE)
#' @return Data frame of motif-based DMR interactions with columns:
#' \itemize{
#'   \item start_chr: Chromosome of the start DMR
#'   \item start: Start position of the start DMR
#'   \item end: End position of the start DMR
#'   \item end_chr: Chromosome of the end DMR
#'   \item end: Start position of the end DMR
#'   \item end_end: End position of the end DMR
#'   \item sim: Similarity value of the motif PWMs between the DMRs
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
#' dmrs_with_motifs <- extractDMRMotifs(dmrs, genome = "hg19", array = "450K")
#' interactions <- computeDMRsInteraction(
#'     dmrs_with_motifs,
#'     genome = "hg19",
#'     array = "450K",
#' )
#' @export
computeDMRsInteraction <- function(dmrs, genome = "hg19", array = "450K", min_sim = 0.7, genomic_locs = NULL, flank_size = 5, 
    find_components=TRUE, query_components_with_jaspar=TRUE, plot.dir = NULL) {
    dmrs <- convertToGRanges(dmrs, genome)
    if (!"pwm" %in% colnames(mcols(dmrs))) {
        .log_info("DMR motifs not precomputed. Extracting motifs...", level = 2)
        dmrs <- extractDMRMotifs(dmrs, genome, array, genomic_locs = genomic_locs, flank_size = flank_size)
    }
    similarity_matrix <- .extractMotifsSimilarity(dmrs, flank_size = flank_size)
    if (!is.null(plot.dir)) {
        dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)
        sim_melt <- reshape2::melt(similarity_matrix, na.rm = TRUE)
        colnames(sim_melt) <- c("DMR1", "DMR2", "Similarity")
        p <- ggplot2::ggplot(sim_melt, ggplot2::aes_string(x = "DMR1", y = "DMR2", fill = "Similarity")) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_viridis_c() +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
            ggplot2::labs(title = "DMR Motif Similarity Heatmap", x = "DMR Index", y = "DMR Index")
        ggplot2::ggsave(filename = file.path(plot.dir, "dmr_motif_similarity_heatmap.png"), plot = p, width = 8, height = 6)
    }
    mask <- !is.na(similarity_matrix) & (abs(similarity_matrix) >= min_sim)
    if (all(!mask, na.rm = TRUE)) {
        .log_warn("No motif-based interactions found between DMRs at similarity >=", min_sim)
        return(NULL)
    }
    components_df <- NULL
    if (find_components) {
        g1 <- igraph::graph_from_adjacency_matrix(mask)
        components <- igraph::components(g1)
        # compute consensus sequence for each connected component
        # create a dataframe with columns component_id, dmrs
        components_df <- data.frame(
            component_id = seq_along(components$csize),
            size = components$csize
        )
        components_df$indices <- lapply(seq_along(components$csize), function(i) {
            which(components$membership == i)
        })
        # Find the average PWM for each component
        components_df$avg_pwm <- lapply(components_df$indices, function(idxs) {
            pwms <- mcols(dmrs)[idxs, "pwm"]
            Reduce("+", pwms) / length(pwms)
        })
        components_df$consensus_sequence <- sapply(components_df$avg_pwm, function(pwm) {
            paste(DNA_BASES[apply(pwm, 2, which.max)], collapse = "")
        })
        # Order by component size
        components_df <- components_df[order(-components_df$size), ]
        if (query_components_with_jaspar){
            # Find similarities to JASPAR motifs
            components_df <- cbind(components_df, comparePWMToJaspar(components_df$avg_pwm, corr_threshold = min_sim))
        }
    }


    rowcol_df <- which(mask, arr.ind = TRUE)
    start_dmrs <- dmrs[rowcol_df[, 1], ]
    end_dmrs <- dmrs[rowcol_df[, 2], ]
    interaction_data_frame <- data.frame(
        chr1 = as.character(GenomeInfoDb::seqnames(start_dmrs)),
        start1 = GenomicRanges::start(start_dmrs),
        end1 = GenomicRanges::end(start_dmrs),
        chr2 = as.character(GenomeInfoDb::seqnames(end_dmrs)),
        start2 = GenomicRanges::start(end_dmrs),
        end2 = GenomicRanges::end(end_dmrs),
        sim = similarity_matrix[rowcol_df]
    )
    list(
        interactions = interaction_data_frame,
        components = components_df
    )
}
