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
#'   \item motif_frequencies: A list column where each entry is a matrix of base frequencies
#'     (rows: positions relative to CpG, columns: bases A, T, G, C)
#' }
#' @examples
#' # Extract motif frequencies for DMRs
#' dmrs <- data.frame(
#'     chr = c("chr16", "chr3"),
#'     start = c(53468112, 37459206),
#'     end = c(53468712, 37493431),
#'     start_cpg = c("cg00000029", "cg00000108"),
#'     end_cpg = c("cg13426503", "cg08730726")
#' )
#' dmrs_with_motifs <- extractDMRMotifs(dmrs, genome = "hg19", array = "450K")
#' # Access motif frequencies for the first DMR
#' motif_freqs_dmr1 <- dmrs_with_motifs$motif_frequencies[[1]]
#' @export
extractDMRMotifs <- function(dmrs, genome, array = "450k", genomic_locs = NULL, flank_size = 5) {
    input_is_df <- is.data.frame(dmrs)
    dmrs <- convert_to_granges(dmrs, genome)
    if (is.null(genomic_locs) || (is.character(genomic_locs) && length(genomic_locs) == 1 && file.exists(genomic_locs))) {
        genomic_locs <- getSortedGenomicLocs(array = array, genome = genome, locations_file = genomic_locs)
    }
    sequences <- getDMRSequences(dmrs, genome, uflank_size = flank_size, dflank_size = flank_size + 1)
    start_inds <- id_to_genomic_locs_index(mcols(dmrs)$start_cpg, genomic_locs)
    end_inds <- id_to_genomic_locs_index(mcols(dmrs)$end_cpg, genomic_locs)
    for (i in seq_along(dmrs)) {
        if (is.na(start_inds[i]) || is.na(end_inds[i])) {
            next
        }
        cpg_range <- start_inds[[i]]:end_inds[[i]]
        sequence <- sequences[[i]]
        start_loc_base <- genomic_locs[start_inds[[i]], "start"]
        seq_cpg_inds <- genomic_locs[cpg_range, "start"] - start_loc_base + 1 + flank_size
        cpg_seqs <- substring(sequence, seq_cpg_inds - flank_size, seq_cpg_inds + flank_size)
        # Apply transpose to get each sequence as a column, and then calculate base frequencies per row
        cpg_seqs <- matrix(unlist(strsplit(cpg_seqs, split = "")), nrow = 2 * flank_size + 1, byrow = FALSE)
        frequencies <- as.matrix(apply(cpg_seqs, 1, function(x) table(factor(toupper(x), levels = BASE_LEVELS)))) # nolint
        mcols(dmrs)$pwm[[i]] <- frequencies / colSums(frequencies) # row: position, column: base
        mcols(dmrs)$consensus_sequence[[i]] <- paste(BASE_LEVELS[apply(frequencies, 2, which.max)], collapse = "")
    }
    if (input_is_df) {
        dmrs <- as.data.frame(dmrs)
        colnames(dmrs)[colnames(dmrs) == "seqnames"] <- "chr"
    }
    invisible(dmrs)
}

# from https://stackoverflow.com/a/13486833/6758862
.corr_significance <- function(x, dfr = nrow(x) - 2, correction="BH") {
    corr_mat <- cor(x)
    above <- row(corr_mat) < col(corr_mat)
    r2 <- corr_mat[above]^2
    fstat <- r2 * dfr / (1 - r2)
    sig_mat <- matrix(NA, nrow = ncol(x), ncol = ncol(x))
    sig_mat[above] <- 1 - pf(fstat, 1, dfr)
    if (!is.null(correction)){
        sig_mat[above] <- p.adjust(sig_mat[above], method=correction)
    }
    sig_mat <- t(sig_mat)
    sig_mat[upper.tri(sig_mat)] <- NA
    list("corr_mat" = corr_mat, "sig_mat" = sig_mat)
}

.extractMotifsSimilaritySignificance <- function(dmrs, correction = "BH", flank_size = 5) {
    if (inherits(dmrs, "GRanges")) {
        pwms <- mcols(dmrs)$pwm
    } else {
        pwms <- dmrs[, "pwm"]
    }
    # Remove CpG position
    pwms <- do.call(cbind, lapply(pwms, function(x) unlist(as.list(x[, -c(flank_size + 1, flank_size + 2)]))))
    .corr_significance(pwms, correction = correction)
}


#' Compute Motif-Based DMR Interactions
#'
#' @description Computes motif-based interactions between DMRs based on their
#' motif correlation. Identifies pairs of DMRs with significant motif correlation
#' and returns a data frame of interactions.
#' @param dmrs Dataframe or GRanges object containing DMR coordinates and motif information
#' @param genome Character. Genome version to use for sequence extraction (e.g., "hg19")
#' @param array Character. Array platform type (e.g., "450K", "EPIC") (default: "450K")
#' @param max_fdr Numeric. Minimum FDR threshold for considering motif correlation significant (default: 0.05)
#' @param genomic_locs Data frame. Optional pre-computed genomic locations. If NULL,
#' locations will be retrieved using getSortedGenomicLocs (default: NULL)
#' @param correction Character. Multiple testing correction method (default: "BH")
#' @param flank_size Integer. Number of base pairs to include as flanking regions around each CpG site (default: 5)
#' @return Data frame of motif-based DMR interactions with columns:
#' \itemize{
#'   \item start_chr: Chromosome of the start DMR
#'   \item start: Start position of the start DMR
#'   \item end: End position of the start DMR
#'   \item end_chr: Chromosome of the end DMR
#'   \item end: Start position of the end DMR
#'   \item end_end: End position of the end DMR
#'   \item corr: Correlation value of the motif between the DMRs
#'   \item fdr: FDR value of the motif correlation between the DMRs
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
#' interactions <- computeMotifBasedDMRsInteraction(
#'     dmrs_with_motifs,
#'     genome = "hg19",
#'     array = "450K",
#'     max_fdr = 0.05
#' )
#' @export
computeMotifBasedDMRsInteraction <- function(dmrs, genome, array, max_fdr = 0.05, min_corr = 0.7, genomic_locs = NULL, correction = "BH", flank_size = 5) {
    dmrs <- convert_to_granges(dmrs, genome)
    if (! "pwm" %in% colnames(mcols(dmrs))) {
        dmrs <- extractDMRMotifs(dmrs, genome, array, genomic_locs = genomic_locs, flank_size = flank_size)
    }
    corr_ret <- .extractMotifsSimilaritySignificance(dmrs,  correction = correction, flank_size = flank_size)
    correlation_matrix <- corr_ret$corr_mat
    significance_matrix <- corr_ret$sig_mat
    mask <- !is.na(significance_matrix) & (significance_matrix <= max_fdr) & (abs(correlation_matrix) >= min_corr)
    if (all(!mask, na.rm = TRUE)) {
        .log_warn("No significant motif-based interactions found between DMRs at FDR <=", max_fdr)
        return(NULL)
    }
    g1 <- igraph::graph_from_adjacency_matrix(mask, mode = "undirected", diag = FALSE)
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
        paste(BASE_LEVELS[apply(pwm, 2, which.max)], collapse = "")
    })
    # Order by component size
    components_df <- components_df[order(-components_df$size), ]


    rowcol_df <- which(mask, arr.ind = TRUE)
    start_dmrs <- dmrs[rowcol_df[, 1], ]
    end_dmrs <- dmrs[rowcol_df[, 2], ]
    interaction_data_frame <- data.frame(
        chr1 = as.character(GenomeInfoDb::seqnames(start_dmrs)),
        start1 = GenomicRanges::start(start_dmrs),
        end1 = GenomicRanges::end(start_dmrs),
        chr2 = as.character(GenomeInfoDb::seqnames(end_dmrs)),
        start2 = GenomicRanges::start(end_dmrs),
        end2  = GenomicRanges::end(end_dmrs),
        corr = correlation_matrix[rowcol_df],
        fdr = significance_matrix[rowcol_df]
    )
    list(
        interactions = interaction_data_frame,
        components = components_df
    )
}