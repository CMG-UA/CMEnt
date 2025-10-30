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
#'' # Access motif frequencies for the first DMR
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
    c <- cor(x)
    above <- row(c) < col(c)
    r2 <- c[above]^2
    fstat <- r2 * dfr / (1 - r2)
    c[above] <- 1 - pf(fstat, 1, dfr)
    if (!is.null(correction)){
        c[above] <- p.adjust(c[above], method=correction)
    }
    corr_mat <- t(c)
    corr_mat[upper.tri(corr_mat)] <- NA
    cor
}

.extractMotifsSimilaritySignificance <- function(dmrs_with_motifs, correction = "BH"){
    if (inherits(dmrs_with_motifs, "GRanges")) {
        pwms <- mcols(dmrs_with_motifs)$pwm
    } else {
        pwms <- dmrs_with_motifs[, "pwm"]
    }
    pwms <- do.call(rbind, lapply(pwms, unlist))
    .corr_significance(pwms, correction = correction)
}


#' Compute Motif-Based DMR Interactions
#'
#' @description Computes motif-based interactions between DMRs based on their
#' motif similarity. Identifies pairs of DMRs with significant motif similarity
#' and returns a data frame of interactions.
#' @param dmrs Dataframe or GRanges object containing DMR coordinates and motif information
#' @param genome Character. Genome version to use for sequence extraction (e.g., "hg19")
#' @param array Character. Array platform type (e.g., "450K", "EPIC") (default: "450K")
#' @param min_fdr Numeric. Minimum FDR threshold for considering motif similarity significant (default: 0.05)
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
#'   \item fdr: FDR value of the motif similarity between the DMRs
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
#'     min_fdr = 0.05
#' )
#' @export
computeMotifBasedDMRsInteraction <- function(dmrs, genome, array, min_fdr = 0.05,  genomic_locs = NULL, correction = "BH", flank_size = 5) {
    dmrs <- convert_to_granges(dmrs, genome)
    if (! "pwm" %in% colnames(mcols(dmrs))) {
        dmrs <- extractDMRMotifs(dmrs, genome, array, genomic_locs = genomic_locs, flank_size = flank_size)
    }
    significance_matrix <- .extractMotifsSimilaritySignificance(dmrs,  correction = correction)
    rowcol_df <- which(!is.na(significance_matrix) && significance_matrix <= min_fdr,  arr.ind = TRUE)
    if (nrow(rowcol_df) == 0) {
        .log_warn("No significant motif-based interactions found between DMRs at FDR <=", min_fdr)
        return(NULL)
    }
    start_dmrs <- dmrs[rowcol_df[, 1], ]
    end_dmrs <- dmrs[rowcol_df[, 2], ]
    interaction_data_frame <- data.frame(
        start_chr = as.character(GenomeInfoDb::seqnames(start_dmrs)),
        start = GenomicRanges::start(start_dmrs),
        end = GenomicRanges::end(start_dmrs),
        end_chr = as.character(GenomeInfoDb::seqnames(end_dmrs)),
        end = GenomicRanges::start(end_dmrs),
        end_end = GenomicRanges::end(end_dmrs),
        fdr = significance_matrix[rowcol_df]
    )
    interaction_data_frame
}