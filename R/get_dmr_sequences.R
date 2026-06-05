#' Extract DNA Sequences for DMRs
#'
#' @description Retrieves the DNA sequences corresponding to genomic regions
#' specified in a GRanges object. This function is useful for extracting the
#' actual DNA sequence of identified DMRs for downstream analyses such as
#' motif finding or sequence composition analysis.
#'
#' @param dmrs GRanges object containing genomic coordinates of DMRs
#' @param genome Character. Genome version to use for sequence extraction, .e.g. "hg38" or "hs1".
#' @param use_online Logical. If TRUE, forces use of online UCSC API instead of
#'   BSgenome packages. If FALSE (default), uses BSgenome packages with online
#'   fallback when packages are unavailable (default: FALSE)
#' @param uflank_size Integer. Number of base pairs to add as flanking regions
#'   upstream of each DMR (default: 0)
#' @param dflank_size Integer. Number of base pairs to add as flanking regions
#'   downstream of each DMR (default: 0)
#' @param batch_size Integer. For online API, number of regions to process per batch (default: 100)
#' @param njobs Integer. For online API, number of cores for parallel processing (default: 1)
#' @return A Character vector containing DNA sequences for each DMR
#'
#' @details
#' The function first attempts to use genome-appropriate BSgenome packages:
#' \itemize{
#'   \item hg19: BSgenome.Hsapiens.UCSC.hg19
#'   \item hg38: BSgenome.Hsapiens.UCSC.hg38
#'   \item hs1: BSgenome.Hsapiens.UCSC.hs1
#'   \item mm10: BSgenome.Mmusculus.UCSC.mm10
#'   \item mm39: BSgenome.Mmusculus.UCSC.mm39
#' }
#'
#' If the required BSgenome package is not installed, the function raises an
#' error with installation instructions. Set `use_online = TRUE` to query
#' sequences from the UCSC Genome Browser REST API instead. The online method
#' processes sequences in batches with optional parallel processing for
#' improved performance with large datasets.
#'
#' For large numbers of DMRs (>10k), consider using parallel processing by setting
#' njobs > 1 when using the online API, or install the appropriate BSgenome package
#' for much faster local sequence retrieval.
#'
#' @examples
#' dmrs <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100000, 100100))
#' \donttest{
#' # Extract sequences for DMRs using BSgenome packages
#' sequences <- getDMRSequences(dmrs, "hg19")
#'
#' # Force use of online UCSC API with parallel processing
#' sequences <- getDMRSequences(dmrs, "hg19", use_online = TRUE, njobs = 4)
#'
#' # Calculate GC content
#' gc_content <- vapply(sequences, function(s) {
#'     (stringr::str_count(s, "G") + stringr::str_count(s, "C")) / nchar(s)
#' }, numeric(1))
#' }
#'
#' @importFrom BSgenome getSeq
#' @importFrom rtracklayer import.chain liftOver
#' @export
getDMRSequences <- function(dmrs, genome, use_online = FALSE, uflank_size = 0, dflank_size = 0,
                            batch_size = 100, njobs = 1) {
    if (!use_online) {
        pkg_name <- .assertBSGenomePackageInstalled(genome, context = "getDMRSequences()")
        use_bsgenome <- TRUE
    } else {
        pkg_name <- NULL
        use_bsgenome <- FALSE
    }
    extended_ret <- .extendGRangesWithFlanks(dmrs, uflank_size, dflank_size)
    dmrs <- extended_ret$granges
    add_na_to_the_start <- extended_ret$start_off_limit
    add_na_to_the_end <- extended_ret$end_off_limit

    if (use_bsgenome) {
        .log_info("Querying sequences using BSgenome package...", level = 2)
        if (!isNamespaceLoaded(pkg_name)) {
            loadNamespace(pkg_name)
        }
        seq_db <- getExportedValue(pkg_name, pkg_name)
        sequences <- Biostrings::getSeq(seq_db, dmrs, as.character = TRUE)
        if (is.list(sequences)) {
            sequences <- vapply(sequences, function(x) paste(x, collapse = ""), character(1))
        }
    } else {
        .log_info("Querying sequences from UCSC Genome Browser API...", level = 2)
        sequences <- .getSequencesFromUCSC(dmrs, genome, batch_size = batch_size, njobs = njobs)
    }
    off_bound_mask <- add_na_to_the_start > 0 | add_na_to_the_end > 0
    if (any(off_bound_mask)) {
        .log_info("  Found ", sum(off_bound_mask), " regions with out-of-bound flanking extensions", level = 2)
        .log_info("Adding 'N' padding for out-of-bound flanking regions...", level = 2)
        sequences[off_bound_mask] <- mapply(function(seq, na_start, na_end) {
            if (is.na(seq)) {
                return(NA_character_)
            }
            seq <- paste0(
                paste(rep("N", na_start), collapse = ""),
                seq,
                paste(rep("N", na_end), collapse = "")
            )
            seq
        }, sequences[off_bound_mask], add_na_to_the_start[off_bound_mask], add_na_to_the_end[off_bound_mask], SIMPLIFY = TRUE)
    }
    sequences
}