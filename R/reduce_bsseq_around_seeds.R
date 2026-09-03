#' Reduce BSseq object around seed locations
#'
#' @description
#' This function reduces a BSseq object to only include methylation data around specified seed locations.
#' It expands the seed locations by a given window size and merges overlapping windows to create a reduced BSseq object.
#'
#' @param bsseq_obj BSseq object containing the methylation data
#' @param seeds Data frame with genomic locations containing 'chr' and 'start' columns.
#' @param expansion_window Integer. DMR maximum stage 2 bi-directional size expansion window (default: 10000).
#'
#' @return BSseq object containing the reduced methylation data
#'
#' @details
#' The reduction is performed by expanding the seed locations by the specified window size to both sides, after division by 2,
#'  and merging overlapping windows.The resulting BSseq object contains only the methylation data within these merged windows.
#'
#' @examples
#'
#' library(bsseq)
#' # Example BSseq object
#' bsseq_obj <- BSseq(
#'     chr = c("chr1", "chr1", "chr2"),
#'     pos = c(100, 200, 300),
#'     M = matrix(c(10, 20, 30), ncol = 1),
#'     Cov = matrix(c(100, 200, 300), ncol = 1)
#' )
#' # Example seeds data frame
#' seeds <- data.frame(chr = c("chr1", "chr2"), start = c(80, 250))
#' reduced_bsseq_obj <- reduceBSseqAroundSeeds(bsseq_obj, seeds, expansion_window = 50)
#'
#' @export
reduceBSseqAroundSeeds <- function(bsseq_obj, seeds, expansion_window) {
    if (is.null(seeds) || nrow(seeds) == 0L || !is.finite(expansion_window) || expansion_window <= 0) {
        .log_warn("No seeds provided or invalid expansion window. Returning the original BSseq object.")
        return(bsseq_obj)
    }
    if (!all(c("chr", "start") %in% colnames(seeds))) {
        stop("Seeds data frame must contain 'chr' and 'start' columns.")
    }
    if (!methods::is(bsseq_obj, "BSseq")) {
        stop("bsseq_obj must be a BSseq object.")
    }
    if (!is.numeric(expansion_window) || length(expansion_window) != 1 || expansion_window <= 0) {
        stop("expansion_window must be a positive numeric value.")
    }
    if (!"end" %in% colnames(seeds)) {
        seeds$end <- seeds$start + 1
    }
    flank <- as.numeric(expansion_window) / 2
    windows <- data.frame(
        chr = as.character(seeds$chr),
        start = pmax(1, as.numeric(seeds$start) - flank),
        end = as.numeric(seeds$end) + flank
    )
    windows <- .mergeGenomicWindows(windows, return_granges = TRUE)
    bsseq_obj <- bsseq::subsetByOverlaps(bsseq_obj, windows)
    .prepareBSseqForBetaHandler(bsseq_obj)
}
