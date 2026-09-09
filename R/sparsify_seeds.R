#' Sparsify seed locations by minimum distance
#'
#' @description
#' Removes seeds inside dense runs where neighboring seeds are closer than a
#' minimum distance on the same chromosome.
#'
#' @param seeds Data frame with genomic seed locations containing `chr` and `start` columns.
#' @param min_distance Non-negative numeric scalar. Seeds with a distance below this
#'   value to both neighboring seeds on the same chromosome are removed.
#'
#' @return A subset of `seeds` with original columns and row names preserved.
#'
#' @details
#' Seeds are compared per chromosome after sorting by `start`. A seed is removed
#' only when both the previous and next seed on that chromosome are closer than
#' `min_distance`. Endpoints of close runs are retained. Seeds at exactly
#' `min_distance` are retained.
#'
#' @examples
#' seeds <- data.frame(
#'     chr = c("chr1", "chr1", "chr1"),
#'     start = c(100, 105, 120)
#' )
#' sparsifySeeds(seeds, min_distance = 10)
#'
#' @export
sparsifySeeds <- function(seeds, min_distance) {
    if (is.null(seeds) || nrow(seeds) == 0L) {
        return(seeds)
    }
    if (!all(c("chr", "start") %in% colnames(seeds))) {
        stop("Seeds data frame must contain 'chr' and 'start' columns.")
    }
    if (!is.numeric(min_distance) ||
        length(min_distance) != 1L ||
        !is.finite(min_distance) ||
        min_distance < 0) {
        stop("min_distance must be a non-negative numeric value.")
    }
    if (min_distance == 0) {
        return(seeds)
    }

    seed_chr <- as.character(seeds$chr)
    seed_start <- suppressWarnings(as.numeric(seeds$start))
    comparable <- !is.na(seed_chr) & nzchar(seed_chr) & is.finite(seed_start)
    if (!any(comparable)) {
        return(seeds)
    }

    keep <- rep(TRUE, nrow(seeds))
    ord <- order(seed_chr[comparable], seed_start[comparable], which(comparable))
    comparable_idx <- which(comparable)[ord]

    for (chr in unique(seed_chr[comparable_idx])) {
        chr_idx <- comparable_idx[seed_chr[comparable_idx] == chr]
        if (length(chr_idx) < 3L) {
            next
        }
        chr_start <- seed_start[chr_idx]
        prev_gap <- c(Inf, diff(chr_start))
        next_gap <- c(diff(chr_start), Inf)
        drop <- prev_gap < min_distance & next_gap < min_distance
        if (any(drop)) {
            keep[chr_idx[drop]] <- FALSE
        }
    }

    seeds[keep, , drop = FALSE]
}
