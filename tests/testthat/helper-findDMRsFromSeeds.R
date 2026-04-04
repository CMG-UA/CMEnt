subsetDenseExampleDmpsChr5And11 <- function(dmps,
                                            target_n = 100L,
                                            max_lookup_dist = 1000L,
                                            min_dense_run = 20L) {
    target_n <- as.integer(target_n)[1]
    max_lookup_dist <- as.integer(max_lookup_dist)[1]
    min_dense_run <- as.integer(min_dense_run)[1]

    stopifnot(target_n > 0L)
    stopifnot(max_lookup_dist >= 0L)
    stopifnot(min_dense_run >= 2L)

    locs <- getSortedGenomicLocs(array = "450k")
    locs <- locs[locs$chr %in% c("chr5", "chr11"), , drop = FALSE]
    locs <- locs[rownames(locs) %in% rownames(dmps), , drop = FALSE]

    dmps_sorted <- dmps[rownames(locs), , drop = FALSE]
    ord <- order(
        as.character(locs[rownames(dmps_sorted), "chr"]),
        locs[rownames(dmps_sorted), "start"]
    )
    dmps_sorted <- dmps_sorted[ord, , drop = FALSE]
    locs_sorted <- locs[rownames(dmps_sorted), , drop = FALSE]

    if (nrow(dmps_sorted) <= target_n) {
        return(dmps_sorted)
    }

    dense_edges <- c(
        FALSE,
        as.character(locs_sorted$chr[-1]) == as.character(locs_sorted$chr[-nrow(locs_sorted)]) &
            (locs_sorted$start[-1] - locs_sorted$start[-nrow(locs_sorted)]) <= max_lookup_dist
    )
    runs <- rle(dense_edges)
    ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L
    dense_run_inds <- which(runs$values & runs$lengths >= (min_dense_run - 1L))

    if (length(dense_run_inds) == 0L) {
        stop("Could not find a dense DMP run for the chr5/chr11 test subset.")
    }

    block_start <- max(1L, starts[dense_run_inds[1]] - 1L)
    block_end <- min(nrow(dmps_sorted), block_start + target_n - 1L)
    dmps_sorted[block_start:block_end, , drop = FALSE]
}
