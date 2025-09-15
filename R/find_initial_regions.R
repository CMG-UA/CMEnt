#' Initial region finding functions
#' @keywords internal

#' Find initial DMR regions from DMPs
#'
#' @param dmps Data frame of DMPs
#' @param max.lookup.dist Maximum distance for DMR lookup
#' @param min.dmps Minimum number of DMPs per region
#' @return Data frame of initial regions
#' @keywords internal
.findInitialRegions <- function(dmps, max.lookup.dist, min.dmps) {
    # Ensure DMPs are sorted
    if (!all(diff(dmps$pos[dmps$chr == dmps$chr[1]]) >= 0)) {
        dmps <- dmps[order(dmps$chr, dmps$pos), ]
    }
    
    # Initialize variables
    regions <- data.frame()
    current_chr <- NULL
    current_start <- NULL
    current_dmps <- c()
    
    # Process each DMP
    for (i in seq_len(nrow(dmps))) {
        dmp <- dmps[i, ]
        
        # Start new region if different chromosome
        if (is.null(current_chr) || dmp$chr != current_chr) {
            if (!is.null(current_chr) && length(current_dmps) >= min.dmps) {
                regions <- rbind(regions, .createRegion(current_dmps, dmps))
            }
            current_chr <- dmp$chr
            current_start <- dmp$pos
            current_dmps <- i
            next
        }
        
        # Check if DMP should be added to current region
        if (dmp$pos - dmps$pos[tail(current_dmps, 1)] <= max.lookup.dist) {
            current_dmps <- c(current_dmps, i)
        } else {
            # Save current region if it meets minimum DMP requirement
            if (length(current_dmps) >= min.dmps) {
                regions <- rbind(regions, .createRegion(current_dmps, dmps))
            }
            # Start new region
            current_start <- dmp$pos
            current_dmps <- i
        }
    }
    
    # Handle last region
    if (length(current_dmps) >= min.dmps) {
        regions <- rbind(regions, .createRegion(current_dmps, dmps))
    }
    
    regions
}

#' Create a region from a set of DMPs
#'
#' @param dmp_indices Indices of DMPs in the region
#' @param dmps Full DMP data frame
#' @return Data frame row representing the region
#' @keywords internal
.createRegion <- function(dmp_indices, dmps) {
    region_dmps <- dmps[dmp_indices, ]
    data.frame(
        chr = region_dmps$chr[1],
        start = min(region_dmps$pos),
        end = max(region_dmps$pos),
        start_dmp = rownames(region_dmps)[1],
        end_dmp = rownames(region_dmps)[nrow(region_dmps)],
        n_dmps = nrow(region_dmps),
        mean_pval = mean(region_dmps$pval_adj),
        stringsAsFactors = FALSE
    )
}