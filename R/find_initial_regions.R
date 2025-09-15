#' Initial Region Finding Functions
#' 
#' Functions for identifying initial DMR candidates from a set of 
#' Differentially Methylated Positions (DMPs).
#'
#' @keywords internal
NULL

#' Find Initial DMR Regions from DMPs
#'
#' This function takes a set of significant DMPs and identifies initial DMR candidates
#' by grouping nearby DMPs based on genomic distance. These initial regions serve as
#' seeds for subsequent DMR expansion.
#'
#' @param dmps Data frame containing DMP information:
#'        \itemize{
#'          \item chr: Chromosome
#'          \item pos: Genomic position
#'          \item dmp: DMP ID
#'          \item pval: Statistical significance
#'          \item delta_beta: Methylation difference
#'        }
#' @param max.lookup.dist Maximum distance (in base pairs) between DMPs to be considered
#'        part of the same initial region
#' @param min.dmps Minimum number of DMPs required to define an initial region
#'
#' @return A data frame containing initial DMR regions with:
#'   \itemize{
#'     \item chr: Chromosome
#'     \item start: Region start position
#'     \item end: Region end position
#'     \item n_dmps: Number of DMPs in the region
#'     \item dmps: Comma-separated list of DMPs
#'   }
#'
#' @details
#' The function:
#' 1. Sorts DMPs by chromosome and position
#' 2. Groups DMPs within max.lookup.dist of each other
#' 3. Filters regions to ensure minimum DMP count
#' 4. Returns regions suitable for expansion
#'
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