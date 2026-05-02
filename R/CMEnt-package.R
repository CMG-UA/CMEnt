#' CMEnt: Characterization of Methylation using positional ENTanglement
#'
#' CMEnt is a comprehensive tool for identifying Differentially Methylated Regions (DMRs)
#' from genomic seeds, commonly Differentially Methylated Positions (DMPs). The package uses a
#' correlation-based approach to expand regions around significant seeds, considering both
#' statistical significance and biological relevance of methylation changes.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{findDMRsFromSeeds}}: Identify DMRs from genomic seeds
#'   \item \code{\link{extractDMRMotifs}}: Extract sequence motifs from DMRs
#'   \item \code{\link{computeDMRsInteraction}}: Compute motif-based DMR interactions
#'   \item \code{\link{plotDMR}}: Visualize individual DMRs with structure, beta values, and motifs
#'   \item \code{\link{plotDMRs}}: Create multi-panel DMR plots
#'   \item \code{\link{plotDMRsCircos}}: Generate circos plots with DMR interactions
#'   \item \code{\link{plotDMRsManhattan}}: Generate genome-wide Manhattan-style plots for DMR scores
#' }
#'
#' @section Input Data:
#' The package accepts multiple types of methylation data input:
#' \itemize{
#'   \item Beta value files: Tab-separated files containing methylation beta values
#'   \item Tabix files: Indexed bed.gz files for efficient access to large datasets
#'   \item BED files: Custom methylation BED files (automatically converted to tabix)
#'   \item Beta matrices: In-memory beta value matrices
#' }
#'
#' @section Platform Support:
#' The package supports methylation data from multiple platforms:
#' \itemize{
#'   \item Human arrays: 450K, EPIC (850k), EPICv2, 27K on hg19, hg38, or hs1
#'   \item Mouse arrays: mm10, mm39 genomes
#'   \item NGS data: Via tabix-indexed files
#'   \item Custom platforms: Via BED file input
#' }
#'
#' @section Workflow:
#' 1. Start with genomic seeds (e.g., DMPs) from any differential analysis
#' 2. Provide methylation data (multiple formats supported)
#' 3. Define parameters for region expansion and connectivity
#' 4. Get DMRs as a GRanges object with comprehensive statistics
#' 5. Optionally extract motifs and compute DMR interactions
#' 6. Visualize with structure plots, heatmaps, and circos plots
#'
#' @keywords internal
"_PACKAGE"

#' @import GenomicRanges
#' @import IRanges
#' @import methods
#' @importFrom stats aggregate as.formula ave p.adjust na.omit plogis predict pt qlogis setNames
#' @importFrom utils capture.output getFromNamespace head read.table tail write.table
#' @importFrom data.table fread
NULL

if (getRversion() >= "2.15.1") {
    utils::globalVariables(
        c(
            ".", ".N", "Beta", "Sample",
            "cases_beta", "cases_beta_sd", "chr", "controls_beta",
            "controls_beta_sd", "dmr_id", "end_bp", "end_score",
            "hover_text", "id", "label", "label_x", "label_y",
            "line_ymax", "line_ymin", "midpoint", "position",
            "rchr", "rend", "right_bp", "rstart", "score_raw",
            "score_smoothed", "start_bp", "start_score", "target_x",
            "target_y", "x", "xend", "xmax", "xmin", "y", "yend",
            "ymax", "ymax_band", "ymin", "ymin_band"
        )
    )
}
