#' DMRsegal: Identify Differentially Methylated Regions from DMPs
#'
#' DMRsegal is a comprehensive tool for identifying Differentially Methylated Regions (DMRs)
#' from pre-computed Differentially Methylated Positions (DMPs). The package uses a
#' correlation-based approach to expand regions around significant DMPs, considering both
#' statistical significance and biological relevance of methylation changes.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{findDMRsFromSeeds}}: The main function for identifying DMRs
#' }
#'
#' @section Input Data:
#' The package accepts two types of methylation data input:
#' \itemize{
#'   \item Beta value files: Tab-separated files containing methylation beta values
#'   \item Tabix files: Indexed bed.gz files for efficient access to large datasets
#' }
#'
#' @section Platform Support:
#' The package supports methylation data from both major Illumina platforms:
#' \itemize{
#'   \item Illumina EPIC array (850k)
#'   \item Illumina 450k array
#'   \item Custom CpG positions (when provided in DMP data)
#' }
#'
#' @section Workflow:
#' 1. Start with pre-computed DMPs from any differential methylation analysis
#' 2. Provide methylation data (beta values or tabix)
#' 3. Define parameters for region expansion
#' 4. Get DMRs as a GRanges object with comprehensive statistics
#'
#' @keywords internal
"_PACKAGE"

#' @import GenomicRanges
#' @import IRanges
#' @import methods
#' @importFrom stats p.adjust na.omit pt setNames
#' @importFrom utils read.table write.table capture.output head install.packages
#' @importFrom data.table fread
NULL
