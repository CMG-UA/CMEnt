#' Example methylation beta values
#'
#' A matrix of simulated DNA methylation beta values containing both
#' differentially methylated regions (DMRs) and background CpG sites.
#'
#' @format A matrix with 1000 rows (CpG sites) and 10 columns (samples):
#' \describe{
#'   \item{Rows}{CpG sites with IDs in the format 'cgXXXXXXX'}
#'   \item{Columns}{Samples, with 5 cases and 5 controls}
#'   \item{Values}{Beta values between 0 and 1, representing methylation levels}
#' }
#'
#' @details The data contains 10 simulated DMRs with correlated methylation
#'   changes between case and control groups. Effect sizes vary between 0.2
#'   and 0.4, with added biological variation.
#'
#' @source Simulated data created to demonstrate DMRSegal functionality
"example_beta"

#' Example DMP analysis results
#'
#' A data frame containing the results of differential methylation analysis
#' at individual CpG sites (DMPs).
#'
#' @format A data frame with 1000 rows and the following columns:
#' \describe{
#'   \item{chr}{Chromosome}
#'   \item{pos}{Genomic position}
#'   \item{dmp}{CpG ID}
#'   \item{pval}{Raw p-value from t-test}
#'   \item{pval_adj}{Benjamini-Hochberg adjusted p-value}
#'   \item{qval}{Q-value}
#'   \item{delta_beta}{Mean methylation difference (Case - Control)}
#'   \item{cases_beta}{Mean beta value in cases}
#'   \item{controls_beta}{Mean beta value in controls}
#'   \item{cases_beta_sd}{Standard deviation of beta values in cases}
#'   \item{controls_beta_sd}{Standard deviation of beta values in controls}
#'   \item{cases_num}{Number of case samples}
#'   \item{controls_num}{Number of control samples}
#' }
#'
#' @source Simulated data created to demonstrate DMRSegal functionality
"example_dmps"

#' Example phenotype data
#'
#' A data frame containing sample metadata for the example dataset.
#'
#' @format A data frame with 10 rows and 4 columns:
#' \describe{
#'   \item{sample_id}{Sample identifier}
#'   \item{group}{Factor with levels "Case" and "Control"}
#'   \item{age}{Simulated age values}
#'   \item{sex}{Factor with levels "M" and "F"}
#' }
#'
#' @source Simulated data created to demonstrate DMRSegal functionality
"example_pheno"

#' True DMR locations in example data
#'
#' A data frame containing the true locations of simulated DMRs in the
#' example dataset.
#'
#' @format A data frame with 10 rows and 3 columns:
#' \describe{
#'   \item{start}{Start position of the DMR}
#'   \item{end}{End position of the DMR}
#'   \item{n_cpgs}{Number of CpG sites in the DMR}
#' }
#'
#' @source Simulated data created to demonstrate DMRSegal functionality
"example_true_dmrs"