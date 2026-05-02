#' Example Beta Matrix
#'
#' A compact methylation beta-value matrix bundled for package examples, tests,
#' and vignettes.
#'
#' @format A numeric matrix with CpG probe IDs as row names and sample IDs as
#' column names.
"beta"

#' Example Phenotype Table
#'
#' Sample-level phenotype metadata corresponding to the bundled example beta
#' matrix.
#'
#' @format A data frame with sample IDs as row names.
"pheno"

#' Example Differentially Methylated Positions
#'
#' Precomputed seed sites used in package examples, tests, and vignettes.
#'
#' @format A data frame keyed by CpG probe IDs.
"dmps"

#' Example Array Type
#'
#' Platform identifier for the bundled example methylation data.
#'
#' @format A character scalar.
"array_type"
