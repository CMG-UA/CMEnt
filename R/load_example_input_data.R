#' Load CMEnt Example Resources
#'
#' @description Load one or more example resources from the
#' \pkg{DMRsegaldata} ExperimentHub package and assign them into the caller's
#' environment using their resource names.
#'
#' @param ... Names of the resources to load, or none to load all available
#' resources. Available resources:
#' \itemize{
#'   \item "beta": Example beta values matrix
#'   \item "pheno": Example phenotype data
#'   \item "dmps": Example differentially methylated positions
#'   \item "array_type": Example array type annotation
#' }
#'
#' @return Invisibly returns the loaded object when a single resource is
#' requested, or a named list of loaded objects when multiple resources are
#' requested.
#'
#' @examples
#'
#' # Load phenotype data into the current environment
#' loadExampleInputData("pheno")
#' head(pheno)
#'
#' # Load multiple resources at once
#' loadExampleInputDataChr5And11("beta", "dmps", "pheno", "array_type")
#' dim(beta)
#'
#' @export
loadExampleInputData <- function(...) {
    .assignExampleInputData(
        resources = .normExInputResources(list(...)),
        envir = parent.frame()
    )
}

#' @rdname loadExampleInputData
#' @description Compatibility wrapper returning the chr5/chr11 example subset.
#' @export
loadExampleInputDataChr5And11 <- function(...) {
    .loadExampleInputDataSubset(
        ...,
        subset = c("chr5", "chr11"),
        envir = parent.frame()
    )
}

#' @rdname loadExampleInputData
#' @description Compatibility wrapper returning the chr21/chr22 example subset.
#' @export
loadExampleInputDataChr21And22 <- function(...) {
    .loadExampleInputDataSubset(
        ...,
        subset = c("chr21", "chr22"),
        envir = parent.frame()
    )
}
