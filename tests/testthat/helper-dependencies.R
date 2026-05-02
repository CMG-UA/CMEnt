is_installed_for_tests <- function(pkg_name) {
    nzchar(system.file(package = pkg_name))
}


skip_if_missing_array_annotation <- function(array = "450K", genome = "hg19") {
    pkg_name <- CMEnt:::.getArrayAnnotationPackage(array = array, genome = genome)
    if (!is_installed_for_tests(pkg_name)) {
        skip(paste0(
            "Package '", pkg_name, "' not installed; required for array/genome annotations in this test."
        ))
    }
}


skip_if_missing_bsgenome <- function(genome = "hg19") {
    pkg_name <- CMEnt:::.resolveBSGenomePackage(genome)
    if (is.null(pkg_name)) {
        skip(paste0("No BSgenome package mapping is available for genome '", genome, "'."))
    }
    if (!is_installed_for_tests(pkg_name)) {
        skip(paste0(
            "Package '", pkg_name, "' not installed; required for sequence-based operations in this test."
        ))
    }
}


skip_if_missing_gene_annotations <- function(genome = "hg38") {
    pkgs <- unname(CMEnt:::.getGeneAnnotationPackages(genome))
    missing <- pkgs[!vapply(pkgs, is_installed_for_tests, logical(1))]
    if (length(missing) > 0L) {
        skip(paste0(
            "Missing gene annotation package(s): ",
            paste(sprintf("'%s'", missing), collapse = ", "),
            "."
        ))
    }
}


skip_if_missing_motif_extraction_deps <- function(array = "450K", genome = "hg19") {
    skip_if_missing_array_annotation(array = array, genome = genome)
    skip_if_missing_bsgenome(genome = genome)
}


skip_if_missing_jaspar <- function(version = getOption("CMEnt.jaspar_version", 2024)) {
    pkg_name <- paste0("JASPAR", version)
    if (!is_installed_for_tests(pkg_name)) {
        skip(paste0(
            "Package '", pkg_name, "' not installed; required for JASPAR lookups in this test."
        ))
    }
}


loadExampleInputDataChr5And11 <- function(resource, use_experiment_hub = TRUE) {
    if (identical(resource, "dmps")) {
        skip_if_missing_array_annotation(array = "450K", genome = "hg19")
    }
    CMEnt::loadExampleInputDataChr5And11(resource, use_experiment_hub = use_experiment_hub)
}


loadExampleInputData <- function(resource, use_experiment_hub = TRUE) {
    if (identical(resource, "dmps")) {
        skip_if_missing_array_annotation(array = "450K", genome = "hg19")
    }
    CMEnt::loadExampleInputData(resource, use_experiment_hub = use_experiment_hub)
}
