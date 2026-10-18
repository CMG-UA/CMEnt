# Container setup only. CMEnt package functions never invoke this installer.
.containerDependencyPackages <- function(command, args) {
    requirements <- switch(
        command,
        buildDMRs = {
            array <- CMEnt:::.normalizeBuildDMRsArray(args$array)
            genome <- CMEnt:::.resolveBuildDMRsGenome(
                args$beta, array = array, genome = args$genome,
                bed_provided = args$bed_provided
            )
            CMEnt:::.buildDMRsDependencyReqs(
                beta = args$beta, array = array, genome = genome,
                annotate_with_genes = args$annotate_with_genes,
                bed_provided = args$bed_provided
            )
        },
        findDMPsArray = {
            if (is.null(args$sorted_locs) || !nzchar(args$sorted_locs) ||
                tolower(args$sorted_locs) == "null") {
                CMEnt:::.arrayAnnotDependencyReqs(args$array, args$genome)
            } else {
                CMEnt:::.makeDependencyRequirements(character(), "")
            }
        },
        findDMPsBSSeq = CMEnt:::.makeDependencyRequirements(
            "DSS", "BSseq DMP calling requires DSS."
        ),
        launchCMEntViewer = .viewerContainerDependencyReqs(args$output_prefix),
        stop("Unsupported CMEnt command '", command, "'.", call. = FALSE)
    )
    unique(requirements$pkg_name)
}

.viewerContainerDependencyReqs <- function(output_prefix) {
    validation <- CMEnt:::.validateOutputPrefix(output_prefix)
    if (!validation$valid) {
        stop(paste(validation$errors, collapse = "\n"), call. = FALSE)
    }
    meta <- readRDS(paste0(output_prefix, ".meta.rds"))
    genome <- CMEnt:::.normalizeBuildDMRsGenome(meta$genome)
    if (is.null(genome)) {
        stop("Viewer metadata file must contain a valid `genome` entry.",
             call. = FALSE)
    }
    array <- CMEnt:::.normalizeBuildDMRsArray(meta$array)
    annotation_genome <- if (genome == "hs1") "hg38" else genome
    has_gene_annotations <- any(grepl(
        annotation_genome, tolower(CMEnt:::supportedOrganisms()$TxDb)
    ))
    CMEnt:::.combineDependencyRequirements(
        if (!is.null(array)) {
            CMEnt:::.arrayAnnotDependencyReqs(array, genome)
        },
        if (has_gene_annotations) CMEnt:::.geneAnnotDependencyReqs(genome),
        CMEnt:::.bsgenomeDependencyReqs(genome),
        CMEnt:::.jasparDependencyRequirements()
    )
}

.installContainerPackages <- function(pkg_names) {
    missing <- pkg_names[!vapply(
        pkg_names, CMEnt:::.isPackageInstalled, logical(1)
    )]
    if (!length(missing)) {
        return(invisible(pkg_names))
    }
    groups <- CMEnt:::.dependencyInstallGroups(missing)
    if (length(groups$cran)) {
        utils::install.packages(groups$cran)
    }
    if (length(groups$bioc)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            utils::install.packages("BiocManager")
        }
        BiocManager::install(groups$bioc, ask = FALSE, update = FALSE)
    }
    if (length(groups$github)) {
        if (!requireNamespace("devtools", quietly = TRUE)) {
            utils::install.packages("devtools")
        }
        for (repo in groups$github) {
            devtools::install_github(repo, upgrade = "never")
        }
    }
    CMEnt:::.assertPackagesInstalled(pkg_names, "CMEnt container setup")
    invisible(pkg_names)
}

.prepareCMEntContainer <- function(command, args) {
    .installContainerPackages(.containerDependencyPackages(command, args))
}
