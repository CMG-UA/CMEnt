.loadContainerSetup <- function() {
    env <- new.env(parent = baseenv())
    source(system.file("docker", "install_dependencies.R", package = "CMEnt"),
           local = env)
    env
}

test_that("container setup resolves the requested genome and input features", {
    setup <- .loadContainerSetup()
    args <- list(beta = "beta.tsv", array = "EPIC", genome = NULL,
                 annotate_with_genes = TRUE, bed_provided = FALSE)
    expect_setequal(setup$.containerDependencyPackages("buildDMRs", args), c(
        "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
        "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db",
        "BSgenome.Hsapiens.UCSC.hg19"
    ))

    args$array <- "EPICv2"
    args$annotate_with_genes <- FALSE
    expect_setequal(setup$.containerDependencyPackages("buildDMRs", args), c(
        "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
        "BSgenome.Hsapiens.UCSC.hg38"
    ))
    args$bed_provided <- TRUE
    expect_equal(setup$.containerDependencyPackages("buildDMRs", args),
                 "BSgenome.Hsapiens.UCSC.hg38")

    args$array <- "Mouse"
    args$bed_provided <- FALSE
    expect_setequal(setup$.containerDependencyPackages("buildDMRs", args), c(
        "IlluminaMouseMethylationanno.12.v1.mm10",
        "BSgenome.Mmusculus.UCSC.mm10"
    ))
    expect_equal(setup$.containerDependencyPackages("findDMPsBSSeq", list()), "DSS")
    args$genome <- "mm10"
    expect_equal(setup$.containerDependencyPackages("findDMPsArray", args),
                 "IlluminaMouseMethylationanno.12.v1.mm10")
    args$sorted_locs <- "locations.rds"
    expect_length(setup$.containerDependencyPackages("findDMPsArray", args), 0L)
})

test_that("viewer container setup uses the saved analysis genome", {
    setup <- .loadContainerSetup()
    prefix <- tempfile()
    paths <- CMEnt:::.viewerOutputFiles(prefix)
    withr::defer(unlink(unlist(paths)))
    file.create(paths$dmrs_file, paths$beta_file)
    saveRDS(list(genome = "hs1", array = "EPICv2"), paths$meta_file)
    cache <- tempfile()
    withr::local_options(CMEnt.jaspar_cache_dir = cache)
    withr::defer(unlink(cache, recursive = TRUE))
    expect_setequal(setup$.containerDependencyPackages(
        "launchCMEntViewer", list(output_prefix = prefix)
    ), c(
        "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
        "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db",
        "BSgenome.Hsapiens.UCSC.hs1", "JASPAR2024"
    ))
    saveRDS(list(genome = "mm39", array = NULL), paths$meta_file)
    expect_setequal(setup$.containerDependencyPackages(
        "launchCMEntViewer", list(output_prefix = prefix)
    ), c("BSgenome.Mmusculus.UCSC.mm39", "JASPAR2024"))
})

test_that("container installation only installs missing selected packages", {
    skip_if_not_installed("BiocManager")
    setup <- .loadContainerSetup()
    state <- new.env(parent = emptyenv())
    state$installed <- "BSgenome.Hsapiens.UCSC.hg38"
    state$calls <- list()
    local_mocked_bindings(
        .isPackageInstalled = function(pkg) pkg %in% state$installed,
        .package = "CMEnt"
    )
    local_mocked_bindings(
        install = function(pkgs, ask, update, ...) {
            state$calls <- c(state$calls, list(list(pkgs, ask, update)))
            state$installed <- union(state$installed, pkgs)
        },
        .package = "BiocManager"
    )
    requested <- c("BSgenome.Hsapiens.UCSC.hg38", "DSS")
    setup$.installContainerPackages(requested)
    expect_equal(state$calls, list(list("DSS", FALSE, FALSE)))
    setup$.installContainerPackages(requested)
    expect_length(state$calls, 1L)
})

test_that("container setup fails before analysis when installation fails", {
    skip_if_not_installed("BiocManager")
    setup <- .loadContainerSetup()
    local_mocked_bindings(
        .isPackageInstalled = function(pkg) FALSE,
        .package = "CMEnt"
    )
    local_mocked_bindings(install = function(...) invisible(NULL),
                          .package = "BiocManager")
    expect_error(setup$.installContainerPackages("DSS"),
                 "CMEnt container setup requires the following package: 'DSS'",
                 fixed = TRUE)
})
