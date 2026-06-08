.onLoad <- function(libname, pkgname) {
    # Set executable permissions for scripts in bin directory
    bin_dir <- system.file("bin", package = pkgname)
    if (dir.exists(bin_dir)) {
        files <- list.files(bin_dir, full.names = TRUE)
        for (f in files) {
            Sys.chmod(f, mode = "0755")
        }
    }
}


.updateCMEntOption <- function(name, default) {
    if (is.null(getOption(name))) {
        options(stats::setNames(list(default), name))
    }
    invisible()
}

.updateCMEntOption("CMEnt.njobs", .defaultNJobs())
.updateCMEntOption("CMEnt.beta_in_mem_threshold_mb", 500)
.updateCMEntOption("CMEnt.verbose", 1)

.updateCMEntOption("CMEnt.use_annotation_cache", TRUE)
.updateCMEntOption(
    "CMEnt.annotation_cache_dir",
    .getOSCacheDir(file.path("R", "CMEnt", "annotation_cache"))
)
.updateCMEntOption(
    "CMEnt.jaspar_cache_dir",
    .getOSCacheDir(file.path("R", "CMEnt", "jaspar_cache"))
)
.updateCMEntOption("CMEnt.scoring_nfold", 5)
.updateCMEntOption("CMEnt.min_motif_similarity", 0.8)
.updateCMEntOption("CMEnt.jaspar_version", 2024)
.updateCMEntOption("CMEnt.jaspar_tax_group", "vertebrates")
.updateCMEntOption("CMEnt.jaspar_corr_threshold", 0.9)
.updateCMEntOption("CMEnt.make_debug_dir", FALSE)
.updateCMEntOption("CMEnt.debug_dir", file.path(tempdir(), "CMEnt_debug"))

options(cli.num_colors = cli::num_ansi_colors())
