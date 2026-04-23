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


update_package_option <- future:::update_package_option
update_package_option("CMEnt.random_seed", mode = "numeric", default = 42)
update_package_option("CMEnt.njobs", mode = "numeric", default = min(8, future::availableCores() - 1))
update_package_option("CMEnt.beta_in_mem_threshold_mb", mode = "numeric", default = 500)
update_package_option("CMEnt.verbose", mode = "numeric", default = 1)
update_package_option("CMEnt.use_bed_cache", mode = "logical", default = FALSE)
update_package_option("CMEnt.bed_cache_dir", mode = "character",
    default =  .getOSCacheDir(file.path("R", "CMEnt", "bed_cache"))
)
update_package_option("CMEnt.use_tabix_cache", mode = "logical", default = TRUE)
update_package_option("CMEnt.tabix_cache_dir", mode = "character",
    default =  .getOSCacheDir(file.path("R", "CMEnt", "tabix_cache"))
)
update_package_option("CMEnt.use_annotation_cache", mode = "logical", default = TRUE)
update_package_option("CMEnt.annotation_cache_dir", mode = "character",
    default = .getOSCacheDir(file.path("R", "CMEnt", "annotation_cache"))
)
update_package_option("CMEnt.h5_cache_dir", mode = "character",
    default = .getOSCacheDir(file.path("R", "CMEnt", "h5_cache"))
)
update_package_option(
    "CMEnt.jaspar_cache_dir", mode = "character",
    default = .getOSCacheDir(file.path("R", "CMEnt", "jaspar_cache"))
)
update_package_option("CMEnt.scoring_nfold", mode = "numeric", default = 5)
update_package_option("CMEnt.min_motif_similarity", mode = "numeric", default = 0.8)
update_package_option("CMEnt.jaspar_version", mode = "numeric", default = 2024)
update_package_option("CMEnt.jaspar_tax_group", mode = "character", default = "vertebrates")
update_package_option("CMEnt.jaspar_corr_threshold", mode = "numeric", default = 0.9)
update_package_option("CMEnt.make_debug_dir", mode = "logical", default = FALSE)
update_package_option("CMEnt.debug_dir", mode = "character", default = file.path(tempdir(), "CMEnt_debug"))