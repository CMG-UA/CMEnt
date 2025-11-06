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
update_package_option("DMRsegal.random_seed", mode = "numeric", default = 42)
update_package_option("DMRsegal.njobs", mode = "numeric", default = min(8, future::availableCores() - 1))
update_package_option("DMRsegal.verbose", mode = "numeric", default = 1)
update_package_option("DMRsegal.use_bed_cache", mode = "logical", default = FALSE)
update_package_option("DMRsegal.bed_cache_dir", mode = "character", default = file.path(path.expand("~"), ".cache", "R", "DMRsegal", "bed_cache"))
update_package_option("DMRsegal.use_tabix_cache", mode = "logical", default = FALSE)
update_package_option("DMRsegal.tabix_cache_dir", mode = "character", default = file.path(path.expand("~"), ".cache", "R", "DMRsegal", "tabix_cache"))
update_package_option("DMRsegal.use_annotation_cache", mode = "logical", default = TRUE)
update_package_option("DMRsegal.annotation_cache_dir", mode = "character", default = file.path(path.expand("~"), ".cache", "R", "DMRsegal", "annotation_cache"))
update_package_option("DMRsegal.jaspar_cache_dir", mode = "character", default = file.path(
    path.expand("~"),
    ".cache", "R", "DMRsegal", "jaspar_cache"
))
update_package_option("DMRsegal.jaspar_version", mode = "numeric", default = 2024)
update_package_option("DMRsegal.jaspar_tax_group", mode = "character", default = "vertebrates")
update_package_option("DMRsegal.jaspar_corr_threshold", mode = "numeric", default = 0.9)
update_package_option("DMRsegal.make_debug_dir", mode = "logical", default = FALSE)
