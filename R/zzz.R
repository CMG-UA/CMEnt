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

utils::globalVariables(c(".Random.seed"))
