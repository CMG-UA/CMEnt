#!/usr/bin/env Rscript

#' Install comb-p Python Package
#'
#' This script helps install comb-p, a Python package for combining p-values
#' with autocorrelation correction, which is used in the benchmarking vignette.
#'
#' @examples
#' \dontrun{
#' source("inst/scripts/install_combp.R")
#' install_combp()
#' }

install_combp <- function(method = "conda", conda = "auto") {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        message("Installing reticulate package...")
        install.packages("reticulate")
    }
    
    library(reticulate)
    
    tryCatch({
        if (py_module_available("comb_p")) {
            message("comb-p is already installed!")
            return(invisible(TRUE))
        }
        
        message("Installing comb-p via conda...")
        message("Note: comb-p requires conda and is only available from bioconda")
        conda_install("combined-pvalues", channel = "bioconda")
        
        if (py_module_available("comb_p")) {
            message("Successfully installed comb-p!")
            return(invisible(TRUE))
        } else {
            warning("comb-p installation completed but module not detected.")
            return(invisible(FALSE))
        }
    }, error = function(e) {
        warning(sprintf("Failed to install comb-p: %s", e$message))
        message("\nManual installation instructions:")
        message("Using conda: conda install -yc bioconda combined-pvalues")
        message("\nNote: comb-p is ONLY available via conda, not pip")
        message("For more information, visit: https://github.com/brentp/combined-pvalues")
        return(invisible(FALSE))
    })
}

verify_combp <- function() {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        message("reticulate package not installed")
        return(FALSE)
    }
    
    library(reticulate)
    
    if (py_module_available("comb_p")) {
        message("comb-p is installed and available!")
        
        tryCatch({
            system("comb-p -h", ignore.stdout = TRUE, ignore.stderr = TRUE)
            message("comb-p command-line tool is also available!")
            return(TRUE)
        }, error = function(e) {
            message("comb-p Python module is available but command-line tool may not be in PATH")
            return(TRUE)
        })
    } else {
        message("comb-p is not installed")
        return(FALSE)
    }
}

if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    
    if (length(args) > 0 && args[1] == "verify") {
        verify_combp()
    } else {
        install_combp()
    }
}
