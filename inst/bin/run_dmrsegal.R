#!/usr/bin/env Rscript

if (!requireNamespace("DMRsegal", quietly = TRUE)) {
    stop("Package 'DMRsegal' must be installed before using this CLI.", call. = FALSE)
}

DMRsegal:::.runDMRsegalCLI(
    cli_args = commandArgs(trailingOnly = TRUE),
    script_name = commandArgs(FALSE)[1]
)
