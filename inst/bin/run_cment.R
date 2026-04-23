#!/usr/bin/env Rscript

if (!requireNamespace("CMEnt", quietly = TRUE)) {
    stop("Package 'CMEnt' must be installed before using this CLI.", call. = FALSE)
}

CMEnt:::.runCMEntCLI(
    cli_args = commandArgs(trailingOnly = TRUE),
    script_name = commandArgs(FALSE)[1]
)
