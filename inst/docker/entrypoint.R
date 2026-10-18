#!/usr/bin/env Rscript

script_file <- sub("^--file=", "", grep(
    "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
))
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "run_cment.R"), local = TRUE)
source(file.path(script_dir, "install_dependencies.R"), local = TRUE)

.runCMEntCLI(
    cli_args = commandArgs(trailingOnly = TRUE),
    script_name = "run_cment.R",
    prepare_dependencies = .prepareCMEntContainer
)
