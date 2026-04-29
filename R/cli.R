#!/usr/bin/env Rscript

.findDMRsFromSeedsCLIOptionList <- function() {
    if (!requireNamespace("optparse", quietly = TRUE)) {
        stop("Package 'optparse' is required for CMEnt CLI support.", call. = FALSE)
    }

    list(
        optparse::make_option(
            "--beta",
            help = "The beta file, with row names the CpGs. Can also be a tabix indexed file or a bed file with at least `bed_chrom_col` and `bed_start_col` columns set, followed by samples with methylation values."
        ),
        optparse::make_option(
            "--seeds_file",
            help = "The seeds tsv file, with row names the seeds."
        ),
        optparse::make_option(
            "--samplesheet",
            help = "The samplesheet (.csv or .tsv) containing samples information"
        ),
        optparse::make_option(
            "--target_col",
            help = "The column in the samplesheet corresponding to the sample name",
            default = "Sample_Name"
        ),
        optparse::make_option(
            "--sample_group_col",
            help = "The column in the samplesheet corresponding to the sample group",
            default = "Sample_Group"
        ),
        optparse::make_option(
            "--sample_group_control",
            help = "Comma-separated values in the sample group column that correspond to the controls, optional",
            default = NULL
        ),
        optparse::make_option(
            "--sample_group_case",
            help = "Comma-separated values in the sample group column that correspond to the cases, optional",
            default = NULL
        ),
        optparse::make_option(
            "--casecontrol_col",
            default = NULL,
            help = "Column in pheno for case/control status"
        ),
        optparse::make_option(
            "--covariates",
            default = NULL,
            help = "Comma-separated columns in samplesheet to be considered as confounders during the analysis"
        ),
        optparse::make_option(
            "--min_seeds",
            default = 2,
            type = "integer",
            help = "The minimum supporting seeds per DMR, minimum 2, defaults to 2"
        ),
        optparse::make_option(
            "--min_adj_seeds",
            default = 2,
            type = "integer",
            help = "The minimum supporting seeds per DMR, after adjusted by array backround CpG content, minimum 2,defaults to 2"
        ),
        optparse::make_option(
            "--min_cpgs",
            default = 50,
            type = "integer",
            help = "The minimum number of the beta file rows (the listed CpGs) per DMR, minimum 2, defaults to 50"
        ),
        optparse::make_option(
            "--max_lookup_dist",
            default = 10000,
            type = "integer",
            help = "The maximum distance from one seed to another to consider belinging in the same DMR, defaults to 10000 (10kb)"
        ),
        optparse::make_option(
            "--min_cpg_delta_beta",
            default = 0.1,
            type = "double",
            help = "The minimum CpG delta beta during DMR expansion, to filter CpGs based on delta beta (default: 0.1)."
        ),
        optparse::make_option(
            "--ignored_sample_groups",
            default = NULL,
            help = "The sample groups to ignore while considering connection and expansion, comma separated. Can also be 'case' or 'control'."
        ),
        optparse::make_option(
            "--array",
            default = "450K",
            help = "Array platform ('450K', '27K', 'EPIC', 'EPICv2', 'Mouse', 'NULL'). Must be 'NULL' if not applicable"
        ),
        optparse::make_option(
            "--genome",
            default = NULL,
            help = "Reference genome identifier. If omitted, CMEnt infers hg19 for 450K, 27K, and EPIC arrays, otherwise hg38."
        ),
        optparse::make_option("--beta_row_names_file", default = NULL),
        optparse::make_option(
            "--max_pval",
            default = 0.05,
            type = "double",
            help = "The maximum p-value to consider positional entanglement, defaults to 0.05"
        ),
        optparse::make_option(
            "--expansion_window",
            default = "auto",
            help = "Stage 2 connectivity is computed only in windows centered on seed-derived Stage 1 DMR neighborhoods, with this total window width in bp. Set <=0 for genome-wide connectivity. Default is -1 for microarrays and 10000 (10 kb) for NGS datasets."
        ),
        optparse::make_option(
            "--max_bridge_seeds_gaps",
            default = 1,
            type = "integer",
            help = "Maximum number of consecutive failed seed edges to bridge during Stage 1 when p-value-driven and flanked by connected edges (default: 1)."
        ),
        optparse::make_option(
            "--max_bridge_extension_gaps",
            default = 1,
            type = "integer",
            help = "Maximum number of consecutive failed CpG edges to extend during Stage 2 when p-value-driven and flanked by connected edges (default: 1)."
        ),
        optparse::make_option(
            "--output_prefix",
            default = NULL,
            help = "Optional prefix for output files"
        ),
        optparse::make_option(
            "--njobs",
            default = 1,
            type = "integer",
            help = "Number of parallel jobs"
        ),
        optparse::make_option(
            "--verbose",
            default = 1,
            type = "integer",
            help = "Level of verbosity for logging messages, from 0 (not verbose) to 3 (very verbose). Default is 1"
        ),
        optparse::make_option(
            "--annotate_with_genes",
            default = TRUE,
            type = "logical",
            help = "Whether to annotate DMRs with gene information (default: TRUE)"
        ),
        optparse::make_option(
            "--bed_provided",
            default = FALSE,
            type = "logical",
            help = "Whether BED file is provided"
        ),
        optparse::make_option(
            "--bed_chrom_col",
            default = "#chrom",
            help = "Column in BED file for chromosome"
        ),
        optparse::make_option(
            "--bed_start_col",
            default = "start",
            help = "Column in BED file for start position"
        ),
        optparse::make_option(
            "--entanglement",
            default = "strong",
            help = "Strategy for testing connectivity between groups: 'strong' (all groups must pass) or 'weak' (at least one group must pass). Default is 'strong'"
        )
    )
}

.launchCMEntViewerCLIOptionList <- function() {
    if (!requireNamespace("optparse", quietly = TRUE)) {
        stop("Package 'optparse' is required for CMEnt CLI support.", call. = FALSE)
    }

    list(
        optparse::make_option(
            "--output_prefix",
            help = "Prefix used when saving DMR analysis results."
        ),
        optparse::make_option(
            "--launch_browser",
            default = TRUE,
            type = "logical",
            help = "Whether to launch the viewer in a browser (default: TRUE)."
        ),
        optparse::make_option(
            "--port",
            default = NULL,
            type = "integer",
            help = "Port number for the Shiny server (default: auto-assigned)."
        ),
        optparse::make_option(
            "--host",
            default = "127.0.0.1",
            help = "Host interface for the Shiny server (default: 127.0.0.1). Use 0.0.0.0 in Docker."
        ),
        optparse::make_option(
            "--diagnostic",
            default = FALSE,
            type = "logical",
            help = "Whether to enable diagnostic viewer panels (default: FALSE)."
        )
    )
}

#' Command Line Interface for Finding DMRs from DMPs
#'
#' @param args Argument list from optparse containing command line parameters
#' @return None, outputs results to files
#'
#' @examples
#' if (!requireNamespace("optparse", quietly = TRUE)) {
#'     install.packages("optparse")
#' }
#' library(optparse)
#'
#' args <- list(
#'     beta = "path/to/beta_values.rds",
#'     seeds_file = "path/to/dmps.tsv",
#'     samplesheet = "path/to/samplesheet.csv",
#'     sample_group_col = "Sample_Group",
#'     casecontrol_col = "casecontrol",
#'     array = "450K",
#'     genome = NULL,
#'     output_prefix = "my_analysis",
#'     njobs = 4,
#'     verbose = 1,
#'     entanglement = "strong"
#' )
#'
#' findDMRsFromSeedsCLI(args)
#'
#' @export
findDMRsFromSeedsCLI <- function(args) {
    options("CMEnt.verbose" = args$verbose)
    pheno <- .processSamplesheet(args)$samplesheet
    array <- args$array
    covariates <- args$covariates
    if (!is.null(covariates)) {
        covariates <- base::strsplit(covariates, ",")[[1]]
    }
    if (!is.null(array) && tolower(array) == "null") {
        array <- NULL
    }
    genome <- args$genome
    if (!is.null(genome) && tolower(genome) == "null") {
        genome <- NULL
    }

    input_args <- list(
        beta = args$beta,
        seeds = args$seeds_file,
        sample_group_col = args$sample_group_col,
        casecontrol_col = args$casecontrol_col,
        covariates = covariates,
        min_cpg_delta_beta = args$min_cpg_delta_beta,
        array = array,
        genome = genome,
        max_pval = args$max_pval,
        max_lookup_dist = args$max_lookup_dist,
        expansion_window = args$expansion_window,
        max_bridge_seeds_gaps = args$max_bridge_seeds_gaps,
        max_bridge_extension_gaps = args$max_bridge_extension_gaps,
        min_seeds = args$min_seeds,
        min_adj_seeds = args$min_adj_seeds,
        min_cpgs = args$min_cpgs,
        ignored_sample_groups = args$ignored_sample_groups,
        output_prefix = args$output_prefix,
        njobs = args$njobs,
        beta_row_names_file = args$beta_row_names_file,
        annotate_with_genes = args$annotate_with_genes,
        pheno = pheno,
        bed_provided = args$bed_provided,
        bed_chrom_col = args$bed_chrom_col,
        bed_start_col = args$bed_start_col,
        entanglement = args$entanglement
    )

    .log_info(
        "Input arguments (without the phenotype): \n\t",
        paste(paste(names(input_args), input_args, sep = ": "), collapse = "\n\t")
    )

    invisible(do.call(findDMRsFromSeeds, input_args))
}

launchCMEntViewerCLI <- function(args) {
    invisible(launchCMEntViewer(
        output_prefix = args$output_prefix,
        launch_browser = args$launch_browser,
        port = args$port,
        host = args$host,
        diagnostic = args$diagnostic
    ))
}

.cmentCLICommandAliases <- function() {
    c(
        findDMRsFromSeeds = "findDMRsFromSeeds",
        find = "findDMRsFromSeeds",
        find_dmrs_from_seeds = "findDMRsFromSeeds",
        launchCMEntViewer = "launchCMEntViewer",
        viewer = "launchCMEntViewer",
        launch_cment_viewer = "launchCMEntViewer"
    )
}

.normalizeCMEntCLICommand <- function(command) {
    if (is.null(command) || length(command) != 1 || !nzchar(command)) {
        return(NULL)
    }

    aliases <- .cmentCLICommandAliases()
    if (!command %in% names(aliases)) {
        return(NULL)
    }

    unname(aliases[[command]])
}

.topLevelCMEntCLIHelp <- function(script_name = NULL) {
    if (is.null(script_name) || !nzchar(script_name)) {
        script_name <- "run_cment.R"
    } else {
        script_name <- basename(script_name)
    }
    if (script_name == "Rscript" || script_name == "R") {
        script_name <- ""
    }

    paste(
        "CMEnt command line interface",
        "",
        "Usage:",
        paste0("  ", script_name, " <command> [options]"),
        "",
        "Commands:",
        "  findDMRsFromSeeds      Identify DMRs from genomic seeds.",
        "  launchCMEntViewer   Launch the interactive viewer for saved outputs.",
        "",
        "Compatibility:",
        "  Passing options without a command defaults to findDMRsFromSeeds.",
        "",
        "Examples:",
        paste0("  ", script_name, " findDMRsFromSeeds --beta beta.tsv.gz --seeds_file seeds.tsv --samplesheet samplesheet.tsv"),
        paste0("  ", script_name, " launchCMEntViewer --output_prefix results/my_analysis"),
        "",
        paste0("Run '", script_name, " help <command>' for command-specific options."),
        sep = "\n"
    )
}

.resolveCMEntCLIInvocation <- function(cli_args) {
    cli_args <- as.character(cli_args)

    if (length(cli_args) == 0) {
        return(list(command = NULL, command_args = character(), top_level_help = TRUE))
    }

    first_arg <- cli_args[[1]]

    if (first_arg %in% c("-h", "--help")) {
        return(list(command = NULL, command_args = character(), top_level_help = TRUE))
    }

    if (identical(first_arg, "help")) {
        if (length(cli_args) == 1) {
            return(list(command = NULL, command_args = character(), top_level_help = TRUE))
        }

        command <- .normalizeCMEntCLICommand(cli_args[[2]])
        if (is.null(command)) {
            stop(
                "Unknown CMEnt command '", cli_args[[2]], "'.\n\n",
                .topLevelCMEntCLIHelp(),
                call. = FALSE
            )
        }

        return(list(
            command = command,
            command_args = "--help",
            top_level_help = FALSE
        ))
    }

    if (startsWith(first_arg, "-")) {
        return(list(
            command = "findDMRsFromSeeds",
            command_args = cli_args,
            top_level_help = FALSE
        ))
    }

    command <- .normalizeCMEntCLICommand(first_arg)
    if (is.null(command)) {
        stop(
            "Unknown CMEnt command '", first_arg, "'.\n\n",
            .topLevelCMEntCLIHelp(),
            call. = FALSE
        )
    }

    list(
        command = command,
        command_args = cli_args[-1],
        top_level_help = FALSE
    )
}

.makeCMEntCLIParser <- function(command, script_name = NULL) {
    if (!requireNamespace("optparse", quietly = TRUE)) {
        stop("Package 'optparse' is required for CMEnt CLI support.", call. = FALSE)
    }

    if (is.null(script_name) || !nzchar(script_name)) {
        script_name <- "run_cment.R"
    } else {
        script_name <- basename(script_name)
    }

    switch(
        command,
        findDMRsFromSeeds = optparse::OptionParser(
            usage = "usage: %prog findDMRsFromSeeds [options]",
            prog = script_name,
            description = "Identify differentially methylated regions from genomic seeds.",
            option_list = .findDMRsFromSeedsCLIOptionList()
        ),
        launchCMEntViewer = optparse::OptionParser(
            usage = "usage: %prog launchCMEntViewer [options]",
            prog = script_name,
            description = "Launch the CMEnt interactive viewer for a saved output prefix.",
            option_list = .launchCMEntViewerCLIOptionList()
        ),
        stop("Unsupported CMEnt command '", command, "'.", call. = FALSE)
    )
}

.isMissingCLIValue <- function(x) {
    is.null(x) || length(x) == 0 || (is.character(x) && !nzchar(x))
}

.requireCLIOptions <- function(args, options) {
    missing_options <- options[vapply(options, function(option_name) {
        .isMissingCLIValue(args[[option_name]])
    }, logical(1))]

    if (length(missing_options) > 0) {
        stop(
            "Missing required option(s): ",
            paste(paste0("--", missing_options), collapse = ", "),
            call. = FALSE
        )
    }
}

.runCMEntCLI <- function(cli_args = commandArgs(trailingOnly = TRUE), script_name = NULL) {
    invocation <- .resolveCMEntCLIInvocation(cli_args)

    if (isTRUE(invocation$top_level_help)) {
        cat(.topLevelCMEntCLIHelp(script_name), "\n")
        return(invisible(0L))
    }

    parser <- .makeCMEntCLIParser(invocation$command, script_name = script_name)
    parsed_args <- optparse::parse_args(
        parser,
        args = invocation$command_args,
        print_help_and_exit = FALSE
    )

    if (isTRUE(parsed_args$help)) {
        optparse::print_help(parser)
        return(invisible(0L))
    }

    switch(
        invocation$command,
        findDMRsFromSeeds = {
            .requireCLIOptions(parsed_args, c("beta", "seeds_file", "samplesheet"))
            invisible(findDMRsFromSeedsCLI(parsed_args))
        },
        launchCMEntViewer = {
            .requireCLIOptions(parsed_args, "output_prefix")
            invisible(launchCMEntViewerCLI(parsed_args))
        }
    )
}
