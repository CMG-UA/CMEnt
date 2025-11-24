#!/usr/bin/env Rscript
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
#'     genome = "hg19",
#'     output_prefix = "my_analysis",
#'     njobs = 4,
#'     verbose = 1,
#'     group_concordance_strategy = "strict"
#' )
#'
#' findDMRsFromSeedsCLI(args)
#'
#' @export
findDMRsFromSeedsCLI <- function(args) {
    options("DMRsegal.verbose" = args$verbose)

    if (!is.null(args$beta) && !is.null(args$tabix)) {
        stop("Either beta or tabix must be provided, not both.")
    }
    if (is.null(args$beta) && is.null(args$tabix)) {
        stop("Either beta or tabix must be provided.")
    }
    pheno <- .processSamplesheet(args)$samplesheet
    .log_info(
        "Head of parsed phenotype:\n\t",
        paste(capture.output(print(head(pheno))), collapse = "\n\t")
    )
    # Prepare arguments for findDMRsFromSeeds
    input_args <- list(
        beta = args$beta,
        seeds = args$seeds_file,
        sample_group_col = args$sample_group_col,
        casecontrol_col = args$casecontrol_col,
        min_cpg_delta_beta = args$min_cpg_delta_beta,
        expansion_step = args$expansion_step,
        array = args$array,
        genome = args$genome,
        max_pval = args$max_pval,
        max_lookup_dist = args$max_lookup_dist,
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
        group_concordance_strategy = args$group_concordance_strategy
    )

    .log_info(
        "Input arguments (without the phenotype): \n\t",
        paste(paste(names(input_args), input_args, sep = ": "), collapse = "\n\t")
    )


    invisible(do.call(findDMRsFromSeeds, input_args))
}
