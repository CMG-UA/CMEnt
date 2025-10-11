#!/usr/bin/env Rscript
#' Command Line Interface for Finding DMRs from DMPs
#'
#' @param args Argument list from optparse containing command line parameters
#' @return None, outputs results to files
#' @export
findDMRsFromSeedsCLI <- function(args) {
    pheno <- .processSamplesheet(args)$samplesheet
    options("DMRSegal.verbose" = args$verbose)
    .log_info(
            "Head of parsed phenotype:\n\t",
            paste(capture.output(print(head(pheno))), collapse = "\n\t")
        )

    if (!is.null(args$beta) && !is.null(args$tabix)) {
        stop("Either beta or tabix must be provided, not both.")
    }
    if (is.null(args$beta) && is.null(args$tabix)) {
        stop("Either beta or tabix must be provided.")
    }

    # Prepare arguments for findDMRsFromSeeds
    input_args <- list(
        beta_file = args$beta,
        tabix_file = args$tabix,
        dmps_file = args$dmps_file,
        pval_col = args$pval_col,
        sample_group_col = args$sample_group_col,
        dmp_group_col = args$dmp_group_col,
        casecontrol_col = args$casecontrol_col,
        min_cpg_delta_beta = args$min_cpg_delta_beta,
        expansion_step = args$expansion_step,
        expansion_relaxation = args$expansion_relaxation,
        array = args$array,
        genome = args$genome,
        max_pval = args$max_pval,
        max_lookup_dist = args$max_lookup_dist,
        min_dmps = args$min_dmps,
        min_adj_dmps = args$min_adj_dmps,
        min_cpgs = args$min_cpgs,
        ignored_sample_groups = args$ignored_sample_groups,
        output_prefix = args$output_prefix,
        njobs = args$njobs,
        verbose = args$verbose,
        beta_row_names_file = args$beta_row_names_file,
        pheno = pheno
    )

    .log_info(
            "Input arguments (without the phenotype): \n\t",
            paste(paste(names(input_args), input_args, sep = ": "), collapse = "\n\t")
        )
    

    invisible(do.call(findDMRsFromSeeds, input_args))
}
