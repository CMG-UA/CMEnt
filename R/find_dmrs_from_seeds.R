#!/usr/bin/env Rscript
#' Command Line Interface for Finding DMRs from DMPs
#'
#' @param args Argument list from optparse containing command line parameters
#' @return None, outputs results to files
#' @export
findDMRsFromSeedsCLI <- function(args) {
    options("DMRsegal.verbose" = args$verbose)

    if (!is.null(args$beta) && !is.null(args$tabix)) {
        stop("Either beta or tabix must be provided, not both.")
    }
    if (is.null(args$beta) && is.null(args$tabix)) {
        stop("Either beta or tabix must be provided.")
    }
    if (!is.null(args$dmp_group_col)) {
        if (is.null(args$dmp_groups_tsv)) {
            stop("dmp_groups_tsv must be provided when dmp_group_col is given.")
        }
        dmps_groups_info <- read.table(args$dmp_groups_tsv,
            header = FALSE,
            sep = "\t",
            stringsAsFactors = FALSE
        )
        samples <- strsplit(dmps_groups_info[, 2], ",")
        dmp_groups_info <- setNames(samples, dmps_groups_info[, 1])
        dmp_groups_info <- as.list(dmp_groups_info)
    }
    pheno <- .processSamplesheet(args)$samplesheet
    .log_info(
        "Head of parsed phenotype:\n\t",
        paste(capture.output(print(head(pheno))), collapse = "\n\t")
    )

    # Prepare arguments for findDMRsFromSeeds
    input_args <- list(
        beta = args$beta,
        tabix_file = args$tabix,
        dmps = args$dmps_file,
        pval_col = args$pval_col,
        sample_group_col = args$sample_group_col,
        dmp_group_col = args$dmp_group_col,
        dmp_groups_info = dmp_groups_info,
        casecontrol_col = args$casecontrol_col,
        min_cpg_delta_beta = args$min_cpg_delta_beta,
        expansion_step = args$expansion_step,
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
        annotate_with_genes = args$annotate_with_genes,
        pheno = pheno
    )

    .log_info(
        "Input arguments (without the phenotype): \n\t",
        paste(paste(names(input_args), input_args, sep = ": "), collapse = "\n\t")
    )


    invisible(do.call(findDMRsFromSeeds, input_args))
}
