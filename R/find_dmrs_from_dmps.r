#!/usr/bin/env Rscript
#' Command Line Interface for Finding DMRs from DMPs
#'
#' @param args Argument list from optparse containing command line parameters
#' @return None, outputs results to files
#' @export
findDMRsFromDMPsCLI <- function(args) {
  pheno <- .processSamplesheet(args)$samplesheet
  
  if (args$verbose){
    message("Head of parsed phenotype:\n\t", 
            paste(capture.output(print(head(pheno))), collapse='\n\t'))
  }

  if (!is.null(args$beta) && !is.null(args$tabix)){
    stop("Either beta or tabix must be provided, not both.")
  }
  if (is.null(args$beta) && is.null(args$tabix)){
    stop("Either beta or tabix must be provided.")
  }

  # Prepare arguments for findDMRsFromDMPs
  input.args <- list(
    beta.file = args$beta,
    tabix.file = args$tabix,
    dmps.tsv.file = args$dmps_tsv_file,
    pval.col = args$pval_col,
    sample_group.col = args$sample_group_col,
    dmp_group.col = args$dmp_group_col,
    casecontrol.col = args$casecontrol_col,
    min.cpg.delta_beta = args$min_cpg_delta_beta,
    expansion.step = args$expansion_step,
    expansion.relaxation = args$expansion_relaxation,
    array = args$array,
    genome = args$genome,
    max.pval = args$max_pval,
    max.lookup.dist = args$max_lookup_dist,
    min.dmps = args$min_dmps,
    min.adj.dmps = args$min_adj_dmps,
    min.cpgs = args$min_cpgs,
    ignored.sample.groups = args$ignored_sample_groups,
    output.prefix = args$output_prefix,
    njobs = args$njobs,
    verbose = args$verbose,
    beta.row.names.file = args$beta_row_names_file,
    dmps.beta.file = args$dmps_beta_file,
    pheno = pheno
  )

  if (args$verbose){
    message("Input arguments (without the phenotype): \n\t", 
            paste(paste(names(input.args), input.args, sep=': '), collapse='\n\t'))
  }
  
  invisible(do.call(findDMRsFromDMPs, input.args))
}


