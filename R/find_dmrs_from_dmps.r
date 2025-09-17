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
    output.dir = args$output_dir,
    pval.col = args$pval_col,
    sample_group.col = args$sample_group_col,
    min.cpg.delta_beta = args$min_cpg_delta_beta,
    expansion.step = args$expansion_step,
    max.pval = args$max_pval,
    max.lookup.dist = args$max_lookup_dist,
    min.dmps = args$min_dmps,
    min.cpgs = args$min_cpgs,
    ignored.sample.groups = args$ignored_sample_groups,
    output.id = args$output_id,
    njobs = args$njobs,
    verbose = args$verbose,
    beta.row.names.file = args$beta_row_names_file,
    pheno = pheno
  )

  if (args$verbose){
    message("Input arguments (without the phenotype): \n\t", 
            paste(paste(names(input.args), input.args, sep=': '), collapse='\n\t'))
  }
  
  invisible(do.call(findDMRsFromDMPs, input.args))
}


