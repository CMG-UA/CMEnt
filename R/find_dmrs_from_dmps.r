#!/usr/bin/env Rscript
stub <- function() {}
thisPath <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    if (length(grep("^-f$", cmdArgs)) > 0) {
        # R console option
        normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
    } else if (length(grep("^--file=", cmdArgs)) > 0) {
        # Rscript/R console option
        normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
    } else if (is.null(attr(stub, "srcref")) == FALSE) {
        # 'source'd via R console
        dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
    } else if (Sys.getenv("RSTUDIO") == "1" && !is.null(rstudioapi::getSourceEditorContext()) && (rstudioapi::getSourceEditorContext()$path != "")) {
        # RStudio
        doc <- rstudioapi::getSourceEditorContext()
        dirname(doc$path)
    } else {
        stop("Cannot find file path")
    }
}

BIN.PATH <- thisPath()
LIB.PATH <- file.path(dirname(BIN.PATH), "lib")
DATA.PATH <- file.path(dirname(BIN.PATH), "data")
if (!file.exists(LIB.PATH)) {
    stop("Cannot find lib path: ", LIB.PATH)
}
if (!file.exists(DATA.PATH)) {
    stop("Cannot find data path: ", DATA.PATH)
}

suppressWarnings(try(source(file.path(LIB.PATH, "diff_meth_base.r")), silent = TRUE))
suppressWarnings(try(source(file.path(LIB.PATH, "dmrs_finder_from_dmps.r")), silent = TRUE))

parser <- getParser(not.add = c("beta", "cont-cov", "cat-cov", "cov-cols", "max-samples-per-case-control"))
parser$add_argument("--beta", help = "The beta file, with row names the CpGs. Required if tabix not provided.", default=NULL)
parser$add_argument("--tabix", help = "The tabix bed.gz file, with the corresponding index in the same directory. Required if beta not provided.", default=NULL)

parser$add_argument("--dmps-tsv-file", help = "The dmps tsv file, with row names the DMPs.")
parser$add_argument("--pval-col", help = "The p-value column in the dmps tsv file, defaults to 'pval_adj'.", default =
                      "pval_adj")

parser$add_argument("--min-dmps",
                    help = "The minimum supporting DMPs per DMR, defaults to 1",
                    default = 1,
                    type = 'numeric')
parser$add_argument("--min-cpgs",
                    help = "The minimum CpGs per DMR, defaults to 50",
                    default = 50,
                    type = 'numeric')
parser$add_argument(
  "--max-lookup-dist",
  help = "The maximum distance from one DMP to another to consider belinging in the same DMR, defaults to 10000 (10kb)",
  default = 10000,
  type = 'numeric'
)
parser$add_argument("--min-cpg-delta-beta", default=0, type='numeric', help="The minimum CpG delta beta during DMR expansion, to filter CpGs based on delta beta, optional.")
parser$add_argument("--ignored-sample-groups", default=NULL, help="The sample groups to ignore when lookig for DMP connectivity and CpG extension, comma separated. Include groups that contain highly heterogeneous samples. Optional.")
parser$add_argument("--expansion-step", default=5, type='integer', help="The DMR expansion step, affects speed/IO trade-off, defaults to 10 cpg sites. Have it small when using tabix file.")
parser$add_argument("--beta-row-names-file", default=NULL)
parser$add_argument("--dmps-beta-file", default=NULL)

parser$add_argument("--max-pval",
                    help = "The maximum p-value to consider positional entanglement, defaults to 0.05",
                    default = 0.05,
                    type = 'numeric')
main <- function(args){
  pheno <- .processSamplesheet(args, has.covariates = F)$samplesheet
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
    beta.row.names.file = args$beta_row_names_file
  )
  if (args$verbose){
    message("Input arguments (without the phenotype): \n\t", paste(paste(names(input.args), input.args, sep=': '), collapse='\n\t'))
    
  }
  input.args[["pheno"]] <- pheno
  invisible(
    do.call(findDMRsFromDMPs, input.args)
  )
}



if (sys.nframe() == 0L) {
  main(parser$parse_args())
}

