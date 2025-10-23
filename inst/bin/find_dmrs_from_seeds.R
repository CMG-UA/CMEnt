#!/usr/bin/env Rscript

library(optparse)
library(DMRsegal)
option_list <- list(
    make_option(c("--beta"), default = NULL, help = "The beta file, with row names the CpGs. Required if tabix not provided"), 
    make_option(c("--tabix"), default = NULL, help = "The tabix bed.gz file, with the corresponding index in the same directory. Required if beta not provided"), 
    make_option("--dmps_file", help = "The dmps tsv file, with row names the DMPs."),
    make_option("--pval_col", default = "pval_adj", help = "The p-value column in the dmps tsv file, defaults to 'pval_adj'"),
    make_option("--min_dmps", default = 1, type = "integer", help = "The minimum supporting DMPs per DMR, defaults to 1"),
    make_option("--min_adj_dmps", default = 1, type = "integer", help = "The minimum adjacent DMPs per DMR, defaults to 1"),
    make_option("--min_cpgs", default = 50, type = "integer", help = "The minimum CpGs per DMR, defaults to 50"),
    make_option("--max_lookup_dist", default = 10000, type = "integer", help = "The maximum distance from one DMP to another to consider belinging in the same DMR, defaults to 10000 (10kb)"),
    make_option("--min_cpg_delta_beta", default = 0, type = "double", help = "The minimum CpG delta beta during DMR expansion, to filter CpGs based on delta beta, optional."),
    make_option("--ignored_sample_groups", default = NULL, help = "The sample groups to ignore when lookig for DMP connectivity and CpG extension, comma separated"),
    make_option("--expansion_step", default = 500, type = "integer", help = "The DMR expansion step in bp, defaults to 500"),
    make_option("--dmp_group_col", default = NULL, help = "Column in DMPs file for grouping DMPs"),
    make_option("--dmp_groups_tsv", default = NULL, help = "TSV file for DMP groups information, required if dmp_group_col is given, two columns: group_id<tab>sample1,sample2,..."),
    make_option("--casecontrol_col", default = "casecontrol", help = "Column in pheno for case/control status"),
    make_option("--array", default = "450K", help = "Array platform: 450K or EPIC"),
    make_option("--genome", default = "hg19", help = "Reference genome: hg19 or hg38"),
    make_option("--beta_row_names_file", default = NULL),
    make_option("--max_pval", default = 0.05, type = "double", help = "The maximum p-value to consider positional entanglement, defaults to 0.05"),
    make_option("--output_prefix", default = NULL, help = "Optional prefix for output files"),
    make_option("--sample_group_col", default = "Sample_Group", help = "Sample group column in samplesheet"),
    make_option("--output_id", default = "dmrs", help = "Output file prefix"),
    make_option("--njobs", default = 1, type = "integer", help = "Number of parallel jobs"),
    make_option("--verbose", default = 1, type = "integer", help = "Level of verbosity for logging messages, from 0 (not verbose) to 3 (very verbose). Default is 1"),
    make_option("--samplesheet", default = NULL, help = "Samplesheet file"),
    make_option("--target_col", default = NULL, help = "Target column in samplesheet"),
    make_option("--sample_group_control", default = NULL, help = "Control group names"),
    make_option("--sample_group_case", default = NULL, help = "Case group names"),
    make_option("--max_missing_cov_ratio", default = 0.3, type = "double", help = "Maximum missing coverage ratio"),
    make_option("--max_samples", default = NULL, type = "integer", help = "Maximum number of samples"),
    make_option("--max_case_samples", default = NULL, type = "integer", help = "Maximum number of case samples"),
    make_option("--max_control_samples", default = NULL, type = "integer", help = "Maximum number of control samples")
)
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)
DMRsegal::findDMRsFromSeedsCLI(args)
