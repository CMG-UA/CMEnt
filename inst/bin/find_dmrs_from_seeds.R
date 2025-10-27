#!/usr/bin/env Rscript

library(optparse)
library(DMRsegal)
option_list <- list(
    make_option("--beta", help = "The beta file, with row names the CpGs. Can also be a tabix indexed file or a bed file with at least `bed_chrom_col` and `bed_start_col` columns set, followed by samples with methylation values."),
    make_option("--dmps_file", help = "The dmps tsv file, with row names the DMPs."),
    make_option("--min_dmps", default = 1, type = "integer", help = "The minimum supporting DMPs per DMR, defaults to 1"),
    make_option("--min_adj_dmps", default = 1, type = "integer", help = "The minimum supporting DMPs per DMR, after adjusted by underlying CpG content, defaults to 1"),
    make_option("--min_cpgs", default = 50, type = "integer", help = "The minimum number of the beta file rows (the listed CpGs) per DMR, defaults to 50"),
    make_option("--max_lookup_dist", default = 10000, type = "integer", help = "The maximum distance from one DMP to another to consider belinging in the same DMR, defaults to 10000 (10kb)"),
    make_option("--min_cpg_delta_beta", default = 0, type = "double", help = "The minimum CpG delta beta during DMR expansion, to filter CpGs based on delta beta, optional."),
    make_option("--ignored_sample_groups", default = NULL, help = "The sample groups to ignore when lookig for DMP connectivity and CpG extension, comma separated"),
    make_option("--expansion_step", default = 500, type = "integer", help = "The DMR expansion step in bp, defaults to 500"),
    make_option("--casecontrol_col", default = "casecontrol", help = "Column in pheno for case/control status"),
    make_option("--array", default = "450K", help = "Array platform: 450K or EPIC"),
    make_option("--genome", default = "hg19", help = "Reference genome: hg19 or hg38"),
    make_option("--beta_row_names_file", default = NULL),
    make_option("--max_pval", default = 0.05, type = "double", help = "The maximum p-value to consider positional entanglement, defaults to 0.05"),
    make_option("--output_prefix", default = NULL, help = "Optional prefix for output files"),
    make_option("--sample_group_col", default = "Sample_Group", help = "Sample group column in samplesheet"),
    make_option("--njobs", default = 1, type = "integer", help = "Number of parallel jobs"),
    make_option("--verbose", default = 1, type = "integer", help = "Level of verbosity for logging messages, from 0 (not verbose) to 3 (very verbose). Default is 1"),
    make_option("--annotate_with_genes", default = TRUE, type = "logical", help = "Whether to annotate DMRs with gene information (default: TRUE)"),
    make_option("--samplesheet", default = NULL, help = "Samplesheet file"),
    make_option("--target_col", default = NULL, help = "Target column in samplesheet"),
    make_option("--sample_group_control", default = NULL, help = "Control group names"),
    make_option("--sample_group_case", default = NULL, help = "Case group names"),
    make_option("--bed_provided", default = FALSE, type = "logical", help = "Whether BED file is provided"),
    make_option("--bed_chrom_col", default = "#chrom", help = "Column in BED file for chromosome"),
    make_option("--bed_start_col", default = "start", help = "Column in BED file for start position"),
    make_option("--bed_end_col", default = "end", help = "Column in BED file for end position. If NULL, will be assumed it does not exist"),
    make_option("--bed_id_col", default = "name", help = "Column in BED file for region ID. If NULL, will be assumed it does not exist"),
    make_option("--bed_score_col", default = "score", help = "Column in BED file for region score. If NULL, will be assumed it does not exist"),
    make_option("--bed_strand_col", default = "strand", help = "Column in BED file for strand. If NULL, will be assumed it does not exist")
)
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)
DMRsegal::findDMRsFromSeedsCLI(args)
