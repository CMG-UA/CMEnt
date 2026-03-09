#!/usr/bin/env Rscript

library(optparse)
library(DMRsegal)
option_list <- list(
    make_option("--beta", help = "The beta file, with row names the CpGs. Can also be a tabix indexed file or a bed file with at least `bed_chrom_col` and `bed_start_col` columns set, followed by samples with methylation values."),
    make_option("--seeds_file", help = "The seeds tsv file, with row names the seeds."),
    make_option("--samplesheet", help = "The samplesheet (.csv or .tsv) containing samples information"),
    make_option("--target_col", help = "The column in the samplesheet corresponding to the sample name", default = "Sample_Name"),
    make_option("--sample_group_col", help = "The column in the samplesheet corresponding to the sample group", default = "Sample_Group"),
    make_option("--sample_group_control", help = "Comma-separated values in the sample group column that correspond to the controls, optional", default = NULL),
    make_option("--sample_group_case", help = "Comma-separated values in the sample group column that correspond to the cases, optional", default = NULL),
    make_option("--casecontrol_col", default = NULL, help = "Column in pheno for case/control status"),
    make_option("--covariates", default = NULL, help = "Comma-separated columns in samplesheet to be considered as confounders during the analysis"),
    make_option("--min_seeds", default = 2, type = "integer", help = "The minimum supporting seeds per DMR, minimum 2, defaults to 2"),
    make_option("--min_adj_seeds", default = 2, type = "integer", help = "The minimum supporting seeds per DMR, after adjusted by array backround CpG content, minimum 2,defaults to 2"),
    make_option("--min_cpgs", default = 50, type = "integer", help = "The minimum number of the beta file rows (the listed CpGs) per DMR, minimum 2, defaults to 50"),
    make_option("--max_lookup_dist", default = 10000, type = "integer", help = "The maximum distance from one seed to another to consider belinging in the same DMR, defaults to 10000 (10kb)"),
    make_option("--min_cpg_delta_beta", default = 0.1, type = "double", help = "The minimum CpG delta beta during DMR expansion, to filter CpGs based on delta beta (default: 0.1)."),
    make_option("--adaptive_min_cpg_delta_beta", default = TRUE, type = "logical", help = "Whether to adaptively increase min_cpg_delta_beta from the seed-level delta-beta distribution (default: TRUE)."),
    make_option("--ignored_sample_groups", default = NULL, help = "The sample groups to ignore while considering connection and expansion, comma separated. Can also be 'case' or 'control'."),
    make_option("--expansion_step", default = 500, type = "integer", help = "The index-specific DMR expansion step, defaults to 500"),
    make_option("--array", default = "450K", help = "Array platform ('450K', '27K', 'EPIC', 'EPICv2', 'Mouse', 'NULL'). Must be 'NULL' if not applicable"),
    make_option("--genome", default = "hg38", help = "Reference genome identifier. For an array-based experiment only hg38, hg19, mm10 and mm39 are supported."),
    make_option("--beta_row_names_file", default = NULL),
    make_option("--max_pval", default = 0.05, type = "double", help = "The maximum p-value to consider positional entanglement, defaults to 0.05"),
    make_option(
        "--expansion_window", default = "auto",
        help = "Stage 2 connectivity is computed only in windows centered on seed-derived Stage 1 DMR neighborhoods, with this total window width in bp. Set <=0 for genome-wide connectivity. Default is -1 for microarrays and 10000 (10 kb) for NGS datasets."
    ),
    make_option("--max_bridge_seeds_gaps", default = 1, type = "integer", help = "Maximum number of consecutive failed seed edges to bridge during Stage 1 when p-value-driven and flanked by connected edges (default: 1)."),
    make_option("--max_bridge_extension_gaps", default = 1, type = "integer", help = "Maximum number of consecutive failed CpG edges to extend during Stage 2 when p-value-driven and flanked by connected edges (default: 1)."),
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
    make_option("--entanglement", default = "strong", help = "Strategy for testing connectivity between groups: 'strong' (all groups must pass) or 'weak' (at least one group must pass). Default is 'strong'")
)
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)
DMRsegal::findDMRsFromSeedsCLI(args)
