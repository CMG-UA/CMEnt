#!/usr/bin/env Rscript

library(argparse)
library(DMRSegal)

parser <- argparse::ArgumentParser()
parser$add_argument("--beta", help = "The beta file, with row names the CpGs. Required if tabix not provided.", default=NULL)
parser$add_argument("--tabix", help = "The tabix bed.gz file, with the corresponding index in the same directory. Required if beta not provided.", default=NULL)
parser$add_argument("--dmps-tsv-file", help = "The dmps tsv file, with row names the DMPs.")
parser$add_argument("--pval-col", help = "The p-value column in the dmps tsv file, defaults to 'pval_adj'.", default = "pval_adj")
parser$add_argument("--min-dmps", help = "The minimum supporting DMPs per DMR, defaults to 1", default = 1, type = 'numeric')
parser$add_argument("--min-cpgs", help = "The minimum CpGs per DMR, defaults to 50", default = 50, type = 'numeric')
parser$add_argument("--max-lookup-dist", help = "The maximum distance from one DMP to another to consider belinging in the same DMR, defaults to 10000 (10kb)", default = 10000, type = 'numeric')
parser$add_argument("--min-cpg-delta-beta", default=0, type='numeric', help="The minimum CpG delta beta during DMR expansion, to filter CpGs based on delta beta, optional.")
parser$add_argument("--ignored-sample-groups", default=NULL, help="The sample groups to ignore when lookig for DMP connectivity and CpG extension, comma separated.")
parser$add_argument("--expansion-step", default=5, type='integer', help="The DMR expansion step, affects speed/IO trade-off, defaults to 10 cpg sites.")
parser$add_argument("--beta-row-names-file", default=NULL)
parser$add_argument("--dmps-beta-file", default=NULL)
parser$add_argument("--max-pval", help = "The maximum p-value to consider positional entanglement, defaults to 0.05", default = 0.05, type = 'numeric')

args <- parser$parse_args()
DMRSegal::findDMRsFromDMPsCLI(args)
