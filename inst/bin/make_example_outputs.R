#!/usr/bin/env Rscript
devtools::load_all()
cases <- c("example_output", "example_outputChr5And11")
files_to_make <- file.path("inst/extdata", paste0(cases, ".rds"))
names(files_to_make) <- cases
# files_to_make <- files_to_make[!file.exists(files_to_make)]
# if (length(files_to_make) == 0) {
#     cat("All example output files already exist. Skipping generation.\n")
#     quit(save = "no")
# }
progressr::handlers("cli")
options("DMRsegal.verbose" = 1)
options("DMRsegal.njobs" = 8)
for (case in names(files_to_make)) {
    dmrsegal_file <- files_to_make[[case]]
    if (case == "example_outputChr5And11") {
        fun <- loadExampleInputDataChr5And11
    } else {
        fun <- loadExampleInputData
    }
    beta <- fun("beta")
    dmps <- fun("dmps")
    pheno <- fun("pheno")
    pheno[, "casecontrol"] <- pheno$Sample_Group == "cancer"
    dmrs_dmrsegal <- progressr::with_progress(DMRsegal::findDMRsFromSeeds(
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        casecontrol_col = "casecontrol",
        covariates = c("Age", "Gender"),
        max_bridge_seeds_gaps  = 1,
        max_bridge_extension_gaps = 1,
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        max_pval = 0.05,
        pval_mode = "parametric",
        entanglement = "weak",
        rank_dmrs = TRUE
    ))
    saveRDS(dmrs_dmrsegal, dmrsegal_file)
}
