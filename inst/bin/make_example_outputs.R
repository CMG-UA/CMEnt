#!/usr/bin/env Rscript
devtools::load_all()
cases <- c("example_outputChr5And11")
files_to_make <- file.path("inst/extdata", paste0(cases, ".rds"))
names(files_to_make) <- cases
# files_to_make <- files_to_make[!file.exists(files_to_make)]
if (length(files_to_make) == 0) {
    cat("All example output files already exist. Skipping generation.\n")
    quit(save = "no")
}
options("CMEnt.verbose" = 1)
options("CMEnt.njobs" = 8)
for (case in names(files_to_make)) {
    cment_file <- files_to_make[[case]]
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    pheno[, "casecontrol"] <- pheno$Sample_Group == "cancer"
    dmrs_cment <- CMEnt::findDMRsFromSeeds(
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        casecontrol_col = "casecontrol",
        covariates = c("Age", "Gender"),
        max_bridge_seeds_gaps  = 0,
        max_bridge_extension_gaps = 3,
        array = "450K",
        genome = "hg19",
        min_seeds = 2,
        min_sites = 3,
        max_lookup_dist = 1000,
        max_pval = 0.05,
        testing_mode = "parametric",
        entanglement = "weak",
        output_prefix = paste0("inst/extdata/", case),
        .score_dmrs = TRUE
    )
    saveRDS(dmrs_cment, cment_file)
}
