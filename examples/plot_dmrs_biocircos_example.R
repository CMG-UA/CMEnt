library(DMRsegal)

if (file.exists("R/plot_dmrs.R")) {
    devtools::load_all(".")
}

if (file.exists("data/beta.rda")) {
    load("data/beta.rda")
    load("data/pheno.rda")
    load("data/array_type.rda")
    dmrs <- readRDS("inst/extdata/example_output.rds")
} else {
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))
    load(system.file("data/array_type.rda", package = "DMRsegal"))
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal"))
}

if (length(dmrs) > 0) {
    dmrs_subset <- dmrs
    
    plotDMRsCircos(
        dmrs = dmrs_subset,
        beta = beta,
        pheno = pheno,
        array = array_type,
        genome = "hg19",
        sample_group_col = "Sample_Group"
    )
    
} else {
    cat("No DMRs found in the example data.\n")
}
