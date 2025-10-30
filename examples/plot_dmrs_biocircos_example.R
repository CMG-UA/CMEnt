library(DMRsegal)

# For development, load the package
if (file.exists("R/plot_dmrs.R")) {
    devtools::load_all(".")
}

# For development, load data directly
# When package is installed, use: system.file("data/beta.rda", package = "DMRsegal")
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
    
    cat("Creating BioCircos plot without interactions...\n")
    plot_obj <- plotDMRsCircos(
        dmrs = dmrs_subset,
        beta = beta,
        pheno = pheno,
        array = array_type,
        genome = "hg19",
        sample_group_col = "Sample_Group",
        plot_interactions = FALSE
    )
    
    print(plot_obj)
    computeMotifBasedDMRsInteraction(dmrs_subset, genome = "hg19", array = array_type, max_fdr = 0.05)
    cat("\n\nCreating BioCircos plot with motif-based interactions...\n")
    plot_obj_with_interactions <- plotDMRsCircos(
        dmrs = dmrs_subset,
        beta = beta,
        pheno = pheno,
        array = array_type,
        genome = "hg19",
        sample_group_col = "Sample_Group",
        plot_interactions = TRUE,
        max_fdr = 0.05
    )
    
    print(plot_obj_with_interactions)
} else {
    cat("No DMRs found in the example data.\n")
}
