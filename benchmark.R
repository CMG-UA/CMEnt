## ----setup, include=FALSE------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = TRUE, message = TRUE)


## ----load_libraries------------------------------------------------------------------------------------------------------------
library(progressr)
# Configure progressr for RStudio
benchmark_output_dir <- "benchmark_data"
dir.create(benchmark_output_dir, showWarnings=FALSE)


## ----load_data-----------------------------------------------------------------------------------------------------------------
library(minfi)
# Get sample type information from MsetEx
sample_types <- pData(minfiData::MsetEx)$Sample_Group
sample_types[sample_types == "GroupA"] <- 'normal'
sample_types[sample_types == "GroupB"] <- 'cancer'
sample_groups <- factor(sample_types, levels = c("normal", "cancer"))

# Subset to a smaller set for the example (3 cancer, 3 normal)
set.seed(123)
samples_per_group <- 3
cancer_idx <- sample(which(sample_groups == "cancer"), samples_per_group)
normal_idx <- sample(which(sample_groups == "normal"), samples_per_group)
selected_samples <- c(cancer_idx, normal_idx)

# Create the final subset - MsetEx is already preprocessed
mset <- minfiData::MsetEx[, selected_samples]
pheno <- data.frame(
    status = sample_groups[selected_samples],
    row.names = colnames(mset)
)
pheno$group <- pheno$status
pheno$casecontrol <- pheno$status == "cancer"


## ----get_dmps------------------------------------------------------------------------------------------------------------------
# Create design matrix for all methods to use consistently
design <- model.matrix(~pheno$status)

# Find DMPs using limma
dmps <- dmpFinder(getBeta(mset), 
                  pheno = pheno$status,
                  type = "categorical")
dmps <- dmps[!is.na(dmps$pval),]
# Add delta beta to dmps results
dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

# Filter significant DMPs
sig_dmps <- dmps[dmps$pval_adj < 0.1, ]


## ----dmrsegal------------------------------------------------------------------------------------------------------------------

dmrsegal_file <- file.path(benchmark_output_dir,"dmrs.dmrsegal.rds")
if (!file.exists(dmrsegal_file)){
  library(DMRSegal)
  # Write beta values to temp file
  beta_file <- tempfile(fileext = ".txt")
  beta_mat <- getBeta(mset)
  
  
  write.table(
      cbind(ID = rownames(beta_mat), beta_mat),
      file = beta_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
  )
  
  beta_file <- DMRSegal::sortBetaFileByCoordinates(beta_file, overwrite=TRUE)
  
  # Write DMPs to temp file
  dmps_file <- tempfile(fileext = ".txt")
  write.table(
      sig_dmps,
      file = dmps_file,
      sep = "\t",
      quote = FALSE,
      row.names = TRUE
  )
  # option(future.debug = TRUE)
  dmrs_segal <- DMRSegal::findDMRsFromDMPs(
        beta_file = beta_file,
        dmps_tsv_file = dmps_file,
        pheno = pheno,
        sample_group_col = "group",
        min_dmps = 1,
        min_cpgs = 3,
        njobs = 8,
        verbose = 2
      )
  saveRDS(dmrs_segal, dmrsegal_file)
  # Remove temporary files
  unlink(c(beta_file, dmps_file))
  detach("package:DMRSegal", unload=TRUE)
}





## ----dmrcate-------------------------------------------------------------------------------------------------------------------
dmrcate_file <- file.path(benchmark_output_dir,"dmrs.dmrcate.rds")
if (!file.exists(dmrcate_file)){
  
  library(DMRcate)
  # DMRcate uses the same mset and design matrix as other methods
  M <- getM(mset)
  
  myannotation<-cpg.annotate("array", M,what="M", arraytype="450K", analysis.type="differential", design=design,coef=2,fdr = 0.05)
  
  # Run DMRcate with the same p-value cutoff we used for DMPs
  
  dmrcate_results <- dmrcate(myannotation, C = 2, lambda = 1000)
  dmrs_dmrcate <- extractRanges(dmrcate_results, genome='hg19')
  saveRDS(dmrs_dmrcate, dmrcate_file)
  detach("package:DMRcate", unload=TRUE)
}


## ----bumphunter----------------------------------------------------------------------------------------------------------------
bumphunter_file <- file.path(benchmark_output_dir,"dmrs.bumphunter.rds")
if (!file.exists(bumphunter_file)){
  library(bumphunter)
  # bumphunter uses the same mset and design matrix as other methods
  # Run bumphunter with the same design matrix
  library(doParallel)
  registerDoParallel(cores = 4)
  gmSet <- preprocessQuantile(mset)
  bumphunter_results <- bumphunter(
      gmSet,
      design = design,
      pickCutoff = TRUE,
      B = 100,  # Number of permutations
      type = "Beta"
  )
  dmrs_bumphunter <- makeGRangesFromDataFrame(bumphunter_results$table)
  saveRDS(dmrs_bumphunter, bumphunter_file)
}


## ----compare_results-----------------------------------------------------------------------------------------------------------
# Function to get DMR stats
get_dmr_stats <- function(dmrs, method) {
    data.frame(
        Method = method,
        Number_of_DMRs = length(dmrs),
        Mean_DMR_Width = mean(width(dmrs)),
        Median_DMR_Width = median(width(dmrs))
    )
}
dmrs_segal <- readRDS(dmrsegal_file)
dmrs_dmrcate <- readRDS(dmrcate_file)
dmrs_bumphunter <- readRDS(bumphunter_file)
# Collect stats
stats_list <- list(
    get_dmr_stats(dmrs_segal, "DMRSegal"),
    get_dmr_stats(dmrs_dmrcate, "DMRcate"),
    get_dmr_stats(dmrs_bumphunter, "bumphunter")
)

# Combine stats
results_comparison <- do.call(rbind, stats_list)
knitr::kable(results_comparison)


## ----overlap-------------------------------------------------------------------------------------------------------------------
# Convert all results to GRanges
gr_segal <- dmrs_segal
gr_dmrcate <- dmrs_dmrcate
gr_bumphunter <- dmrs_bumphunter

# Find overlaps
overlap_segal_dmrcate <- length(subsetByOverlaps(gr_segal, gr_dmrcate))
overlap_segal_bumphunter <- length(subsetByOverlaps(gr_segal, gr_bumphunter))
overlap_dmrcate_bumphunter <- length(subsetByOverlaps(gr_dmrcate, gr_bumphunter))

# Create overlap matrix
overlap_matrix <- matrix(
    c(
        1, overlap_segal_dmrcate / length(gr_segal), overlap_segal_bumphunter / length(gr_segal),
        overlap_segal_dmrcate/length(gr_dmrcate), 1, overlap_dmrcate_bumphunter / length(gr_dmrcate),
        overlap_segal_bumphunter/length(gr_bumphunter), overlap_dmrcate_bumphunter / length(gr_bumphunter), 1
    ),
    nrow = 3,
    dimnames = list(
        c("DMRSegal", "DMRcate", "bumphunter"),
        c("DMRSegal", "DMRcate", "bumphunter")
    )
)

# Plot heatmap
pheatmap::pheatmap(
    overlap_matrix,
    main = "DMR Overlap Between Methods",
    display_numbers = TRUE, cluster_rows=F, cluster_cols=F
)

