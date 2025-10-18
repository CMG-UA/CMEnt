## ----setup, include=FALSE------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = TRUE, message = TRUE)


## ----load_libraries------------------------------------------------------------------------------------------------------------
library(progressr)
progressr::handlers(global = TRUE)
benchmark_output_dir <- "benchmark_data"
dir.create(benchmark_output_dir, showWarnings=FALSE)


## ----load_data-----------------------------------------------------------------------------------------------------------------
library(minfi)
# Get sample type information from MsetEx
sample_types <- pData(minfiData::MsetEx)$Sample_Group
sample_types[sample_types == "GroupA"] <- 'normal'
sample_types[sample_types == "GroupB"] <- 'cancer'
sample_groups <- factor(sample_types, levels = c("normal", "cancer"))
selected_samples <- seq_along(sample_groups)

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
dmps <- dmpFinder(mset,
                  pheno = pheno$status,
                  type = "categorical",
                  shrinkVar = TRUE)
dmps <- dmps[!is.na(dmps$pval), ]
# Add delta beta to dmps results
dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")

# Filter significant DMPs
sig_dmps <- dmps[dmps$pval_adj < 0.05, ]


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
  library(profvis)
profvis({
  dmrs_segal <- DMRSegal::findDMRsFromSeeds(
        beta_file = beta_file,
        dmps_file = dmps_file,
        pheno = pheno,
        sample_group_col = "group",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        expansion_relaxation=1,
        njobs = 1,
        verbose = 2
      )
})
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
  grset <- preprocessQuantile(mset)
  
  bumphunter_results <- bumphunter(
      grset,
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

intersection <- function(gr1, gr2) {
    sum(width(GenomicRanges::intersect(gr1, gr2)))
}

smallest_width <- function(gr1, gr2) {
    pmin(sum(width(gr1)), sum(width(gr2)))
}

# IoU metric
iou <- function(gr1, gr2) {
    intersection <- sum(width(GenomicRanges::intersect(gr1, gr2)))
    union <- sum(width(GenomicRanges::union(gr1, gr2)))
    if (union == 0) return(0)
    intersection / union
}

# Find overlaps
overlap_segal_dmrcate <- length(subsetByOverlaps(gr_segal, gr_dmrcate))
overlap_segal_bumphunter <- length(subsetByOverlaps(gr_segal, gr_bumphunter))
overlap_dmrcate_bumphunter <- length(subsetByOverlaps(gr_dmrcate, gr_bumphunter))
# Calculate intersection over smallest width
intersection_segal_dmrcate <- intersection(gr_segal, gr_dmrcate) / smallest_width(gr_segal, gr_dmrcate)
intersection_segal_bumphunter <- intersection(gr_segal, gr_bumphunter) / smallest_width(gr_segal, gr_bumphunter)
intersection_dmrcate_bumphunter <- intersection(gr_dmrcate, gr_bumphunter) / smallest_width(gr_dmrcate, gr_bumphunter)


iou_segal_dmrcate <- iou(gr_segal, gr_dmrcate)
iou_segal_bumphunter <- iou(gr_segal, gr_bumphunter)
iou_dmrcate_bumphunter <- iou(gr_dmrcate, gr_bumphunter)

# Create overlap matrix
overlap_matrix <- t(matrix(
    c(
        length(gr_segal), overlap_segal_dmrcate, overlap_segal_bumphunter,
        iou_segal_dmrcate, length(gr_dmrcate), overlap_dmrcate_bumphunter,
        iou_segal_bumphunter, iou_dmrcate_bumphunter, length(gr_bumphunter)
    ),
    nrow = 3,
    dimnames = list(
        c("DMRSegal", "DMRcate", "bumphunter"),
        c("DMRSegal", "DMRcate", "bumphunter")
    )
))

library(ggplot2)
library(tidyverse)
library(ggnewscale)

mat_long <- as.data.frame(overlap_matrix) %>%
  rownames_to_column("row") %>%
  pivot_longer(-row, names_to = "col", values_to = "value")
mat_long <- mat_long %>%
  mutate(row = factor(row, levels = rownames(overlap_matrix)),
         col = factor(col, levels = colnames(overlap_matrix)))

plot_data = mat_long %>%
  mutate(triangle =
    case_when(as.numeric(col) < as.numeric(row) ~ 'Lower',
     as.numeric(col) > as.numeric(row) ~ 'Upper',
        TRUE ~ 'Diag')) %>%
  group_split(triangle)


ggplot() +
  geom_tile(data=plot_data[[1]], aes(x=col, y=row, fill=value)) +
  scale_fill_gradient(low = "white", high = "blue", name = "Counts") +

  new_scale_fill() +  # Part of ggnewscale, allows multiple fill scales
  # Lower triangle heatmap with a green color scale
  geom_tile(data = plot_data[[2]], aes(x=col, y=row, fill=value)) +
  scale_fill_gradient(low = "white", high = "green", name = "Iou") +
    new_scale_fill() +
    # Upper triangle heatmap with a red color scale
    geom_tile(data = plot_data[[3]], aes(x=col, y=row, fill=value)) +
    scale_fill_gradient(low = "white", high = "red", name = "Overlap") +
    scale_x_discrete(limits=rev) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(benchmark_output_dir,"dmr_overlap_heatmap.png"), width=6, height=5)