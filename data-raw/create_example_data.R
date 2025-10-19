# Script to create example data for DMRSegal package
# This script uses REAL CpG IDs from 450K array annotation to ensure
# compatibility with the package's annotation checking
set.seed(42)

library(DMRSegal)

# Get real CpG annotations from 450K array
sorted_locs <- getSortedGenomicLocs(array = "450K", genome = "hg19")

# Sample CpGs from chromosome 1 for faster processing
chr1_cpgs <- sorted_locs[sorted_locs$chr == "chr1", ]
chr1_cpgs <- chr1_cpgs[order(chr1_cpgs$pos), ]

# Create simulated methylation data with realistic patterns using REAL CpG IDs
createMethylationData <- function(available_cpgs, n_cpgs = 1000, n_samples = 10, n_dmrs = 10) {
    # Sample real CpGs from available annotations
    if (n_cpgs > nrow(available_cpgs)) {
        n_cpgs <- nrow(available_cpgs)
        warning("Requested more CpGs than available, using ", n_cpgs)
    }
    
    # Sample CpGs ensuring they're spread across the chromosome
    cpg_indices <- sort(sample(seq_len(nrow(available_cpgs)), n_cpgs))
    selected_cpgs <- available_cpgs[cpg_indices, ]
    cpg_ids <- rownames(selected_cpgs)
    pos <- selected_cpgs$pos
    chrom <- as.character(selected_cpgs$chr[1])
    
    # Create clusters of CpGs to represent potential DMRs
    dmr_starts <- sample(seq_len(n_cpgs - 20), n_dmrs)
    dmr_lengths <- sample(5:15, n_dmrs, replace = TRUE)
    
    # Create sample groups
    sample_groups <- factor(rep(c("Control", "Case"), each = n_samples/2))
    sample_names <- paste0("sample", 1:n_samples)
    
    # Create beta values with biological variation
    beta <- matrix(
        rbeta(n_cpgs * n_samples, shape1 = 5, shape2 = 5),
        nrow = n_cpgs,
        ncol = n_samples
    )
    
    # Add differential methylation in DMR regions
    for (i in 1:n_dmrs) {
        start <- dmr_starts[i]
        end <- min(dmr_starts[i] + dmr_lengths[i], n_cpgs)
        dmr_idx <- start:end
        effect_size <- runif(1, 0.2, 0.4) * sample(c(-1, 1), 1)
        
        # Add correlated changes within DMRs
        case_idx <- which(sample_groups == "Case")
        beta[dmr_idx, case_idx] <- beta[dmr_idx, case_idx] + 
            effect_size + rnorm(length(dmr_idx) * length(case_idx), 0, 0.05)
        
        # Ensure beta values stay in [0,1]
        beta[dmr_idx, case_idx] <- pmin(pmax(beta[dmr_idx, case_idx], 0), 1)
    }
    
    # Add row and column names (using REAL CpG IDs)
    rownames(beta) <- cpg_ids
    colnames(beta) <- sample_names
    
    # Create sample metadata
    pheno <- data.frame(
        sample_id = sample_names,
        group = sample_groups,
        age = round(rnorm(n_samples, mean = 50, sd = 10)),
        sex = factor(sample(c("M", "F"), n_samples, replace = TRUE)),
        row.names = sample_names
    )
    
    # Calculate t-tests for each CpG
    pvals <- sapply(1:n_cpgs, function(i) {
        tryCatch({
            t.test(beta[i, sample_groups == "Case"],
                   beta[i, sample_groups == "Control"])$p.value
        }, error = function(e) 1.0)
    })
    
    # Create DMPs data frame with REAL CpG IDs as rownames
    dmps <- data.frame(
        chr = chrom,
        pos = pos,
        pval = pvals,
        pval_adj = p.adjust(pvals, method = "BH"),
        mean_diff = rowMeans(beta[, sample_groups == "Case"]) - 
                   rowMeans(beta[, sample_groups == "Control"]),
        row.names = cpg_ids,
        stringsAsFactors = FALSE
    )
    
    # Create true DMRs reference
    true_dmrs <- data.frame(
        dmr_id = seq_len(n_dmrs),
        start_pos = sapply(seq_len(n_dmrs), function(i) {
            dmr_start <- dmr_starts[i]
            pos[dmr_start]
        }),
        end_pos = sapply(seq_len(n_dmrs), function(i) {
            dmr_start <- dmr_starts[i]
            dmr_end <- min(dmr_starts[i] + dmr_lengths[i], n_cpgs)
            pos[dmr_end]
        }),
        n_cpgs = sapply(seq_len(n_dmrs), function(i) {
            dmr_start <- dmr_starts[i]
            dmr_end <- min(dmr_starts[i] + dmr_lengths[i], n_cpgs)
            dmr_end - dmr_start + 1
        })
    )
    
    # Return all created objects
    list(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        true_dmrs = true_dmrs
    )
}

# Create example dataset using REAL 450K CpG IDs
cat("Creating example data using real 450K CpG annotations...\n")
example_data <- createMethylationData(
    available_cpgs = chr1_cpgs,
    n_cpgs = 1000,
    n_samples = 10,
    n_dmrs = 10
)

cat("Generated data with", nrow(example_data$beta), "CpGs and", 
    ncol(example_data$beta), "samples\n")
cat("First 5 CpG IDs:", paste(head(rownames(example_data$beta), 5), collapse = ", "), "\n")

# Create output directory if needed
if (!dir.exists("inst/extdata")) {
    dir.create("inst/extdata", recursive = TRUE)
}

# Create beta matrix file
beta_file <- "inst/extdata/example_beta.txt"
write.table(
    cbind(ID = rownames(example_data$beta), example_data$beta),
    file = beta_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
cat("Wrote beta file to", beta_file, "\n")

# Create gzipped version
if (file.exists(paste0(beta_file, ".gz"))) {
    file.remove(paste0(beta_file, ".gz"))
}
system(paste("gzip -c", beta_file, ">", paste0(beta_file, ".gz")))
cat("Created gzipped version:", paste0(beta_file, ".gz"), "\n")

# Create data directory if needed
if (!dir.exists("data")) {
    dir.create("data")
}

# Save R objects
example_beta <- example_data$beta
example_dmps <- example_data$dmps
example_pheno <- example_data$pheno
example_true_dmrs <- example_data$true_dmrs

save(example_beta, file = "data/example_beta.rda", compress = "xz")
save(example_dmps, file = "data/example_dmps.rda", compress = "xz")
save(example_pheno, file = "data/example_pheno.rda", compress = "xz")
save(example_true_dmrs, file = "data/example_true_dmrs.rda", compress = "xz")

cat("\nSummary of generated data:\n")
cat("- Beta matrix:", nrow(example_beta), "CpGs x", ncol(example_beta), "samples\n")
cat("- DMPs:", nrow(example_dmps), "positions\n")
cat("- Significant DMPs (adj.p < 0.05):", sum(example_dmps$pval_adj < 0.05), "\n")
cat("- Phenotype data:", nrow(example_pheno), "samples\n")
cat("- True DMRs:", nrow(example_true_dmrs), "regions\n")
cat("\nAll example data files created successfully!\n")
