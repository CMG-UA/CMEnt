# Script to create example data for DMRSegal package
set.seed(42)

# Create simulated methylation data with realistic patterns
createMethylationData <- function(n_cpgs, n_samples, n_dmrs = 10) {
    # Create chromosome positions with realistic spacing
    chrom <- "chr1"
    pos <- sort(sample(1:1000000, n_cpgs))
    
    # Create clusters of CpGs to represent potential DMRs
    dmr_starts <- sample(1:(n_cpgs - 20), n_dmrs)
    dmr_lengths <- sample(5:15, n_dmrs, replace = TRUE)
    dmr_cpgs <- unlist(lapply(1:n_dmrs, function(i) {
        dmr_starts[i]:(dmr_starts[i] + dmr_lengths[i])
    }))
    
    # Create sample groups
    sample_groups <- factor(rep(c("Control", "Case"), each = n_samples/2))
    
    # Create beta values with biological variation
    beta <- matrix(
        rbeta(n_cpgs * n_samples, shape1 = 5, shape2 = 5),
        nrow = n_cpgs,
        ncol = n_samples
    )
    
    # Add differential methylation in DMR regions
    for (i in 1:n_dmrs) {
        dmr_idx <- dmr_starts[i]:(dmr_starts[i] + dmr_lengths[i])
        effect_size <- runif(1, 0.2, 0.4) * sample(c(-1, 1), 1)
        
        # Add correlated changes within DMRs
        case_idx <- which(sample_groups == "Case")
        beta[dmr_idx, case_idx] <- beta[dmr_idx, case_idx] + 
            effect_size + rnorm(length(dmr_idx) * length(case_idx), 0, 0.05)
        
        # Ensure beta values stay in [0,1]
        beta[dmr_idx, case_idx] <- pmin(pmax(beta[dmr_idx, case_idx], 0), 1)
    }
    
    # Create mix of 450k and EPIC CpG IDs and sample names
    n_450k <- floor(n_cpgs * 0.6)  # 60% 450k probes
    n_epic <- n_cpgs - n_450k      # 40% EPIC-specific probes
    
    cpg_ids <- c(
        paste0("cg", sprintf("%07d", 1:n_450k)),              # 450k style
        paste0("cg", sprintf("%08d", (n_450k + 1):n_cpgs))    # EPIC style
    )
    sample_names <- paste0("sample", 1:n_samples)
    
    # Add row and column names
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
        t.test(beta[i, sample_groups == "Case"],
               beta[i, sample_groups == "Control"])$p.value
    })
    
    # Create DMPs data frame
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
    
    # Return all created objects
    list(
        beta = beta,
        dmps = dmps,
        pheno = pheno,
        true_dmrs = data.frame(
            start = pos[dmr_starts],
            end = pos[dmr_starts + dmr_lengths],
            n_cpgs = dmr_lengths + 1
        )
    )
}

# Create example dataset
example_data <- createMethylationData(
    n_cpgs = 1000,
    n_samples = 10,
    n_dmrs = 10
)

# Create beta matrix file
beta_file <- "inst/extdata/example_beta.txt"
write.table(
    cbind(ID = rownames(example_data$beta), example_data$beta),
    file = beta_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# Create gzipped version
system(paste("gzip -c", beta_file, ">", paste0(beta_file, ".gz")))

# Save R objects
example_beta <- example_data$beta
example_dmps <- example_data$dmps
example_pheno <- example_data$pheno
example_true_dmrs <- example_data$true_dmrs

save(example_beta, file = "data/example_beta.rda", compress = "xz")
save(example_dmps, file = "data/example_dmps.rda", compress = "xz")
save(example_pheno, file = "data/example_pheno.rda", compress = "xz")
save(example_true_dmrs, file = "data/example_true_dmrs.rda", compress = "xz")