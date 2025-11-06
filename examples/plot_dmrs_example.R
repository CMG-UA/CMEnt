#!/usr/bin/env Rscript
# Example script demonstrating DMR visualization
#
# This script shows how to use the plotting functions to visualize
# DMRs identified by findDMRsFromSeeds

library(DMRsegal)

# Load example DMR results from benchmark data
dmrs <- readRDS("benchmark_data/dmrs.DMRsegal.rds")

cat("Loaded", length(dmrs), "DMRs\n")
cat("First DMR summary:\n")
print(dmrs[1])

# ============================================================================
# Example 1: Simple DMR structure plot
# ============================================================================
cat("\n=== Example 1: Simple DMR structure plot ===\n")

# Plot the first DMR
p1 <- plotDMR(dmrs, dmr_index = 1)
print(p1)

# Save to file
ggsave("dmr_plot_example_1.png", p1, width = 10, height = 4, dpi = 300)
cat("Saved plot to: dmr_plot_example_1.png\n")

# ============================================================================
# Example 2: Plot multiple DMRs in a grid
# ============================================================================
cat("\n=== Example 2: Multiple DMRs in a grid ===\n")

    
# Plot individually
for (i in 1:3) {
    p <- plotDMR(dmrs, dmr_index = i)
    fname <- sprintf("dmr_plot_example_2_%d.png", i)
    ggsave(fname, p, width = 10, height = 4, dpi = 300)
    cat("Saved plot to:", fname, "\n")
}


# ============================================================================
# Example 3: Plot DMRs with different characteristics
# ============================================================================
cat("\n=== Example 3: Highlighting specific DMR features ===\n")

# Find DMRs with high number of DMPs
dmr_data <- as.data.frame(S4Vectors::mcols(dmrs))
high_dmp_idx <- which(dmr_data$seeds_num >= 3)[1:3]

cat("Plotting DMRs with high number of DMPs:\n")
for (idx in high_dmp_idx) {
    if (!is.na(idx)) {
        p <- plotDMR(dmrs, dmr_index = idx)
        fname <- sprintf("dmr_high_dmps_%d.png", idx)
        ggsave(fname, p, width = 10, height = 4, dpi = 300)
        cat(sprintf("  DMR #%d: %d DMPs, %d CpGs - saved to %s\n",
                    idx, dmr_data$seeds_num[idx], dmr_data$cpgs_num[idx], fname))
    }
}

# ============================================================================
# Example 4: Create publication-quality plot
# ============================================================================
cat("\n=== Example 4: Publication-quality plot ===\n")

# Select an interesting DMR (one with good stats)
interesting_idx <- which(
    dmr_data$seeds_num >= 2 & 
    dmr_data$cpgs_num >= 10 & 
    abs(dmr_data$delta_beta) > 0.2
)[1]

if (!is.na(interesting_idx)) {
    p4 <- plotDMR(
        dmrs, 
        dmr_index = interesting_idx,
        show_delta_beta = TRUE,
        show_pvalues = TRUE
    )
    
    # Enhance for publication
    p4 <- p4 + ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13),
        axis.title = ggplot2::element_text(size = 12),
        axis.text = ggplot2::element_text(size = 10)
    )
    
    print(p4)
    ggsave("dmr_publication_quality.png", p4, width = 12, height = 5, dpi = 600)
    ggsave("dmr_publication_quality.pdf", p4, width = 12, height = 5)
    cat("Saved publication-quality plots\n")
} else {
    cat("No DMRs meeting the criteria for publication plot\n")
}

# ============================================================================
# Example 5: Summary statistics plot
# ============================================================================
cat("\n=== Example 5: Summary statistics ===\n")

# Create a summary plot showing DMR characteristics
library(ggplot2)

summary_data <- data.frame(
    dmr_id = seq_along(dmrs),
    n_seeds = dmr_data$seeds_num,
    n_cpgs = dmr_data$cpgs_num,
    delta_beta = dmr_data$delta_beta,
    width = GenomicRanges::width(dmrs)
)

# DMPs vs CpGs scatter plot
p5a <- ggplot(summary_data, aes(x = n_seeds, y = n_cpgs)) +
    geom_point(aes(color = abs(delta_beta)), size = 2, alpha = 0.7) +
    scale_color_gradient2(low = "blue", mid = "yellow", high = "red",
                         midpoint = 0.3, name = "|Delta Beta|") +
    labs(
        title = "DMR Composition",
        x = "Number of DMPs",
        y = "Number of CpGs"
    ) +
    theme_minimal()

print(p5a)
ggsave("dmr_summary_composition.png", p5a, width = 8, height = 6, dpi = 300)

# Volcano-like plot
p5b <- ggplot(summary_data, aes(x = delta_beta, y = neg_log10_pval)) +
    geom_point(aes(size = n_cpgs, color = n_seeds), alpha = 0.6) +
    scale_color_gradient(low = "lightblue", high = "darkblue", name = "DMPs") +
    scale_size_continuous(name = "CpGs", range = c(1, 8)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(
        title = "DMR Effect Sizes and Significance",
        x = "Delta Beta",
        y = "-log10(adjusted p-value)"
    ) +
    theme_minimal()

print(p5b)
ggsave("dmr_summary_volcano.png", p5b, width = 8, height = 6, dpi = 300)

cat("\n=== All examples completed ===\n")
cat("Generated plots:\n")
cat("  - dmr_plot_example_1.png: Basic DMR structure\n")
cat("  - dmr_plot_example_2.png: Multiple DMRs grid\n")
cat("  - dmr_high_dmps_*.png: DMRs with many DMPs\n")
cat("  - dmr_publication_quality.png/pdf: High-quality plot\n")
cat("  - dmr_summary_composition.png: DMR composition overview\n")
cat("  - dmr_summary_volcano.png: Effect sizes and significance\n")
