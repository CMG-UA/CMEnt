test_that("Circos PCA excludes constant sites while retaining their heatmap values", {
    beta <- rbind(
        cg1 = rep(0.5, 6),
        cg2 = c(0.1, 0.2, 0.5, 0.6, 0.8, 0.9),
        cg3 = c(0.9, 0.6, 0.5, 0.2, 0.3, 0.1)
    )
    colnames(beta) <- paste0("S", seq_len(ncol(beta)))
    sorted_locs <- data.frame(
        chr = "chr1", start = c(100L, 200L, 300L),
        row.names = rownames(beta)
    )
    beta_handler <- getBetaHandler(beta, sorted_locs = sorted_locs)
    pheno <- data.frame(group = rep("case", ncol(beta)), row.names = colnames(beta))
    dmrs <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 300L))
    dmrs$sites <- "cg2,cg3"

    # Constant sites must not change which samples represent the variable sites.
    set.seed(42)
    expect_no_warning(expected <- CMEnt:::.prepareCircosHeatmapData(
        dmrs, beta_handler, pheno, "group", max_num_samples_per_group = 2
    ))
    dmrs$sites <- "cg1,cg2,cg3"
    set.seed(42)
    expect_warning(result <- CMEnt:::.prepareCircosHeatmapData(
        dmrs, beta_handler, pheno, "group", max_num_samples_per_group = 2
    ), "Excluding 1 zero-variance sites from PCA for sample group 'case'", fixed = TRUE)

    expect_identical(result$reduced_pheno, expected$reduced_pheno)
    samples <- rownames(result$reduced_pheno)
    expect_length(samples, 2L)
    expect_identical(colnames(result$heatmap_df)[-(1:3)], samples)
    expect_equal(as.matrix(result$heatmap_df[, samples]), beta[, samples, drop = FALSE])
    expect_equal(rownames(result$heatmap_df), rownames(beta))
})

test_that("Circos PCA warns and uses the first samples when all sites are constant", {
    for (site_values in list(0.5, c(0.2, 0.8))) {
        beta <- matrix(
            site_values, nrow = length(site_values), ncol = 6,
            dimnames = list(paste0("cg", seq_along(site_values)), paste0("S", 1:6))
        )
        sorted_locs <- data.frame(
            chr = "chr1", start = seq_along(site_values) * 100L,
            row.names = rownames(beta)
        )
        beta_handler <- getBetaHandler(beta, sorted_locs = sorted_locs)
        pheno <- data.frame(
            group = c("case", "control", "case", "case", "control", "case"),
            row.names = c("S4", "S5", "S2", "S1", "S6", "S3")
        )
        dmrs <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 300L))
        dmrs$sites <- paste(rownames(beta), collapse = ",")

        expect_warning(result <- CMEnt:::.prepareCircosHeatmapData(
            dmrs, beta_handler, pheno, "group", max_num_samples_per_group = 2
        ), paste0(
            "All sites have zero variance within sample group 'case'. ",
            "Skipping PCA and using the first 2 samples for the Circos heatmap."
        ), fixed = TRUE)

        samples <- c("S4", "S2", "S5", "S6")
        expect_identical(result$reduced_pheno, pheno[samples, , drop = FALSE])
        expect_identical(colnames(result$heatmap_df)[-(1:3)], samples)
        expect_equal(as.matrix(result$heatmap_df[, samples]), beta[, samples, drop = FALSE])
    }
})

test_that("Circos constant sites need no warning when sample reduction is unnecessary", {
    beta <- matrix(0.5, nrow = 1, ncol = 3, dimnames = list("cg1", c("S1", "S2", "S3")))
    beta_handler <- getBetaHandler(
        beta, sorted_locs = data.frame(chr = "chr1", start = 100L, row.names = "cg1")
    )
    pheno <- data.frame(group = rep("case", ncol(beta)), row.names = colnames(beta))
    dmrs <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 100L))
    dmrs$sites <- "cg1"

    for (limit in c(0, 3)) {
        expect_no_warning(result <- CMEnt:::.prepareCircosHeatmapData(
            dmrs, beta_handler, pheno, "group", max_num_samples_per_group = limit
        ))
        expect_identical(result$reduced_pheno, pheno)
        expect_equal(as.matrix(result$heatmap_df[, colnames(beta)]), beta)
    }
})
