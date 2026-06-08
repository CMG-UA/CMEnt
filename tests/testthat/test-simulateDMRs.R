suppressPackageStartupMessages({
    library(testthat)
    library(CMEnt)
    library(bsseq)
})
options("CMEnt.verbose" = 0)

create_simulation_bsseq <- function(n_sites = 80, n_samples = 6, seed = 11) {
    set.seed(seed)
    chr <- rep("chr1", n_sites)
    cluster_id <- rep(seq_len(ceiling(n_sites / 8)), each = 8)[seq_len(n_sites)]
    within_cluster <- ave(seq_along(cluster_id), cluster_id, FUN = seq_along) - 1L
    pos <- 2000L * cluster_id + 80L * within_cluster
    cov <- matrix(
        rpois(n_sites * n_samples, lambda = 20) + 1L,
        nrow = n_sites,
        ncol = n_samples
    )
    latent_sample <- rnorm(n_samples, 0, 0.5)
    eta <- matrix(qlogis(seq(0.25, 0.75, length.out = n_sites)), nrow = n_sites, ncol = n_samples)
    eta <- eta + outer(sin(seq(0, 2 * pi, length.out = n_sites)), latent_sample)
    prob <- plogis(eta)
    meth <- matrix(
        rbinom(n_sites * n_samples, size = c(cov), prob = c(prob)),
        nrow = n_sites,
        ncol = n_samples
    )
    BSseq(
        chr = chr,
        pos = pos,
        M = meth,
        Cov = cov,
        sampleNames = paste0("sample_", seq_len(n_samples))
    )
}

create_simulation_microarray <- function(n_sites = 90, n_samples = 6, seed = 17) {
    set.seed(seed)
    chr <- rep("chr1", n_sites)
    cluster_id <- rep(seq_len(ceiling(n_sites / 6)), each = 6)[seq_len(n_sites)]
    within_cluster <- ave(seq_along(cluster_id), cluster_id, FUN = seq_along) - 1L
    pos <- 50000L * cluster_id + 300L * within_cluster
    beta <- matrix(
        stats::plogis(
            matrix(stats::rnorm(n_sites * n_samples, sd = 0.9), nrow = n_sites, ncol = n_samples) +
                outer(seq(-1, 1, length.out = n_sites), stats::rnorm(n_samples, sd = 0.2))
        ),
        nrow = n_sites,
        ncol = n_samples
    )
    beta <- pmin(pmax(beta, 0.01), 0.99)
    rownames(beta) <- paste0("cg", seq_len(n_sites))
    colnames(beta) <- paste0("sample_", seq_len(n_samples))
    sorted_locs <- data.frame(
        chr = chr,
        start = pos,
        end = pos + 1L,
        row.names = rownames(beta),
        stringsAsFactors = FALSE
    )
    list(beta = beta, sorted_locs = sorted_locs)
}

test_that("simulateDMRs returns dmrseq-like outputs for BSseq input", {
    bs <- create_simulation_bsseq()
    set.seed(123)
    sim <- simulateDMRs(
        beta = bs,
        num_dmrs = 4,
        delta_max0 = 0.25,
        min_sites = 5,
        max_sites = 20
    )

    expect_equal(sim$assay, "BSseq")
    expect_s4_class(sim$simulated, "BSseq")
    expect_s4_class(sim$gr.dmrs, "GRanges")
    expect_equal(length(sim$gr.dmrs), 4)
    expect_equal(ncol(sim$simulated), ncol(bs))
    expect_equal(colnames(sim$simulated), colnames(bs))
    expect_equal(rownames(sim$pheno), colnames(bs))
    expect_equal(colnames(sim$pheno), "Sample_Group")
    expect_equal(sim$pheno$Sample_Group, c(rep("control", 3), rep("case", 3)))
    expect_equal(
        as.data.frame(SummarizedExperiment::colData(sim$simulated)),
        sim$pheno
    )
    expect_equal(nrow(sim$truth), 4)
    expect_true(all(c("seqnames", "start", "end", "delta_beta", "num_sites") %in% colnames(sim$truth)))
    expect_true(all(bsseq::getCoverage(sim$simulated, type = "M") <= bsseq::getCoverage(sim$simulated, type = "Cov")))
})

test_that("simulateDMRs fits and stores correlation calibration metadata", {
    bs <- create_simulation_bsseq()
    set.seed(456)
    sim <- simulateDMRs(
        beta = bs,
        num_dmrs = 4,
        delta_max0 = 0.25,
        min_sites = 5,
        max_sites = 20
    )

    expect_true(all(c(
        "background_corr_target", "expected_correlation",
        "corr_target", "corr_sd_used", "corr_metric_estimate", "neighbor_window"
    ) %in% colnames(sim$truth)))
    expect_true(all(c(
        "background_corr_target", "expected_correlation",
        "corr_target", "corr_sd_used", "corr_metric_estimate", "neighbor_window"
    ) %in% names(GenomicRanges::mcols(sim$gr.dmrs))))
    expect_true(any(is.finite(sim$truth$corr_target)))
    expect_true(all(is.finite(sim$truth$corr_sd_used)))
    expect_true(all(is.finite(sim$truth$corr_metric_estimate)))
    expect_true(all(sim$truth$expected_correlation == 0.7))
    expect_true(all(sim$truth$corr_target >= sim$truth$expected_correlation))
    finite_bg <- is.finite(sim$truth$background_corr_target)
    expect_true(all(sim$truth$corr_target[finite_bg] >= sim$truth$background_corr_target[finite_bg]))
    expect_true(all(sim$truth$neighbor_window == 5L))
})

test_that("simulateDMRs reports full touched regions by default", {
    bs <- create_simulation_bsseq()
    set.seed(456)
    sim <- simulateDMRs(
        beta = bs,
        num_dmrs = 4,
        delta_max0 = 0.25,
        min_sites = 5,
        max_sites = 20
    )

    expect_equal(as.character(GenomicRanges::seqnames(sim$gr.dmrs)), as.character(GenomicRanges::seqnames(sim$selected_regions)))
    expect_equal(sim$truth$start, GenomicRanges::start(sim$selected_regions))
    expect_equal(sim$truth$end, GenomicRanges::end(sim$selected_regions))
})

test_that("simulateDMRs uses expected correlation as a local-background floor", {
    bs <- create_simulation_bsseq()
    set.seed(456)
    sim_low_target <- simulateDMRs(
        beta = bs,
        num_dmrs = 4,
        delta_max0 = 0.25,
        min_sites = 5,
        max_sites = 20,
        expected_correlation = 0.2
    )
    set.seed(456)
    sim_high_target <- simulateDMRs(
        beta = bs,
        num_dmrs = 4,
        delta_max0 = 0.25,
        min_sites = 5,
        max_sites = 20,
        expected_correlation = 0.6
    )

    expect_equal(sim_low_target$truth$background_corr_target, sim_high_target$truth$background_corr_target)
    expect_equal(sim_low_target$truth$corr_target, pmax(0.2, sim_low_target$truth$background_corr_target))
    expect_equal(sim_high_target$truth$corr_target, pmax(0.6, sim_high_target$truth$background_corr_target))
    expect_true(all(sim_high_target$truth$corr_target >= sim_low_target$truth$corr_target))
    expect_true(all(sim_high_target$truth$corr_metric_estimate >= sim_low_target$truth$corr_metric_estimate))
})

test_that("simulateDMRs exposes the correlated-neighbor window", {
    bs <- create_simulation_bsseq()
    set.seed(456)
    sim <- simulateDMRs(
        beta = bs,
        num_dmrs = 4,
        delta_max0 = 0.25,
        min_sites = 5,
        max_sites = 20,
        neighbor_window = 7L
    )

    expect_equal(sim$neighbor_window, 7L)
    expect_true(all(sim$truth$neighbor_window == 7L))
    expect_true(all(GenomicRanges::mcols(sim$gr.dmrs)$neighbor_window == 7L))
})

test_that("simulateDMRs is reproducible with external seed for BSseq input", {
    bs <- create_simulation_bsseq()
    set.seed(42)
    sim1 <- simulateDMRs(beta = bs, num_dmrs = 3, min_sites = 5, max_sites = 20)
    set.seed(42)
    sim2 <- simulateDMRs(beta = bs, num_dmrs = 3, min_sites = 5, max_sites = 20)

    expect_equal(bsseq::getCoverage(sim1$simulated, type = "M"), bsseq::getCoverage(sim2$simulated, type = "M"))
    expect_equal(sim1$truth, sim2$truth)
})

test_that("simulateDMRs collapses duplicate input loci before simulation", {
    bs <- create_simulation_bsseq()
    bs_dup <- bs[c(1L, seq_len(nrow(bs))), ]

    expect_warning(
        {
            set.seed(42)
            sim <- simulateDMRs(beta = bs_dup, num_dmrs = 3, min_sites = 5, max_sites = 20)
        },
        NA
    )

    loc_key <- paste0(as.character(GenomicRanges::seqnames(sim$simulated)), ":", GenomicRanges::start(sim$simulated))
    expect_equal(anyDuplicated(loc_key), 0L)
    expect_equal(sim$duplicate_loci_collapsed, 1L)
})

test_that("simulateDMRs reorders collapsed loci before segmenting", {
    meth <- cov <- matrix(1L, nrow = 9L, ncol = 2L)
    chr <- rep("chr2", 9L)
    pos <- c(1000L, 1100L, 1200L, 100L, 200L, 300L, 100L, 200L, 300L)

    collapsed <- CMEnt:::.collapseSimulationDuplicateLoci(meth, cov, chr, pos)
    raw_segments <- CMEnt:::.findContiguousSegments(collapsed$chr, collapsed$pos, max_gap = 500L)
    ordered <- CMEnt:::.orderSimulationLoci(collapsed)
    segments <- CMEnt:::.findContiguousSegments(ordered$chr, ordered$pos, max_gap = 500L)

    expect_gt(length(unique(raw_segments)), 1L)
    expect_equal(ordered$pos, c(100L, 200L, 300L, 1000L, 1100L, 1200L))
    expect_equal(lengths(split(ordered$pos, segments)), c("1" = 3L, "2" = 3L))
})

test_that("simulateDMRs keeps sample names and returns case/control phenotype", {
    bs <- create_simulation_bsseq()
    groups <- c("untreated", "treated", "untreated", "treated", "untreated", "treated")
    set.seed(42)
    sim <- simulateDMRs(
        beta = bs,
        groups = groups,
        case_group = "treated",
        num_dmrs = 3,
        min_sites = 5,
        max_sites = 20
    )

    expect_equal(colnames(sim$simulated), colnames(bs))
    expect_equal(
        as.character(SummarizedExperiment::colData(sim$simulated)$Sample_Group),
        c("control", "case", "control", "case", "control", "case")
    )
    expect_equal(sim$pheno$Sample_Group, c("control", "case", "control", "case", "control", "case"))
    expect_equal(sim$case_group, "case")
    expect_equal(sim$input_case_group, "treated")
    expect_equal(unname(sim$input_groups), groups)
})

test_that("simulateDMRs supports microarray beta input", {
    array_input <- create_simulation_microarray()
    set.seed(321)
    sim <- simulateDMRs(
        beta = array_input$beta,
        sorted_locs = array_input$sorted_locs,
        num_dmrs = 4,
        delta_max0 = 0.25,
        min_sites = 4,
        max_sites = 20,
        max_gap = 500L
    )

    expect_equal(sim$assay, "microarray")
    expect_true(is.matrix(sim$simulated))
    expect_true(is.data.frame(sim$beta_locs))
    expect_equal(ncol(sim$simulated), ncol(array_input$beta))
    expect_equal(colnames(sim$simulated), colnames(array_input$beta))
    expect_equal(rownames(sim$pheno), colnames(array_input$beta))
    expect_equal(colnames(sim$pheno), "Sample_Group")
    expect_equal(length(sim$gr.dmrs), 4)
    expect_equal(rownames(sim$beta_locs), rownames(sim$simulated))
    expect_true(all(sim$simulated >= 0 & sim$simulated <= 1, na.rm = TRUE))
})

test_that("simulateDMRs extends an existing samplesheet with case/control groups", {
    array_input <- create_simulation_microarray()
    samplesheet <- data.frame(
        age = seq(40, 45),
        batch = rep(c("A", "B"), 3),
        row.names = rev(colnames(array_input$beta)),
        stringsAsFactors = FALSE
    )

    set.seed(321)
    sim <- simulateDMRs(
        beta = array_input$beta,
        sorted_locs = array_input$sorted_locs,
        samplesheet = samplesheet,
        sample_group_col = "condition",
        num_dmrs = 4,
        delta_max0 = 0.25,
        min_sites = 4,
        max_sites = 20,
        max_gap = 500L
    )

    expect_equal(rownames(sim$pheno), colnames(array_input$beta))
    expect_equal(colnames(sim$pheno), c("age", "batch", "condition"))
    expect_equal(sim$pheno$age, samplesheet[colnames(array_input$beta), "age"])
    expect_equal(sim$pheno$condition, c(rep("control", 3), rep("case", 3)))
})

test_that("simulateDMRs residualizes covariates before simulation", {
    array_input <- create_simulation_microarray(n_sites = 90, n_samples = 8)
    batch <- rep(c("A", "A", "B", "B"), 2)
    groups <- rep(c("control", "case"), 4)
    batch_effect <- ifelse(batch == "B", 1.5, 0)
    beta <- plogis(qlogis(array_input$beta) + matrix(
        batch_effect,
        nrow = nrow(array_input$beta),
        ncol = length(batch),
        byrow = TRUE
    ))
    colnames(beta) <- colnames(array_input$beta)
    rownames(beta) <- rownames(array_input$beta)
    samplesheet <- data.frame(
        batch = batch,
        row.names = colnames(beta),
        stringsAsFactors = FALSE
    )

    set.seed(321)
    sim <- simulateDMRs(
        beta = beta,
        sorted_locs = array_input$sorted_locs,
        groups = groups,
        case_group = "case",
        samplesheet = samplesheet,
        covariates = "batch",
        num_dmrs = 4,
        delta_max0 = 0.25,
        min_sites = 4,
        max_sites = 20,
        max_gap = 500L,
        resample_counts = FALSE
    )

    before <- qlogis(CMEnt:::.clampSimulationBeta(beta))
    after <- qlogis(CMEnt:::.clampSimulationBeta(sim$simulated))
    before_batch_delta <- mean(before[, batch == "B"]) - mean(before[, batch == "A"])
    after_batch_delta <- mean(after[, batch == "B"]) - mean(after[, batch == "A"])
    expect_lt(abs(after_batch_delta), abs(before_batch_delta) * 0.25)
    expect_equal(sim$pheno$batch, samplesheet[colnames(beta), "batch"])
    expect_equal(sim$covariates, "batch")
})

test_that("simulateDMRs supports BSseq input provided as an rds path", {
    bs <- create_simulation_bsseq()
    bs_rds <- tempfile(fileext = ".rds")
    saveRDS(bs, bs_rds)
    on.exit(unlink(bs_rds), add = TRUE)

    set.seed(777)
    sim <- simulateDMRs(
        beta = bs_rds,
        num_dmrs = 3,
        delta_max0 = 0.25,
        min_sites = 5,
        max_sites = 20
    )

    expect_equal(sim$assay, "BSseq")
    expect_s4_class(sim$simulated, "BSseq")
    expect_equal(ncol(sim$simulated), ncol(bs))
})

test_that("simulateDMRs restores smooth autocorrelated site profiles within DMRs", {
    pos <- c(100L, 160L, 240L, 330L, 430L, 540L, 660L, 790L)
    cov <- matrix(400L, nrow = length(pos), ncol = 6)
    eta <- matrix(
        qlogis(seq(0.25, 0.75, length.out = length(pos))),
        nrow = length(pos),
        ncol = 6
    )
    meth <- matrix(
        round(plogis(eta) * cov),
        nrow = nrow(cov),
        ncol = ncol(cov)
    )
    bs <- BSseq(
        chr = rep("chr1", length(pos)),
        pos = pos,
        M = meth,
        Cov = cov,
        sampleNames = paste0("sample_", seq_len(ncol(cov)))
    )

    set.seed(2025)
    sim <- simulateDMRs(
        beta = bs,
        num_dmrs = 1,
        delta_max0 = 0.3,
        min_sites = length(pos),
        max_sites = length(pos),
        truth_min_delta_beta = 0,
        resample_counts = FALSE
    )

    case_samples <- which(sim$groups == sim$case_group)
    orig_eta <- CMEnt:::.convertCountsToLogits(
        bsseq::getCoverage(bs, type = "M")[, case_samples, drop = FALSE],
        bsseq::getCoverage(bs, type = "Cov")[, case_samples, drop = FALSE]
    )
    sim_eta <- CMEnt:::.convertCountsToLogits(
        bsseq::getCoverage(sim$simulated, type = "M")[, case_samples, drop = FALSE],
        bsseq::getCoverage(sim$simulated, type = "Cov")[, case_samples, drop = FALSE]
    )
    observed_profile <- rowMeans(sim_eta - orig_eta)
    observed_profile <- observed_profile / max(abs(observed_profile))

    base_profile <- CMEnt:::.simulationEffectProfile(
        pos = pos,
        degree = 4L,
        flank_fraction = 0.2
    )

    expect_gt(stats::cor(observed_profile[-length(observed_profile)], observed_profile[-1L]), 0.3)
    expect_lt(stats::cor(abs(observed_profile), base_profile), 0.995)
})
