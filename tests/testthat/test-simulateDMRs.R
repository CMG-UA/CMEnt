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
    expect_equal(
        colnames(sim$simulated),
        c(paste0("Condition1_Rep", seq_len(3)), paste0("Condition2_Rep", seq_len(3)))
    )
    expect_equal(nrow(sim$truth), 4)
    expect_true(all(c("seqnames", "start", "end", "delta_beta", "num_sites") %in% colnames(sim$truth)))
    expect_true(all(bsseq::getCoverage(sim$simulated, type = "M") <= bsseq::getCoverage(sim$simulated, type = "Cov")))
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

test_that("simulateDMRs uses simDMRs sample names for custom groups", {
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

    expect_equal(
        colnames(sim$simulated),
        c(paste0("Condition1_Rep", seq_len(3)), paste0("Condition2_Rep", seq_len(3)))
    )
    expect_equal(as.character(SummarizedExperiment::colData(sim$simulated)$Sample_Group), rep(c("Condition1", "Condition2"), each = 3))
    expect_equal(sim$case_group, "Condition2")
    expect_equal(sim$input_case_group, "treated")
    expect_equal(unname(sim$input_groups), c("untreated", "untreated", "untreated", "treated", "treated", "treated"))
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
    expect_equal(length(sim$gr.dmrs), 4)
    expect_equal(rownames(sim$beta_locs), rownames(sim$simulated))
    expect_true(all(sim$simulated >= 0 & sim$simulated <= 1, na.rm = TRUE))
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
