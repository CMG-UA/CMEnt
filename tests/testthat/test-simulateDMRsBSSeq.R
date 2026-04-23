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

test_that("simulateDMRsBSSeq returns dmrseq-like outputs", {
    bs <- create_simulation_bsseq()
    sim <- simulateDMRsBSSeq(
        bs,
        num_dmrs = 4,
        delta_max0 = 0.25,
        min_cpgs = 5,
        max_cpgs = 20,
        seed = 123
    )

    expect_s4_class(sim$bs, "BSseq")
    expect_s4_class(sim$gr.dmrs, "GRanges")
    expect_equal(length(sim$gr.dmrs), 4)
    expect_equal(ncol(sim$bs), ncol(bs))
    expect_equal(
        colnames(sim$bs),
        c(paste0("Condition1_Rep", seq_len(3)), paste0("Condition2_Rep", seq_len(3)))
    )
    expect_equal(nrow(sim$truth), 4)
    expect_true(all(c("seqnames", "start", "end", "delta_beta", "num_cpgs") %in% colnames(sim$truth)))
    expect_true(all(bsseq::getCoverage(sim$bs, type = "M") <= bsseq::getCoverage(sim$bs, type = "Cov")))
})

test_that("simulateDMRsBSSeq is reproducible with seed", {
    bs <- create_simulation_bsseq()
    sim1 <- simulateDMRsBSSeq(bs, num_dmrs = 3, seed = 42, min_cpgs = 5, max_cpgs = 20)
    sim2 <- simulateDMRsBSSeq(bs, num_dmrs = 3, seed = 42, min_cpgs = 5, max_cpgs = 20)

    expect_equal(bsseq::getCoverage(sim1$bs, type = "M"), bsseq::getCoverage(sim2$bs, type = "M"))
    expect_equal(sim1$truth, sim2$truth)
})

test_that("simulateDMRsBSSeq collapses duplicate input loci before simulation", {
    bs <- create_simulation_bsseq()
    bs_dup <- bs[c(1L, seq_len(nrow(bs))), ]

    expect_warning(
        sim <- simulateDMRsBSSeq(bs_dup, num_dmrs = 3, seed = 42, min_cpgs = 5, max_cpgs = 20),
        NA
    )

    loc_key <- paste0(as.character(GenomicRanges::seqnames(sim$bs)), ":", GenomicRanges::start(sim$bs))
    expect_equal(anyDuplicated(loc_key), 0L)
    expect_equal(sim$duplicate_loci_collapsed, 1L)
})

test_that("simulateDMRsBSSeq uses simDMRs sample names for custom groups", {
    bs <- create_simulation_bsseq()
    groups <- c("untreated", "treated", "untreated", "treated", "untreated", "treated")
    sim <- simulateDMRsBSSeq(
        bs,
        groups = groups,
        case_group = "treated",
        num_dmrs = 3,
        seed = 42,
        min_cpgs = 5,
        max_cpgs = 20
    )

    expect_equal(
        colnames(sim$bs),
        c(paste0("Condition1_Rep", seq_len(3)), paste0("Condition2_Rep", seq_len(3)))
    )
    expect_equal(as.character(SummarizedExperiment::colData(sim$bs)$Sample_Group), rep(c("Condition1", "Condition2"), each = 3))
    expect_equal(sim$case_group, "Condition2")
    expect_equal(sim$input_case_group, "treated")
    expect_equal(unname(sim$input_groups), c("untreated", "untreated", "untreated", "treated", "treated", "treated"))
})
