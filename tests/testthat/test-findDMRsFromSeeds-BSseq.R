suppressPackageStartupMessages({
    library(testthat)
    library(DMRsegal)
    library(bsseq)
    library(GenomicRanges)
})

test_that("findDMRsFromSeeds works with BSseq input", {
    set.seed(42)
    n_loci <- 200
    n_cases <- 5
    n_controls <- 5
    n_samples <- n_cases + n_controls
    cov <- matrix(rpois(n_loci * n_samples, lambda = 20), ncol = n_samples)
    base_meth <- runif(n_loci, 0.3, 0.7)
    m_controls <- matrix(rbinom(n_loci * n_controls, size = cov[, 1:n_controls], prob = base_meth), ncol = n_controls)
    dmr_region_idx <- 50:70
    base_meth_dmr <- base_meth
    base_meth_dmr[dmr_region_idx] <- base_meth_dmr[dmr_region_idx] + 0.3
    base_meth_dmr[base_meth_dmr > 1] <- 0.95
    m_cases <- matrix(rbinom(n_loci * n_cases, size = cov[, (n_controls + 1):n_samples], prob = base_meth_dmr), ncol = n_cases)
    met <- cbind(m_controls, m_cases)
    gr <- GRanges(
        seqnames = rep("chr1", n_loci),
        ranges = IRanges(start = seq(10000, by = 50, length.out = n_loci), width = 1)
    )
    names(gr) <- paste0("cg", seq_len(n_loci))
    sample_names <- c(paste0("Control", seq_len(n_controls)), paste0("Case", seq_len(n_cases)))
    bsseq_obj <- BSseq(
        M = met, Cov = cov, gr = gr,
        sampleNames = sample_names
    )
    pheno <- data.frame(
        Sample_ID = sample_names,
        Group = c(rep("Control", n_controls), rep("Case", n_cases)),
        stringsAsFactors = FALSE
    )
    rownames(pheno) <- sample_names
    seeds <- data.frame(
        cpg_id = paste0("cg", dmr_region_idx),
        pval = rep(0.001, length(dmr_region_idx))
    )
    expect_error({
        dmrs <- findDMRsFromSeeds(
            beta = bsseq_obj,
            seeds = seeds,
            pheno = pheno,
            sample_group_col = "Group",
            min_cpgs = 5,
            max_lookup_dist = 200,
            njobs = 1,
            annotate_with_genes = FALSE
        )
    }, NA)
})
