test_that("findDMPsBSSeq handles HDF5-backed BSseq methylation matrices", {
    skip_if_not(is_installed_for_tests("DSS"))
    skip_if_not(is_installed_for_tests("bsseq"))
    skip_if_not(is_installed_for_tests("HDF5Array"))

    set.seed(4)
    n_sites <- 12L
    n_samples <- 6L
    cov <- matrix(rpois(n_sites * n_samples, lambda = 20) + 1L, nrow = n_sites)
    prob <- matrix(0.4, nrow = n_sites, ncol = n_samples)
    prob[, 4:6] <- 0.75
    meth <- matrix(
        rbinom(n_sites * n_samples, size = as.vector(cov), prob = as.vector(prob)),
        nrow = n_sites
    )
    bs <- bsseq::BSseq(
        chr = rep("chr1", n_sites),
        pos = seq_len(n_sites) * 100L,
        M = meth,
        Cov = cov,
        sampleNames = paste0("s", seq_len(n_samples))
    )
    bs_hdf5 <- HDF5Array::saveHDF5SummarizedExperiment(
        bs,
        dir = tempfile("bsseq_hdf5_"),
        replace = TRUE
    )
    pheno <- data.frame(
        Sample_ID = paste0("s", seq_len(n_samples)),
        Sample_Group = rep(c("control", "case"), each = 3L)
    )

    output_file <- tempfile(fileext = ".tsv.gz")
    dmps <- CMEnt::findDMPsBSSeq(
        bsseq = bs_hdf5,
        samplesheet = pheno,
        sample_group_col = "Sample_Group",
        id_col = "Sample_ID",
        chr = "all",
        case_group = "case",
        output_file = output_file,
        njobs = 1
    )

    expect_equal(nrow(dmps), n_sites)
    expect_true(all(is.finite(dmps$delta_beta)))
    written <- utils::read.table(output_file, header = TRUE, sep = "\t", check.names = FALSE)
    expect_equal(colnames(written), colnames(dmps))
    expect_equal(nrow(written), nrow(dmps))
})
