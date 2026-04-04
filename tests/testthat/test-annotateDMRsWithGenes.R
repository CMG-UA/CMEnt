options("DMRsegal.verbose" = 0)
test_that("annotateDMRsWithGenes matches between sequential and parallel execution", {
    skip_on_cran()
    skip_if_offline()
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("org.Hs.eg.db")

    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr2"),
        ranges = IRanges::IRanges(
            start = c(1000000, 150000000, 2000000),
            width = c(1000, 2000, 1500)
        ),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg38")
    )

    sequential <- annotateDMRsWithGenes(dmrs, genome = "hg38", njobs = 1)
    parallel <- annotateDMRsWithGenes(dmrs, genome = "hg38", njobs = 2)

    expect_s4_class(sequential, "GRanges")
    expect_s4_class(parallel, "GRanges")
    expect_identical(as.data.frame(sequential), as.data.frame(parallel))
})
