test_that("getSortedGenomicLocs supports hs1 for 450K arrays", {
    skip_on_cran()
    skip_if_offline()
    if (!nzchar(system.file(package = "IlluminaHumanMethylation450kanno.ilmn12.hg19"))) {
        skip("Package 'IlluminaHumanMethylation450kanno.ilmn12.hg19' not installed")
    }

    locs <- getSortedGenomicLocs(array = "450K", genome = "hs1")

    expect_true(is.data.frame(locs))
    expect_true(nrow(locs) > 0)
    expect_true(all(c("chr", "start", "end") %in% colnames(locs)))
})

test_that("annotateDMRsWithGenes supports hs1 via lifted hg38 gene models", {
    skip_on_cran()
    skip_if_offline()
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("org.Hs.eg.db")

    dmrs <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 1000000, end = 1001000),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hs1")
    )

    annotated <- annotateDMRsWithGenes(dmrs, genome = "hs1")

    expect_s4_class(annotated, "GRanges")
    expect_true(all(c("in_promoter_of", "in_gene_body_of") %in% colnames(S4Vectors::mcols(annotated))))
})
