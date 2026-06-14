options("CMEnt.verbose" = 0)

test_that("gene bodies exclude promoter intervals", {
    genes <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(c(100L, 500L, 900L), c(300L, 800L, 1000L)),
        strand = c("+", "-", "+")
    )
    names(genes) <- c("gene1", "gene2", "gene3")

    promoters <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(c(50L, 750L, 250L), c(120L, 850L, 550L)),
        strand = c("+", "-", "+")
    )
    S4Vectors::mcols(promoters)$name <- c("gene1", "gene2", "gene3")

    gene_bodies <- CMEnt:::.trimGeneBodiesByPromoters(genes, promoters)

    expect_equal(as.character(GenomicRanges::seqnames(gene_bodies)), rep("chr1", 3))
    expect_equal(GenomicRanges::start(gene_bodies), c(121L, 551L, 900L))
    expect_equal(GenomicRanges::end(gene_bodies), c(249L, 749L, 1000L))
    expect_equal(names(gene_bodies), c("gene1", "gene2", "gene3"))
})

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
    expect_true(all(
        c("delta_beta_promoter", "delta_beta_gene_body") %in%
            colnames(S4Vectors::mcols(sequential))
    ))
    expect_identical(as.data.frame(sequential), as.data.frame(parallel))
})

test_that("feature-specific delta beta uses DMR sites within annotated regions", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(100, 500), end = c(300, 700)),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg38")
    )
    S4Vectors::mcols(dmrs)$sites <- c("cg1,cg2,cg3", "cg4")

    site_locs <- data.frame(
        chr = "chr1",
        start = c(110L, 160L, 250L, 600L),
        end = c(110L, 160L, 250L, 600L),
        row.names = paste0("cg", 1:4)
    )
    annotation_specs <- list(
        list(
            delta_column = "delta_beta_promoter",
            features = GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 175))
        ),
        list(
            delta_column = "delta_beta_gene_body",
            features = GenomicRanges::GRanges(
                "chr1",
                IRanges::IRanges(c(200, 590), c(300, 610))
            )
        )
    )
    delta_beta <- c(cg1 = 0.2, cg2 = 0.4, cg3 = -0.1, cg4 = 0.8)

    annotated_delta <- CMEnt:::.annotateDMRSiteDeltaBetaByFeature(
        dmrs = dmrs,
        annotation_specs = annotation_specs,
        site_locs = site_locs,
        site_delta_beta = delta_beta,
        aggfun = mean,
        genome = "hg38"
    )

    expect_equal(annotated_delta$delta_beta_promoter, c(0.3, NA_real_))
    expect_equal(annotated_delta$delta_beta_gene_body, c(-0.1, 0.8))
})
