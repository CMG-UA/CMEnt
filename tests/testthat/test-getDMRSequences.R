test_that("getDMRSequences works with BSgenome packages when available", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

    dmrs <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 10000, end = 10050)
    )

    sequences <- getDMRSequences(dmrs, genome = "hg19")

    expect_type(sequences, "character")
    expect_length(sequences, 1)
    expect_equal(nchar(sequences[1]), 51)
})

test_that("getDMRSequences falls back to online API when BSgenome not available", {
    skip_on_cran()
    skip_if_offline()

    dmrs <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 10000, end = 10050)
    )

    sequences <- getDMRSequences(dmrs, genome = "hg19", use_online = TRUE)

    expect_type(sequences, "character")
    expect_length(sequences, 1)
    expect_equal(nchar(sequences[1]), 51)
    expect_false(is.na(sequences[1]))
})

test_that("getDMRSequences handles multiple regions", {
    skip_on_cran()
    skip_if_offline()

    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(10000, 20000), end = c(10050, 20050))
    )

    sequences <- getDMRSequences(dmrs, genome = "hg19", use_online = TRUE)

    expect_type(sequences, "character")
    expect_length(sequences, 2)
    expect_equal(nchar(sequences[1]), 51)
    expect_equal(nchar(sequences[2]), 51)
})

test_that("getDMRSequences works with different genomes", {
    skip_on_cran()
    skip_if_offline()

    dmrs <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 10000, end = 10050)
    )

    seq_hg38 <- getDMRSequences(dmrs, genome = "hg38", use_online = TRUE)
    expect_type(seq_hg38, "character")
    expect_length(seq_hg38, 1)
})
