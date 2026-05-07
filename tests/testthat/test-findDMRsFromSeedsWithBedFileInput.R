options("CMEnt.verbose" = 0)

test_that("findDMRsFromSeeds matches in-memory results for a full synthetic BED input", {
    fixture <- makeSyntheticFindDMRsFixture()
    baseline <- runSyntheticFindDMRs(fixture)
    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))

    writeSyntheticBedFile(
        fixture,
        bed_file,
        end_col = "end",
        id_col = "id",
        include_score = TRUE,
        include_strand = TRUE
    )
    bed_dmrs <- runSyntheticFindDMRs(
        fixture,
        beta = bed_file,
        seeds = makeSyntheticSeedsWithIds(fixture),
        seeds_id_col = "ID",
        bed_provided = TRUE,
        bed_chrom_col = "chrom",
        bed_start_col = "start"
    )

    expect_equal(as.character(GenomicRanges::seqnames(bed_dmrs)), as.character(GenomicRanges::seqnames(baseline)))
    expect_equal(GenomicRanges::start(bed_dmrs), GenomicRanges::start(baseline))
    expect_equal(GenomicRanges::end(bed_dmrs), GenomicRanges::end(baseline))
    expect_equal(as.integer(S4Vectors::mcols(bed_dmrs)$seeds_num), as.integer(S4Vectors::mcols(baseline)$seeds_num))
    expect_match(S4Vectors::mcols(bed_dmrs)$id, "^chr1:chr1:100-chr1:300$")
})

test_that("findDMRsFromSeeds detects BED files by extension", {
    fixture <- makeSyntheticFindDMRsFixture()
    baseline <- runSyntheticFindDMRs(fixture)
    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))

    writeSyntheticBedFile(fixture, bed_file)
    bed_dmrs <- runSyntheticFindDMRs(
        fixture,
        beta = bed_file,
        seeds = makeSyntheticSeedsWithIds(fixture),
        seeds_id_col = "ID",
        bed_chrom_col = "chrom",
        bed_start_col = "start"
    )

    expect_equal(GenomicRanges::start(bed_dmrs), GenomicRanges::start(baseline))
    expect_equal(GenomicRanges::end(bed_dmrs), GenomicRanges::end(baseline))
    expect_equal(as.integer(S4Vectors::mcols(bed_dmrs)$seeds_num), as.integer(S4Vectors::mcols(baseline)$seeds_num))
})

test_that("findDMRsFromSeeds rejects BED seeds that are not in chr:pos format", {
    fixture <- makeSyntheticFindDMRsFixture()
    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))

    writeSyntheticBedFile(fixture, bed_file)

    expect_error(
        runSyntheticFindDMRs(
            fixture,
            beta = bed_file,
            bed_provided = TRUE,
            bed_chrom_col = "chrom",
            bed_start_col = "start"
        ),
        "must be in 'chr:pos' format"
    )
})

test_that("findDMRsFromSeeds handles BED files without a chr prefix", {
    fixture <- makeSyntheticFindDMRsFixture()
    baseline <- runSyntheticFindDMRs(fixture)
    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))

    writeSyntheticBedFile(fixture, bed_file, drop_chr_prefix = TRUE)
    bed_dmrs <- runSyntheticFindDMRs(
        fixture,
        beta = bed_file,
        seeds = makeSyntheticSeedsWithIds(fixture, include_chr_prefix = FALSE),
        seeds_id_col = "ID",
        bed_provided = TRUE,
        bed_chrom_col = "chrom",
        bed_start_col = "start"
    )

    expect_equal(GenomicRanges::start(bed_dmrs), GenomicRanges::start(baseline))
    expect_equal(GenomicRanges::end(bed_dmrs), GenomicRanges::end(baseline))
    expect_equal(as.integer(S4Vectors::mcols(bed_dmrs)$seeds_num), as.integer(S4Vectors::mcols(baseline)$seeds_num))
    expect_match(S4Vectors::mcols(bed_dmrs)$id, "^1:1:100-1:300$")
})

test_that("findDMRsFromSeeds respects custom BED coordinate column names", {
    fixture <- makeSyntheticFindDMRsFixture()
    baseline <- runSyntheticFindDMRs(fixture)
    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))

    writeSyntheticBedFile(
        fixture,
        bed_file,
        chrom_col = "my_chr",
        start_col = "my_pos",
        end_col = "my_end"
    )
    bed_dmrs <- runSyntheticFindDMRs(
        fixture,
        beta = bed_file,
        seeds = makeSyntheticSeedsWithIds(fixture),
        seeds_id_col = "ID",
        bed_provided = TRUE,
        bed_chrom_col = "my_chr",
        bed_start_col = "my_pos"
    )

    expect_equal(GenomicRanges::start(bed_dmrs), GenomicRanges::start(baseline))
    expect_equal(GenomicRanges::end(bed_dmrs), GenomicRanges::end(baseline))
    expect_equal(as.integer(S4Vectors::mcols(bed_dmrs)$seeds_num), as.integer(S4Vectors::mcols(baseline)$seeds_num))
})
