options("CMEnt.verbose" = 0)

test_that("buildDMRs matches in-memory results for a full synthetic BED input", {
    skip_on_ci()
    fixture <- makeSyntheticBuildDMRsFixture()
    baseline <- runSyntheticBuildDMRs(fixture)
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
    bed_dmrs <- runSyntheticBuildDMRs(
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

test_that("buildDMRs detects BED files by extension", {
    skip_on_ci()
    fixture <- makeSyntheticBuildDMRsFixture()
    baseline <- runSyntheticBuildDMRs(fixture)
    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))

    writeSyntheticBedFile(fixture, bed_file)
    bed_dmrs <- runSyntheticBuildDMRs(
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

test_that("buildDMRs rejects BED seeds that are not in chr:pos format", {
    fixture <- makeSyntheticBuildDMRsFixture()
    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))

    writeSyntheticBedFile(fixture, bed_file)

    expect_error(
        runSyntheticBuildDMRs(
            fixture,
            beta = bed_file,
            bed_provided = TRUE,
            bed_chrom_col = "chrom",
            bed_start_col = "start"
        ),
        "must be in 'chr:pos' format"
    )
})

test_that("custom BED preprocessing returns locations without reparsing tabix output", {
    fixture <- makeSyntheticBuildDMRsFixture()
    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))

    writeSyntheticBedFile(fixture, bed_file)

    local_mocked_bindings(
        genomicLocsFromTabix = function(...) stop("tabix output should not be reparsed"),
        .package = "CMEnt"
    )
    ret <- readCustomMethylationBedData(
        bed_file,
        pheno = fixture$pheno,
        chrom_col = "chrom",
        start_col = "start",
        njobs = 2
    )

    expect_true(file.exists(ret$tabix_file))
    locations <- as.data.frame(ret$locations)
    expect_equal(rownames(locations), paste0(fixture$locs$chr, ":", fixture$locs$start))
    expect_equal(locations$chr, fixture$locs$chr)
    expect_equal(as.integer(locations$start), fixture$locs$start)
})

test_that("custom BED preprocessing falls back when tabix tools are unavailable", {
    fixture <- makeSyntheticBuildDMRsFixture()
    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))

    writeSyntheticBedFile(fixture, bed_file)

    local_mocked_bindings(
        .tabixToolsAvailable = function() FALSE,
        convertBetaToTabix = function(...) stop("tabix conversion should not run"),
        .package = "CMEnt"
    )
    ret <- suppressWarnings(readCustomMethylationBedData(
        bed_file,
        pheno = fixture$pheno,
        chrom_col = "chrom",
        start_col = "start"
    ))

    expect_null(ret$tabix_file)
    expect_true(file.exists(ret$beta_file))
    expect_equal(readLines(ret$beta_file, n = 1L), paste(c("site_id", rownames(fixture$pheno)), collapse = "\t"))
    expect_s4_class(ret$locations, "DelayedDataFrame")
})

test_that("buildDMRs handles BED files without a chr prefix", {
    skip_on_ci()
    fixture <- makeSyntheticBuildDMRsFixture()
    baseline <- runSyntheticBuildDMRs(fixture)
    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))

    writeSyntheticBedFile(fixture, bed_file, drop_chr_prefix = TRUE)
    bed_dmrs <- runSyntheticBuildDMRs(
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

test_that("buildDMRs respects custom BED coordinate column names", {
    skip_on_ci()
    fixture <- makeSyntheticBuildDMRsFixture()
    baseline <- runSyntheticBuildDMRs(fixture)
    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))

    writeSyntheticBedFile(
        fixture,
        bed_file,
        chrom_col = "my_chr",
        start_col = "my_pos",
        end_col = "my_end"
    )
    bed_dmrs <- runSyntheticBuildDMRs(
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
