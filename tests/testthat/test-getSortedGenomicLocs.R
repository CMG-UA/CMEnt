test_that("getSortedGenomicLocs preserves array probe IDs and bypasses legacy caches", {
    calls <- new.env(parent = emptyenv())
    calls$cache_reads <- character()

    probe_locs <- GenomicRanges::GRanges(
        seqnames = c("chr2", "chr1"),
        ranges = IRanges::IRanges(start = c(200L, 100L), width = 1L)
    )
    names(probe_locs) <- c("cg00000002", "cg00000001")

    local_mocked_bindings(
        .readBiocFileCacheRDS = function(cache_dir, rname) {
            calls$cache_reads <- c(calls$cache_reads, rname)
            NULL
        },
        .assertArrayAnnotPkgInstalled = function(array, genome, context) "mock.annotation",
        .loadAnnotationLocations = function(pkg_name, source_genome) probe_locs,
        .saveBiocFileCacheRDS = function(object, cache_dir, rname) invisible(NULL),
        .package = "CMEnt"
    )
    withr::local_options(list(
        CMEnt.annotation_cache_dir = tempfile("cment-anno-cache-")
    ))

    locs <- getSortedGenomicLocs(array = "EPIC", genome = "hg19")

    expect_identical(calls$cache_reads, "epic_hg19_locations_v2")
    expect_identical(rownames(locs), c("cg00000001", "cg00000002"))
    expect_identical(locs$name, rownames(locs))
})
