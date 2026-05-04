test_that("BiocFileCache RDS helpers always resolve paths inside the cache directory", {
    cache_dir <- tempfile("cment-bfc-cache-")
    work_dir <- tempfile("cment-bfc-work-")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    withr::local_dir(work_dir)

    CMEnt:::.saveBiocFileCacheRDS(list(value = 1L), cache_dir, "test_entry")
    cache_file <- CMEnt:::.getBiocFileCachePath(
        CMEnt:::.getBiocFileCache(cache_dir),
        rname = "test_entry",
        ext = ".rds",
        create = FALSE
    )

    expect_true(file.exists(cache_file))
    expect_false(file.exists(file.path(work_dir, basename(cache_file))))
    expect_equal(CMEnt:::.readBiocFileCacheRDS(cache_dir, "test_entry")$value, 1L)

    CMEnt:::.saveBiocFileCacheRDS(list(value = 2L), cache_dir, "test_entry")

    expect_true(file.exists(cache_file))
    expect_false(file.exists(file.path(work_dir, basename(cache_file))))
    expect_equal(readRDS(cache_file)$value, 2L)
})

test_that("BiocFileCache path helper returns absolute paths for existing non-RDS entries", {
    cache_dir <- tempfile("cment-bfc-cache-")
    work_dir <- tempfile("cment-bfc-work-")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

    bfc <- CMEnt:::.getBiocFileCache(cache_dir)
    created_path <- CMEnt:::.getBiocFileCachePath(
        bfc,
        rname = "test_chain",
        ext = ".over.chain"
    )
    writeLines("chain", created_path)

    withr::local_dir(work_dir)
    resolved_path <- CMEnt:::.getBiocFileCachePath(
        bfc,
        rname = "test_chain",
        ext = ".over.chain",
        create = FALSE
    )

    expect_identical(normalizePath(resolved_path), normalizePath(created_path))
    expect_true(file.exists(resolved_path))
    expect_false(file.exists(file.path(work_dir, basename(created_path))))
})
