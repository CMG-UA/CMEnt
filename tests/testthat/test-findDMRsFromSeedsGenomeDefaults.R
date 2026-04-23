options("CMEnt.verbose" = 0)

test_that("findDMRsFromSeeds infers hg19 for legacy human arrays", {
    beta <- matrix(0.5, nrow = 2, ncol = 2)

    expect_equal(CMEnt:::.resolveFindDMRsGenome(beta = beta, array = "450K", genome = NULL), "hg19")
    expect_equal(CMEnt:::.resolveFindDMRsGenome(beta = beta, array = "27K", genome = NULL), "hg19")
    expect_equal(CMEnt:::.resolveFindDMRsGenome(beta = beta, array = "EPIC", genome = NULL), "hg19")
})

test_that("findDMRsFromSeeds infers hg38 for non-legacy inputs", {
    beta <- matrix(0.5, nrow = 2, ncol = 2)

    expect_equal(CMEnt:::.resolveFindDMRsGenome(beta = beta, array = "EPICv2", genome = NULL), "hg38")
    expect_equal(CMEnt:::.resolveFindDMRsGenome(beta = beta, array = "NULL", genome = NULL), "hg38")

    bed_file <- tempfile(fileext = ".bed")
    withr::defer(unlink(bed_file))
    writeLines("chrom\tstart\tSample1", bed_file)

    expect_equal(
        CMEnt:::.resolveFindDMRsGenome(beta = bed_file, array = "450K", genome = NULL, bed_provided = TRUE),
        "hg38"
    )
})

test_that("findDMRsFromSeeds keeps explicit genome overrides", {
    beta <- matrix(0.5, nrow = 2, ncol = 2)

    expect_equal(CMEnt:::.resolveFindDMRsGenome(beta = beta, array = "450K", genome = "hg38"), "hg38")
})
