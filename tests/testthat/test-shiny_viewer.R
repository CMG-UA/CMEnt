test_that("validateOutputPrefix detects missing files", {
    temp_dir <- tempdir()
    fake_prefix <- file.path(temp_dir, "nonexistent_analysis")
    result <- DMRsegal:::.validateOutputPrefix(fake_prefix)
    expect_false(result$valid)
    expect_true(length(result$errors) >= 2)
    expect_true(any(grepl("DMRs file not found", result$errors)))
    expect_true(any(grepl("Beta values file not found", result$errors)))
})

test_that("validateOutputPrefix warns about missing optional files", {
    temp_dir <- tempdir()
    prefix <- file.path(temp_dir, "test_shiny_viewer")

    dmrs_file <- paste0(prefix, ".dmrs.tsv.gz")
    gz <- gzfile(dmrs_file, "w")
    write.table(
        data.frame(
            chr = "chr1", start = 1000, end = 2000, width = 1000,
            cpgs_num = 5, seeds_num = 2, delta_beta = 0.3,
            cases_beta = 0.7, controls_beta = 0.4, cpgs = "cg1,cg2,cg3,cg4,cg5",
            seeds = "cg2,cg3", upstream_cpgs = "cg1", downstream_cpgs = "cg4,cg5",
            start_seed = "cg2", end_seed = "cg3", start_seed_pos = 1200, end_seed_pos = 1800,
            start_cpg = "cg1", end_cpg = "cg5", id = "chr1:cg1-cg5"
        ),
        gz, sep = "\t", row.names = FALSE, quote = FALSE
    )
    close(gz)

    beta_file <- paste0(prefix, ".seeds_beta.tsv.gz")
    gz <- gzfile(beta_file, "w")
    write.table(
        data.frame(cpg = c("cg1", "cg2", "cg3"), S1 = c(0.5, 0.6, 0.7), S2 = c(0.4, 0.5, 0.6)),
        gz, sep = "\t", row.names = FALSE, quote = FALSE
    )
    close(gz)

    result <- DMRsegal:::.validateOutputPrefix(prefix)
    expect_true(result$valid)
    expect_equal(length(result$errors), 0)
    expect_true(length(result$warnings) >= 2)
    expect_true(any(grepl("Interactions file not found", result$warnings)))
    expect_true(any(grepl("Components file not found", result$warnings)))

    unlink(dmrs_file)
    unlink(beta_file)
})

test_that("Shiny module UI functions return tagList", {
    expect_s3_class(DMRsegal:::.overviewUI("test"), "shiny.tag.list")
    expect_s3_class(DMRsegal:::.plotDMRUI("test"), "shiny.tag.list")
    expect_s3_class(DMRsegal:::.plotDMRsUI("test"), "shiny.tag.list")
    expect_s3_class(DMRsegal:::.manhattanUI("test"), "shiny.tag.list")
    expect_s3_class(DMRsegal:::.blockFormationUI("test"), "shiny.tag.list")
    expect_s3_class(DMRsegal:::.circosUI("test"), "shiny.tag.list")
})

test_that("createViewerUI returns a shiny UI object", {
    ui <- DMRsegal:::.createViewerUI()
    expect_s3_class(ui, "bslib_page")
})
