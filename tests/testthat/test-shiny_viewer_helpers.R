suppressPackageStartupMessages({
    library(testthat)
    library(CMEnt)
})
options("CMEnt.verbose" = 0)

test_that("viewer file paths and circos cache reload are deterministic", {
    prefix <- file.path(tempdir(), "viewer-helper")
    files <- CMEnt:::.viewerOutputFiles(prefix)

    expect_equal(files$dmrs_file, paste0(prefix, ".dmrs.tsv.gz"))
    expect_equal(files$components_file, paste0(prefix, ".dmr_components.tsv"))
    expect_equal(CMEnt:::.reloadViewerCircosCache(list())$output_prefix, NULL)

    write.table(
        data.frame(component_id = 1L, indices = "1,2"),
        files$components_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    withr::defer(unlink(files$components_file), teardown_env())

    reloaded <- CMEnt:::.reloadViewerCircosCache(list(output_prefix = prefix))
    expect_null(reloaded$interactions)
    expect_s3_class(reloaded$components, "data.frame")
    expect_equal(reloaded$components$component_id, 1L)
})

test_that("viewer tabular cache handles missing, empty, valid, and invalid files", {
    missing_file <- tempfile()
    expect_null(CMEnt:::.readViewerTabularCache(missing_file, "missing"))

    empty_file <- tempfile()
    writeLines("", empty_file)
    withr::defer(unlink(empty_file), teardown_env())
    expect_equal(nrow(CMEnt:::.readViewerTabularCache(empty_file, "empty")), 0)

    valid_file <- tempfile()
    writeLines(c("a\tb", "1\t2"), valid_file)
    withr::defer(unlink(valid_file), teardown_env())
    expect_equal(CMEnt:::.readViewerTabularCache(valid_file, "valid")$a, 1L)

    invalid_file <- tempfile()
    dir.create(invalid_file)
    withr::defer(unlink(invalid_file, recursive = TRUE), teardown_env())
    expect_warning(
        expect_null(CMEnt:::.readViewerTabularCache(invalid_file, "invalid")),
        "Failed to load invalid"
    )
})

test_that("viewer metadata validation resolves sample group columns", {
    pheno <- data.frame(
        Sample_Group = c("case", "control"),
        other = c("x", "y"),
        row.names = c("S1", "S2")
    )

    expect_equal(CMEnt:::.resolveViewerSampleGroupCol(pheno), "Sample_Group")
    expect_equal(CMEnt:::.resolveViewerSampleGroupCol(pheno, "other"), "other")
    expect_equal(
        CMEnt:::.resolveViewerSampleGroupCol(data.frame(batch = "A", row.names = "S1")),
        "batch"
    )
    only_casecontrol <- setNames(data.frame(1, row.names = "S1"), "__casecontrol__")
    expect_error(
        CMEnt:::.resolveViewerSampleGroupCol(only_casecontrol),
        "usable sample group"
    )
})

test_that("viewer metadata loader validates shape and normalizes metadata", {
    prefix <- tempfile("viewer-meta-")
    meta_file <- paste0(prefix, ".meta.rds")
    withr::defer(unlink(meta_file), teardown_env())

    saveRDS(list(
        pheno = data.frame(group = c("case", "control"), row.names = c("S1", "S2")),
        genome = "HG19",
        array = "450k",
        sample_group_col = "missing"
    ), meta_file)
    meta <- CMEnt:::.loadCMEntViewerMeta(prefix)
    expect_equal(meta$genome, "hg19")
    expect_equal(meta$array, "450K")
    expect_equal(meta$sample_group_col, "group")

    saveRDS(1L, meta_file)
    expect_error(CMEnt:::.loadCMEntViewerMeta(prefix), "must contain a list")

    pheno <- data.frame(group = "case")
    rownames(pheno) <- ""
    saveRDS(list(pheno = pheno, genome = "hg19"), meta_file)
    expect_error(CMEnt:::.loadCMEntViewerMeta(prefix), "non-empty row names")
})

test_that("viewer asset and task helpers expose expected fallbacks", {
    expect_equal(CMEnt:::.viewerAssetHref("logo.png"), "cment-viewer-assets/logo.png")
    expect_match(as.character(CMEnt:::.viewerBrandUI(NULL)), "cment-viewer-wordmark")

    expect_equal(CMEnt:::.viewerTaskMessage("single_dmr_plot")$message, "Generating DMR plot...")
    expect_equal(CMEnt:::.viewerTaskMessage("unknown")$message, "Processing...")
    expect_true(is.null(CMEnt:::.viewerDevPackagePath()) || dir.exists(CMEnt:::.viewerDevPackagePath()))
})

test_that("viewer plotly panel range supports multiple ggplot build layouts", {
    direct <- list(layout = list(panel_params = list(list(
        `x.range` = c(1, 3),
        `y.range` = c(2, 4)
    ))))
    expect_equal(CMEnt:::.viewerPlotlyPanelRange(direct, "x"), c(1, 3))
    expect_equal(CMEnt:::.viewerPlotlyPanelRange(direct, "y"), c(2, 4))

    nested <- list(layout = list(panel_params = list(list(
        x = list(range = c(10, 20)),
        y = list(continuous_range = c(30, 40))
    ))))
    expect_equal(CMEnt:::.viewerPlotlyPanelRange(nested, "x"), c(10, 20))
    expect_equal(CMEnt:::.viewerPlotlyPanelRange(nested, "y"), c(30, 40))

    expect_null(CMEnt:::.viewerPlotlyPanelRange(list(layout = list(panel_params = list(list()))), "x"))
})
