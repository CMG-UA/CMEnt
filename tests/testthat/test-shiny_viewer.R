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

    meta_file <- paste0(prefix, ".meta.rds")
    saveRDS(list(
        pheno = data.frame(Sample_Name = c("S1", "S2"), Sample_Group = c("case", "control")),
        genome = "hg19",
        array = "450K"
    ), meta_file)

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

test_that("Manhattan UI exposes an optional plotting region input", {
    ui_html <- as.character(DMRsegal:::.manhattanUI("test"))

    expect_match(ui_html, "Plotting Region (optional)", fixed = TRUE)
    expect_match(ui_html, "chr7:100000-2000000", fixed = TRUE)
})

test_that("Viewer UI uses shinycssloaders instead of the custom busy overlay", {
    skip_if_not_installed("shinycssloaders")

    ui <- DMRsegal:::.createViewerUI()
    ui_html <- as.character(ui)

    expect_true(grepl("shiny-spinner-output-container", ui_html, fixed = TRUE))
    expect_false(grepl("dmrsegal-viewer-busy-overlay", ui_html, fixed = TRUE))
    expect_false(grepl("busy-overlay.js", ui_html, fixed = TRUE))
})

test_that("Viewer spinner helpers delegate to shinycssloaders page spinners", {
    show_helper <- paste(deparse(body(DMRsegal:::.viewerShowPageSpinner)), collapse = "\n")
    hide_helper <- paste(deparse(body(DMRsegal:::.viewerHidePageSpinner)), collapse = "\n")

    expect_match(show_helper, "showPageSpinner", fixed = TRUE)
    expect_match(hide_helper, "hidePageSpinner", fixed = TRUE)
})

test_that("Viewer no longer ships the custom busy overlay script", {
    script_path <- test_path("..", "..", "inst", "shiny", "DMRsegalViewer", "www", "busy-overlay.js")
    expect_false(file.exists(script_path))
})

test_that("Circos viewer falls back to lightweight interactive settings without precomputed interactions", {
    params <- list(mode = "auto", method = "components")

    sanitized <- DMRsegal:::.sanitizeViewerCircosParams(
        params,
        has_precomputed_interactions = FALSE
    )

    expect_false(DMRsegal:::.viewerHasPrecomputedCircosInteractions(list()))
    expect_equal(DMRsegal:::.circosViewerAutoMethodChoices(FALSE), c("blocks" = "blocks", "quick" = "quick"))
    expect_equal(sanitized$method, "blocks")
    expect_false(sanitized$query_components_with_jaspar)
    expect_false(sanitized$use_precomputed_interactions)
    expect_true(sanitized$method_fallback)
})

test_that("Circos viewer recognizes usable precomputed interaction caches", {
    precomputed <- list(
        components = data.frame(
            component_id = 1L,
            indices = I(list(1:2))
        ),
        interactions = data.frame(component_id = 1L, index1 = 1L, index2 = 2L)
    )
    empty_cache <- list(
        components = data.frame(),
        interactions = data.frame()
    )
    broken_components <- list(
        components = data.frame(component_id = 1L, indices = NA_character_),
        interactions = data.frame(component_id = 1L, index1 = 1L, index2 = 2L)
    )

    expect_true(DMRsegal:::.viewerHasUsableCircosComponents(precomputed$components))
    expect_false(DMRsegal:::.viewerHasUsableCircosComponents(broken_components$components))
    expect_true(DMRsegal:::.viewerHasPrecomputedCircosCache(precomputed))
    expect_true(DMRsegal:::.viewerHasPrecomputedCircosCache(empty_cache))
    expect_false(DMRsegal:::.viewerHasPrecomputedCircosCache(list()))
    expect_true(DMRsegal:::.viewerHasPrecomputedCircosInteractions(precomputed))
    expect_false(DMRsegal:::.viewerHasPrecomputedCircosInteractions(empty_cache))
    expect_false(DMRsegal:::.viewerHasPrecomputedCircosInteractions(list()))
    expect_false(DMRsegal:::.viewerHasPrecomputedCircosInteractions(broken_components))

    interaction_args <- DMRsegal:::.viewerCircosInteractionArgs(
        list(),
        has_precomputed_interactions = FALSE
    )

    expect_s3_class(interaction_args$components, "data.frame")
    expect_s3_class(interaction_args$interactions, "data.frame")
    expect_equal(nrow(interaction_args$components), 0)
    expect_equal(nrow(interaction_args$interactions), 0)
})

test_that("Circos viewer status UI reflects cache availability and actions", {
    missing_status <- DMRsegal:::.circosInteractionStatusUI(
        ns = shiny::NS("test"),
        has_precomputed_cache = FALSE,
        has_precomputed_interactions = FALSE,
        is_computing_interactions = FALSE
    )
    empty_status <- DMRsegal:::.circosInteractionStatusUI(
        ns = shiny::NS("test"),
        has_precomputed_cache = TRUE,
        has_precomputed_interactions = FALSE,
        is_computing_interactions = FALSE
    )
    ready_status <- DMRsegal:::.circosInteractionStatusUI(
        ns = shiny::NS("test"),
        has_precomputed_cache = TRUE,
        has_precomputed_interactions = TRUE,
        is_computing_interactions = FALSE
    )

    expect_match(as.character(missing_status), "Precomputed Circos cache not found or is incomplete.", fixed = TRUE)
    expect_match(as.character(missing_status), "Compute Interactions", fixed = TRUE)
    expect_match(as.character(empty_status), "Precomputed Circos cache loaded, but no interaction links are available.", fixed = TRUE)
    expect_match(as.character(empty_status), "Recompute Interactions", fixed = TRUE)
    expect_match(as.character(ready_status), "Precomputed Circos cache available.", fixed = TRUE)
})

test_that("Circos UI renders the static in-app viewer", {
    ui <- DMRsegal:::.circosUI("test")
    ui_html <- as.character(ui)

    expect_match(ui_html, "This viewer uses the static circlize plot in the app.", fixed = TRUE)
    expect_match(ui_html, "shiny-plot-output", fixed = TRUE)
    expect_false(grepl("plotly", ui_html, fixed = TRUE))
})

test_that("plotly viewer preserves Manhattan chromosome tick labels", {
    skip_if_not_installed("plotly")

    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(100, 200, 100), width = c(40, 40, 40)),
        seqinfo = GenomeInfoDb::Seqinfo(
            seqnames = c("chr1", "chr2"),
            seqlengths = c(1000, 1500),
            genome = "hg19"
        )
    )
    S4Vectors::mcols(dmrs)$score <- c(0.62, 0.71, 0.64)
    S4Vectors::mcols(dmrs)$in_promoter_of <- c("GENE1", NA, NA)
    S4Vectors::mcols(dmrs)$in_gene_body_of <- c(NA, "GENE2", NA)
    S4Vectors::mcols(dmrs)$block_id <- c("chr1_block1", "chr1_block1", NA)

    widget <- DMRsegal:::.configureViewerPlotly(
        plotDMRsManhattan(dmrs, genome = "hg19")
    )
    built_widget <- plotly::plotly_build(widget)

    expect_equal(unname(built_widget$x$layout$xaxis$ticktext), c("chr1", "chr2"))
})

test_that("viewer background Manhattan task forwards region scope", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(100, 200), width = c(40, 40)),
        seqinfo = GenomeInfoDb::Seqinfo(
            seqnames = "chr1",
            seqlengths = 1000,
            genome = "hg19"
        )
    )
    S4Vectors::mcols(dmrs)$score <- c(0.62, 0.71)
    S4Vectors::mcols(dmrs)$in_promoter_of <- c("GENE1", NA)
    S4Vectors::mcols(dmrs)$in_gene_body_of <- c(NA, "GENE2")

    expect_error(
        DMRsegal:::.viewerRunBackgroundTaskFromData(
            task_type = "manhattan_plot",
            data = list(dmrs = dmrs, genome = "hg19"),
            params = list(
                region = "chr1:800-900",
                point_size = 1.1,
                point_alpha = 0.75
            )
        ),
        "No DMRs overlap the requested plotting region"
    )
})

test_that("plotly viewer converts full-height block rectangles to finite ranges", {
    skip_if_not_installed("plotly")

    x <- c(
        seq(1e6, by = 5e4, length.out = 10),
        seq(2e7, by = 5e4, length.out = 10),
        seq(4e7, by = 5e4, length.out = 10)
    )
    y <- c(
        seq(0.58, 0.74, length.out = 15),
        seq(0.74, 0.58, length.out = 15)
    )
    dmrs <- GenomicRanges::GRanges(
        seqnames = rep("chr1", length(x)),
        ranges = IRanges::IRanges(start = x, width = 100),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
    )
    S4Vectors::mcols(dmrs)$score <- y

    widget <- DMRsegal:::.configureViewerPlotly(
        plotDMRBlockFormation(
            dmrs = dmrs,
            chromosome = "chr1",
            block_gap_mode = "fixed",
            block_gap_fixed_bp = 1e6,
            genome = "hg19"
        )
    )
    built_widget <- plotly::plotly_build(widget)
    rect_traces <- Filter(function(trace) {
        identical(trace$type, "scatter") && identical(trace$fill, "toself")
    }, built_widget$x$data)

    expect_true(length(rect_traces) > 0)
    expect_true(all(vapply(rect_traces, function(trace) all(is.finite(trace$y)), logical(1))))
})

test_that("plotly viewer keeps block-formation line traces ordered", {
    skip_if_not_installed("plotly")

    x <- c(
        seq(1e6, by = 5e4, length.out = 10),
        seq(2e7, by = 5e4, length.out = 10),
        seq(4e7, by = 5e4, length.out = 10)
    )
    y <- c(
        seq(0.58, 0.74, length.out = 15),
        seq(0.74, 0.58, length.out = 15)
    )
    dmrs <- GenomicRanges::GRanges(
        seqnames = rep("chr1", length(x)),
        ranges = IRanges::IRanges(start = x, width = 100),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
    )
    S4Vectors::mcols(dmrs)$score <- y

    widget <- DMRsegal:::.configureViewerPlotly(
        plotDMRBlockFormation(
            dmrs = dmrs,
            chromosome = "chr1",
            block_gap_mode = "fixed",
            block_gap_fixed_bp = 1e6,
            genome = "hg19"
        )
    )
    built_widget <- plotly::plotly_build(widget)
    line_traces <- Filter(function(trace) {
        identical(trace$type, "scatter") &&
            identical(trace$mode, "lines") &&
            !identical(trace$fill, "toself")
    }, built_widget$x$data)

    expect_true(length(line_traces) > 0)
    expect_true(all(vapply(line_traces, function(trace) {
        x_vals <- suppressWarnings(as.numeric(trace$x))
        sum(diff(x_vals[is.finite(x_vals)]) < 0, na.rm = TRUE) == 0
    }, logical(1))))
})

test_that("loadViewerCircosCache treats blank cache files as empty tables", {
    temp_dir <- tempdir()
    prefix <- file.path(temp_dir, paste0("blank-circos-cache-", as.integer(Sys.time())))
    interactions_file <- paste0(prefix, ".dmr_interactions.tsv")
    components_file <- paste0(prefix, ".dmr_components.tsv")

    writeLines("", interactions_file)
    writeLines("", components_file)

    cache <- DMRsegal:::.loadViewerCircosCache(prefix)

    expect_s3_class(cache$interactions, "data.frame")
    expect_s3_class(cache$components, "data.frame")
    expect_equal(nrow(cache$interactions), 0)
    expect_equal(nrow(cache$components), 0)
})

test_that("Circos heatmap prep skips DMR CpGs missing from viewer beta values", {
    beta <- matrix(
        c(
            0.5, 0.4,
            0.6, 0.5,
            0.7, 0.6
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(c("cg1", "cg2", "cg3"), c("S1", "S2"))
    )
    sorted_locs <- data.frame(
        chr = c("chr1", "chr1", "chr1"),
        start = c(100L, 200L, 300L),
        end = c(100L, 200L, 300L),
        row.names = c("cg1", "cg2", "cg3")
    )
    beta_handler <- DMRsegal::getBetaHandler(beta = beta, sorted_locs = sorted_locs)
    pheno <- data.frame(
        Sample_Group = c("case", "control"),
        row.names = c("S1", "S2")
    )
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(100L, 200L), end = c(350L, 450L)),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
    )
    S4Vectors::mcols(dmrs)$cpgs <- c("cg1,cg2,cg4", "cg2,cg3,cg5")

    heatmap_data <- NULL
    expect_warning(
        heatmap_data <- DMRsegal:::.prepareCircosHeatmapData(
            dmrs = dmrs,
            beta_handler = beta_handler,
            pheno = pheno,
            sample_group_col = "Sample_Group",
            max_num_samples_per_group = 0
        ),
        "Skipping 2 supporting CpGs absent from beta values"
    )

    expect_s3_class(heatmap_data$heatmap_df, "data.frame")
    expect_equal(rownames(heatmap_data$heatmap_df), c("cg1", "cg2", "cg3"))
    expect_true(all(stats::complete.cases(heatmap_data$heatmap_df[, c("chr", "start", "end")])))
})

test_that("plotDMR restores showtext auto hooks after drawing", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")
    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    orig_showtext_auto <- DMRsegal:::.isShowtextAutoEnabled()
    orig_device <- grDevices::dev.cur()
    on.exit({
        showtext::showtext_auto(enable = orig_showtext_auto)
        while (grDevices::dev.cur() > orig_device && grDevices::dev.cur() > 1) {
            grDevices::dev.off()
        }
    }, add = TRUE)

    showtext::showtext_auto(enable = FALSE)
    expect_false(DMRsegal:::.isShowtextAutoEnabled())

    expect_no_error(
        plotDMR(
            dmrs = dmrs[seq_len(min(2, length(dmrs)))],
            dmr_index = 1,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            max_cpgs = 20,
            max_samples_per_group = 4
        )
    )

    expect_false(DMRsegal:::.isShowtextAutoEnabled())
})
