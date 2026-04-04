options("DMRsegal.verbose" = 0)
test_that("plotDMRsManhattan returns a ggplot object", {
    skip_if_not_installed("ggplot2")

    dmrs_path <- system.file("extdata/example_output.rds", package = "DMRsegal")
    if (!nzchar(dmrs_path) || !file.exists(dmrs_path)) {
        skip("Benchmark DMRs not available")
    }
    dmrs <- readRDS(dmrs_path)
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("Benchmark DMRs not available")
    }

    p <- plotDMRsManhattan(dmrs, genome = "hg19")
    expect_s3_class(p, "ggplot")
})


test_that("plotDMRsManhattan draws block rectangles when block ids are present", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(100, 180, 260, 100), width = c(40, 40, 40, 40)),
        seqinfo = GenomeInfoDb::Seqinfo(
            seqnames = c("chr1", "chr2"),
            seqlengths = c(1000, 1500),
            genome = "hg19"
        )
    )
    S4Vectors::mcols(dmrs)$score <- c(0.62, 0.71, 0.66, 0.64)
    S4Vectors::mcols(dmrs)$in_promoter_of <- c("GENE1", NA, NA, NA)
    S4Vectors::mcols(dmrs)$in_gene_body_of <- c(NA, "GENE2", NA, NA)
    S4Vectors::mcols(dmrs)$block_id <- c("chr1_block1", "chr1_block1", NA, NA)

    p <- plotDMRsManhattan(dmrs, genome = "hg19")
    geom_classes <- vapply(p$layers, function(layer) class(layer$geom)[1], character(1))

    expect_true("GeomRect" %in% geom_classes)
})

test_that("Manhattan block rectangles are full-height overview bands", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(100, 180, 260), width = c(40, 40, 40)),
        seqinfo = GenomeInfoDb::Seqinfo(
            seqnames = "chr1",
            seqlengths = 1000,
            genome = "hg19"
        )
    )
    S4Vectors::mcols(dmrs)$score <- c(0.62, 0.71, 0.66)
    S4Vectors::mcols(dmrs)$in_promoter_of <- c("GENE1", NA, NA)
    S4Vectors::mcols(dmrs)$block_id <- c("chr1_block1", "chr1_block1", NA)

    built <- ggplot2::ggplot_build(plotDMRsManhattan(dmrs, genome = "hg19"))
    rect_layer <- built$data[[1]]

    expect_true(all(is.infinite(rect_layer$ymin)))
    expect_true(all(is.infinite(rect_layer$ymax)))
    expect_true(all((rect_layer$xmax - rect_layer$xmin) > 0))
})

test_that("Manhattan plotting spans use full chromosome lengths", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(100, 220, 120), width = c(40, 40, 40)),
        seqinfo = GenomeInfoDb::Seqinfo(
            seqnames = c("chr1", "chr2"),
            seqlengths = c(1000, 2000),
            genome = "hg19"
        )
    )
    S4Vectors::mcols(dmrs)$score <- c(0.62, 0.71, 0.64)
    S4Vectors::mcols(dmrs)$in_promoter_of <- c("GENE1", NA, NA)
    S4Vectors::mcols(dmrs)$in_gene_body_of <- c(NA, "GENE2", NA)

    spans <- DMRsegal:::.resolveManhattanPlotSpans(dmrs, genome = "hg19")

    expect_equal(spans$chr, c("chr1", "chr2"))
    expect_equal(spans$plot_start, c(1, 1))
    expect_equal(spans$plot_end, c(1000, 2000))
    expect_equal(spans$offset, c(0, 1000))
})

test_that("Manhattan plotting region controls the plotted scope", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(180, 620, 120), width = c(40, 40, 40)),
        seqinfo = GenomeInfoDb::Seqinfo(
            seqnames = c("chr1", "chr2"),
            seqlengths = c(1000, 2000),
            genome = "hg19"
        )
    )
    S4Vectors::mcols(dmrs)$score <- c(0.62, 0.71, 0.64)
    S4Vectors::mcols(dmrs)$in_promoter_of <- c("GENE1", NA, NA)
    S4Vectors::mcols(dmrs)$in_gene_body_of <- c(NA, "GENE2", NA)

    spans <- DMRsegal:::.resolveManhattanPlotSpans(
        dmrs,
        genome = "hg19",
        region = "chr1:200-800"
    )
    expect_equal(spans$chr, "chr1")
    expect_equal(spans$plot_start, 200)
    expect_equal(spans$plot_end, 800)

    built <- ggplot2::ggplot_build(plotDMRsManhattan(dmrs, genome = "hg19", region = "chr1:200-800"))
    point_layer <- built$data[[length(built$data)]]

    expect_equal(length(unique(point_layer$x)), 2)
    expect_true(all(point_layer$x >= 1))
    expect_true(all(point_layer$x <= 601))
})

test_that(".normalizeCircosRegion parses region strings", {
    region_df <- DMRsegal:::.normalizeCircosRegion("chr7:1000-2500")
    expect_equal(nrow(region_df), 1)
    expect_equal(region_df$chr[1], "chr7")
    expect_equal(region_df$start[1], 1000)
    expect_equal(region_df$end[1], 2500)
})

test_that(".selectCircosInteractions prioritizes JASPAR-matched interactions", {
    link_data <- data.frame(
        component_id = c(1, 1, 2, 2, 3),
        sim = c(0.91, 0.90, 0.89, 0.88, 0.87),
        has_jaspar_match = c(TRUE, FALSE, TRUE, FALSE, FALSE),
        component_best_score = c(1, 1, 2, 2, 3),
        stringsAsFactors = FALSE
    )

    selected <- DMRsegal:::.selectCircosInteractions(link_data, max_components = 3)
    expect_equal(nrow(selected), 3)
    expect_equal(sum(selected$has_jaspar_match), 2)
    expect_true(all(selected$has_jaspar_match[seq_len(2)]))
})

test_that("plotDMRBlockFormation returns a ggplot diagnostics object", {
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

    p <- plotDMRBlockFormation(
        dmrs = dmrs,
        chromosome = "chr1",
        block_gap_mode = "fixed",
        block_gap_fixed_bp = 1e6,
        genome = "hg19"
    )
    expect_s3_class(p, "ggplot")

    geom_classes <- vapply(p$layers, function(layer) class(layer$geom)[1], character(1))
    expect_true("GeomPoint" %in% geom_classes)
    expect_true("GeomLine" %in% geom_classes)
    expect_true("GeomSegment" %in% geom_classes)

    built <- ggplot2::ggplot_build(p)
    line_layer <- built$data[[which(geom_classes == "GeomLine")[1]]]
    expect_equal(sum(diff(line_layer$x[is.finite(line_layer$x)]) < 0, na.rm = TRUE), 0)

    rect_layers <- built$data[which(geom_classes == "GeomRect")]
    expect_true(all(vapply(rect_layers, function(layer) {
        all(is.finite(layer$ymin)) && all(is.finite(layer$ymax))
    }, logical(1))))
})

test_that("plotDMRBlockFormation validates chromosome availability", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 100, width = 40),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
    )
    S4Vectors::mcols(dmrs)$score <- 0.7

    expect_error(
        plotDMRBlockFormation(dmrs = dmrs, chromosome = "chr2", genome = "hg19"),
        "not found"
    )
})
