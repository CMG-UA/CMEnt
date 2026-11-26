test_that("heatmaps retain samples named site when selecting samples", {
    sample_ids <- c("site", "site.1", "ID", "Sample A")
    beta <- matrix(
        seq(0.1, 0.8, length.out = 8), nrow = 2,
        dimnames = list(c("cg1", "cg2"), sample_ids)
    )
    locs <- data.frame(chr = "chr1", start = c(100L, 200L), row.names = rownames(beta))
    pheno <- data.frame(group = c("case", "case", "control", "control"), row.names = sample_ids)

    for (metadata in list(NULL, pheno)) {
        for (limit in c(1L, 10L)) {
            plot <- CMEnt:::.plotBetaHeatmap(
                dmr_data = data.frame(seeds = "cg1,cg2"), beta_data = beta,
                total_shown_positions = locs, pheno = metadata,
                sample_group_col = "group", max_samples_per_group = limit
            )
            retained_samples <- unique(as.character(plot$data$Sample))
            expect_true("site" %in% retained_samples)
            if (limit == 10L) {
                expect_setequal(retained_samples, sample_ids)
            }
            expected <- beta[cbind(
                match(plot$data$site, rownames(beta)),
                match(plot$data$Sample, sample_ids)
            )]
            expect_equal(plot$data$Beta, unname(expected))
            expect_false(anyNA(plot$data$SampleIndex))
            rendered <- ggplot2::ggplot_build(plot)$data[[1]]
            expect_equal(nrow(rendered), nrow(plot$data))
            expect_false(anyNA(rendered$y))
        }
    }
})

test_that("Tabix samples can share names with BED fields", {
    skip_if_not(CMEnt:::.tabixToolsAvailable())
    sample_ids <- c(
        "ID", "Name", "SCORE", "id", "name", "score", "strand",
        "chr", "start", "end", "stop", "#chrom", "site", "Sample A",
        ".CMEnt_sample_1", ".CMEnt_site_id"
    )
    beta <- matrix(
        seq(0.1, 0.9, length.out = 2 * length(sample_ids)), nrow = 2,
        dimnames = list(c("cg1", "cg2"), sample_ids)
    )
    locs <- data.frame(
        chr = "chr1", start = c(100L, 200L), end = c(101L, 201L),
        row.names = rownames(beta)
    )
    query_locs <- locs
    rownames(query_locs) <- paste(query_locs$chr, query_locs$start, sep = ":")
    files <- character()
    withr::defer(unlink(files))

    for (layout in c("BED3", "BED6", "BED6-name")) {
        input_file <- tempfile(fileext = ".tsv")
        tabix_file <- tempfile(fileext = ".bed.gz")
        files <- c(files, input_file, tabix_file, paste0(tabix_file, ".tbi"))
        if (layout != "BED6") {
            input <- data.frame(
                "#chrom" = locs$chr, start = locs$start, end = locs$end, check.names = FALSE
            )
            if (layout == "BED6-name") {
                input$name <- rownames(beta)
                input$score <- 0
                input$strand <- "*"
            }
            input <- data.frame(input, beta, check.names = FALSE)
            write.table(input, input_file, sep = "\t", row.names = FALSE, quote = FALSE)
            convertBetaToTabix(input_file, .bed_file = input_file, output_file = tabix_file)
        } else {
            write.table(beta, input_file, sep = "\t", col.names = NA, quote = FALSE)
            convertBetaToTabix(input_file, sorted_locs = locs, output_file = tabix_file)
        }
        handler <- getBetaHandler(tabix_file, sorted_locs = query_locs)
        expect_identical(handler$getBetaColNames(), sample_ids, info = layout)

        values <- handler$getBeta(
            row_names = rev(rownames(query_locs)), col_names = rev(sample_ids)
        )
        expect_identical(colnames(values), rev(sample_ids), info = layout)
        expect_equal(
            unname(as.matrix(values)),
            unname(beta[2:1, rev(seq_along(sample_ids)), drop = FALSE]),
            info = layout
        )
        values <- handler$getBeta(row_names = rownames(query_locs))
        expect_identical(colnames(values), sample_ids, info = layout)
        expect_equal(unname(as.matrix(values)), unname(beta), info = layout)
    }
})
