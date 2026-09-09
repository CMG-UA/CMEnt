options("CMEnt.verbose" = 0)

plot_fixture <- makeSyntheticPlotFixture()

test_that("plotDMR creates a gtable object from synthetic DMR metadata", {
    skip_if_not_installed("ggplot2")

    p <- suppressWarnings(plotDMR(
        plot_fixture$dmrs,
        dmr_index = 1,
        plot_motif = FALSE,
        beta_locs = plot_fixture$locs,
        array = NULL,
        genome = "hg38"
    ))

    expect_s3_class(p, "gtable")
    expect_true(inherits(p, "gTree"))
    expect_true(inherits(p, "grob"))
})

test_that("plotDMR handles invalid dmr_index", {
    skip_if_not_installed("ggplot2")

    expect_error(
        plotDMR(plot_fixture$dmrs, dmr_index = 0, plot_motif = FALSE, beta_locs = plot_fixture$locs, array = NULL),
        "is out of bounds"
    )
    expect_error(
        plotDMR(plot_fixture$dmrs, dmr_index = length(plot_fixture$dmrs) + 1, plot_motif = FALSE, beta_locs = plot_fixture$locs, array = NULL),
        "is out of bounds"
    )
})

test_that("plotDMR and plotDMRs validate PDF output paths", {
    skip_if_not_installed("ggplot2")

    output_file <- tempfile(fileext = ".png")
    expect_silent(suppressWarnings(plotDMR(plot_fixture$dmrs, 1, plot_motif = FALSE, output_file = output_file, beta_locs = plot_fixture$locs, array = NULL)))
    expect_true(file.exists(output_file))
    expect_error(plotDMRs(plot_fixture$dmrs, 1:2, plot_motif = FALSE, output_file = "dmrs.png"), ".pdf")
})

test_that("plotDMRs returns one plot per requested DMR", {
    skip_if_not_installed("ggplot2")

    p <- suppressWarnings(plotDMRs(
        plot_fixture$dmrs,
        dmr_indices = 1:2,
        ncol = 2,
        plot_motif = FALSE,
        beta_locs = plot_fixture$locs,
        array = NULL,
        genome = "hg38"
    ))

    expect_length(p, 2L)
    expect_true(all(vapply(p, inherits, logical(1), what = "gtable")))
})

test_that("plotDMRs uses top_n when dmr_indices is NULL", {
    skip_if_not_installed("ggplot2")

    p <- suppressWarnings(plotDMRs(
        plot_fixture$dmrs,
        dmr_indices = NULL,
        top_n = 1,
        plot_motif = FALSE,
        beta_locs = plot_fixture$locs,
        array = NULL,
        genome = "hg38"
    ))

    expect_length(p, 1L)
    expect_true(inherits(p[[1]], "gtable"))
})

test_that("plotDMR works without a title", {
    skip_if_not_installed("ggplot2")

    p <- suppressWarnings(plotDMR(
        plot_fixture$dmrs,
        dmr_index = 1,
        plot_title = FALSE,
        plot_motif = FALSE,
        beta_locs = plot_fixture$locs,
        array = NULL,
        genome = "hg38"
    ))

    expect_s3_class(p, "gtable")
})

test_that(".plotDMRStructure retains only seed sites when no extension sites are present", {
    skip_if_not_installed("ggplot2")

    ret <- CMEnt:::.plotDMRStructure(
        dmrs = plot_fixture$dmrs,
        dmr_index = 1,
        beta_locs = plot_fixture$locs,
        plot_title = FALSE,
        .ret_details = TRUE
    )

    expect_setequal(rownames(ret$total_locs), c("cgA", "cgB", "cgC"))
})

test_that(".plotDMRStructure accepts character-encoded genomic positions", {
    skip_if_not_installed("ggplot2")

    dmrs <- plot_fixture$dmrs
    S4Vectors::mcols(dmrs)$start_seed_pos <- as.character(S4Vectors::mcols(dmrs)$start_seed_pos)
    S4Vectors::mcols(dmrs)$end_seed_pos <- as.character(S4Vectors::mcols(dmrs)$end_seed_pos)
    locs <- plot_fixture$locs
    locs$start <- as.character(locs$start)
    locs$end <- as.character(locs$end)

    ret <- CMEnt:::.plotDMRStructure(
        dmrs = dmrs,
        dmr_index = 1,
        beta_locs = locs,
        plot_title = FALSE,
        .ret_details = TRUE
    )

    expect_no_error(ggplot2::ggplot_build(ret$structure_plot))
    expect_type(ret$total_locs$start, "integer")
})

test_that("plotDMR accepts data.frame DMR input", {
    skip_if_not_installed("ggplot2")

    p <- suppressWarnings(plotDMR(
        as.data.frame(plot_fixture$dmrs),
        dmr_index = 1,
        plot_motif = FALSE,
        beta_locs = plot_fixture$locs,
        array = NULL,
        genome = "hg38"
    ))

    expect_s3_class(p, "gtable")
})

test_that("plotDMR preserves overlapping extension site IDs without rowname mangling", {
    skip_if_not_installed("ggplot2")

    site_ids <- sprintf("cg%08d", 1:5)
    beta <- matrix(
        c(
            0.1, 0.9,
            0.2, 0.8,
            0.3, 0.7,
            0.4, 0.6,
            0.5, 0.5
        ),
        nrow = length(site_ids),
        byrow = TRUE,
        dimnames = list(site_ids, c("S1", "S2"))
    )
    sorted_locs <- data.frame(
        chr = rep("chr1", length(site_ids)),
        start = seq(100L, 140L, by = 10L),
        end = seq(100L, 140L, by = 10L),
        row.names = site_ids,
        stringsAsFactors = FALSE
    )
    beta_handler <- getBetaHandler(beta = beta, sorted_locs = sorted_locs)
    pheno <- data.frame(
        Sample_Group = c("case", "control"),
        row.names = c("S1", "S2"),
        stringsAsFactors = FALSE
    )

    dmrs <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 100L, end = 140L),
        seqinfo = GenomeInfoDb::Seqinfo(seqnames = "chr1", genome = "hg19")
    )
    S4Vectors::mcols(dmrs)$seeds <- "cg00000002,cg00000004"
    S4Vectors::mcols(dmrs)$sites <- paste(site_ids, collapse = ",")
    S4Vectors::mcols(dmrs)$upstream_sites <- "cg00000001,cg00000002,cg00000003"
    S4Vectors::mcols(dmrs)$downstream_sites <- "cg00000003,cg00000004,cg00000005"
    S4Vectors::mcols(dmrs)$start_seed_pos <- 110L
    S4Vectors::mcols(dmrs)$end_seed_pos <- 130L
    S4Vectors::mcols(dmrs)$score <- 0.75
    S4Vectors::mcols(dmrs)$cv_accuracy <- 0.8
    S4Vectors::mcols(dmrs)$seeds_num <- 2L
    S4Vectors::mcols(dmrs)$delta_beta <- 0.3
    S4Vectors::mcols(dmrs)$in_promoter_of <- NA_character_
    S4Vectors::mcols(dmrs)$in_gene_body_of <- NA_character_

    ret <- CMEnt:::.plotDMRStructure(
        dmrs = dmrs,
        dmr_index = 1,
        beta_locs = sorted_locs,
        plot_title = FALSE,
        .ret_details = TRUE
    )

    expect_true("cg00000003" %in% rownames(ret$total_locs))
    expect_false("cg000000031" %in% rownames(ret$total_locs))
    expect_setequal(rownames(ret$total_locs), site_ids)

    expect_no_error(
        suppressWarnings(plotDMR(
            dmrs = dmrs,
            dmr_index = 1,
            beta = beta_handler,
            pheno = pheno,
            genome = "hg19",
            array = NULL,
            sample_group_col = "Sample_Group",
            plot_motif = FALSE,
            plot_title = FALSE
        ))
    )
})

test_that("plotDMR drops unresolved extension site IDs from total_locs", {
    skip_if_not_installed("ggplot2")

    dmrs <- plot_fixture$dmrs[1]
    S4Vectors::mcols(dmrs)$upstream_sites <- "cgA,cgMissingUpstream"
    S4Vectors::mcols(dmrs)$downstream_sites <- "cgMissingDownstream,cgC"

    ret <- NULL
    expect_warning(
        ret <- CMEnt:::.plotDMRStructure(
            dmrs = dmrs,
            dmr_index = 1,
            beta_locs = plot_fixture$locs,
            plot_title = FALSE,
            .ret_details = TRUE
        ),
        "dropping 2 extension site ID"
    )

    expect_setequal(rownames(ret$total_locs), c("cgA", "cgB", "cgC"))
    expect_false(anyNA(ret$total_locs$start))
    expect_false(any(is.na(rownames(ret$total_locs))))
})

test_that("plotDMR derives plot span from all resolved merged metadata IDs", {
    skip_if_not_installed("ggplot2")

    dmrs <- plot_fixture$dmrs[1]
    S4Vectors::mcols(dmrs)$seeds <- "cgA,cgC"
    S4Vectors::mcols(dmrs)$sites <- "cgA,cgB,cgC"
    S4Vectors::mcols(dmrs)$upstream_sites <- "cgC"
    S4Vectors::mcols(dmrs)$downstream_sites <- "cgA"

    ret <- CMEnt:::.plotDMRStructure(
        dmrs = dmrs,
        dmr_index = 1,
        beta_locs = plot_fixture$locs,
        plot_title = FALSE,
        .ret_details = TRUE
    )

    expect_setequal(rownames(ret$total_locs), c("cgA", "cgC"))
    expect_false(anyNA(ret$total_locs$start))
    expect_true(all(c(100L, 300L) %in% ret$breaks))
})

test_that("plotDMR ignores inside-seed extension metadata for extension labels", {
    skip_if_not_installed("ggplot2")

    dmrs <- plot_fixture$dmrs[1]
    S4Vectors::mcols(dmrs)$seeds <- "cgA,cgC"
    S4Vectors::mcols(dmrs)$sites <- "cgA,cgB,cgC"
    S4Vectors::mcols(dmrs)$upstream_sites <- ""
    S4Vectors::mcols(dmrs)$downstream_sites <- "cgB"

    expect_no_error(
        CMEnt:::.plotDMRStructure(
            dmrs = dmrs,
            dmr_index = 1,
            beta_locs = plot_fixture$locs,
            plot_title = FALSE,
            .ret_details = TRUE
        )
    )
})

test_that("plotDMR plot structure contains expected components", {
    skip_if_not_installed("ggplot2")

    p <- suppressWarnings(plotDMR(
        plot_fixture$dmrs,
        dmr_index = 1,
        plot_motif = FALSE,
        beta_locs = plot_fixture$locs,
        array = NULL,
        genome = "hg38"
    ))

    expect_true(inherits(p, "gtable"))
    expect_true(inherits(p, "gTree"))
    expect_true(inherits(p, "grob"))
})

test_that(".plotDMRStructure includes existing gene overlap metadata in title", {
    skip_if_not_installed("ggplot2")

    dmrs <- plot_fixture$dmrs
    S4Vectors::mcols(dmrs)$in_promoter_of[1] <- "GENE1"
    S4Vectors::mcols(dmrs)$in_gene_body_of[1] <- "GENE1,GENE2"

    ret <- suppressWarnings(CMEnt:::.plotDMRStructure(
        dmrs,
        dmr_index = 1,
        beta_locs = plot_fixture$locs,
        array = NULL,
        genome = "hg38",
        .ret_details = TRUE
    ))

    expect_match(ret$structure_plot$labels$title, "Overlapping Promoters: GENE1", fixed = TRUE)
    expect_match(ret$structure_plot$labels$title, "Overlapping Gene Bodies: GENE1, GENE2", fixed = TRUE)
})

test_that(".plotPWM clarifies when a motif logo is consensus-only", {
    skip_if_not_installed("ggplot2")

    dmr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = 1, width = 12))
    S4Vectors::mcols(dmr)$pwm <- list(matrix(
        rep(c(1, 0, 0, 0), 12),
        nrow = 4,
        dimnames = list(Biostrings::DNA_BASES, NULL)
    ))
    S4Vectors::mcols(dmr)$consensus_seq <- "AAAAAAAAAAAA"
    S4Vectors::mcols(dmr)$seeds <- "cg00000001"

    p <- CMEnt:::.plotPWM(dmr, genome = "hg19", array = NULL, beta_locs = NULL)

    expect_s3_class(p, "ggplot")
    expect_match(p$labels$subtitle, "Consensus-only logo", fixed = TRUE)
    expect_identical(p$labels$y, "Relative base weight")
})

test_that("plotDMR with beta and pheno accepts precomputed PWM metadata", {
    skip_if_not_installed("ggplot2")

    dmrs <- plot_fixture$dmrs
    S4Vectors::mcols(dmrs)$pwm <- replicate(
        length(dmrs),
        matrix(
            rep(c(0.7, 0.1, 0.1, 0.1), 12),
            nrow = 4,
            dimnames = list(Biostrings::DNA_BASES, NULL)
        ),
        simplify = FALSE
    )
    S4Vectors::mcols(dmrs)$consensus_seq <- c("AAAAAAAAAAAA", "CCCCCCCCCCCC")

    p <- suppressWarnings(plotDMR(
        dmrs = dmrs,
        dmr_index = 1,
        beta = plot_fixture$beta_handler,
        pheno = plot_fixture$pheno,
        genome = "hg38",
        array = NULL,
        sample_group_col = "Sample_Group"
    ))

    expect_s3_class(p, "gtable")
})

test_that(".plotBetaHeatmap omits missing tiles and prioritizes covered seed samples", {
    skip_if_not_installed("ggplot2")

    site_ids <- c("cgA", "cgB", "cgC")
    samples <- paste0(rep(c("A", "B"), each = 4), seq_len(4))
    beta <- matrix(
        NA_real_,
        nrow = length(site_ids),
        ncol = length(samples),
        dimnames = list(site_ids, samples)
    )
    beta[c("cgA", "cgB"), c("A1", "A2", "B1", "B2")] <- c(0.1, 0.2, 0.8, 0.9)
    beta["cgC", ] <- seq(0.2, 0.9, length.out = length(samples))
    pheno <- data.frame(
        group = rep(c("A", "B"), each = 4),
        row.names = samples,
        stringsAsFactors = FALSE
    )
    total_locs <- data.frame(
        chr = rep("chr1", length(site_ids)),
        start = c(100L, 200L, 300L),
        row.names = site_ids,
        stringsAsFactors = FALSE
    )

    p <- CMEnt:::.plotBetaHeatmap(
        dmr_data = data.frame(seeds = "cgA,cgB"),
        beta_data = beta,
        total_shown_positions = total_locs,
        pheno = pheno,
        max_samples_per_group = 2,
        sample_group_col = "group"
    )

    expect_false(any(!is.finite(p$data$Beta)))
    expect_setequal(as.character(unique(p$data$Sample)), c("A1", "A2", "B1", "B2"))
    expect_s3_class(p$theme$panel.grid.major.x, "element_blank")
})
