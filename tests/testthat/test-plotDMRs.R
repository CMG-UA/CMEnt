options("CMEnt.verbose" = 0)
test_that("plotDMR creates a gtable object", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- suppressWarnings(plotDMR(dmrs, dmr_index = 2, output_file = "Rplots.pdf"))

    expect_s3_class(p, "gtable")
    expect_true(inherits(p, "gTree"))
    expect_true(inherits(p, "grob"))
})

test_that("plotDMR handles invalid dmr_index", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    expect_error(
        plotDMR(dmrs, dmr_index = 0),
        "is out of bounds"
    )

    expect_error(
        plotDMR(dmrs, dmr_index = length(dmrs) + 1),
        "is out of bounds"
    )
})


test_that("plotDMR works without a title", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    custom_title <- "Test DMR Title"
    p <- suppressWarnings(plotDMR(dmrs, dmr_index = 1, plot_title = FALSE))

    expect_s3_class(p, "gtable")
})

test_that("plotDMRs creates a combined plot", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    n_dmrs <- min(4, length(dmrs))
    p <- suppressWarnings(plotDMRs(dmrs, dmr_indices = 1:n_dmrs, ncol = 2))

    expect_true(!is.null(p))
    expect_true(inherits(p[[1]], "gtable"))
})

test_that("plotDMRs handles NULL dmr_indices", {
    skip_if_not_installed("ggplot2")
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- suppressWarnings(plotDMRs(dmrs, dmr_indices = NULL))

    expect_true(!is.null(p))
    expect_true(inherits(p[[1]], "gtable"))
})

test_that("plotDMRs respects ncol parameter", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p1 <- suppressWarnings(plotDMRs(dmrs, dmr_indices = 1:4, ncol = 2))
    expect_true(!is.null(p1))

    p2 <- suppressWarnings(plotDMRs(dmrs, dmr_indices = 1:4, ncol = 4))
    expect_true(!is.null(p2))
})


test_that("plotDMR handles DMRs with no extended CpGs", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    dmr_data <- as.data.frame(S4Vectors::mcols(dmrs))
    no_extended_idx <- which(dmr_data$start_seed == dmr_data$start_cpg & dmr_data$end_seed == dmr_data$end_cpg)

    if (length(no_extended_idx) > 0) {
        p <- suppressWarnings(plotDMR(dmrs, dmr_index = no_extended_idx[1]))
        expect_s3_class(p, "gtable")
    } else {
        skip("No DMRs without extended CpGs found")
    }
})

test_that("plotDMR handles DMRs with multiple seeds", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    dmr_data <- as.data.frame(S4Vectors::mcols(dmrs))
    multi_seed_idx <- which(dmr_data$seeds_num >= 3)

    if (length(multi_seed_idx) > 0) {
        p <- suppressWarnings(plotDMR(dmrs, dmr_index = multi_seed_idx[1]))
        expect_s3_class(p, "gtable")
        expect_true(inherits(p, "gTree"))
    } else {
        skip("No DMRs with multiple seeds found")
    }
})

test_that("plotDMR preserves overlapping extension CpG IDs without rowname mangling", {
    skip_if_not_installed("ggplot2")

    cpg_ids <- sprintf("cg%08d", 1:5)
    beta <- matrix(
        c(
            0.1, 0.9,
            0.2, 0.8,
            0.3, 0.7,
            0.4, 0.6,
            0.5, 0.5
        ),
        nrow = length(cpg_ids),
        byrow = TRUE,
        dimnames = list(cpg_ids, c("S1", "S2"))
    )
    sorted_locs <- data.frame(
        chr = rep("chr1", length(cpg_ids)),
        start = seq(100L, 140L, by = 10L),
        end = seq(100L, 140L, by = 10L),
        row.names = cpg_ids,
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
    S4Vectors::mcols(dmrs)$cpgs <- paste(cpg_ids, collapse = ",")
    S4Vectors::mcols(dmrs)$upstream_cpgs <- "cg00000001,cg00000002,cg00000003"
    S4Vectors::mcols(dmrs)$downstream_cpgs <- "cg00000003,cg00000004,cg00000005"
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
    expect_setequal(rownames(ret$total_locs), cpg_ids)

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

test_that("plotDMR plot structure contains expected components", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- suppressWarnings(plotDMR(dmrs, dmr_index = 1))

    expect_true(inherits(p, "gtable"))
    expect_true(inherits(p, "gTree"))
    expect_true(inherits(p, "grob"))
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

test_that("plotDMR with beta and pheno includes PWM plot", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }
    pheno <- loadExampleInputData("pheno")
    beta <- loadExampleInputData("beta")

    if (is.character(pheno) && length(pheno) == 1 && file.exists(pheno)) {
        pheno_env <- new.env(parent = emptyenv())
        load(pheno, envir = pheno_env)
        pheno <- pheno_env$pheno
    }
    if (is.character(beta) && length(beta) == 1 && file.exists(beta)) {
        beta_env <- new.env(parent = emptyenv())
        load(beta, envir = beta_env)
        beta <- beta_env$beta
    }

    if (!is.data.frame(pheno) || !(is.matrix(beta) || is.data.frame(beta))) {
        skip("Example beta/pheno data not available in a plottable format")
    }

    p <- suppressWarnings(plotDMR(
        dmrs = dmrs,
        dmr_index = 1,
        beta = beta,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        genome = "hg19",
        array = "450K"
    ))

    expect_true(inherits(p, "gtable"))

    if ("pwm" %in% colnames(S4Vectors::mcols(dmrs))) {
        expect_true(nrow(p) >= 3)
    }
})
