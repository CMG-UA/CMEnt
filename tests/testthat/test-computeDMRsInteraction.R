library(testthat)

.loadMotifsAndReduce <- function(dmrs) {
    dmrs <- dmrs[seq_len(min(200, length(dmrs)))]
    if (! "pwm" %in% names(mcols(dmrs))) {
        mcols(dmrs)$pwm <- lapply(seq_along(dmrs), function(i) {
            matrix(runif(4 * 12), nrow = 4, dimnames = list(c("A", "T", "G", "C"), NULL))
        })
        mcols(dmrs)$consensus_seq <- lapply(seq_along(dmrs), function(i) {
            paste(sample(c("A", "T", "G", "C"), 10, replace = TRUE), collapse = "")
        })
    }
    dmrs
}

test_that("computeDMRsInteraction returns correct structure with valid input", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)


    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.7,
        query_components_with_jaspar = FALSE
    ))

    expect_type(result, "list")
    expect_true("interactions" %in% names(result))
    expect_true("components" %in% names(result))
    expect_true("dmrs" %in% names(result))

    if (!is.null(result$interactions)) {
        expect_s3_class(result$interactions, "data.frame")
        expect_true(
            all(
                c("chr1", "start1", "end1", "chr2", "start2", "end2", "sim") %in%
                    colnames(result$interactions)
            )
        )
    }

    expect_s3_class(result$components, "data.frame")
    expect_s4_class(result$dmrs, "GRanges")
    expect_true("component_ids" %in% colnames(mcols(result$dmrs)))
    expect_true(
        all(
            c("component_id", "size", "indices", "avg_pwm", "consensus_seq") %in%
                colnames(result$components)
        )
    )
})

test_that("computeDMRsInteraction handles GRanges input", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.7,
        query_components_with_jaspar = FALSE
    ))

    expect_type(result, "list")
    expect_true("interactions" %in% names(result))
    expect_true(nrow(result$interactions) > 0)
    expect_true("components" %in% names(result))
    expect_true(nrow(result$components) > 0)
})

test_that("computeDMRsInteraction works with not precomputed motifs", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))

    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.7
    ))

    expect_type(result, "list")
    expect_true("interactions" %in% names(result))
    expect_true("components" %in% names(result))
})

test_that("computeDMRsInteraction handles different similarity thresholds", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    result_high <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.9,
    ))

    result_low <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.5
    ))

    expect_type(result_high, "list")
    expect_type(result_low, "list")

    if (!is.null(result_high$interactions) && !is.null(result_low$interactions)) {
        expect_true(nrow(result_low$interactions) >= nrow(result_high$interactions))
    }
})

test_that("computeDMRsInteraction returns NULL when no interactions found", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }
    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.99
    ))

    if (is.null(result$interactions)) {
        expect_null(result$interactions)
    } else {
        expect_true(nrow(result$interactions) == 0 || is.null(result$interactions))
    }
})

test_that("computeDMRsInteraction handles custom flank_size", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    result_default <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        flank_size = 5
    ))

    result_custom <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        flank_size = 10
    ))

    expect_type(result_default, "list")
    expect_type(result_custom, "list")
})

test_that("computeDMRsInteraction does not collapse into a giant component at strict threshold", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.8,
        min_component_size = 2,
        query_components_with_jaspar = FALSE
    ))
    expect_true("components" %in% names(result))
    if (nrow(result$components) > 0) {
        expect_true(max(result$components$size) < (0.9 * length(dmrs)))
    }
})

test_that("computeDMRsInteraction validates similarity values", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.7
    ))

    if (!is.null(result$interactions) && nrow(result$interactions) > 0) {
        expect_true(all(result$interactions$sim >= -1 & result$interactions$sim <= 1))
        expect_true(all(result$interactions$sim >= 0.7))
    }
})


test_that("computeDMRsInteraction assigns contiguous positive component IDs when ranks exist", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    mcols(dmrs)$rank <- seq_along(dmrs)

    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.7,
        min_component_size = 1,
        query_components_with_jaspar = FALSE
    ))

    expect_true("component_id" %in% colnames(result$components))
    if (nrow(result$components) > 0) {
        expect_true(all(result$components$component_id >= 1))
        expect_false(any(result$components$component_id == 0))
        expect_equal(result$components$component_id, seq_len(nrow(result$components)))
        expect_true(all(result$components$size > 0))
    }
    if (nrow(result$interactions) > 0) {
        expect_true("component_id" %in% colnames(result$interactions))
    }
})

test_that("computeDMRsInteraction components are ordered by size", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.5,
        query_components_with_jaspar = FALSE
    ))

    expect_true(
        is.unsorted(result$components$size, strictly = FALSE) == FALSE ||
            all(diff(result$components$size) <= 0)
    )
})

test_that("computeDMRsInteraction consensus sequences are valid", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.7,
        query_components_with_jaspar = FALSE
    ))

    if (nrow(result$components) > 0) {
        expect_true(all(sapply(result$components$consensus_seq, function(x) {
            all(strsplit(x, "")[[1]] %in% c("A", "T", "G", "C"))
        })))
    }
})

test_that("computeDMRsInteraction works with different array types", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)


    result_450k <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.7,
        query_components_with_jaspar = FALSE
    ))

    result_epic <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "EPIC",
        min_similarity = 0.7,
        query_components_with_jaspar = FALSE
    ))

    expect_type(result_450k, "list")
    expect_type(result_epic, "list")
})

test_that("computeDMRsInteraction creates plot when plot_dir is specified", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)


    temp_dir <- tempdir()
    plot_dir <- file.path(temp_dir, "test_interaction_plots")

    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.7,
        plot_dir = plot_dir,
        query_components_with_jaspar = FALSE
    ))

    expect_true(dir.exists(plot_dir))

    if (!is.null(result$interactions) && nrow(result$interactions) > 0) {
        expect_true(file.exists(file.path(plot_dir, "dmr_motif_similarity_heatmap.png")))
    }

    unlink(plot_dir, recursive = TRUE)
})

test_that("computeDMRsInteraction avg_pwm has correct dimensions", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)


    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.7,
        flank_size = 5,
        query_components_with_jaspar = FALSE
    ))

    if (nrow(result$components) > 0) {
        for (i in seq_len(nrow(result$components))) {
            pwm <- result$components$avg_pwm[[i]]
            expect_true(is.matrix(pwm))
            expect_equal(nrow(pwm), 4)
            expect_equal(ncol(pwm), 2 * 5 + 2)
        }
    }
})

test_that("computeDMRsInteraction annotates returned DMRs with component_ids", {
    pwm_a <- matrix(
        rep(c(1, 0, 0, 0), 12),
        nrow = 4,
        dimnames = list(c("A", "T", "G", "C"), NULL)
    )
    pwm_t <- matrix(
        rep(c(0, 1, 0, 0), 12),
        nrow = 4,
        dimnames = list(c("A", "T", "G", "C"), NULL)
    )
    pwm_uniform <- matrix(
        0.25,
        nrow = 4,
        ncol = 12,
        dimnames = list(c("A", "T", "G", "C"), NULL)
    )
    dmrs <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 6),
        ranges = IRanges::IRanges(start = seq(1, 600, by = 100), width = 50),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg38")
    )
    mcols(dmrs)$pwm <- list(pwm_a, pwm_a, pwm_a, pwm_t, pwm_t, pwm_uniform)
    mcols(dmrs)$consensus_seq <- rep("AAAAAAAAAAAA", length(dmrs))

    result <- computeDMRsInteraction(
        dmrs,
        genome = "hg38",
        array = NULL,
        min_similarity = 0.8,
        min_component_size = 2,
        query_components_with_jaspar = FALSE
    )

    expect_true("dmrs" %in% names(result))
    expect_s4_class(result$dmrs, "GRanges")
    expect_true("component_ids" %in% colnames(mcols(result$dmrs)))
    expect_equal(
        as.character(mcols(result$dmrs)$component_ids),
        c("1", "1", "1", "2", "2", NA)
    )
})
