.reduceAndClearMotifs <- function(dmrs, size = 200) {
    if (length(dmrs) > size) {
        dmrs <- dmrs[seq_len(size)]
    }
    if ("pwm" %in% names(mcols(dmrs))) {
        mcols(dmrs)$pwm <- NULL
        mcols(dmrs)$consensus_seq <- NULL
    }
    dmrs
}

.loadMotifsAndReduce <- function(dmrs) {
    dmrs <- .reduceAndClearMotifs(dmrs)
    mcols(dmrs)$pwm <- lapply(seq_along(dmrs), function(i) {
        matrix(runif(4 * 12), nrow = 4, dimnames = list(c("A", "T", "G", "C"), NULL))
    })
    mcols(dmrs)$consensus_seq <- lapply(seq_along(dmrs), function(i) {
        paste(sample(c("A", "T", "G", "C"), 10, replace = TRUE), collapse = "")
    })
    dmrs
}

test_that("computeDMRsInteraction returns correct structure with valid input", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
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
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
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
    skip_if_missing_motif_extraction_deps(array = "450K", genome = "hg19")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    dmrs <- .reduceAndClearMotifs(dmrs)
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
})

test_that("computeDMRsInteraction handles different similarity thresholds", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    result_high <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.9,
        query_components_with_jaspar = FALSE
    ))

    result_low <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.5,
        query_components_with_jaspar = FALSE
    ))

    expect_type(result_high, "list")
    expect_type(result_low, "list")

    if (!is.null(result_high$interactions) && !is.null(result_low$interactions)) {
        expect_true(nrow(result_low$interactions) >= nrow(result_high$interactions))
    }
})

test_that("computeDMRsInteraction returns NULL when no interactions found", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }
    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.99,
        query_components_with_jaspar = FALSE
    ))

    if (is.null(result$interactions)) {
        expect_null(result$interactions)
    } else {
        expect_true(nrow(result$interactions) == 0 || is.null(result$interactions))
    }
})

test_that("computeDMRsInteraction handles custom motif_site_flank_size", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    result_default <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        motif_site_flank_size = 5,
        query_components_with_jaspar = FALSE
    ))

    result_custom <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        motif_site_flank_size = 10,
        query_components_with_jaspar = FALSE
    ))

    expect_type(result_default, "list")
    expect_type(result_custom, "list")
})

test_that("computeDMRsInteraction does not collapse into a giant component at strict threshold", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
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
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.7,
        query_components_with_jaspar = FALSE
    ))

    if (!is.null(result$interactions) && nrow(result$interactions) > 0) {
        expect_true(all(result$interactions$sim >= -1 & result$interactions$sim <= 1))
        expect_true(all(result$interactions$sim >= 0.7))
    }
})

test_that("computeDMRsInteraction assigns contiguous positive component IDs when scores exist", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)

    mcols(dmrs)$score <- seq_along(dmrs)

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
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
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
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
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
            all(base::strsplit(x, "")[[1]] %in% c("A", "T", "G", "C"))
        })))
    }
})

test_that("computeDMRsInteraction works with different array types", {
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
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
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
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
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "CMEnt", mustWork = FALSE))
    dmrs <- .loadMotifsAndReduce(dmrs)


    result <- suppressWarnings(computeDMRsInteraction(
        dmrs,
        genome = "hg19",
        array = "450K",
        min_similarity = 0.7,
        motif_site_flank_size = 5,
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

test_that("extractDMRMotifs keeps seeds grouped per DMR", {
    skip_if_missing_motif_extraction_deps(array = "450K", genome = "hg38")

    locs <- getSortedGenomicLocs(array = "450K", genome = "hg38")
    chr1_inds <- which(locs$chr == "chr1")
    skip_if(length(chr1_inds) < 6, "Not enough chr1 sites in annotation")

    site_ids <- rownames(locs)[chr1_inds[1:6]]
    dmrs <- data.frame(
        chr = c("chr1", "chr1"),
        start = as.integer(locs[site_ids[c(1, 4)], "start"]),
        end = as.integer(locs[site_ids[c(3, 6)], "start"]),
        start_site = site_ids[c(1, 4)],
        end_site = site_ids[c(3, 6)],
        start_seed = site_ids[c(1, 4)],
        end_seed = site_ids[c(3, 6)],
        seeds = c(
            paste(site_ids[1:3], collapse = ","),
            paste(site_ids[4:6], collapse = ",")
        ),
        stringsAsFactors = FALSE
    )

    out <- suppressWarnings(extractDMRMotifs(dmrs, genome = "hg38", array = "450K"))
    expect_equal(nrow(out), 2)
    expect_true("pwm" %in% colnames(out))
    expect_true(all(vapply(out$pwm, is.matrix, logical(1))))
    expect_true(all(vapply(out$pwm, ncol, integer(1)) == 12))
})

test_that("extractDMRMotifs centers seed windows when DMR starts upstream of first seed", {
    motif_site_flank_size <- 5L
    dmrs <- data.frame(
        chr = "chr1",
        start = 100L,
        end = 140L,
        start_site = "cg_upstream",
        end_site = "cg_seed2",
        start_seed = "cg_seed1",
        end_seed = "cg_seed2",
        seeds = "cg_seed1,cg_seed2",
        stringsAsFactors = FALSE
    )
    beta_locs <- data.frame(
        chr = c("chr1", "chr1"),
        start = c(110L, 130L),
        end = c(111L, 131L),
        row.names = c("cg_seed1", "cg_seed2"),
        stringsAsFactors = FALSE
    )

    seq_len <- dmrs$end - dmrs$start + 1L + motif_site_flank_size + motif_site_flank_size + 1L
    seq_chars <- rep("A", seq_len)
    site_positions <- beta_locs$start - dmrs$start + 1L + motif_site_flank_size
    seq_chars[site_positions] <- "C"
    seq_chars[site_positions + 1L] <- "G"
    sequence <- paste(seq_chars, collapse = "")

    testthat::local_mocked_bindings(
        getDMRSequences = function(...) sequence,
        .package = "CMEnt"
    )

    out <- extractDMRMotifs(
        dmrs,
        genome = "hg38",
        array = NULL,
        beta_locs = beta_locs,
        motif_site_flank_size = motif_site_flank_size
    )

    consensus_seq <- as.character(out$consensus_seq[[1]])
    expect_equal(consensus_seq, "AAAAACGAAAAA")
    expect_equal(substr(consensus_seq, motif_site_flank_size + 1L, motif_site_flank_size + 2L), "CG")
})

test_that("extractDMRMotifs ignores seed windows not centered on C", {
    motif_site_flank_size <- 5L
    dmrs <- data.frame(
        chr = "chr1",
        start = 100L,
        end = 140L,
        start_site = "cg_seed1",
        end_site = "cg_seed3",
        start_seed = "cg_seed1",
        end_seed = "cg_seed3",
        seeds = "cg_seed1,cg_seed2,cg_seed3",
        stringsAsFactors = FALSE
    )
    beta_locs <- data.frame(
        chr = rep("chr1", 3),
        start = c(110L, 120L, 130L),
        end = c(111L, 121L, 131L),
        row.names = c("cg_seed1", "cg_seed2", "cg_seed3"),
        stringsAsFactors = FALSE
    )

    seq_len <- dmrs$end - dmrs$start + 1L + motif_site_flank_size + motif_site_flank_size + 1L
    seq_chars <- rep("A", seq_len)
    center_positions <- beta_locs$start - dmrs$start + 1L + motif_site_flank_size
    seq_chars[center_positions[3]] <- "C"
    sequence <- paste(seq_chars, collapse = "")

    testthat::local_mocked_bindings(
        getDMRSequences = function(...) sequence,
        .package = "CMEnt"
    )

    out <- suppressWarnings(extractDMRMotifs(
        dmrs,
        genome = "hg38",
        array = NULL,
        beta_locs = beta_locs,
        motif_site_flank_size = motif_site_flank_size
    ))

    consensus_seq <- as.character(out$consensus_seq[[1]])
    expect_equal(consensus_seq, "AAAAACAAAAAA")
    expect_equal(substr(consensus_seq, motif_site_flank_size + 1L, motif_site_flank_size + 1L), "C")
})

test_that("motif site windows are normalized when extracted on the complementary orientation", {
    normalized <- CMEnt:::.normalizeMotifSiteSequences(
        c("TTTTTGCTTTTT", "AAAATTGCCCAA"),
        motif_site_flank_size = 5L
    )

    expect_equal(normalized[[1]], "AAAAACGAAAAA")
    expect_equal(substr(normalized[[2]], 6L, 6L), "C")
})

test_that("motif similarity tolerates DMRs without valid PWMs", {
    pwm <- matrix(
        rep(c(0, 1, 0, 0), 12),
        nrow = 4,
        dimnames = list(Biostrings::DNA_BASES, NULL)
    )
    dmrs <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 2),
        ranges = IRanges::IRanges(start = c(100, 200), width = 20),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg38")
    )
    mcols(dmrs)$pwm <- list(pwm, NULL)

    sim <- CMEnt:::.extractMotifsSimilarity(dmrs, motif_site_flank_size = 5)

    expect_equal(dim(sim), c(2L, 2L))
    expect_equal(sim[1, 2], 0)
    expect_equal(sim[2, 1], 0)
})

test_that("getBackgroundArrayMotif uses start-anchored site windows for array probes", {
    skip_if_missing_motif_extraction_deps(array = "450K", genome = "hg19")

    cache_dir <- tempfile("cment-bg-cache-")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    withr::local_options(list(CMEnt.annotation_cache_dir = cache_dir))

    locs <- getSortedGenomicLocs(array = "450K", genome = "hg19")
    skip_if(nrow(locs) == 0, "Annotation locations not available")
    expect_true(all((locs$end - locs$start + 1L) == 2L))

    bg_pwm <- suppressWarnings(getBackgroundArrayMotif(
        genome = "hg19", array = "450K",
        motif_site_flank_size = 5,
        .sorted_locs = locs[seq_len(1000), , drop = FALSE]
    ))

    expect_true(is.matrix(bg_pwm))
    expect_equal(dim(bg_pwm), c(4, 12))
    expect_true(all(is.finite(bg_pwm)))
    expect_true(all(abs(colSums(bg_pwm) - 1) < 1e-8))
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

test_that("computeDMRsInteraction writes complete tabular outputs", {
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

    output_prefix <- file.path(tempdir(), paste0("dmr-interaction-", as.integer(Sys.time())))
    interactions_file <- paste0(output_prefix, ".dmr_interactions.tsv")
    components_file <- paste0(output_prefix, ".dmr_components.tsv")

    result <- computeDMRsInteraction(
        dmrs,
        genome = "hg38",
        array = NULL,
        min_similarity = 0.8,
        min_component_size = 2,
        query_components_with_jaspar = FALSE,
        output_prefix = output_prefix
    )

    expect_true(file.exists(interactions_file))
    expect_true(file.exists(components_file))

    saved_interactions <- read.delim(
        interactions_file,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    saved_components <- read.delim(
        components_file,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    expect_equal(nrow(saved_interactions), nrow(result$interactions))
    expect_equal(nrow(saved_components), nrow(result$components))
    expect_true(all(nzchar(saved_components$indices)))

    unlink(c(interactions_file, components_file))
})
