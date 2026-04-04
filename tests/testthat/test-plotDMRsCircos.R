options("DMRsegal.verbose" = 0)

test_that("plotDMRsCircos creates a circos plot", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[seq_len(min(5, length(dmrs))), drop = FALSE]
    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            min_similarity = 0.8
        )
    )
})

test_that("plotDMRsCircos works with interactions", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[as.character(dmrs@seqnames) %in% c("chr5", "chr11"), , drop = FALSE]

    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            min_similarity = 0.8
        )
    )
})

test_that("plotDMRsCircos handles BetaHandler input", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }
    
    dmrs_subset <- dmrs[seq_len(min(3, length(dmrs)))]

    beta_handler <- getBetaHandler(beta, array = array_type, genome = "hg19")

    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta_handler,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group"
        )
    )
})

test_that("plotDMRsCircos handles data frame DMRs input", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[seq_len(min(3, length(dmrs)))]
    dmrs_df <- convertToDataFrame(dmrs_subset)

    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs_df,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group"
        )
    )
})

test_that("plotDMRsCircos validates inputs", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))

    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[seq_len(min(3, length(dmrs)))]

    expect_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = "nonexistent_file.txt",
            pheno = pheno,
            genome = "hg19"
        ),
        "does not exist"
    )

    bad_pheno <- pheno
    colnames(bad_pheno) <- c("col1", "col2", "col3")

    expect_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = bad_pheno,
            genome = "hg19",
            sample_group_col = "Sample_Group"
        ),
        "not found in pheno"
    )
})

test_that("plotDMRsCircos supports chromosome and region filters", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[seq_len(min(12, length(dmrs)))]
    target_chr <- as.character(GenomicRanges::seqnames(dmrs_subset)[1])
    target_start <- GenomicRanges::start(dmrs_subset)[1]
    target_end <- GenomicRanges::end(dmrs_subset)[1]
    target_region <- sprintf("%s:%d-%d", target_chr, target_start, target_end)

    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            chromosomes = target_chr,
            region = target_region,
            max_components = 5
        )
    )

    expect_error(
        plotDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            region = "chr5-100-200"
        ),
        "chr:start-end"
    )
})

test_that(".filterDMRsByScopeForCircos preserves GRanges when some regions are empty", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(100L, 1000L), width = c(80L, 80L))
    )

    scoped <- DMRsegal:::.filterDMRsByScopeForCircos(
        dmrs,
        region_df = data.frame(
            chr = c("chr3", "chr1"),
            start = c(1L, 90L),
            end = c(100L, 200L),
            stringsAsFactors = FALSE
        )
    )

    expect_s4_class(scoped, "GRanges")
    expect_length(scoped, 1L)
    expect_equal(as.character(GenomicRanges::seqnames(scoped)), "chr1:90-200")
    expect_equal(GenomicRanges::start(scoped), 100L)
    expect_equal(GenomicRanges::end(scoped), 179L)
})

test_that("plotDMRsCircos extracts motifs only for scoped DMRs", {
    skip_if_not_installed("mockery")
    library(mockery)

    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(100L, 1000L), width = c(80L, 80L)),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
    )
    S4Vectors::mcols(dmrs)$delta_beta <- c(0.3, -0.2)
    S4Vectors::mcols(dmrs)$cpgs <- c("cg1,cg2", "cg3,cg4")
    S4Vectors::mcols(dmrs)$seeds <- c("cg1,cg2", "cg3,cg4")

    beta <- matrix(
        c(0.4, 0.5, 0.6, 0.7),
        ncol = 1,
        dimnames = list(c("cg1", "cg2", "cg3", "cg4"), "S1")
    )
    sorted_locs <- data.frame(
        chr = c("chr1", "chr1", "chr2", "chr2"),
        start = c(100L, 140L, 1000L, 1040L),
        end = c(100L, 140L, 1000L, 1040L),
        row.names = c("cg1", "cg2", "cg3", "cg4")
    )
    beta_handler <- getBetaHandler(beta = beta, sorted_locs = sorted_locs)
    pheno <- data.frame(Sample_Group = "case", row.names = "S1")

    extracted_n <- NA_integer_
    stub(
        plotDMRsCircos,
        "extractDMRMotifs",
        function(dmrs, ...) {
            extracted_n <<- length(dmrs)
            S4Vectors::mcols(dmrs)$pwm <- replicate(
                length(dmrs),
                matrix(0.25, nrow = 4, ncol = 10),
                simplify = FALSE
            )
            dmrs
        }
    )
    stub(plotDMRsCircos, ".prepareCircosLinkData", function(...) NULL)
    stub(
        plotDMRsCircos,
        ".getCytobandData",
        function(genome) {
            data.frame(
                V1 = c("chr1", "chr2"),
                V2 = c(1L, 1L),
                V3 = c(1e6L, 1e6L),
                V4 = c("p", "p"),
                V5 = c("gneg", "gneg")
            )
        }
    )

    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs,
            beta = beta_handler,
            pheno = pheno,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            region = "chr1:50-250"
        )
    )
    expect_equal(extracted_n, 1L)
})

test_that("plotDMRsCircos skips non-drawable gneg-only ideograms without warning", {
    skip_if_not_installed("mockery")
    library(mockery)

    dmrs <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 100L, width = 80L),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
    )
    S4Vectors::mcols(dmrs)$delta_beta <- 0.3
    S4Vectors::mcols(dmrs)$cpgs <- "cg1,cg2"
    S4Vectors::mcols(dmrs)$seeds <- "cg1,cg2"

    beta <- matrix(
        c(0.4, 0.5),
        ncol = 1,
        dimnames = list(c("cg1", "cg2"), "S1")
    )
    sorted_locs <- data.frame(
        chr = c("chr1", "chr1"),
        start = c(100L, 140L),
        end = c(100L, 140L),
        row.names = c("cg1", "cg2")
    )
    beta_handler <- getBetaHandler(beta = beta, sorted_locs = sorted_locs)
    pheno <- data.frame(Sample_Group = "case", row.names = "S1")

    stub(
        plotDMRsCircos,
        "extractDMRMotifs",
        function(dmrs, ...) {
            S4Vectors::mcols(dmrs)$pwm <- list(matrix(0.25, nrow = 4, ncol = 10))
            dmrs
        }
    )
    stub(plotDMRsCircos, ".prepareCircosLinkData", function(...) NULL)
    stub(
        plotDMRsCircos,
        ".getCytobandData",
        function(genome) {
            data.frame(
                V1 = "chr1",
                V2 = 1L,
                V3 = 1e6L,
                V4 = "p",
                V5 = "gneg"
            )
        }
    )

    expect_no_warning(
        plotDMRsCircos(
            dmrs = dmrs,
            beta = beta_handler,
            pheno = pheno,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            region = "chr1:50-250"
        )
    )
})

test_that("plotDMRsCircos reuses precomputed interactions without extracting motifs again", {
    skip_if_not_installed("mockery")
    library(mockery)

    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(100L, 1000L), width = c(80L, 80L)),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
    )
    S4Vectors::mcols(dmrs)$delta_beta <- c(0.3, -0.2)
    S4Vectors::mcols(dmrs)$cpgs <- c("cg1,cg2", "cg3,cg4")
    S4Vectors::mcols(dmrs)$seeds <- c("cg1,cg2", "cg3,cg4")

    beta <- matrix(
        c(0.4, 0.5, 0.6, 0.7),
        ncol = 1,
        dimnames = list(c("cg1", "cg2", "cg3", "cg4"), "S1")
    )
    sorted_locs <- data.frame(
        chr = c("chr1", "chr1", "chr1", "chr1"),
        start = c(100L, 140L, 1000L, 1040L),
        end = c(100L, 140L, 1000L, 1040L),
        row.names = c("cg1", "cg2", "cg3", "cg4")
    )
    beta_handler <- getBetaHandler(beta = beta, sorted_locs = sorted_locs)
    pheno <- data.frame(Sample_Group = "case", row.names = "S1")

    extracted_motifs <- FALSE
    stub(
        plotDMRsCircos,
        "extractDMRMotifs",
        function(...) {
            extracted_motifs <<- TRUE
            stop("extractDMRMotifs should not run when precomputed interactions are supplied")
        }
    )
    stub(
        plotDMRsCircos,
        ".getCytobandData",
        function(genome) {
            data.frame(
                V1 = "chr1",
                V2 = 1L,
                V3 = 2e6L,
                V4 = "p",
                V5 = "gneg"
            )
        }
    )

    components <- data.frame(
        component_id = 1L,
        size = 2L,
        indices = "1,2",
        stringsAsFactors = FALSE
    )
    interactions <- data.frame(
        index1 = 1L,
        chr1 = "chr1",
        start1 = 100L,
        end1 = 180L,
        index2 = 2L,
        chr2 = "chr1",
        start2 = 1000L,
        end2 = 1080L,
        sim = 0.95,
        component_id = 1L,
        stringsAsFactors = FALSE
    )

    expect_no_error(
        plotDMRsCircos(
            dmrs = dmrs,
            beta = beta_handler,
            pheno = pheno,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            components = components,
            interactions = interactions,
            query_components_with_jaspar = FALSE
        )
    )
    expect_false(extracted_motifs)
})

test_that(".selectCircosRegions respects region caps and block priority", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr1", "chr2", "chr2", "chr3"),
        ranges = IRanges::IRanges(
            start = c(100, 300, 9000000, 100, 500, 100),
            width = c(50, 50, 50, 60, 60, 40)
        ),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
    )
    S4Vectors::mcols(dmrs)$score <- c(1, 2, 6, 3, 4, 5)
    S4Vectors::mcols(dmrs)$delta_beta <- c(0.4, 0.3, 0.1, 0.35, 0.34, 0.2)
    S4Vectors::mcols(dmrs)$block_id <- c("chr1_b1", "chr1_b1", "chr1_b2", "chr2_b1", "chr2_b1", NA)

    selected <- DMRsegal:::.selectCircosRegions(
        dmrs = dmrs,
        method = "blocks",
        n_regions = 3,
        region_flank_bp = 100,
        max_regions_per_chr = 1,
        min_inter_region_bp = 1000
    )

    expect_s3_class(selected, "data.frame")
    expect_lte(nrow(selected), 3)
    expect_true(all(table(selected$chr) <= 1))
    expect_equal(selected$chr[1], "chr1")
})

test_that(".selectCircosRegions supports component and hybrid candidate selection", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr2", "chr2", "chr3"),
        ranges = IRanges::IRanges(
            start = c(100, 400, 150, 600, 200),
            width = c(50, 50, 60, 60, 40)
        ),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
    )
    S4Vectors::mcols(dmrs)$score <- c(1, 2, 4, 5, 3)
    S4Vectors::mcols(dmrs)$delta_beta <- c(0.45, 0.41, 0.3, 0.28, 0.34)
    S4Vectors::mcols(dmrs)$block_id <- c("chr1_b1", "chr1_b1", NA, NA, "chr3_b1")

    components <- data.frame(
        component_id = c(1, 2),
        size = c(3, 2),
        jaspar_names = c("TFAP2A", NA),
        stringsAsFactors = FALSE
    )
    components$indices <- list(c(1, 3, 5), c(2, 4))

    selected_components <- DMRsegal:::.selectCircosRegions(
        dmrs = dmrs,
        method = "components",
        n_regions = 3,
        region_flank_bp = 100,
        max_regions_per_chr = 2,
        min_inter_region_bp = 1000,
        components = components
    )
    selected_hybrid <- DMRsegal:::.selectCircosRegions(
        dmrs = dmrs,
        method = "hybrid",
        n_regions = 3,
        region_flank_bp = 100,
        max_regions_per_chr = 2,
        min_inter_region_bp = 1000,
        components = components
    )

    expect_s3_class(selected_components, "data.frame")
    expect_true(nrow(selected_components) >= 1)
    expect_true(all(c("chr", "start", "end") %in% colnames(selected_components)))
    expect_s3_class(selected_hybrid, "data.frame")
    expect_true(nrow(selected_hybrid) >= 1)
    expect_true(any(selected_hybrid$chr %in% c("chr1", "chr3")))
})

test_that(".selectCircosRegions accepts serialized component index strings", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr2"),
        ranges = IRanges::IRanges(
            start = c(100, 400, 150),
            width = c(50, 50, 60)
        ),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
    )
    S4Vectors::mcols(dmrs)$score <- c(1, 2, 3)
    S4Vectors::mcols(dmrs)$delta_beta <- c(0.45, 0.41, 0.3)

    components <- data.frame(
        component_id = 1L,
        size = 3L,
        indices = "1,3",
        stringsAsFactors = FALSE
    )

    selected <- DMRsegal:::.selectCircosRegions(
        dmrs = dmrs,
        method = "components",
        n_regions = 2,
        region_flank_bp = 100,
        max_regions_per_chr = 2,
        min_inter_region_bp = 1000,
        components = components
    )

    expect_s3_class(selected, "data.frame")
    expect_true(nrow(selected) >= 1)
    expect_true(any(selected$chr %in% c("chr1", "chr2")))
})

test_that("plotAutoDMRsCircos returns selected regions invisibly", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    dmrs_subset <- dmrs[seq_len(min(12, length(dmrs)))]
    selected <- NULL
    expect_no_error(
        selected <- plotAutoDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            method = "quick",
            n_regions = 2,
            max_regions_per_chr = 1,
            query_components_with_jaspar = FALSE
        )
    )

    expect_s3_class(selected, "data.frame")
    expect_true(all(c("chr", "start", "end") %in% colnames(selected)))
    expect_lte(nrow(selected), 2)
})

test_that("plotAutoDMRsCircos forwards plot arguments through dots", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }

    expect_false("max_dmrs_per_chr" %in% names(formals(plotAutoDMRsCircos)))

    dmrs_subset <- dmrs[seq_len(min(12, length(dmrs)))]
    expect_no_error(
        selected <- plotAutoDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            sample_group_col = "Sample_Group",
            method = "quick",
            n_regions = 1,
            max_regions_per_chr = 1,
            query_components_with_jaspar = FALSE,
            max_dmrs_per_chr = 1,
            max_cpgs_per_dmr = 1,
            max_num_samples_per_group = 2
        )
    )
    expect_s3_class(selected, "data.frame")
    expect_error(
        plotAutoDMRsCircos(
            dmrs = dmrs_subset,
            beta = beta,
            pheno = pheno,
            region = data.frame(chr = "chr5", start = 1L, end = 10L)
        ),
        "manages `region` directly",
        fixed = TRUE
    )
})

test_that("plotAutoDMRsCircos supports components and hybrid selection", {
    beta <- loadExampleInputDataChr5And11("beta")
    pheno <- loadExampleInputDataChr5And11("pheno")
    array_type <- loadExampleInputDataChr5And11("array_type")

    dmrs <- readRDS(system.file("extdata/example_outputChr5And11.rds", package = "DMRsegal"))
    if (is.null(dmrs) || length(dmrs) == 0) {
        skip("No DMRs available for testing")
    }
    dmrs_subset <- dmrs[as.character(GenomicRanges::seqnames(dmrs)) %in% c("chr5", "chr11"), ]

    interaction_state <- suppressWarnings(computeDMRsInteraction(
        dmrs = dmrs_subset,
        genome = "hg19",
        array = array_type,
        min_similarity = 0.8,
        query_components_with_jaspar = FALSE
    ))
    if (is.null(interaction_state$components) || nrow(interaction_state$components) == 0) {
        skip("No interaction components available for testing")
    }

    selected_components <- NULL
    selected_hybrid <- NULL
    expect_no_error(
        selected_components <- plotAutoDMRsCircos(
            dmrs = interaction_state$dmrs,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            components = interaction_state$components,
            interactions = interaction_state$interactions,
            sample_group_col = "Sample_Group",
            method = "components",
            n_regions = 2,
            max_regions_per_chr = 1,
            query_components_with_jaspar = FALSE
        )
    )
    expect_no_error(
        selected_hybrid <- plotAutoDMRsCircos(
            dmrs = interaction_state$dmrs,
            beta = beta,
            pheno = pheno,
            array = array_type,
            genome = "hg19",
            components = interaction_state$components,
            interactions = interaction_state$interactions,
            sample_group_col = "Sample_Group",
            method = "hybrid",
            n_regions = 2,
            max_regions_per_chr = 1,
            query_components_with_jaspar = FALSE
        )
    )

    expect_s3_class(selected_components, "data.frame")
    expect_s3_class(selected_hybrid, "data.frame")
    expect_true(all(c("chr", "start", "end") %in% colnames(selected_components)))
    expect_true(all(c("chr", "start", "end") %in% colnames(selected_hybrid)))
})
