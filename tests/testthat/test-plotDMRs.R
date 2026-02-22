library(testthat)

test_that("plotDMR creates a gtable object", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- suppressWarnings(plotDMR(dmrs, dmr_index = 2, output_file = "Rplots.pdf"))

    expect_s3_class(p, "gtable")
    expect_true(inherits(p, "gTree"))
    expect_true(inherits(p, "grob"))
})

test_that("plotDMR handles invalid dmr_index", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
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

test_that("plotDMR works with different array types", {
    skip_if_not_installed("ggplot2")

    dmrs_450k <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs_450k) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }
    p1 <- suppressWarnings(plotDMR(dmrs_450k, dmr_index = 1, array = "450K", genome = "hg19"))
    expect_s3_class(p1, "gtable")

    dmrs_epic <- remapDMRsArray(dmrs_450k, from_array = "450K", to_array = "EPIC", from_genome = "hg19", to_genome = "hg19")
    p2 <- suppressWarnings(plotDMR(dmrs_epic, dmr_index = 1, array = "EPIC", genome = "hg19"))
    expect_s3_class(p2, "gtable")
})

test_that("plotDMR works with different genome versions", {
    skip_if_not_installed("ggplot2")

    dmrs_hg19 <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs_hg19) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }
    options("DMRsegal.verbose" = 2)
    p1 <- suppressWarnings(plotDMR(dmrs_hg19, dmr_index = 1, array = "450K", genome = "hg19"))
    expect_s3_class(p1, "gtable")

    dmrs_hg38 <- remapDMRsArray(dmrs_hg19, from_array = "450K", to_array = "450K", from_genome = "hg19", to_genome = "hg38")

    p2 <- suppressWarnings(plotDMR(dmrs_hg38, dmr_index = 1, array = "450K", genome = "hg38"))
    expect_s3_class(p2, "gtable")
})

test_that("plotDMR works without a title", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    custom_title <- "Test DMR Title"
    p <- suppressWarnings(plotDMR(dmrs, dmr_index = 1, plot_title = FALSE))

    expect_s3_class(p, "gtable")
})

test_that("plotDMRs creates a combined plot", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    n_dmrs <- min(4, length(dmrs))
    p <- suppressWarnings(plotDMRs(dmrs, dmr_indices = 1:n_dmrs, ncol = 2))

    expect_true(!is.null(p))
    expect_true(inherits(p[[1]], "gtable"))
})

test_that("plotDMRs handles NULL dmr_indices", {
    skip_if_not_installed("ggplot2")
    options("DMRsegal.verbose" = 3)
    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- suppressWarnings(plotDMRs(dmrs, dmr_indices = NULL))

    expect_true(!is.null(p))
    expect_true(inherits(p[[1]], "gtable"))
})

test_that("plotDMRs respects ncol parameter", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p1 <- suppressWarnings(plotDMRs(dmrs, dmr_indices = 1:4, ncol = 2))
    expect_true(!is.null(p1))

    p2 <- suppressWarnings(plotDMRs(dmrs, dmr_indices = 1:4, ncol = 4))
    expect_true(!is.null(p2))
})


test_that("plotDMR handles DMRs with no extended CpGs", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
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

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
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

test_that("plotDMR plot structure contains expected components", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- suppressWarnings(plotDMR(dmrs, dmr_index = 1))

    expect_true(inherits(p, "gtable"))
    expect_true(inherits(p, "gTree"))
    expect_true(inherits(p, "grob"))
})

test_that("plotDMR with beta and pheno includes PWM plot", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    pheno_file <- system.file("data/pheno.rda", package = "DMRsegal")
    beta_file <- system.file("data/beta.rda", package = "DMRsegal")

    if (!file.exists(pheno_file) || !file.exists(beta_file)) {
        pheno_file <- "../../data/pheno.rda"
        beta_file <- "../../data/beta.rda"
    }

    if (!file.exists(pheno_file) || !file.exists(beta_file)) {
        skip("Data files not available")
    }

    load(pheno_file)
    load(beta_file)

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
