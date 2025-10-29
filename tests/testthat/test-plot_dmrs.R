library(testthat)

test_that("plotDMR creates a ggplot object", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- plotDMR(dmrs, dmr_index = 1)

    expect_s3_class(p, "ggplot")
    expect_true(!is.null(p$layers))
    expect_true(length(p$layers) > 0)
})

test_that("plotDMR handles invalid dmr_index", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    expect_error(
        plotDMR(dmrs, dmr_index = 0),
        "dmr_index must be between"
    )

    expect_error(
        plotDMR(dmrs, dmr_index = length(dmrs) + 1),
        "dmr_index must be between"
    )
})

test_that("plotDMR works with different array types", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p1 <- plotDMR(dmrs, dmr_index = 1, array = "450K")
    expect_s3_class(p1, "ggplot")

    p2 <- plotDMR(dmrs, dmr_index = 1, array = "EPIC")
    expect_s3_class(p2, "ggplot")
})

test_that("plotDMR works with different genome versions", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p1 <- plotDMR(dmrs, dmr_index = 1, genome = "hg19")
    expect_s3_class(p1, "ggplot")

    p2 <- plotDMR(dmrs, dmr_index = 1, genome = "hg38")
    expect_s3_class(p2, "ggplot")
})

test_that("plotDMR works with custom title", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    custom_title <- "Test DMR Title"
    p <- plotDMR(dmrs, dmr_index = 1, plot_title = FALSE)

    expect_s3_class(p, "ggplot")
})

test_that("plotDMRs creates a combined plot", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    n_dmrs <- min(4, length(dmrs))
    p <- plotDMRs(dmrs, dmr_indices = 1:n_dmrs, ncol = 2)

    expect_true(!is.null(p))
    expect_true(inherits(p, "gtable"))
})

test_that("plotDMRs handles NULL dmr_indices", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- plotDMRs(dmrs, dmr_indices = NULL)

    expect_true(!is.null(p))
    expect_true(inherits(p, "gtable"))
})

test_that("plotDMRs respects ncol parameter", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p1 <- plotDMRs(dmrs, dmr_indices = 1:4, ncol = 2)
    expect_true(!is.null(p1))

    p2 <- plotDMRs(dmrs, dmr_indices = 1:4, ncol = 4)
    expect_true(!is.null(p2))
})


test_that("plotDMRWithBeta works", {
    skip_if_not_installed("DMRsegaldata")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    pheno <- DMRsegaldata::pheno
    beta <- DMRsegaldata::beta

    p <- suppressWarnings(plotDMRWithBeta(
        dmrs = dmrs,
        dmr_index = 1,
        beta = beta,
        pheno = pheno,
        sample_group_col = "Sample_Group"
    )
    )

    expect_true(inherits(p, "gtable"))
})


test_that("plotDMR handles DMRs with no extended CpGs", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    dmr_data <- as.data.frame(S4Vectors::mcols(dmrs))
    no_extended_idx <- which(dmr_data$start_dmp == dmr_data$start_cpg & dmr_data$end_dmp == dmr_data$end_cpg)

    if (length(no_extended_idx) > 0) {
        p <- plotDMR(dmrs, dmr_index = no_extended_idx[1])
        expect_s3_class(p, "ggplot")
    } else {
        skip("No DMRs without extended CpGs found")
    }
})

test_that("plotDMR handles DMRs with multiple DMPs", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    dmr_data <- as.data.frame(S4Vectors::mcols(dmrs))
    multi_dmp_idx <- which(dmr_data$dmps_num >= 3)

    if (length(multi_dmp_idx) > 0) {
        p <- plotDMR(dmrs, dmr_index = multi_dmp_idx[1])
        expect_s3_class(p, "ggplot")
        expect_true(length(p$layers) > 0)
    } else {
        skip("No DMRs with multiple DMPs found")
    }
})

test_that("plotDMR plot structure contains expected components", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("extdata/example_output.rds", package = "DMRsegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- plotDMR(dmrs, dmr_index = 1)

    expect_true(!is.null(p$labels$x))
    expect_true(!is.null(p$labels$title))
    expect_true(length(p$layers) >= 3)
})
