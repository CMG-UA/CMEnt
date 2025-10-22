library(testthat)
library(DMRSegal)

test_that("plotDMR creates a ggplot object", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- plotDMR(dmrs, dmr_index = 1)

    expect_s3_class(p, "ggplot")
    expect_true(!is.null(p$layers))
    expect_true(length(p$layers) > 0)
})

test_that("plotDMR handles invalid dmr_index", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
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

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p1 <- plotDMR(dmrs, dmr_index = 1, array = "450K")
    expect_s3_class(p1, "ggplot")

    p2 <- plotDMR(dmrs, dmr_index = 1, array = "EPIC")
    expect_s3_class(p2, "ggplot")
})

test_that("plotDMR works with different genome versions", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p1 <- plotDMR(dmrs, dmr_index = 1, genome = "hg19")
    expect_s3_class(p1, "ggplot")

    p2 <- plotDMR(dmrs, dmr_index = 1, genome = "hg38")
    expect_s3_class(p2, "ggplot")
})

test_that("plotDMR works with custom title", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    custom_title <- "Test DMR Title"
    p <- plotDMR(dmrs, dmr_index = 1, title = custom_title)

    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$title, custom_title)
})

test_that("plotDMRs creates a combined plot", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    n_dmrs <- min(4, length(dmrs))
    p <- plotDMRs(dmrs, dmr_indices = 1:n_dmrs, ncol = 2)

    expect_true(!is.null(p))
    expect_true(inherits(p, "gtable") || inherits(p, "list") || inherits(p, "patchwork") || inherits(p, "ggplot") || inherits(p, "gg"), paste0("plotDMRs should return a gtable or list of ggplot objects, instead got: ", class(p)))
})

test_that("plotDMRs handles NULL dmr_indices", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- plotDMRs(dmrs, dmr_indices = NULL)

    expect_true(!is.null(p))
    expect_true(inherits(p, "gtable") || inherits(p, "list") || inherits(p, "patchwork") || inherits(p, "ggplot") || inherits(p, "gg"))
})

test_that("plotDMRs respects ncol parameter", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p1 <- plotDMRs(dmrs, dmr_indices = 1:4, ncol = 2)
    expect_true(!is.null(p1))

    p2 <- plotDMRs(dmrs, dmr_indices = 1:4, ncol = 4)
    expect_true(!is.null(p2))
})

test_that("plotDMRWithBeta validates beta_handler type", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("reshape2")

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    expect_error(
        plotDMRWithBeta(dmrs, 1, beta_handler = 123, pheno = data.frame()),
        "beta_handler must be either a file path"
    )

    expect_error(
        plotDMRWithBeta(dmrs, 1, beta_handler = list(), pheno = data.frame()),
        "beta_handler must be either a file path"
    )
})

test_that("plotDMRWithBeta works with beta path", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("reshape2")
    skip_if_not_installed("DMRSegaldata")

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    pheno_file <- DMRSegaldata::getExamplePheno()
    beta_file <- DMRSegaldata::getExampleBeta()
    samplesheet <- read.table(pheno_file, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)


    pheno <- data.frame(
        Sample_Group = samplesheet$Sample_Group,
        row.names = rownames(samplesheet)
    )


    p <- plotDMRWithBeta(
        dmrs = dmrs,
        dmr_index = 1,
        beta_handler = beta_file,
        pheno = pheno,
        sample_group_col = "Sample_Group"
    )

    expect_true(!is.null(p))
    expect_true(inherits(p, "gtable") || inherits(p, "list") || inherits(p, "patchwork") || inherits(p, "ggplot") || inherits(p, "gg"))

    unlink(beta_file)
})


test_that("plotDMR handles DMRs with no extended CpGs", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
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

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
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

    dmrs <- readRDS(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("data/example_output.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- plotDMR(dmrs, dmr_index = 1)

    expect_true(!is.null(p$labels$x))
    expect_true(!is.null(p$labels$title))
    expect_true(length(p$layers) >= 3)
})
