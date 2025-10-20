library(testthat)
library(DMRSegal)

test_that("plotDMR creates a ggplot object", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- plotDMR(dmrs, dmr_index = 1)

    expect_s3_class(p, "ggplot")
    expect_true(!is.null(p$layers))
    expect_true(length(p$layers) > 0)
})

test_that("plotDMR handles invalid dmr_index", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
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

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p1 <- plotDMR(dmrs, dmr_index = 1, array = "450K")
    expect_s3_class(p1, "ggplot")

    p2 <- plotDMR(dmrs, dmr_index = 1, array = "EPIC")
    expect_s3_class(p2, "ggplot")
})

test_that("plotDMR works with different genome versions", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p1 <- plotDMR(dmrs, dmr_index = 1, genome = "hg19")
    expect_s3_class(p1, "ggplot")

    p2 <- plotDMR(dmrs, dmr_index = 1, genome = "hg38")
    expect_s3_class(p2, "ggplot")
})

test_that("plotDMR works with custom title", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    custom_title <- "Test DMR Title"
    p <- plotDMR(dmrs, dmr_index = 1, title = custom_title)

    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$title, custom_title)
})

test_that("plotDMR handles show_pvalues parameter", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p1 <- plotDMR(dmrs, dmr_index = 1, show_pvalues = TRUE)
    expect_s3_class(p1, "ggplot")

    p2 <- plotDMR(dmrs, dmr_index = 1, show_pvalues = FALSE)
    expect_s3_class(p2, "ggplot")
})

test_that("plotDMRs creates a combined plot", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    n_dmrs <- min(4, length(dmrs))
    p <- plotDMRs(dmrs, dmr_indices = 1:n_dmrs, ncol = 2)

    expect_true(!is.null(p))
    expect_true(inherits(p, "gtable") || inherits(p, "list"))
})

test_that("plotDMRs handles NULL dmr_indices", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- plotDMRs(dmrs, dmr_indices = NULL)

    expect_true(!is.null(p))
    expect_true(inherits(p, "gtable") || inherits(p, "list"))
})

test_that("plotDMRs respects ncol parameter", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
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
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
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

test_that("plotDMRWithBeta works with file path", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("reshape2")
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    library(minfi)
    mset <- minfiData::MsetEx[1:1000, ]
    beta_mat <- getBeta(mset)

    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        row.names = colnames(mset)
    )

    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    p <- plotDMRWithBeta(
        dmrs = dmrs,
        dmr_index = 1,
        beta_handler = beta_file,
        pheno = pheno,
        sample_group_col = "Sample_Group"
    )

    expect_true(!is.null(p))
    expect_true(inherits(p, "gtable") || inherits(p, "list"))

    unlink(beta_file)
})

test_that("plotDMRWithBeta works with BetaFileHandler object", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("reshape2")
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    library(minfi)
    mset <- minfiData::MsetEx[1:1000, ]
    beta_mat <- getBeta(mset)

    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        row.names = colnames(mset)
    )

    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    beta_handler <- BetaFileHandler$new(
        beta_file = beta_file,
        array = "450K",
        genome = "hg19",
        verbose = 0
    )

    p <- plotDMRWithBeta(
        dmrs = dmrs,
        dmr_index = 1,
        beta_handler = beta_handler,
        pheno = pheno,
        sample_group_col = "Sample_Group"
    )

    expect_true(!is.null(p))
    expect_true(inherits(p, "gtable") || inherits(p, "list"))

    unlink(beta_file)
})

test_that("plotDMRWithBeta respects max_cpgs parameter", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("reshape2")
    skip_if_not_installed("minfi")
    skip_if_not_installed("minfiData")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    library(minfi)
    mset <- minfiData::MsetEx[1:1000, ]
    beta_mat <- getBeta(mset)

    pheno <- data.frame(
        Sample_Group = pData(mset)$Sample_Group,
        row.names = colnames(mset)
    )

    beta_file <- tempfile(fileext = ".txt")
    write.table(
        cbind(ID = rownames(beta_mat), beta_mat),
        file = beta_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    beta_file <- sortBetaFileByCoordinates(beta_file, overwrite = TRUE)

    p1 <- plotDMRWithBeta(
        dmrs = dmrs,
        dmr_index = 1,
        beta_handler = beta_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        max_cpgs = 50
    )
    expect_true(!is.null(p1))

    p2 <- plotDMRWithBeta(
        dmrs = dmrs,
        dmr_index = 1,
        beta_handler = beta_file,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        max_cpgs = 200
    )
    expect_true(!is.null(p2))

    unlink(beta_file)
})

test_that("plotDMR handles DMRs with no extended CpGs", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    dmr_data <- as.data.frame(S4Vectors::mcols(dmrs))
    no_extended_idx <- which(dmr_data$sup_cpgs_num == 0)

    if (length(no_extended_idx) > 0) {
        p <- plotDMR(dmrs, dmr_index = no_extended_idx[1])
        expect_s3_class(p, "ggplot")
    } else {
        skip("No DMRs without extended CpGs found")
    }
})

test_that("plotDMR handles DMRs with multiple DMPs", {
    skip_if_not_installed("ggplot2")

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
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

    dmrs <- readRDS(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))
    if (length(dmrs) == 0 || !file.exists(system.file("../benchmark_data/dmrs.dmrsegal.rds", package = "DMRSegal", mustWork = FALSE))) {
        skip("Benchmark DMRs not available")
    }

    p <- plotDMR(dmrs, dmr_index = 1)

    expect_true(!is.null(p$labels$x))
    expect_true(!is.null(p$labels$title))
    expect_true(length(p$layers) >= 3)
})
