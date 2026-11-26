test_that("plotAutoDMRsCircos validates sample_group_col before preparing DMRs", {
    expect_true("sample_group_col" %in% names(formals(plotAutoDMRsCircos)))

    error_message <- "plotAutoDMRsCircos() requires 'sample_group_col' when phenotype data is provided."
    expect_error(
        plotAutoDMRsCircos(
            dmrs = stop("DMRs must not be evaluated before validation"),
            beta = matrix(0.5),
            pheno = data.frame(group = "control")
        ),
        error_message,
        fixed = TRUE
    )
    for (group_col in list(NULL, "", NA_character_, character(), c("group", "other"))) {
        expect_error(
            plotAutoDMRsCircos(
                dmrs = stop("DMRs must not be evaluated before validation"),
                beta = matrix(0.5),
                pheno = data.frame(group = "control"),
                sample_group_col = group_col,
                method = "components"
            ),
            error_message,
            fixed = TRUE
        )
    }
})

test_that("plotAutoDMRsCircos forwards sample_group_col and allows omitted beta tracks", {
    plot_args <- NULL
    local_mocked_bindings(
        plotDMRsCircos = function(...) plot_args <<- list(...),
        .package = "CMEnt"
    )
    dmrs <- GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(start = 100L, width = 50L),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
    )

    plotAutoDMRsCircos(
        dmrs, beta = matrix(0.5), pheno = data.frame(group = "control"),
        sample_group_col = "group", method = "quick"
    )
    expect_identical(plot_args$sample_group_col, "group")

    for (inputs in list(
        list(beta = NULL, pheno = NULL),
        list(beta = NULL, pheno = data.frame(group = "control")),
        list(beta = matrix(0.5), pheno = NULL)
    )) {
        expect_no_error(do.call(plotAutoDMRsCircos, c(list(dmrs = dmrs, method = "quick"), inputs)))
        expect_null(plot_args$sample_group_col)
    }
})
