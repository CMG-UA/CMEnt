options("CMEnt.verbose" = 0)

makeWeakEntanglementInput <- function() {
    sites_beta <- rbind(
        c(0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.31, 0.36, 0.41, 0.46, 0.51, 0.56),
        c(0.30, 0.35, 0.40, 0.50, 0.55, 0.45, 0.31, 0.36, 0.41, 0.51, 0.56, 0.46)
    )
    pheno <- data.frame(sample = seq_len(ncol(sites_beta)))
    pheno[["__casecontrol__"]] <- c(rep(0L, 6), rep(1L, 6))

    list(
        sites_beta = sites_beta,
        pheno = pheno,
        group_inds = list(g1 = 1:6, g2 = 7:12),
        testing_mode_per_group = c(g1 = "parametric", g2 = "parametric"),
        empirical_strategy_per_group = c(g1 = "auto", g2 = "auto")
    )
}

test_that("weak entanglement keeps edges that fail only after strong-mode correction", {
    input <- makeWeakEntanglementInput()

    strong <- CMEnt:::.testConnectivityBatch(
        sites_beta = input$sites_beta,
        group_inds = input$group_inds,
        pheno = input$pheno,
        testing_mode_per_group = input$testing_mode_per_group,
        empirical_strategy_per_group = input$empirical_strategy_per_group,
        max_pval = 0.05,
        entanglement = "strong",
        aggfun = mean,
        ntries = 10,
        mid_p = TRUE
    )
    weak <- CMEnt:::.testConnectivityBatch(
        sites_beta = input$sites_beta,
        group_inds = input$group_inds,
        pheno = input$pheno,
        testing_mode_per_group = input$testing_mode_per_group,
        empirical_strategy_per_group = input$empirical_strategy_per_group,
        max_pval = 0.05,
        entanglement = "weak",
        aggfun = mean,
        ntries = 10,
        mid_p = TRUE
    )

    expect_false(strong$connected[[1]])
    expect_identical(strong$reason[[1]], "pval>max_pval")
    expect_true(weak$connected[[1]])
    expect_identical(weak$reason[[1]], "")
})

test_that("weak entanglement reports the best passing group p-value", {
    g1x <- seq(0.2, 0.8, length.out = 8)
    g1y <- g1x + c(0.01, -0.01, 0, 0.01, -0.01, 0, 0.01, -0.01)
    g2x <- seq(0.2, 0.8, length.out = 8)
    g2y <- c(0.5, 0.1, 0.8, 0.3, 0.7, 0.2, 0.4, 0.6)
    sites_beta <- rbind(c(g1x, g2x), c(g1y, g2y))
    pheno <- data.frame(sample = seq_len(ncol(sites_beta)))
    pheno[["__casecontrol__"]] <- c(rep(0L, 8), rep(1L, 8))
    group_inds <- list(g1 = 1:8, g2 = 9:16)
    testing_mode_per_group <- c(g1 = "parametric", g2 = "parametric")
    empirical_strategy_per_group <- c(g1 = "auto", g2 = "auto")

    weak <- CMEnt:::.testConnectivityBatch(
        sites_beta = sites_beta,
        group_inds = group_inds,
        pheno = pheno,
        testing_mode_per_group = testing_mode_per_group,
        empirical_strategy_per_group = empirical_strategy_per_group,
        max_pval = 0.05,
        entanglement = "weak",
        aggfun = mean
    )

    expect_true(weak$connected[[1]])
    expect_lt(weak$pval[[1]], 0.05)
})

test_that("weak entanglement reports the best p-value that caused disconnection", {
    g1x <- seq(0.2, 0.8, length.out = 8)
    g1y <- c(0.1646, 0.3944, 0.3580, 0.2304, 0.8687, 0.8215, 0.9774, 0.4475)
    g2x <- seq(0.2, 0.8, length.out = 8)
    g2y <- c(0.5, 0.1, 0.8, 0.3, 0.7, 0.2, 0.4, 0.6)
    sites_beta <- rbind(c(g1x, g2x), c(g1y, g2y))
    pheno <- data.frame(sample = seq_len(ncol(sites_beta)))
    pheno[["__casecontrol__"]] <- c(rep(0L, 8), rep(1L, 8))
    group_inds <- list(g1 = 1:8, g2 = 9:16)
    testing_mode_per_group <- c(g1 = "parametric", g2 = "parametric")
    empirical_strategy_per_group <- c(g1 = "auto", g2 = "auto")

    weak <- CMEnt:::.testConnectivityBatch(
        sites_beta = sites_beta,
        group_inds = group_inds,
        pheno = pheno,
        testing_mode_per_group = testing_mode_per_group,
        empirical_strategy_per_group = empirical_strategy_per_group,
        max_pval = 0.05,
        entanglement = "weak",
        aggfun = mean
    )

    expected <- min(vapply(group_inds, function(idx) {
        group_m <- CMEnt:::.transformBeta(
            sites_beta[, idx, drop = FALSE],
            pheno = pheno[idx, , drop = FALSE],
            covariates = NULL
        )
        r <- stats::cor(group_m[1, ], group_m[2, ])
        df <- length(idx) - 2L
        tstat <- r * sqrt(df / (1 - r * r))
        -2 * expm1(stats::pt(abs(tstat), df = df, log.p = TRUE))
    }, numeric(1)))

    expect_false(weak$connected[[1]])
    expect_identical(weak$reason[[1]], "pval>max_pval;pval>max_pval")
    expect_equal(weak$pval[[1]], expected)
    expect_gt(weak$pval[[1]], 0.05)
})

test_that("buildDMRs accepts weak entanglement on deterministic synthetic data", {
    fixture <- makeSyntheticBuildDMRsFixture()

    strong <- runSyntheticBuildDMRs(fixture, entanglement = "strong")
    weak <- runSyntheticBuildDMRs(fixture, entanglement = "weak")

    expect_s4_class(weak, "GRanges")
    expect_identical(S4Vectors::mcols(weak)$id, S4Vectors::mcols(strong)$id)
})

test_that("entanglement parameter validates correctly", {
    fixture <- makeSyntheticBuildDMRsFixture()

    expect_error(
        runSyntheticBuildDMRs(fixture, entanglement = "invalid"),
        "is not a prefix"
    )
})
