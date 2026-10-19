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
            sites_beta,
            pheno = pheno,
            covariates = NULL,
            cols = idx
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

test_that("undersized groups cannot supply a passing result in weak mode", {
    x <- seq(0.2, 0.8, length.out = 8)
    y <- c(0.5, 0.1, 0.8, 0.3, 0.7, 0.2, 0.4, 0.6)
    sites_beta <- rbind(c(0.2, 0.8, x), c(0.2, 0.8, y), c(0.2, 0.8, y))
    pheno <- data.frame(sample = seq_len(ncol(sites_beta)))

    ret <- CMEnt:::.testConnectivityBatch(
        sites_beta = sites_beta,
        group_inds = list(small = 1:2, eligible = 3:10),
        pheno = pheno,
        testing_mode_per_group = c(small = "parametric", eligible = "parametric"),
        empirical_strategy_per_group = c(small = "auto", eligible = "auto"),
        max_pval = 0.05,
        entanglement = "weak"
    )

    expect_identical(ret$connected, c(FALSE, TRUE))
    expect_identical(ret$reason, c("df<1;pval>max_pval", ""))
    expect_identical(ret$failing_groups, c("small;eligible", ""))
    expect_gt(ret$pval[[1]], 0.05)
    expect_lt(ret$pval[[2]], 0.05)
})

test_that("weak mode rejects correlations when all groups are undersized", {
    sites_beta <- rbind(c(0.2, 0.4, 0.6), c(0.3, 0.5, 0.7))
    pheno <- data.frame(sample = seq_len(ncol(sites_beta)))

    ret <- CMEnt:::.testConnectivityBatch(
        sites_beta = sites_beta,
        group_inds = list(empty = integer(0), single = 1L, pair = 2:3),
        pheno = pheno,
        testing_mode_per_group = c(empty = "parametric", single = "parametric", pair = "parametric"),
        empirical_strategy_per_group = c(empty = "auto", single = "auto", pair = "auto"),
        max_pval = 0.05,
        entanglement = "weak"
    )

    expect_false(ret$connected[[1]])
    expect_true(is.na(ret$pval[[1]]))
    expect_identical(ret$reason[[1]], "df<1;df<1;df<1")
    expect_identical(ret$failing_groups[[1]], "empty;single;pair")
})

test_that("weak mode cannot use constant sites as evidence of correlation", {
    sites_beta <- rbind(rep(0.2, 6), seq(0.2, 0.7, length.out = 6), rep(0.7, 6), rep(0.8, 6))
    pheno <- data.frame(sample = seq_len(ncol(sites_beta)))

    for (testing_mode in c("parametric", "empirical")) {
        ret <- CMEnt:::.testConnectivityBatch(
            sites_beta = sites_beta,
            group_inds = list(g1 = 1:6),
            pheno = pheno,
            testing_mode_per_group = c(g1 = testing_mode),
            empirical_strategy_per_group = c(g1 = "permutations"),
            max_pval = 0.05,
            entanglement = "weak",
            ntries = 10
        )

        expect_identical(ret$connected, rep(FALSE, 3))
        expect_true(all(is.na(ret$pval)))
        expect_identical(ret$reason, rep("na r", 3))
        expect_identical(ret$failing_groups, rep("g1", 3))
    }
})

test_that("constant groups cannot connect weak-mode pairs or lower another group's p-value", {
    x <- seq(0.2, 0.7, length.out = 6)
    y <- c(0.5, 0.1, 0.8, 0.3, 0.7, 0.2)
    sites_beta <- rbind(c(rep(0.2, 6), x), c(rep(0.3, 6), y), c(rep(0.4, 6), y))
    pheno <- data.frame(sample = seq_len(ncol(sites_beta)))

    for (testing_mode in c("parametric", "empirical")) {
        reference <- CMEnt:::.testConnectivityBatch(
            sites_beta = sites_beta,
            group_inds = list(variable = 7:12),
            pheno = pheno,
            testing_mode_per_group = c(variable = testing_mode),
            empirical_strategy_per_group = c(variable = "permutations"),
            max_pval = 0.05,
            entanglement = "weak",
            ntries = factorial(6)
        )
        ret <- CMEnt:::.testConnectivityBatch(
            sites_beta = sites_beta,
            group_inds = list(constant = 1:6, variable = 7:12),
            pheno = pheno,
            testing_mode_per_group = c(constant = testing_mode, variable = testing_mode),
            empirical_strategy_per_group = c(constant = "permutations", variable = "permutations"),
            max_pval = 0.05,
            entanglement = "weak",
            ntries = factorial(6)
        )

        expect_identical(ret$connected, c(FALSE, TRUE))
        expect_equal(ret$pval, reference$pval)
        expect_gt(ret$pval[[2]], 0)
        expect_identical(ret$reason, c("na r;pval>max_pval", ""))
        expect_identical(ret$failing_groups, c("constant;variable", ""))
    }
})

test_that("buildDMRs accepts weak entanglement on deterministic synthetic data", {
    fixture <- makeSyntheticBuildDMRsFixture()

    expect_message(
        strong <- runSyntheticBuildDMRs(fixture, entanglement = "strong"),
        "No DMRs remain after filtering based on min_seeds.",
        fixed = TRUE
    )
    expect_message(
        weak <- runSyntheticBuildDMRs(fixture, entanglement = "weak"),
        "No DMRs remain after filtering based on min_seeds.",
        fixed = TRUE
    )

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
