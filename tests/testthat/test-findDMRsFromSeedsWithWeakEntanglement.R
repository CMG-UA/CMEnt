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

test_that("findDMRsFromSeeds accepts weak entanglement on deterministic synthetic data", {
    fixture <- makeSyntheticFindDMRsFixture()

    strong <- runSyntheticFindDMRs(fixture, entanglement = "strong")
    weak <- runSyntheticFindDMRs(fixture, entanglement = "weak")

    expect_s4_class(weak, "GRanges")
    expect_identical(S4Vectors::mcols(weak)$id, S4Vectors::mcols(strong)$id)
})

test_that("entanglement parameter validates correctly", {
    fixture <- makeSyntheticFindDMRsFixture()

    expect_error(
        runSyntheticFindDMRs(fixture, entanglement = "invalid"),
        "is not a prefix"
    )
})
