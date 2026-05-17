test_that("combinePvalues returns input p-values when loci have no neighbors", {
    pvals <- c(a = 0.01, b = 0.2, c = 0.8)
    positions <- c(1, 1000, 2000)

    observed <- combinePvalues(
        pvals = pvals,
        positions = positions,
        max_dist = 10,
        lag_step = 10
    )

    expect_equal(as.numeric(observed), as.numeric(pvals))
    expect_equal(names(observed), names(pvals))
})

test_that("combinePvalues applies comb-p z-score SLK combination", {
    pvals <- c(0.5, 0.1, 0.1)
    positions <- c(10, 20, 30)
    acf <- data.frame(
        lag_min = 0,
        lag_max = 100,
        correlation = 0
    )

    observed <- combinePvalues(
        pvals = pvals,
        positions = positions,
        acf = acf
    )
    expected <- stats::pnorm(
        sum(stats::qnorm(pvals, lower.tail = FALSE)) / sqrt(length(pvals)),
        lower.tail = FALSE
    )

    expect_equal(as.numeric(observed), rep(expected, length(pvals)))
})

test_that("combinePvalues preserves original order and chromosome boundaries", {
    pvals <- c(first = 0.01, second = 0.2, third = 0.03)
    positions <- c(100, 10, 120)
    chr <- c("chr1", "chr2", "chr1")
    acf <- data.frame(
        lag_min = 0,
        lag_max = 50,
        correlation = 0
    )

    observed <- combinePvalues(
        pvals = pvals,
        positions = positions,
        chr = chr,
        acf = acf
    )
    expected_chr1 <- stats::pnorm(
        sum(stats::qnorm(c(0.01, 0.03), lower.tail = FALSE)) / sqrt(2),
        lower.tail = FALSE
    )

    expect_equal(as.numeric(observed), c(expected_chr1, 0.2, expected_chr1))
    expect_equal(names(observed), names(pvals))
})

test_that("combinePvalues can adjust SLK p-values with p.adjust methods", {
    pvals <- c(0.01, 0.02, 0.5)
    positions <- c(10, 20, 1000)
    acf <- data.frame(
        lag_min = 0,
        lag_max = 100,
        correlation = 0
    )

    slk <- combinePvalues(pvals, positions, acf = acf)
    adjusted <- combinePvalues(pvals, positions, acf = acf, method = "fdr")

    expect_equal(
        as.numeric(adjusted),
        as.numeric(stats::p.adjust(slk, method = "fdr"))
    )
    expect_equal(attr(adjusted, "slk_pvalues"), slk)
})
