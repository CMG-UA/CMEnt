library(testthat)

test_that("bridgeConnectivityGaps bridges only short p-value-driven gaps", {
    corr_ret <- data.frame(
        connected = c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE),
        pval = c(0.01, 0.20, 0.01, 0.20, 0.21, 0.01),
        reason = c("", "pval>max_pval", "", "pval>max_pval", "pval>max_pval", "")
    )

    bridged <- DMRsegal:::.bridgeConnectivityGaps(
        corr_ret = corr_ret,
        max_bridge_seeds_gaps = 1L
    )
    expect_equal(bridged$connected, c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE))
    expect_equal(attr(bridged, "bridged_edges"), 1L)
})

test_that("bridgeConnectivityGaps does not bridge structural failures", {
    corr_ret <- data.frame(
        connected = c(TRUE, FALSE, TRUE),
        pval = c(0.01, NA_real_, 0.01),
        reason = c("", "exceeded max distance", "")
    )

    bridged <- DMRsegal:::.bridgeConnectivityGaps(
        corr_ret = corr_ret,
        max_bridge_seeds_gaps = 1L
    )
    expect_equal(bridged$connected, corr_ret$connected)
    expect_equal(attr(bridged, "bridged_edges"), 0L)
})
