options("CMEnt.verbose" = 0)

test_that(".rowPairCorrelationStats matches centered matrix calculations", {
    mat <- matrix(
        c(
            0.10, 0.20, NA, 0.40,
            0.30, NA, 0.50, 0.60,
            0.40, 0.45, 0.55, 0.65,
            NA, NA, NA, NA,
            0.20, 0.30, 0.40, 0.50
        ),
        nrow = 5,
        byrow = TRUE
    )
    start_inds <- c(1L, 2L, 4L)
    end_inds <- c(2L, 3L, 5L)
    x_mat <- mat[start_inds, , drop = FALSE]
    y_mat <- mat[end_inds, , drop = FALSE]

    x_means <- rowMeans(x_mat, na.rm = TRUE)
    y_means <- rowMeans(y_mat, na.rm = TRUE)
    x_centered <- x_mat - x_means
    y_centered <- y_mat - y_means
    valid_pairs <- !is.na(x_mat) & !is.na(y_mat)

    ret <- CMEnt:::.rowPairCorrelationStats(mat, start_inds, end_inds)

    expect_equal(ret$x_means, x_means)
    expect_equal(ret$y_means, y_means)
    expect_equal(ret$sum_xy, rowSums(x_centered * y_centered, na.rm = TRUE))
    expect_equal(ret$sum_x2, rowSums(x_centered^2, na.rm = TRUE))
    expect_equal(ret$sum_y2, rowSums(y_centered^2, na.rm = TRUE))
    expect_equal(ret$n_valid, rowSums(valid_pairs))
})

test_that(".testConnectivityBatch can use supplied n_sites without forcing sites_beta dimensions", {
    dim.lazy_beta <- function(x) {
        stop("sites_beta dimensions were forced")
    }
    sites_beta <- structure(list(), class = "lazy_beta")
    pheno <- data.frame(sample = integer(0))

    ret <- CMEnt:::.testConnectivityBatch(
        sites_beta = sites_beta,
        group_inds = list(),
        pheno = pheno,
        max_pval = 0.05,
        entanglement = "strong",
        aggfun = mean,
        testing_mode_per_group = character(0),
        empirical_strategy_per_group = character(0),
        n_sites = 1L
    )

    expect_equal(nrow(ret), 0L)
})

test_that("empirical Monte Carlo connectivity is reproducible under the same RNG seed", {
    sites_beta <- rbind(
        c(0.30, 0.35, 0.40, 0.45, 0.50),
        c(0.32, 0.36, 0.39, 0.44, 0.48)
    )
    pheno <- data.frame(sample = seq_len(5))
    pheno[["__casecontrol__"]] <- c(0L, 0L, 1L, 1L, 1L)
    group_inds <- list(g1 = seq_len(5))

    set.seed(42)
    ret1 <- CMEnt:::.testConnectivityBatch(
        sites_beta = sites_beta,
        group_inds = group_inds,
        pheno = pheno,
        max_pval = 0.05,
        entanglement = "strong",
        aggfun = mean,
        testing_mode_per_group = c(g1 = "empirical"),
        empirical_strategy_per_group = c(g1 = "montecarlo"),
        ntries = 25,
        mid_p = TRUE
    )
    set.seed(42)
    ret2 <- CMEnt:::.testConnectivityBatch(
        sites_beta = sites_beta,
        group_inds = group_inds,
        pheno = pheno,
        max_pval = 0.05,
        entanglement = "strong",
        aggfun = mean,
        testing_mode_per_group = c(g1 = "empirical"),
        empirical_strategy_per_group = c(g1 = "montecarlo"),
        ntries = 25,
        mid_p = TRUE
    )

    expect_equal(ret1, ret2)
})

test_that("ntries = 0 enumerates the full permutation set when it is small enough", {
    sites_beta <- rbind(
        c(0.25, 0.35, 0.45, 0.55),
        c(0.30, 0.40, 0.50, 0.45)
    )
    pheno <- data.frame(sample = seq_len(4))
    pheno[["__casecontrol__"]] <- c(0L, 0L, 1L, 1L)
    group_inds <- list(g1 = seq_len(4))

    ret_default <- CMEnt:::.testConnectivityBatch(
        sites_beta = sites_beta,
        group_inds = group_inds,
        pheno = pheno,
        max_pval = 0.5,
        entanglement = "strong",
        aggfun = mean,
        testing_mode_per_group = c(g1 = "empirical"),
        empirical_strategy_per_group = c(g1 = "permutations"),
        ntries = 0,
        mid_p = FALSE
    )
    ret_full <- CMEnt:::.testConnectivityBatch(
        sites_beta = sites_beta,
        group_inds = group_inds,
        pheno = pheno,
        max_pval = 0.5,
        entanglement = "strong",
        aggfun = mean,
        testing_mode_per_group = c(g1 = "empirical"),
        empirical_strategy_per_group = c(g1 = "permutations"),
        ntries = factorial(4),
        mid_p = FALSE
    )

    expect_equal(ret_default, ret_full)
})

test_that(".testConnectivityBatch treats zero-variance correlations as connected with p-value zero", {
    sites_beta <- rbind(
        c(0.20, 0.30, 0.40, 0.50, 0.60, 0.70),
        rep(0.10, 6)
    )
    pheno <- data.frame(sample = seq_len(6))
    pheno[["__casecontrol__"]] <- c(0L, 0L, 0L, 1L, 1L, 1L)
    group_inds <- list(g1 = seq_len(6))

    for (testing_mode in c("parametric", "empirical")) {
        ret <- CMEnt:::.testConnectivityBatch(
            sites_beta = sites_beta,
            group_inds = group_inds,
            pheno = pheno,
            max_pval = 0.05,
            entanglement = "strong",
            aggfun = mean,
            testing_mode_per_group = c(g1 = testing_mode),
            empirical_strategy_per_group = c(g1 = "montecarlo"),
            ntries = 10,
            mid_p = TRUE
        )

        expect_true(ret$connected[[1]])
        expect_identical(ret$reason[[1]], "")
        expect_equal(ret$pval[[1]], 0)
    }
})

test_that(".testConnectivityBatch marks edges as failing when empirical permutation p-values cannot reach threshold", {
    set.seed(1)
    sites_beta <- matrix(runif(5 * 12, min = 0.05, max = 0.95), nrow = 5, ncol = 12)
    pheno <- data.frame(dummy = seq_len(12))
    pheno[["__casecontrol__"]] <- c(rep(0L, 6), rep(1L, 6))
    group_inds <- list(g1 = 1:6, g2 = 7:12)

    expect_message(
        ret_strong <- CMEnt:::.testConnectivityBatch(
            sites_beta = sites_beta,
            group_inds = group_inds,
            pheno = pheno,
            max_pval = 1e-5,
            entanglement = "strong",
            aggfun = median,
            testing_mode_per_group = c(g1 = "empirical", g2 = "empirical"),
            empirical_strategy_per_group = c(g1 = "permutations", g2 = "permutations"),
            ntries = 50,
            mid_p = FALSE
        ),
        "sufficient small empirical p-value"
    )

    expect_true(all(!ret_strong$connected))
    expect_true(all(ret_strong$reason == "pval>max_pval"))
    expect_true(all(ret_strong$pval == 1))

    ret_weak <- suppressWarnings(CMEnt:::.testConnectivityBatch(
        sites_beta = sites_beta,
        group_inds = group_inds,
        pheno = pheno,
        max_pval = 1e-5,
        entanglement = "weak",
        aggfun = median,
        testing_mode_per_group = c(g1 = "empirical", g2 = "empirical"),
        empirical_strategy_per_group = c(g1 = "permutations", g2 = "permutations"),
        ntries = 50,
        mid_p = FALSE
    ))

    expect_true(all(!ret_weak$connected))
    expect_true(all(grepl("pval>max_pval", ret_weak$reason, fixed = TRUE)))
})
