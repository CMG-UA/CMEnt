options("CMEnt.verbose" = 0)

test_that("findDMRsFromSeeds Stage 2 single pass connectivity array outputs with ugap", {
    set.seed(1)
    site_ids <- paste0("cg", 1:10)
    beta <- matrix(runif(60), nrow = 10)
    rownames(beta) <- site_ids
    colnames(beta) <- paste0("S", 1:6)
    locs <- data.frame(chr = rep("chr1", 10), start = seq(100, 1000, 100), end = seq(101, 1001, 100), row.names = site_ids, stringsAsFactors = FALSE)
    bh <- getBetaHandler(beta = beta, sorted_locs = locs)
    pheno <- data.frame(g = c("A", "A", "A", "B", "B", "B"), row.names = colnames(beta), stringsAsFactors = FALSE)
    gi <- list(A = 1:3, B = 4:6)
    conn <- data.frame(
        connected = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
        pval = NA_real_, reason = rep("", 10),
        stringsAsFactors = FALSE
    )
    conn$reason[10] <- "end-of-input"
    splits <- matrix(c(1, 8), ncol = 2)
    expect_no_error(
        CMEnt:::.buildConnectivityArraySinglePass(
            beta_handler = bh, beta_locs = locs, pheno = pheno, group_inds = gi,
            testing_mode_per_group = c(A = "parametric", B = "parametric"),
            empirical_strategy_per_group = c(A = "auto", B = "auto"), max_pval = 0.05,
            max_lookup_dist = 1000, connectivity_array = conn, ugap = 1L, dgap = 0L,
            splits = splits, njobs = 1
        )
    )
})


test_that("bridge recheck follows runs containing newly bridged edges", {
    site_ids <- paste0("cg", seq_len(6))
    beta <- matrix(
        rep(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), each = length(site_ids)),
        nrow = length(site_ids),
        byrow = FALSE
    )
    rownames(beta) <- site_ids
    colnames(beta) <- paste0("S", seq_len(6))
    locs <- data.frame(
        chr = rep("chr1", length(site_ids)),
        start = seq(100, 600, 100),
        end = seq(101, 601, 100),
        row.names = site_ids,
        stringsAsFactors = FALSE
    )
    bh <- getBetaHandler(beta = beta, sorted_locs = locs)
    pheno <- data.frame(
        g = c("A", "A", "A", "B", "B", "B"),
        row.names = colnames(beta),
        stringsAsFactors = FALSE
    )
    gi <- list(A = 1:3, B = 4:6)
    conn <- data.frame(
        connected = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
        pval = c(0.01, rep(NA_real_, 5)),
        reason = c("", rep("pval>max_pval", 4), "end-of-input"),
        stringsAsFactors = FALSE
    )
    splits <- matrix(c(1, 5), ncol = 2)

    ret1 <- CMEnt:::.buildConnectivityArraySinglePass(
        beta_handler = bh, beta_locs = locs, pheno = pheno,
        group_inds = gi, testing_mode_per_group = c(A = "parametric", B = "parametric"),
        empirical_strategy_per_group = c(A = "auto", B = "auto"), max_pval = 0.05,
        max_lookup_dist = 1000, connectivity_array = conn, ugap = 0L, dgap = 1L,
        splits = splits, njobs = 1
    )
    expect_true(ret1$connectivity_array$connected[[2]])
    expect_equal(ret1$recheck, 2L)

    ret2 <- CMEnt:::.buildConnectivityArraySinglePass(
        beta_handler = bh, beta_locs = locs, pheno = pheno,
        group_inds = gi, testing_mode_per_group = c(A = "parametric", B = "parametric"),
        empirical_strategy_per_group = c(A = "auto", B = "auto"), max_pval = 0.05,
        max_lookup_dist = 1000, connectivity_array = ret1$connectivity_array,
        recheck = ret1$recheck, ugap = 0L, dgap = 1L,
        splits = splits, njobs = 1
    )
    expect_true(ret2$connectivity_array$connected[[3]])
    expect_equal(ret2$recheck, 3L)
})

test_that("bridge recheck keeps the minimum p-value observed for bridged edges", {
    site_ids <- paste0("cg", seq_len(5))
    x <- seq(0.2, 0.8, length.out = 8)
    y <- x + c(0.01, -0.01, 0, 0.01, -0.01, 0, 0.01, -0.01)
    beta <- rbind(
        rep(0.5, length(x)),
        x,
        rep(0.3, length(x)),
        y,
        rep(0.4, length(x))
    )
    rownames(beta) <- site_ids
    colnames(beta) <- paste0("S", seq_along(x))
    locs <- data.frame(
        chr = rep("chr1", length(site_ids)),
        start = seq(100, 500, 100),
        end = seq(101, 501, 100),
        row.names = site_ids,
        stringsAsFactors = FALSE
    )
    bh <- getBetaHandler(beta = beta, sorted_locs = locs)
    pheno <- data.frame(
        g = rep("A", length(x)),
        row.names = colnames(beta),
        stringsAsFactors = FALSE
    )
    pheno[["__casecontrol__"]] <- rep(c(0L, 1L), each = 4L)
    conn <- data.frame(
        connected = c(TRUE, FALSE, FALSE, FALSE, FALSE),
        pval = c(0.01, 1e-12, NA_real_, NA_real_, NA_real_),
        reason = c("", "pval>max_pval", "pval>max_pval", "pval>max_pval", "end-of-input"),
        stringsAsFactors = FALSE
    )

    ret <- CMEnt:::.buildConnectivityArraySinglePass(
        beta_handler = bh, beta_locs = locs, pheno = pheno,
        group_inds = list(A = seq_along(x)),
        testing_mode_per_group = c(A = "parametric"),
        empirical_strategy_per_group = c(A = "auto"),
        max_pval = 0.05,
        max_lookup_dist = 1000,
        connectivity_array = conn,
        ugap = 0L,
        dgap = 2L,
        splits = matrix(c(1, 4), ncol = 2),
        njobs = 1
    )

    expect_true(ret$connectivity_array$connected[[2]])
    expect_true(ret$connectivity_array$connected[[3]])
    expect_equal(ret$connectivity_array$pval[[2]], 1e-12)
    expect_false(is.na(ret$connectivity_array$pval[[3]]))
    expect_lt(ret$connectivity_array$pval[[3]], 0.05)
})

test_that("downstream bridge recheck shifts back when the forward candidate is too distant", {
    site_ids <- paste0("cg", seq_len(5))
    x <- seq(0.2, 0.8, length.out = 8)
    y <- x + c(0.01, -0.01, 0, 0.01, -0.01, 0, 0.01, -0.01)
    beta <- rbind(
        x,
        rep(0.5, length(x)),
        y,
        rep(0.4, length(x)),
        rev(x)
    )
    rownames(beta) <- site_ids
    colnames(beta) <- paste0("S", seq_along(x))
    locs <- data.frame(
        chr = rep("chr1", length(site_ids)),
        start = c(100L, 200L, 300L, 400L, 10000L),
        end = c(101L, 201L, 301L, 401L, 10001L),
        row.names = site_ids,
        stringsAsFactors = FALSE
    )
    bh <- getBetaHandler(beta = beta, sorted_locs = locs)
    pheno <- data.frame(
        g = rep("A", length(x)),
        row.names = colnames(beta),
        stringsAsFactors = FALSE
    )
    pheno[["__casecontrol__"]] <- rep(c(0L, 1L), each = 4L)
    conn <- data.frame(
        connected = c(FALSE, TRUE, FALSE, FALSE, FALSE),
        pval = c(NA_real_, 0.01, NA_real_, NA_real_, NA_real_),
        reason = c("pval>max_pval", "", "pval>max_pval", "pval>max_pval", "end-of-input"),
        stringsAsFactors = FALSE
    )

    ret <- CMEnt:::.buildConnectivityArraySinglePass(
        beta_handler = bh, beta_locs = locs, pheno = pheno,
        group_inds = list(A = seq_along(x)),
        testing_mode_per_group = c(A = "parametric"),
        empirical_strategy_per_group = c(A = "auto"),
        max_pval = 0.05,
        max_lookup_dist = 500,
        connectivity_array = conn,
        ugap = 0L,
        dgap = 2L,
        splits = matrix(c(1, 4), ncol = 2),
        njobs = 1
    )

    expect_true(ret$connectivity_array$connected[[1]])
    expect_true(ret$connectivity_array$connected[[2]])
    expect_false(ret$connectivity_array$connected[[3]])
    expect_false(ret$connectivity_array$connected[[4]])
    expect_equal(ret$recheck, 1L)
    expect_lt(ret$connectivity_array$pval[[1]], 0.05)
    expect_equal(ret$connectivity_array$reason[[1]], "bridged")
    expect_equal(ret$connectivity_array$reason[[2]], "")
    expect_equal(ret$connectivity_array$pval[[2]], 0.01)
})


test_that("connectivity chunk size is derived from available RAM", {
    chunk_size <- CMEnt:::.connectivityChunkSize(
        n_samples = 6,
        njobs = 2,
        n_pairs = 5000,
        available_ram_bytes = 1024^2
    )
    expect_equal(chunk_size, floor(0.9 * 1024^2 / (2 * 6 * 8 * 12)))
    expect_equal(
        CMEnt:::.connectivityChunkSize(6, 2, 10, available_ram_bytes = 1024^2),
        10L
    )
})
