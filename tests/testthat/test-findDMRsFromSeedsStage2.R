options("CMEnt.verbose" = 0)

test_that("findDMRsFromSeeds Stage 2 single pass connectivity array outputs with ugap", {
    set.seed(1)
    cpg_ids <- paste0("cg", 1:10)
    beta <- matrix(runif(60), nrow = 10)
    rownames(beta) <- cpg_ids
    colnames(beta) <- paste0("S", 1:6)
    locs <- data.frame(chr = rep("chr1", 10), start = seq(100, 1000, 100), end = seq(101, 1001, 100), row.names = cpg_ids, stringsAsFactors = FALSE)
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
            pval_mode_per_group = c(A = "parametric", B = "parametric"),
            empirical_strategy_per_group = c(A = "auto", B = "auto"), max_pval = 0.05,
            max_lookup_dist = 1000, connectivity_array = conn, ugap = 1L, dgap = 0L,
            splits = splits, njobs = 1, chunk_size = 100
        )
    )
})


test_that("end-of-input is not bridged", {
    set.seed(1)
    cpg_ids <- c("cg1", "cg2", "cg3", "cg4")
    beta <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.71, 0.61, 0.51, 0.41, 0.31, 0.21), nrow = 4, byrow = TRUE)
    rownames(beta) <- cpg_ids
    colnames(beta) <- paste0("S", 1:6)
    locs <- data.frame(
        chr = c("chr1", "chr1", "chr2", "chr2"),
        start = c(100, 200, 100, 200),
        end = c(101, 201, 101, 201),
        row.names = cpg_ids,
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
        connected = c(TRUE, FALSE, TRUE, FALSE),
        pval = c(0.01, NA, 0.01, NA),
        reason = c("", "end-of-input", "", "end-of-input"),
        stringsAsFactors = FALSE
    )
    splits <- matrix(c(1, 3), ncol = 2)
    ret <- CMEnt:::.buildConnectivityArraySinglePass(
        beta_handler = bh, beta_locs = locs, pheno = pheno,
        group_inds = gi, pval_mode_per_group = c(A = "parametric", B = "parametric"),
        empirical_strategy_per_group = c(A = "auto", B = "auto"), max_pval = 0.05,
        max_lookup_dist = 1000, connectivity_array = conn, ugap = 0L, dgap = 1L,
        splits = splits, njobs = 1, chunk_size = 100
    )
    expect_equal(ret$connectivity_array$connected, c(TRUE, FALSE, TRUE, FALSE))
})


test_that("upstream gap recheck skips out-of-bounds runs without warnings", {
    set.seed(1)
    cpg_ids <- c("cg1", "cg2", "cg3", "cg4")
    beta <- matrix(c(
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
        0.11, 0.21, 0.31, 0.41, 0.51, 0.61,
        0.7, 0.6, 0.5, 0.4, 0.3, 0.2,
        0.71, 0.61, 0.51, 0.41, 0.31, 0.21
    ), nrow = 4, byrow = TRUE)
    rownames(beta) <- cpg_ids
    colnames(beta) <- paste0("S", 1:6)
    locs <- data.frame(
        chr = c("chr1", "chr1", "chr2", "chr2"),
        start = c(100, 200, 100, 200),
        end = c(101, 201, 101, 201),
        row.names = cpg_ids,
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
        connected = c(TRUE, FALSE, TRUE, FALSE),
        pval = c(0.01, NA, 0.01, NA),
        reason = c("", "end-of-input", "", "end-of-input"),
        stringsAsFactors = FALSE
    )
    splits <- matrix(c(1, 3), ncol = 2)
    ret <- expect_no_warning(
        CMEnt:::.buildConnectivityArraySinglePass(
            beta_handler = bh, beta_locs = locs, pheno = pheno,
            group_inds = gi, pval_mode_per_group = c(A = "parametric", B = "parametric"),
            empirical_strategy_per_group = c(A = "auto", B = "auto"), max_pval = 0.05,
            max_lookup_dist = 1000, connectivity_array = conn, ugap = 2L, dgap = 0L,
            splits = splits, njobs = 1, chunk_size = 100
        )
    )
    expect_equal(ret$connectivity_array$connected, c(TRUE, FALSE, TRUE, FALSE))
})


test_that("downstream gap recheck skips out-of-bounds runs", {
    set.seed(1)
    cpg_ids <- c("cg1", "cg2", "cg3", "cg4")
    beta <- matrix(c(
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
        0.11, 0.21, 0.31, 0.41, 0.51, 0.61,
        0.7, 0.6, 0.5, 0.4, 0.3, 0.2,
        0.71, 0.61, 0.51, 0.41, 0.31, 0.21
    ), nrow = 4, byrow = TRUE)
    rownames(beta) <- cpg_ids
    colnames(beta) <- paste0("S", 1:6)
    locs <- data.frame(
        chr = c("chr1", "chr1", "chr2", "chr2"),
        start = c(100, 200, 100, 200),
        end = c(101, 201, 101, 201),
        row.names = cpg_ids,
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
        connected = c(FALSE, TRUE, FALSE, TRUE),
        pval = c(NA, 0.01, NA, 0.01),
        reason = c("", "", "end-of-input", "end-of-input"),
        stringsAsFactors = FALSE
    )
    splits <- matrix(c(1, 3), ncol = 2)
    ret <- expect_no_warning(
        CMEnt:::.buildConnectivityArraySinglePass(
            beta_handler = bh, beta_locs = locs, pheno = pheno,
            group_inds = gi, pval_mode_per_group = c(A = "parametric", B = "parametric"),
            empirical_strategy_per_group = c(A = "auto", B = "auto"), max_pval = 0.05,
            max_lookup_dist = 1000, connectivity_array = conn, ugap = 0L, dgap = 2L,
            splits = splits, njobs = 1, chunk_size = 100
        )
    )
    expect_equal(ret$connectivity_array$connected, c(FALSE, TRUE, FALSE, TRUE))
})


test_that("chunk_size is capped by memory budget option", {
    set.seed(1)
    cpg_ids <- paste0("cg", seq_len(2000))
    beta <- matrix(runif(2000 * 6), nrow = 2000)
    rownames(beta) <- cpg_ids
    colnames(beta) <- paste0("S", 1:6)
    locs <- data.frame(
        chr = rep("chr1", 2000),
        start = seq(100, by = 10, length.out = 2000),
        end = seq(101, by = 10, length.out = 2000),
        row.names = cpg_ids,
        stringsAsFactors = FALSE
    )

    bh <- getBetaHandler(beta = beta, sorted_locs = locs)
    pheno <- data.frame(
        g = c("A", "A", "A", "B", "B", "B"),
        row.names = colnames(beta),
        stringsAsFactors = FALSE
    )
    gi <- list(A = 1:3, B = 4:6)

    withr::local_options(list(
        CMEnt.max_chunk_memory_mb = 1,
        CMEnt.chunk_memory_multiplier = 12,
        CMEnt.verbose = 0
    ))

    ret <- CMEnt:::.buildConnectivityArraySinglePass(
        beta_handler = bh,
        beta_locs = locs,
        pheno = pheno,
        group_inds = gi,
        pval_mode_per_group = c(A = "parametric", B = "parametric"),
        empirical_strategy_per_group = c(A = "auto", B = "auto"),
        max_pval = 0.05,
        max_lookup_dist = 1000,
        chunk_size = 5000,
        njobs = 1
    )

    # With a 1 MB memory budget and 6 columns, chunk_size must be capped below 5000.
    expect_true(nrow(ret$splits) > 1)
})
