
library(testthat)

test_that("findDMRsFromSeeds with expansion_window and max_bridge_seeds_gaps parameters", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    # Test with expansion_window and max_bridge_seeds_gaps
    expect_message(
        dmrs_expanded <- findDMRsFromSeeds(
            .score_dmrs = FALSE,
            beta = beta,
            seeds = dmps,
            pheno = pheno,
            sample_group_col = "Sample_Group",
            min_seeds = 2,
            min_cpgs = 3,
            max_lookup_dist = 1000,
            expansion_window = 1, # Expand DMRs by 1bp
            max_bridge_seeds_gaps = 2, # Allow bridging up to 2 seeds apart
            annotate_with_genes = FALSE,
            verbose = 2
        ), "Stage 2 connectivity restricted to 264 seed-derived windows" # Expect the windows to be the same number as the DMRs at that point
    )

    # Assertions
    expect_s4_class(dmrs_expanded, "GRanges")
    if (!is.null(dmrs_expanded)) {
        expect_equal(length(dmrs_expanded), 114L)
        expect_true(all(c("cpgs_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_expanded))))
        dmr_df <- as.data.frame(dmrs_expanded)
        # Expansion windows are hard thresholds: each final DMR stays inside its seed-derived window.
        win_df <- DMRsegal:::.buildConnectivityWindowsFromDMRs(
            dmrs = data.frame(
                chr = as.character(dmr_df$seqnames),
                start_seed_pos = dmr_df$start_seed_pos,
                end_seed_pos = dmr_df$end_seed_pos,
                stringsAsFactors = FALSE
            ),
            expansion_window = 1
        )
        inside_window <- vapply(seq_len(nrow(dmr_df)), function(i) {
            chr_wins <- win_df[win_df$chr == as.character(dmr_df$seqnames[i]), , drop = FALSE]
            any(dmr_df$start[i] >= chr_wins$start & dmr_df$end[i] <= chr_wins$end)
        }, logical(1))
        expect_true(all(inside_window))

        # Projection from subset indices must preserve original CpG IDs.
        all_beta_ids <- rownames(beta)
        used_cpgs <- unique(unlist(strsplit(dmr_df$cpgs, ",", fixed = TRUE), use.names = FALSE))
        used_cpgs <- used_cpgs[nzchar(used_cpgs)]
        expect_true(all(used_cpgs %in% all_beta_ids))
    }
})

test_that("findDMRsFromSeeds handles min_cpg_delta_beta filtering", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    # Test with no delta beta filtering
    dmrs_no_filter <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        min_cpg_delta_beta = 0,
        max_lookup_dist = 1000,
        annotate_with_genes = FALSE
    )

    # Test with delta beta filtering
    dmrs_with_filter <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        min_cpg_delta_beta = 0.1, # Filter out small changes
        max_lookup_dist = 1000,
        annotate_with_genes = FALSE
    )


    # Assertions
    expect_true(is.null(dmrs_no_filter) || inherits(dmrs_no_filter, "GRanges"))
    expect_true(is.null(dmrs_with_filter) || inherits(dmrs_with_filter, "GRanges"))

    # Filtered results should have fewer or equal DMRs
    if (!is.null(dmrs_no_filter) && !is.null(dmrs_with_filter)) {
        expect_true(length(dmrs_with_filter) <= length(dmrs_no_filter))
    }
})

test_that("findDMRsFromSeeds does not bridge across chromosome boundaries", {
    cpg_ids <- c("cgA", "cgB", "cgC", "cgD")
    beta <- matrix(
        c(
            0.10, 0.20, 0.30, 0.15, 0.25, 0.35,
            0.11, 0.21, 0.31, 0.14, 0.24, 0.34,
            0.60, 0.50, 0.40, 0.65, 0.55, 0.45,
            0.59, 0.49, 0.39, 0.66, 0.56, 0.46
        ),
        nrow = 4,
        byrow = TRUE
    )
    rownames(beta) <- cpg_ids
    colnames(beta) <- paste0("S", seq_len(6))

    locs <- data.frame(
        chr = c("chr1", "chr1", "chr2", "chr2"),
        start = c(100, 200, 100, 200),
        end = c(101, 201, 101, 201),
        row.names = cpg_ids,
        stringsAsFactors = FALSE
    )

    beta_handler <- getBetaHandler(beta = beta, sorted_locs = locs)
    pheno <- data.frame(
        Sample_Group = c("A", "A", "A", "B", "B", "B"),
        casecontrol = c(0, 0, 0, 1, 1, 1),
        row.names = colnames(beta),
        stringsAsFactors = FALSE
    )
    seeds <- data.frame(pval = rep(1e-6, length(cpg_ids)), row.names = cpg_ids)

    dmrs <- expect_no_error(findDMRsFromSeeds(
        .score_dmrs = FALSE,
        beta = beta_handler,
        seeds = seeds,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        casecontrol_col = "casecontrol",
        min_seeds = 2,
        min_cpgs = 2,
        min_cpg_delta_beta = 0,
        adaptive_min_cpg_delta_beta = FALSE,
        max_lookup_dist = 1000,
        max_pval = 0.05,
        pval_mode = "parametric",
        max_bridge_seeds_gaps = 1,
        max_bridge_extension_gaps = 0,
        annotate_with_genes = FALSE,
        njobs = 1,
        verbose = 0
    ))

    expect_s4_class(dmrs, "GRanges")
    expect_true(length(dmrs) > 0)

    dmr_df <- as.data.frame(dmrs)
    cpg_chr <- setNames(locs$chr, rownames(locs))
    expect_true(all(cpg_chr[dmr_df$start_seed] == cpg_chr[dmr_df$end_seed]))
})

test_that("findDMRsFromSeeds stores all seed IDs including the terminal seed", {
    cpg_ids <- c("cgA", "cgB", "cgC")
    beta <- matrix(
        c(
            0.10, 0.20, 0.30, 0.15, 0.25, 0.35,
            0.11, 0.21, 0.31, 0.14, 0.24, 0.34,
            0.12, 0.22, 0.32, 0.13, 0.23, 0.33
        ),
        nrow = 3,
        byrow = TRUE
    )
    rownames(beta) <- cpg_ids
    colnames(beta) <- paste0("S", seq_len(6))

    locs <- data.frame(
        chr = rep("chr1", length(cpg_ids)),
        start = c(100, 200, 300),
        end = c(101, 201, 301),
        row.names = cpg_ids,
        stringsAsFactors = FALSE
    )

    beta_handler <- getBetaHandler(beta = beta, sorted_locs = locs)
    pheno <- data.frame(
        Sample_Group = c("A", "A", "A", "B", "B", "B"),
        casecontrol = c(0, 0, 0, 1, 1, 1),
        row.names = colnames(beta),
        stringsAsFactors = FALSE
    )
    seeds <- data.frame(pval = rep(1e-6, length(cpg_ids)), row.names = cpg_ids)

    dmrs <- expect_no_error(findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        beta = beta_handler,
        seeds = seeds,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        casecontrol_col = "casecontrol",
        min_seeds = 2,
        min_cpgs = 2,
        min_cpg_delta_beta = 0,
        adaptive_min_cpg_delta_beta = FALSE,
        max_lookup_dist = 1000,
        max_pval = 0.05,
        pval_mode = "parametric",
        max_bridge_seeds_gaps = 1,
        max_bridge_extension_gaps = 0,
        annotate_with_genes = FALSE,
        njobs = 1,
        verbose = 0
    ))

    expect_s4_class(dmrs, "GRanges")
    expect_equal(length(dmrs), 1L)

    dmr_df <- as.data.frame(dmrs)
    expect_identical(dmr_df$seeds[[1]], paste(cpg_ids, collapse = ","))
    expect_identical(dmr_df$start_seed[[1]], cpg_ids[[1]])
    expect_identical(dmr_df$end_seed[[1]], cpg_ids[[length(cpg_ids)]])
    expect_equal(dmr_df$seeds_num[[1]], length(cpg_ids))
    expect_equal(length(strsplit(dmr_df$seeds[[1]], ",", fixed = TRUE)[[1]]), dmr_df$seeds_num[[1]])
})

test_that("findDMRsFromSeeds Stage 2 expansion matches between sequential and chunked parallel execution", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")

    withr::local_options(list(DMRsegal.parallel_dmr_chunk_size = 1L))

    args <- list(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        expansion_window = 1,
        max_bridge_seeds_gaps = 2,
        annotate_with_genes = FALSE,
        verbose = 0
    )

    dmrs_seq <- do.call(findDMRsFromSeeds, c(args, list(njobs = 1)))
    dmrs_par <- do.call(findDMRsFromSeeds, c(args, list(njobs = 2)))

    expect_s4_class(dmrs_seq, "GRanges")
    expect_s4_class(dmrs_par, "GRanges")

    seq_df <- as.data.frame(dmrs_seq)
    par_df <- as.data.frame(dmrs_par)

    ord_cols <- c("seqnames", "start", "end", "start_seed", "end_seed", "seeds")
    seq_df <- seq_df[do.call(order, seq_df[ord_cols]), , drop = FALSE]
    par_df <- par_df[do.call(order, par_df[ord_cols]), , drop = FALSE]
    rownames(seq_df) <- NULL
    rownames(par_df) <- NULL

    expect_identical(seq_df, par_df)
})
