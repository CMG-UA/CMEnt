options("CMEnt.verbose" = 0)

test_that("findDMRsFromSeeds with expansion_window and max_bridge_seeds_gaps parameters", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    dmps <- subsetDenseExampleDmpsChr5And11(dmps)
    pheno <- loadExampleInputDataChr5And11("pheno")

    # Test with expansion_window and max_bridge_seeds_gaps
    dmrs_expanded <- findDMRsFromSeeds(
            .score_dmrs = FALSE,
            extract_motifs = FALSE,
            annotate_with_genes = FALSE,
            beta = beta,
            seeds = dmps,
            pheno = pheno,
            sample_group_col = "Sample_Group",
            min_seeds = 2,
            min_sites = 3,
            max_lookup_dist = 1000,
            expansion_window = 1, # Expand DMRs by 1bp
            max_bridge_seeds_gaps = 2 # Allow bridging up to 2 seeds apart
        )
    # Assertions
    expect_s4_class(dmrs_expanded, "GRanges")
    if (!is.null(dmrs_expanded)) {
        expect_gt(length(dmrs_expanded), 0L)
        expect_true(all(c("sites_num", "seeds_num", "delta_beta") %in% names(mcols(dmrs_expanded))))
        dmr_df <- as.data.frame(dmrs_expanded)
        # Expansion windows are hard thresholds: each final DMR stays inside its seed-derived window.
        win_df <- CMEnt:::.buildConnectivityWindowsFromDMRs(
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

        # Projection from subset indices must preserve original site IDs.
        all_beta_ids <- rownames(beta)
        used_sites <- unique(unlist(base::strsplit(dmr_df$sites, ",", fixed = TRUE), use.names = FALSE))
        used_sites <- used_sites[nzchar(used_sites)]
        expect_true(all(used_sites %in% all_beta_ids))
    }
})

test_that("findDMRsFromSeeds handles ext_site_delta_beta filtering", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    dmps <- subsetDenseExampleDmpsChr5And11(dmps)

    # Test with no delta beta filtering
    dmrs_no_filter <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_sites = 3,
        ext_site_delta_beta = NA_real_,
        max_lookup_dist = 1000
    )

    # Test with delta beta extension
    dmrs_with_db_ext <- findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_sites = 3,
        ext_site_delta_beta = 0.1, # Extend DMRs to include sites with delta beta >= 0.1
        max_lookup_dist = 1000
    )


    # Assertions
    expect_true(is.null(dmrs_no_filter) || inherits(dmrs_no_filter, "GRanges"))
    expect_true(is.null(dmrs_with_db_ext) || inherits(dmrs_with_db_ext, "GRanges"))

    # Extended results by delta beta comparison should have more or equal DMRs
    if (!is.null(dmrs_no_filter) && !is.null(dmrs_with_db_ext)) {
        expect_true(length(dmrs_with_db_ext) >= length(dmrs_no_filter))
    }
})

test_that("ext_site_delta_beta uses NA as the off switch and 0 as an active threshold", {
    sites_beta <- matrix(
        c(
            0.10, 0.85, 0.40, 0.20, 0.75, 0.35,
            0.80, 0.20, 0.60, 0.70, 0.10, 0.90
        ),
        nrow = 2,
        byrow = TRUE
    )
    pheno <- data.frame(
        Sample_Group = c("grp1", "grp1", "grp1", "grp2", "grp2", "grp2"),
        casecontrol = c(0, 0, 0, 1, 1, 1),
        stringsAsFactors = FALSE
    )
    pheno[, "__casecontrol__"] <- pheno$casecontrol
    group_inds <- split(seq_len(nrow(pheno)), pheno$Sample_Group)
    pval_mode_per_group <- stats::setNames(rep("parametric", length(group_inds)), names(group_inds))
    empirical_strategy_per_group <- stats::setNames(rep("permutations", length(group_inds)), names(group_inds))

    no_force_connect <- CMEnt:::.testConnectivityBatch(
        sites_beta = sites_beta,
        group_inds = group_inds,
        pheno = pheno,
        pval_mode_per_group = pval_mode_per_group,
        empirical_strategy_per_group = empirical_strategy_per_group,
        max_pval = 0.05,
        force_connect_delta_beta = NA_real_,
        max_lookup_dist = 1000,
        site_starts = c(100L, 200L),
        entanglement = "strong",
        aggfun = stats::mean,
        ntries = 10,
        mid_p = TRUE
    )

    zero_threshold_force_connect <- CMEnt:::.testConnectivityBatch(
        sites_beta = sites_beta,
        group_inds = group_inds,
        pheno = pheno,
        pval_mode_per_group = pval_mode_per_group,
        empirical_strategy_per_group = empirical_strategy_per_group,
        max_pval = 0.05,
        force_connect_delta_beta = 0,
        max_lookup_dist = 1000,
        site_starts = c(100L, 200L),
        entanglement = "strong",
        aggfun = stats::mean,
        ntries = 10,
        mid_p = TRUE
    )

    expect_false(no_force_connect$connected[[1]])
    expect_true(zero_threshold_force_connect$connected[[1]])
    expect_identical(
        zero_threshold_force_connect$reason[[1]],
        "abs(delta_beta)>=force_connect_delta_beta"
    )
})

test_that("findDMRsFromSeeds handles adjusted seeds filtering for array data", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    dmps <- subsetDenseExampleDmpsChr5And11(dmps)

    dmrs_adj <- expect_no_error(suppressWarnings(findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        array = "450K",
        genome = "hg19",
        min_seeds = 2,
        min_adj_seeds = 3,
        min_sites = 3,
        ext_site_delta_beta = NA_real_,
        max_lookup_dist = 1000,
        njobs = 1
    )))

    expect_true(is.null(dmrs_adj) || inherits(dmrs_adj, "GRanges"))

    if (!is.null(dmrs_adj) && length(dmrs_adj) > 0L) {
        dmr_df <- as.data.frame(dmrs_adj)
        expect_true(all(c("sites_num_bg", "seeds_num_adj") %in% names(dmr_df)))
        expect_true(all(is.finite(dmr_df$sites_num_bg)))
        expect_true(all(dmr_df$sites_num_bg >= 1))
        expect_true(all(dmr_df$seeds_num_adj >= 3))
    }
})

test_that("findDMRsFromSeeds does not bridge across chromosome boundaries", {
    site_ids <- c("cgA", "cgB", "cgC", "cgD")
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
    rownames(beta) <- site_ids
    colnames(beta) <- paste0("S", seq_len(6))

    locs <- data.frame(
        chr = c("chr1", "chr1", "chr2", "chr2"),
        start = c(100, 200, 100, 200),
        end = c(101, 201, 101, 201),
        row.names = site_ids,
        stringsAsFactors = FALSE
    )

    beta_handler <- getBetaHandler(beta = beta, sorted_locs = locs)
    pheno <- data.frame(
        Sample_Group = c("A", "A", "A", "B", "B", "B"),
        casecontrol = c(0, 0, 0, 1, 1, 1),
        row.names = colnames(beta),
        stringsAsFactors = FALSE
    )
    seeds <- data.frame(pval = rep(1e-6, length(site_ids)), row.names = site_ids)

    dmrs <- expect_no_error(findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta_handler,
        seeds = seeds,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        casecontrol_col = "casecontrol",
        min_seeds = 2,
        min_sites = 2,
        ext_site_delta_beta = NA_real_,
        max_lookup_dist = 1000,
        max_pval = 0.05,
        pval_mode = "parametric",
        max_bridge_seeds_gaps = 1,
        max_bridge_extension_gaps = 0,
        njobs = 1
    ))

    expect_s4_class(dmrs, "GRanges")
    expect_true(length(dmrs) > 0)

    dmr_df <- as.data.frame(dmrs)
    site_chr <- setNames(locs$chr, rownames(locs))
    expect_true(all(site_chr[dmr_df$start_seed] == site_chr[dmr_df$end_seed]))
})

test_that("findDMRsFromSeeds stores all seed IDs including the terminal seed", {
    site_ids <- c("cgA", "cgB", "cgC")
    beta <- matrix(
        c(
            0.10, 0.20, 0.30, 0.15, 0.25, 0.35,
            0.11, 0.21, 0.31, 0.14, 0.24, 0.34,
            0.12, 0.22, 0.32, 0.13, 0.23, 0.33
        ),
        nrow = 3,
        byrow = TRUE
    )
    rownames(beta) <- site_ids
    colnames(beta) <- paste0("S", seq_len(6))

    locs <- data.frame(
        chr = rep("chr1", length(site_ids)),
        start = c(100, 200, 300),
        end = c(101, 201, 301),
        row.names = site_ids,
        stringsAsFactors = FALSE
    )

    beta_handler <- getBetaHandler(beta = beta, sorted_locs = locs)
    pheno <- data.frame(
        Sample_Group = c("A", "A", "A", "B", "B", "B"),
        casecontrol = c(0, 0, 0, 1, 1, 1),
        row.names = colnames(beta),
        stringsAsFactors = FALSE
    )
    seeds <- data.frame(pval = rep(1e-6, length(site_ids)), row.names = site_ids)

    dmrs <- expect_no_error(findDMRsFromSeeds(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta_handler,
        seeds = seeds,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        casecontrol_col = "casecontrol",
        min_seeds = 2,
        min_sites = 2,
        ext_site_delta_beta = NA_real_,
        max_lookup_dist = 1000,
        max_pval = 0.05,
        pval_mode = "parametric",
        max_bridge_seeds_gaps = 1,
        max_bridge_extension_gaps = 0,
        njobs = 1
    ))

    expect_s4_class(dmrs, "GRanges")
    expect_equal(length(dmrs), 1L)

    dmr_df <- as.data.frame(dmrs)
    expect_identical(dmr_df$seeds[[1]], paste(site_ids, collapse = ","))
    expect_identical(dmr_df$start_seed[[1]], site_ids[[1]])
    expect_identical(dmr_df$end_seed[[1]], site_ids[[length(site_ids)]])
    expect_equal(dmr_df$seeds_num[[1]], length(site_ids))
    expect_equal(length(base::strsplit(dmr_df$seeds[[1]], ",", fixed = TRUE)[[1]]), dmr_df$seeds_num[[1]])
})

test_that("findDMRsFromSeeds Stage 2 expansion matches between sequential and chunked parallel execution", {
    beta <- loadExampleInputDataChr5And11("beta")
    dmps <- loadExampleInputDataChr5And11("dmps")
    pheno <- loadExampleInputDataChr5And11("pheno")
    dmps <- subsetDenseExampleDmpsChr5And11(dmps)

    args <- list(
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        beta = beta,
        seeds = dmps,
        pheno = pheno,
        sample_group_col = "Sample_Group",
        min_seeds = 2,
        min_sites = 3,
        max_lookup_dist = 1000,
        expansion_window = 1,
        max_bridge_seeds_gaps = 2
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
