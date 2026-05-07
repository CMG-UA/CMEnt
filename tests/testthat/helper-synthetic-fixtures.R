makeSyntheticFindDMRsFixture <- function() {
    site_ids <- paste0("cg", LETTERS[1:6])
    beta <- matrix(
        c(
            0.10, 0.12, 0.14, 0.78, 0.80, 0.82,
            0.11, 0.13, 0.15, 0.77, 0.79, 0.81,
            0.12, 0.14, 0.16, 0.76, 0.78, 0.80,
            0.50, 0.51, 0.49, 0.52, 0.50, 0.51,
            0.20, 0.22, 0.21, 0.23, 0.24, 0.22,
            0.80, 0.78, 0.79, 0.81, 0.82, 0.80
        ),
        nrow = 6,
        byrow = TRUE,
        dimnames = list(site_ids, paste0("S", seq_len(6)))
    )
    locs <- data.frame(
        chr = c(rep("chr1", 4), rep("chr2", 2)),
        start = c(100L, 200L, 300L, 400L, 100L, 200L),
        end = c(101L, 201L, 301L, 401L, 101L, 201L),
        row.names = site_ids,
        stringsAsFactors = FALSE
    )
    pheno <- data.frame(
        Sample_Group = rep("cohort", 6),
        casecontrol = c(0L, 0L, 0L, 1L, 1L, 1L),
        Age = c(50L, 54L, 58L, 50L, 54L, 58L),
        Gender = c("F", "M", "F", "F", "M", "F"),
        row.names = colnames(beta),
        stringsAsFactors = FALSE
    )
    seeds <- data.frame(
        pval = c(1e-6, 1e-6, 1e-6, 0.5, 1e-6, 1e-6),
        row.names = site_ids,
        stringsAsFactors = FALSE
    )

    list(
        beta = beta,
        beta_handler = getBetaHandler(beta = beta, sorted_locs = locs),
        locs = locs,
        pheno = pheno,
        seeds = seeds
    )
}

runSyntheticFindDMRs <- function(fixture = makeSyntheticFindDMRsFixture(),
                                 beta = fixture$beta_handler,
                                 seeds = fixture$seeds,
                                 ...) {
    args <- list(
        beta = beta,
        seeds = seeds,
        pheno = fixture$pheno,
        sample_group_col = "Sample_Group",
        casecontrol_col = "casecontrol",
        .score_dmrs = FALSE,
        extract_motifs = FALSE,
        annotate_with_genes = FALSE,
        min_seeds = 2,
        min_sites = 2,
        max_lookup_dist = 150,
        max_pval = 0.05,
        testing_mode = "parametric",
        ext_site_delta_beta = NA_real_,
        max_bridge_seeds_gaps = 0L,
        njobs = 1
    )
    args <- utils::modifyList(args, list(...))
    suppressWarnings(do.call(findDMRsFromSeeds, args))
}

makeSyntheticSeedsWithIds <- function(fixture = makeSyntheticFindDMRsFixture(),
                                      include_chr_prefix = TRUE,
                                      id_col = "ID") {
    ids <- paste0(
        if (include_chr_prefix) fixture$locs$chr else sub("^chr", "", fixture$locs$chr),
        ":",
        fixture$locs$start
    )
    seeds <- fixture$seeds
    seeds[[id_col]] <- ids[match(rownames(seeds), rownames(fixture$locs))]
    seeds
}

writeSyntheticBedFile <- function(fixture,
                                  path,
                                  chrom_col = "chrom",
                                  start_col = "start",
                                  end_col = NULL,
                                  id_col = NULL,
                                  include_score = FALSE,
                                  include_strand = FALSE,
                                  drop_chr_prefix = FALSE) {
    beta_df <- as.data.frame(fixture$beta, stringsAsFactors = FALSE)
    bed_df <- data.frame(row.names = seq_len(nrow(fixture$locs)))
    chrom <- fixture$locs$chr
    if (drop_chr_prefix) {
        chrom <- sub("^chr", "", chrom)
    }
    bed_df[[chrom_col]] <- chrom
    bed_df[[start_col]] <- fixture$locs$start
    if (!is.null(end_col)) {
        bed_df[[end_col]] <- fixture$locs$end
    }
    if (!is.null(id_col)) {
        bed_df[[id_col]] <- rownames(fixture$locs)
    }
    if (include_score) {
        bed_df[["score"]] <- 0
    }
    if (include_strand) {
        bed_df[["strand"]] <- "*"
    }
    for (sample_id in colnames(beta_df)) {
        bed_df[[sample_id]] <- beta_df[[sample_id]]
    }
    write.table(
        bed_df,
        file = path,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE,
        col.names = TRUE
    )
    invisible(bed_df)
}

makeSyntheticPlotFixture <- function() {
    fixture <- makeSyntheticFindDMRsFixture()
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(100L, 100L), end = c(300L, 200L)),
        seqinfo = GenomeInfoDb::Seqinfo(
            seqnames = c("chr1", "chr2"),
            genome = "hg38"
        )
    )

    S4Vectors::mcols(dmrs)$id <- c("chr1:cgA-cgC", "chr2:cgE-cgF")
    S4Vectors::mcols(dmrs)$start_seed <- c("cgA", "cgE")
    S4Vectors::mcols(dmrs)$end_seed <- c("cgC", "cgF")
    S4Vectors::mcols(dmrs)$start_seed_pos <- c(100L, 100L)
    S4Vectors::mcols(dmrs)$end_seed_pos <- c(300L, 200L)
    S4Vectors::mcols(dmrs)$start_site <- c("cgA", "cgE")
    S4Vectors::mcols(dmrs)$end_site <- c("cgC", "cgF")
    S4Vectors::mcols(dmrs)$seeds <- c("cgA,cgB,cgC", "cgE,cgF")
    S4Vectors::mcols(dmrs)$sites <- c("cgA,cgB,cgC", "cgE,cgF")
    S4Vectors::mcols(dmrs)$upstream_sites <- c("", "")
    S4Vectors::mcols(dmrs)$downstream_sites <- c("", "")
    S4Vectors::mcols(dmrs)$upstream_expansion_length <- c(0L, 0L)
    S4Vectors::mcols(dmrs)$downstream_expansion_length <- c(0L, 0L)
    S4Vectors::mcols(dmrs)$score <- c(0.90, 0.35)
    S4Vectors::mcols(dmrs)$cv_accuracy <- c(0.85, 0.62)
    S4Vectors::mcols(dmrs)$seeds_num <- c(3L, 2L)
    S4Vectors::mcols(dmrs)$sites_num <- c(3L, 2L)
    S4Vectors::mcols(dmrs)$delta_beta <- c(0.66, 0.02)
    S4Vectors::mcols(dmrs)$pval <- c(0.0023, 0.0131)
    S4Vectors::mcols(dmrs)$in_promoter_of <- NA_character_
    S4Vectors::mcols(dmrs)$in_gene_body_of <- NA_character_

    c(fixture, list(dmrs = dmrs))
}

makeSyntheticDMRBetaStats <- function() {
    data.frame(
        dmr_id = c(1L, 1L, 1L, 2L, 2L, 2L),
        cases_beta = c(0.90, 0.82, 0.18, -0.60, -0.55, -0.10),
        controls_beta = c(0.12, 0.10, 0.40, -0.15, -0.18, -0.30),
        cases_beta_sd = c(0.03, 0.04, 0.08, 0.05, 0.04, 0.03),
        controls_beta_sd = c(0.02, 0.02, 0.03, 0.02, 0.03, 0.02),
        stringsAsFactors = FALSE
    )
}
