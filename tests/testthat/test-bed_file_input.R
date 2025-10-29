library(testthat)

create_dmps_with_chr_pos <- function(dmps, beta_mat, locs) {
    dmp_row_names <- rownames(dmps)
    dmp_indices <- match(dmp_row_names, rownames(beta_mat))
    valid_indices <- !is.na(dmp_indices)

    dmps_subset <- dmps[valid_indices, , drop = FALSE]
    dmps_subset$ID <- paste0(as.character(locs$chr[dmp_indices[valid_indices]]), ":", locs$start[dmp_indices[valid_indices]])
    dmps_subset <- dmps_subset[!grepl("NA", dmps_subset$ID), , drop = FALSE]

    dmps_subset
}

create_dmps_without_chr_prefix <- function(dmps, beta_mat, locs) {
    dmp_row_names <- rownames(dmps)
    dmp_indices <- match(dmp_row_names, rownames(beta_mat))
    valid_indices <- !is.na(dmp_indices)

    dmps_subset <- dmps[valid_indices, , drop = FALSE]
    chr_without_prefix <- gsub("^chr", "", as.character(locs$chr[dmp_indices[valid_indices]]))
    dmps_subset$ID <- paste0(chr_without_prefix, ":", locs$start[dmp_indices[valid_indices]])
    dmps_subset <- dmps_subset[!grepl("NA", dmps_subset$ID), , drop = FALSE]

    dmps_subset
}


test_that("findDMRsFromSeeds works with minimal bed file", {
    
    skip_on_ci()
    load(system.file("data/beta.rda", package = "DMRsegal"))
    load(system.file("data/dmps.rda", package = "DMRsegal"))
    load(system.file("data/pheno.rda", package = "DMRsegal"))
    load(system.file("data/array_type.rda", package = "DMRsegal"))
    beta_handler <- getBetaHandler(beta, array = array_type, genome = "hg19")
    beta_mat <- as.matrix(beta_handler$getBeta())
    locs <- beta_handler$getBetaLocs()

    bed_file <- tempfile(fileext = ".bed")
    sample_cols <- rownames(pheno)

    bed_data <- data.frame(
        chrom = as.character(locs$chr),
        start = locs$start,
        stringsAsFactors = FALSE
    )
    for (sample in sample_cols) {
        bed_data[[sample]] <- beta_mat[, sample]
    }

    write.table(bed_data, file = bed_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

    dmps_with_chr_pos <- create_dmps_with_chr_pos(dmps, beta_mat, locs)
    dmrs <- findDMRsFromSeeds(
        beta = bed_file,
        dmps = dmps_with_chr_pos,
        dmps_tsv_id_col = "ID",
        pheno = pheno,
        sample_group_col = "Sample_Group",
        bed_provided = TRUE,
        bed_chrom_col = "chrom",
        bed_start_col = "start",
        min_dmps = 2,
        min_cpgs = 3,
        max_lookup_dist = 1000,
        njobs = 1,
        memory_threshold_mb = 500,
        verbose = 3
    )

    unlink(bed_file)

    expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
    if (!is.null(dmrs) && length(dmrs) > 0) {
        expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs))))
    }
})

# test_that("findDMRsFromSeeds works with full bed file including all optional columns", {
#     skip_on_ci()

#     load(system.file("data/beta.rda", package = "DMRsegal"))
#     load(system.file("data/dmps.rda", package = "DMRsegal"))
#     load(system.file("data/pheno.rda", package = "DMRsegal"))
#     load(system.file("data/array_type.rda", package = "DMRsegal"))

#     beta_handler <- getBetaHandler(beta, array = array_type, genome = "hg19")
#     beta_mat <- as.matrix(beta_handler$getBeta())
#     locs <- beta_handler$getBetaLocs()

#     bed_file <- tempfile(fileext = ".bed")
#     sample_cols <- rownames(pheno)

#     bed_data <- data.frame(
#         chrom = as.character(locs$chr),
#         start = locs$start,
#         end = locs$start + 1,
#         id = rownames(beta_mat),
#         score = rep(0, nrow(beta_mat)),
#         strand = rep("*", nrow(beta_mat)),
#         stringsAsFactors = FALSE
#     )
#     for (sample in sample_cols) {
#         bed_data[[sample]] <- beta_mat[, sample]
#     }

#     write.table(bed_data, file = bed_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#     dmps_with_chr_pos <- create_dmps_with_chr_pos(dmps, beta_mat, locs)

#     dmrs <- findDMRsFromSeeds(
#         beta = bed_file,
#         dmps = dmps_with_chr_pos,
#         dmps_tsv_id_col = "ID",
#         pheno = pheno,
#         sample_group_col = "Sample_Group",
#         bed_provided = TRUE,
#         bed_chrom_col = "chrom",
#         bed_start_col = "start",
#         min_dmps = 2,
#         min_cpgs = 3,
#         max_lookup_dist = 1000,
#         njobs = 1,
#         memory_threshold_mb = 500,
#         verbose = 2
#     )

#     unlink(bed_file)

#     expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
#     if (!is.null(dmrs) && length(dmrs) > 0) {
#         expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs))))
#     }
# })

# test_that("findDMRsFromSeeds detects bed file by extension", {
    
#     skip_on_ci()

#     load(system.file("data/beta.rda", package = "DMRsegal"))
#     load(system.file("data/dmps.rda", package = "DMRsegal"))
#     load(system.file("data/pheno.rda", package = "DMRsegal"))
#     load(system.file("data/array_type.rda", package = "DMRsegal"))

#     beta_handler <- getBetaHandler(beta, array = array_type, genome = "hg19")
#     beta_mat <- as.matrix(beta_handler$getBeta())
#     locs <- beta_handler$getBetaLocs()

#     bed_file <- tempfile(fileext = ".bed")
#     sample_cols <- rownames(pheno)

#     bed_data <- data.frame(
#         chrom = as.character(locs$chr),
#         start = locs$start,
#         stringsAsFactors = FALSE
#     )
#     for (sample in sample_cols) {
#         bed_data[[sample]] <- beta_mat[, sample]
#     }

#     write.table(bed_data, file = bed_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#     dmps_with_chr_pos <- create_dmps_with_chr_pos(dmps, beta_mat, locs)

#     dmrs <- findDMRsFromSeeds(
#         beta = bed_file,
#         dmps = dmps_with_chr_pos,
#         dmps_tsv_id_col = "ID",
#         pheno = pheno,
#         sample_group_col = "Sample_Group",
#         bed_chrom_col = "chrom",
#         bed_start_col = "start",
#         min_dmps = 2,
#         min_cpgs = 3,
#         max_lookup_dist = 1000,
#         njobs = 1,
#         memory_threshold_mb = 500,
#         verbose = 2
#     )

#     unlink(bed_file)

#     expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
#     if (!is.null(dmrs) && length(dmrs) > 0) {
#         expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs))))
#     }
# })

# test_that("findDMRsFromSeeds throws error when DMP IDs are not in chr:pos format with bed file", {
    
#     skip_on_ci()

#     load(system.file("data/beta.rda", package = "DMRsegal"))
#     load(system.file("data/dmps.rda", package = "DMRsegal"))
#     load(system.file("data/pheno.rda", package = "DMRsegal"))
#     load(system.file("data/array_type.rda", package = "DMRsegal"))

#     beta_handler <- getBetaHandler(beta, array = array_type, genome = "hg19")
#     beta_mat <- as.matrix(beta_handler$getBeta())
#     locs <- beta_handler$getBetaLocs()

#     bed_file <- tempfile(fileext = ".bed")
#     sample_cols <- rownames(pheno)

#     bed_data <- data.frame(
#         chrom = as.character(locs$chr),
#         start = locs$start,
#         stringsAsFactors = FALSE
#     )
#     for (sample in sample_cols) {
#         bed_data[[sample]] <- beta_mat[, sample]
#     }

#     write.table(bed_data, file = bed_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#     expect_error(
#         findDMRsFromSeeds(
#             beta = bed_file,
#             dmps = dmps,
#             pheno = pheno,
#             sample_group_col = "Sample_Group",
#             bed_provided = TRUE,
#             bed_chrom_col = "chrom",
#             bed_start_col = "start",
#             min_dmps = 2,
#             min_cpgs = 3,
#             njobs = 1,
#             verbose = 2
#         ),
#         "must be in 'chr:pos' format"
#     )

#     unlink(bed_file)
# })

# test_that("findDMRsFromSeeds works with bed file without chr prefix in chromosome names", {
    
#     skip_on_ci()

#     load(system.file("data/beta.rda", package = "DMRsegal"))
#     load(system.file("data/dmps.rda", package = "DMRsegal"))
#     load(system.file("data/pheno.rda", package = "DMRsegal"))
#     load(system.file("data/array_type.rda", package = "DMRsegal"))

#     beta_handler <- getBetaHandler(beta, array = array_type, genome = "hg19")
#     beta_mat <- as.matrix(beta_handler$getBeta())
#     locs <- beta_handler$getBetaLocs()

#     bed_file <- tempfile(fileext = ".bed")
#     sample_cols <- rownames(pheno)

#     bed_data <- data.frame(
#         chrom = gsub("^chr", "", as.character(locs$chr)),
#         start = locs$start,
#         stringsAsFactors = FALSE
#     )
#     for (sample in sample_cols) {
#         bed_data[[sample]] <- beta_mat[, sample]
#     }

#     write.table(bed_data, file = bed_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#     dmps_with_chr_pos <- create_dmps_without_chr_prefix(dmps, beta_mat, locs)

#     dmrs <- findDMRsFromSeeds(
#         beta = bed_file,
#         dmps = dmps_with_chr_pos,
#         dmps_tsv_id_col = "ID",
#         pheno = pheno,
#         sample_group_col = "Sample_Group",
#         bed_provided = TRUE,
#         bed_chrom_col = "chrom",
#         bed_start_col = "start",
#         min_dmps = 2,
#         min_cpgs = 3,
#         max_lookup_dist = 1000,
#         njobs = 1,
#         memory_threshold_mb = 500,
#         verbose = 2
#     )

#     unlink(bed_file)

#     expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
#     if (!is.null(dmrs) && length(dmrs) > 0) {
#         expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs))))
#     }
# })

# test_that("findDMRsFromSeeds works with bed file and custom column names", {
    
#     skip_on_ci()

#     load(system.file("data/beta.rda", package = "DMRsegal"))
#     load(system.file("data/dmps.rda", package = "DMRsegal"))
#     load(system.file("data/pheno.rda", package = "DMRsegal"))
#     load(system.file("data/array_type.rda", package = "DMRsegal"))

#     beta_handler <- getBetaHandler(beta, array = array_type, genome = "hg19")
#     beta_mat <- as.matrix(beta_handler$getBeta())
#     locs <- beta_handler$getBetaLocs()

#     bed_file <- tempfile(fileext = ".bed")
#     sample_cols <- rownames(pheno)

#     bed_data <- data.frame(
#         my_chr = as.character(locs$chr),
#         my_pos = locs$start,
#         my_end = locs$start + 1,
#         stringsAsFactors = FALSE
#     )
#     for (sample in sample_cols) {
#         bed_data[[sample]] <- beta_mat[, sample]
#     }

#     write.table(bed_data, file = bed_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#     dmps_with_chr_pos <- create_dmps_with_chr_pos(dmps, beta_mat, locs)

#     dmrs <- findDMRsFromSeeds(
#         beta = bed_file,
#         dmps = dmps_with_chr_pos,
#         dmps_tsv_id_col = "ID",
#         pheno = pheno,
#         sample_group_col = "Sample_Group",
#         bed_provided = TRUE,
#         bed_chrom_col = "my_chr",
#         bed_start_col = "my_pos",
#         min_dmps = 2,
#         min_cpgs = 3,
#         max_lookup_dist = 1000,
#         njobs = 1,
#         memory_threshold_mb = 500,
#         verbose = 2
#     )

#     unlink(bed_file)

#     expect_true(is.null(dmrs) || inherits(dmrs, "GRanges"))
#     if (!is.null(dmrs) && length(dmrs) > 0) {
#         expect_true(all(c("cpgs_num", "dmps_num", "delta_beta") %in% names(mcols(dmrs))))
#     }
# })
