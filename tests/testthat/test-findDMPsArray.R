make_array_dmp_fixture <- function() {
    site_ids <- paste0("cg", seq_len(6L))
    beta <- matrix(
        c(
            0.10, 0.12, 0.11, 0.82, 0.80, 0.81,
            0.74, 0.72, 0.73, 0.20, 0.22, 0.21,
            0.40, 0.42, 0.41, 0.43, 0.44, 0.45,
            0.15, 0.17, 0.16, 0.78, 0.76, 0.77,
            0.65, 0.64, 0.66, 0.35, 0.36, 0.34,
            0.49, 0.50, 0.51, 0.48, 0.47, 0.46
        ),
        nrow = length(site_ids),
        byrow = TRUE,
        dimnames = list(site_ids, paste0("s", seq_len(6L)))
    )
    locs <- data.frame(
        chr = c("chr2", "chr1", "chr1", "chr2", "chr1", "chr2"),
        start = c(300L, 100L, 250L, 200L, 150L, 100L),
        row.names = site_ids,
        stringsAsFactors = FALSE
    )
    sample_ids <- rev(colnames(beta))
    pheno <- data.frame(
        Sample_ID = sample_ids,
        Sample_Group = ifelse(sample_ids %in% paste0("s", 4:6), "case", "control"),
        Age = rev(rep(c(50L, 60L, 70L), 2L)),
        stringsAsFactors = FALSE
    )
    list(beta = beta, locs = locs, pheno = pheno)
}

test_that("findDMPsArray returns BSseq-like DMP columns sorted by location", {
    fixture <- make_array_dmp_fixture()

    dmps <- findDMPsArray(
        beta = fixture$beta,
        samplesheet = fixture$pheno,
        sorted_locs = fixture$locs,
        case_group = "case",
        covariates = "Age",
        chr = "all"
    )

    expect_equal(
        colnames(dmps),
        c("chr", "start", "end", "site_id", "pval", "qval", "delta_beta", "score")
    )
    expect_equal(nrow(dmps), nrow(fixture$beta))
    expect_equal(dmps$end, dmps$start)
    expect_equal(paste(dmps$chr, dmps$start), sort(paste(dmps$chr, dmps$start)))
    expect_equal(
        dmps$delta_beta[match("cg1", dmps$site_id)],
        mean(fixture$beta["cg1", 4:6]) - mean(fixture$beta["cg1", 1:3]),
        tolerance = 1e-8
    )
    expect_true(all(is.finite(dmps$pval)))
    expect_true(all(is.finite(dmps$qval)))
    expect_true(all(dmps$score >= 0 & dmps$score <= 1))
})

test_that("findDMPsArray writes the same output shape, including gzipped files", {
    fixture <- make_array_dmp_fixture()
    output_file <- tempfile(fileext = ".tsv.gz")

    dmps <- findDMPsArray(
        beta = fixture$beta,
        samplesheet = fixture$pheno,
        sorted_locs = fixture$locs,
        case_group = "case",
        chr = "chr1",
        output_file = output_file
    )
    written <- utils::read.table(output_file, header = TRUE, sep = "\t", check.names = FALSE)

    expect_equal(colnames(written), colnames(dmps))
    expect_equal(nrow(written), nrow(dmps))
    expect_true(all(written$chr == "chr1"))
})


test_that("findDMPsArray drops collinear covariate columns without dropping samples", {
    fixture <- make_array_dmp_fixture()
    fixture$pheno$Batch <- ifelse(fixture$pheno$Sample_Group == "case", "case_batch", "control_batch")

    unadjusted <- findDMPsArray(
        beta = fixture$beta,
        samplesheet = fixture$pheno,
        sorted_locs = fixture$locs,
        case_group = "case",
        chr = "all"
    )
    adjusted <- NULL
    expect_warning(
        adjusted <- findDMPsArray(
            beta = fixture$beta,
            samplesheet = fixture$pheno,
            sorted_locs = fixture$locs,
            case_group = "case",
            covariates = "Batch",
            chr = "all"
        ),
        "design columns are collinear"
    )

    expect_equal(adjusted$site_id, unadjusted$site_id)
    expect_equal(adjusted$pval, unadjusted$pval, tolerance = 1e-12)
    expect_equal(adjusted$delta_beta, unadjusted$delta_beta, tolerance = 1e-12)
})
