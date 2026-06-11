test_that("example input loaders assign into the caller environment", {
    skip_if_experimenthub_unavailable()

    direct_env <- new.env(parent = baseenv())
    evalq({
        CMEnt::loadExampleInputData("pheno")
    }, envir = direct_env)

    expect_true(exists("pheno", envir = direct_env, inherits = FALSE))
    expect_s3_class(get("pheno", envir = direct_env, inherits = FALSE), "data.frame")

    subset_env <- new.env(parent = baseenv())
    evalq({
        CMEnt::loadExampleInputDataChr5And11("pheno")
    }, envir = subset_env)

    expect_true(exists("pheno", envir = subset_env, inherits = FALSE))
    expect_s3_class(get("pheno", envir = subset_env, inherits = FALSE), "data.frame")
})

test_that("chromosome subset loader filters beta without full annotation validation", {
    beta <- matrix(
        seq_len(8),
        nrow = 4,
        dimnames = list(c("cg_chr22", "cg_missing", "cg_chr21", "cg_chr1"), c("S1", "S2"))
    )
    dmps <- data.frame(
        pval = c(0.01, 0.02, 0.03),
        row.names = c("cg_chr22", "cg_chr21", "cg_missing")
    )
    locs <- data.frame(
        chr = c("chr1", "chr21", "chr22"),
        start = c(1L, 2L, 3L),
        end = c(2L, 3L, 4L)
    )
    rownames(locs) <- c("cg_chr1", "cg_chr21", "cg_chr22")

    local_mocked_bindings(
        .fetchExampleInputData = function(resources) {
            list(beta = beta, dmps = dmps, pheno = data.frame(), array_type = "450K")[resources]
        },
        getSortedGenomicLocs = function(...) locs,
        getBetaHandler = function(...) stop("full beta validation should not run"),
        .package = "CMEnt"
    )

    target_env <- new.env(parent = baseenv())
    CMEnt:::.loadExampleInputDataSubset(
        "beta",
        "dmps",
        subset = c("chr21", "chr22"),
        envir = target_env
    )

    expect_equal(rownames(target_env$beta), c("cg_chr21", "cg_chr22"))
    expect_equal(rownames(target_env$dmps), c("cg_chr21", "cg_chr22"))
})
