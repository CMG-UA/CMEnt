test_that("example input loaders assign into the caller environment", {
    skip_if_not_installed("ExperimentHub")

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
