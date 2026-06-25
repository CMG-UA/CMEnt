suppressPackageStartupMessages({
    library(testthat)
    library(CMEnt)
})
options("CMEnt.verbose" = 0)

test_that("beta file helpers normalize headers and detect mismatched rows", {
    beta_file <- tempfile(fileext = ".tsv")
    withr::defer(unlink(beta_file), teardown_env())
    writeLines(c(
        "S1\tS2",
        "cg1\t0.1\t0.2",
        "cg2\t0.3\t0.4"
    ), beta_file)

    info <- CMEnt:::.getBetaFileHeaderInfo(beta_file)
    expect_true(info$missing_row_id_header)
    expect_equal(info$file_beta_col_names, c(".row_id", "S1", "S2"))
    expect_equal(info$sample_col_names, c("S1", "S2"))

    beta_data <- CMEnt:::.readBetaFileData(beta_file)
    expect_equal(colnames(beta_data), c(".row_id", "S1", "S2"))
    expect_equal(beta_data$.row_id, c("cg1", "cg2"))

    bad_file <- tempfile(fileext = ".tsv")
    withr::defer(unlink(bad_file), teardown_env())
    writeLines(c("id\tS1", "cg1\t0.1\tunexpected\textra"), bad_file)
    expect_error(
        CMEnt:::.readBetaFileData(bad_file),
        "column count mismatch"
    )
})

test_that("tabix beta sample detection excludes structural BED columns", {
    names <- c("#chrom", "start", "end", "id", "score", "strand", "S1", "S2")

    expect_equal(
        CMEnt:::.getBetaSampleColInds(names, is_tabix = TRUE),
        c(7L, 8L)
    )
    expect_equal(CMEnt:::.getBetaSampleColInds(character()), integer(0))
    expect_equal(CMEnt:::.getBetaSampleColInds("only_id"), integer(0))
})

test_that("phenotype list columns must contain flat scalar values", {
    expect_equal(
        CMEnt:::.coercePhenoColumn(list("case", "control"), "group"),
        c("case", "control")
    )
    expect_error(
        CMEnt:::.coercePhenoColumn(list("case", c("control", "extra")), "group"),
        "length != 1"
    )
    expect_error(
        CMEnt:::.coercePhenoColumn(list(list("case")), "group"),
        "nested lists"
    )
})

test_that("readSamplesheet filters case and control labels", {
    samplesheet <- tempfile(fileext = ".csv")
    withr::defer(unlink(samplesheet), teardown_env())
    writeLines(c(
        "sample,group,batch",
        "S10,control,A",
        "S2,case,B",
        "S1,control,A",
        "S3,case,B"
    ), samplesheet)

    pheno <- CMEnt:::.readSamplesheet(
        samplesheet_file = samplesheet,
        samplesheet_file_sep = ",",
        sample_group_col = "group",
        sample_group_case = "case",
        sample_group_control = "control",
        target_col = "sample",
        subset = c("S1", "S2", "S3")
    )

    expect_equal(rownames(pheno), c("S1", "S2", "S3"))
    expect_equal(pheno$casecontrol, c(0, 1, 1))
    expect_error(
        CMEnt:::.readSamplesheet(samplesheet, ",", "missing", "case", "control", "sample"),
        "Sample group column"
    )
    expect_error(
        CMEnt:::.readSamplesheet(samplesheet, ",", "group", "other", "control", "sample"),
        "case label other"
    )
    expect_error(
        CMEnt:::.readSamplesheet(samplesheet, ",", "group", NULL, NULL, "sample"),
        "Both sample_group_case and sample_group_control"
    )
})

test_that("processSamplesheet infers csv and tsv separators", {
    samplesheet <- tempfile(fileext = ".tsv")
    withr::defer(unlink(samplesheet), teardown_env())
    writeLines(c(
        "sample\tgroup",
        "S1\tcontrol",
        "S2\tcase"
    ), samplesheet)

    args <- list(
        samplesheet = samplesheet,
        target_col = "sample",
        sample_group_col = "group",
        sample_group_control = "control",
        sample_group_case = "case"
    )
    pheno <- CMEnt:::.processSamplesheet(args)$samplesheet

    expect_equal(rownames(pheno), c("S1", "S2"))
    expect_equal(pheno$casecontrol, c(0, 1))
})

test_that("logging and parallel helpers handle quiet and fallback branches", {
    withr::local_options(list(CMEnt.verbose = 0, CMEnt.biocparallel_backend = "unknown"))

    expect_equal(CMEnt:::.fmt_dur(NULL), "")
    expect_match(CMEnt:::.formatLogOutput("hello world", lead = ">", level = 2), ">")
    expect_warning(CMEnt:::.log_warn("careful"), "careful")
    expect_error(CMEnt:::.log_error("boom"), "boom")
    expect_equal(CMEnt:::.node_size(), if (8L * .Machine$sizeof.pointer == 32L) 28L else 56L)

    expect_s4_class(CMEnt:::.makeBiocParallelParam(1), "SerialParam")
    expect_s4_class(CMEnt:::.makeBiocParallelParam(2, n_tasks = 1), "SerialParam")
    expect_error(CMEnt:::.makeBiocParallelParam(0), "positive integer")
    expect_warning(
        param <- CMEnt:::.makeBiocParallelParam(2, parallel_backend = "unknown"),
        "Unsupported CMEnt BiocParallel backend"
    )
    expect_true(inherits(param, "BiocParallelParam"))

    withr::local_envvar(`_R_CHECK_LIMIT_CORES_` = "TRUE")
    expect_lte(CMEnt:::.defaultNJobs(), 2L)
    expect_lte(BiocParallel::bpworkers(CMEnt:::.makeBiocParallelParam(8, parallel_backend = "auto")), 2L)
})

test_that("registry post-processing selects, renames, derives, and indexes columns", {
    registry <- data.frame(
        chrom = c("chr1", "chr2"),
        pos = c(10L, 20L),
        sample = c("S1", "S2"),
        stringsAsFactors = FALSE
    )

    processed <- CMEnt:::.postProcessRegistry(
        registry,
        select_columns = c("chrom", "pos"),
        rename = c(chrom = "chr", pos = "start"),
        derive = list(
            locus = list(cols = c("chr", "start"), fun = function(chr, start) paste0(chr, ":", start))
        ),
        indices = "locus"
    )

    expect_equal(colnames(processed), c("chr", "start", "locus"))
    expect_equal(rownames(processed), c("chr1:10", "chr2:20"))
    expect_error(
        CMEnt:::.postProcessRegistry(registry, derive = list(missing = list(cols = "absent", fun = identity))),
        "Cannot derive column"
    )
    expect_error(
        CMEnt:::.postProcessRegistry(registry, indices = "absent"),
        "Cannot set indices"
    )
})

test_that("covariate transformation residualizes nonsingular covariates", {
    pheno <- data.frame(
        batch = c("A", "A", "B", "B"),
        age = c(20, 30, 40, 50)
    )
    beta <- matrix(
        c(0.2, 0.3, 0.7, 0.8, 0.4, 0.5, 0.6, 0.7),
        nrow = 2,
        byrow = TRUE
    )

    model <- CMEnt:::.prepareCovariateModel(pheno, covariates = c("batch", "age"))
    expect_false(model$is_singular)
    expect_equal(dim(model$covariate_matrix), c(4, 3))
    expect_true("batchB" %in% colnames(model$covariate_matrix))
    expect_error(
        CMEnt:::.prepareCovariateModel(pheno, covariates = "missing"),
        "Not all covariates"
    )
    expect_error(
        CMEnt:::.prepareCovariateModel(data.frame(batch = c("A", NA)), covariates = "batch"),
        "missing values"
    )

    transformed <- CMEnt:::.transformBeta(beta, pheno, covariate_model = model)
    expect_equal(dim(transformed), dim(beta))
    expect_false(anyNA(transformed))

    constant <- CMEnt:::.prepareCovariateModel(
        data.frame(batch = c("A", "A", "A")),
        covariates = "batch"
    )
    expect_false(constant$is_singular)
    expect_equal(dim(constant$covariate_matrix), c(3, 1))
    expect_equal(constant$dropped_columns, "batch")

    aliased <- NULL
    aliased <- CMEnt:::.prepareCovariateModel(
        data.frame(batch = c("A", "B", "C", "A"), duplicate_batch = c("A", "B", "C", "A")),
        covariates = c("batch", "duplicate_batch")
    )
    expect_equal(nrow(aliased$covariate_matrix), 4)
    expect_false(aliased$is_singular)
    expect_true(length(aliased$dropped_columns) > 0)
    expect_equal(
        CMEnt:::.remove_confounder_effect(beta, covariate_matrix = matrix(, nrow = 4, ncol = 0)),
        beta
    )
})

test_that("small pure utility helpers handle edge cases", {
    expect_equal(CMEnt:::.splitCsvValues(" a, ,b,NA "), c("a", "b", "NA"))
    expect_equal(CMEnt:::.splitCsvValues(NA_character_), character(0))
    expect_equal(CMEnt:::.splitCsvIndices("1, two, 3"), c(1L, 3L))
    expect_equal(CMEnt:::.downsampleFlankIndices(1:10, 3), c(1L, 5L, 9L))
    expect_equal(CMEnt:::.downsampleFlankIndices(1:3, 10), 1:3)
    expect_null(CMEnt:::.getDerivedOutputPath(NULL, ".txt"))
    expect_equal(CMEnt:::.getDerivedOutputPath("prefix", ".txt"), "prefix.txt")
    expect_equal(CMEnt:::.resolveBSGenomePackage("mm39"), "BSgenome.Mmusculus.UCSC.mm39")
    expect_null(CMEnt:::.resolveBSGenomePackage("unknown"))
})

test_that("convertToGRanges fills missing GRanges genome without replacing seqlevels", {
    dmrs <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(10L, 20L), width = 5L)
    )
    old_seqlevels <- GenomeInfoDb::seqlevels(dmrs)

    expect_warning(
        converted <- CMEnt:::.convertToGRanges(dmrs, genome = "hg19"),
        "no genome information"
    )

    expect_equal(GenomeInfoDb::seqlevels(converted), old_seqlevels)
    expect_true(all(GenomeInfoDb::genome(converted) == "hg19"))
})
