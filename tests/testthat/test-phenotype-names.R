test_that("covariate residualization handles escaped phenotype column names", {
    fixture <- makeSyntheticBuildDMRsFixture()
    expected <- CMEnt:::.transformBeta(fixture$beta, fixture$pheno, covariates = "Age")

    for (column_name in c("Age (years)", "age`years", "path\\score")) {
        pheno <- fixture$pheno
        colnames(pheno)[colnames(pheno) == "Age"] <- column_name
        expect_equal(
            CMEnt:::.transformBeta(fixture$beta, pheno, covariates = column_name),
            expected,
            info = column_name
        )
    }
})

test_that("array DMP calling preserves colliding and escaped covariate names", {
    fixture <- makeSyntheticBuildDMRsFixture()
    pheno <- fixture$pheno
    pheno$Sample_ID <- rownames(pheno)
    pheno$Sample_Group <- ifelse(pheno$casecontrol == 1L, "case", "control")
    pheno$condition <- pheno$Age
    pheno$condition.1 <- "existing"
    args <- list(
        beta = fixture$beta, samplesheet = pheno, sample_group_col = "Sample_Group",
        array = "450K", genome = "hg38", sorted_locs = fixture$locs,
        case_group = "case", chr = "all"
    )
    expected <- do.call(findDMPsArray, c(args, list(covariates = "Age")))
    unadjusted <- do.call(findDMPsArray, args)
    expect_false(isTRUE(all.equal(expected$pval, unadjusted$pval)))

    for (column_name in c("condition", "Age (years)", "age`years", "path\\score")) {
        args$samplesheet[[column_name]] <- pheno$Age
        expect_equal(
            do.call(findDMPsArray, c(args, list(covariates = column_name))),
            expected,
            info = column_name
        )
    }
})

test_that("BSseq DMP calling preserves colliding and escaped covariate names", {
    skip_if_not_installed("DSS")
    skip_if_not_installed("bsseq")
    withr::local_seed(42)
    coverage <- matrix(40L, nrow = 12, ncol = 6)
    prob <- matrix(rep(c(0.25, 0.35, 0.45, 0.55, 0.65, 0.75), each = 12), nrow = 12)
    meth <- matrix(
        rbinom(length(coverage), size = as.vector(coverage), prob = as.vector(prob)),
        nrow = 12
    )
    bs <- bsseq::BSseq(
        chr = rep("chr1", 12), pos = seq_len(12) * 100L,
        M = meth, Cov = coverage, sampleNames = paste0("S", 1:6)
    )
    pheno <- data.frame(
        Sample_ID = colnames(bs), Sample_Group = rep(c("control", "case"), each = 3),
        Age = rep(c(50, 54, 58), 2)
    )
    pheno$condition <- pheno$Age
    pheno$condition.1 <- "existing"
    args <- list(
        bsseq = bs, samplesheet = pheno, sample_group_col = "Sample_Group",
        case_group = "case", chr = "all", njobs = 1
    )
    expected <- do.call(findDMPsBSSeq, c(args, list(covariates = "Age")))
    unadjusted <- do.call(findDMPsBSSeq, args)
    expect_false(isTRUE(all.equal(expected$pval, unadjusted$pval)))

    for (column_name in c("condition", "Age (years)", "age`years", "path\\score")) {
        args$samplesheet[[column_name]] <- pheno$Age
        expect_equal(
            do.call(findDMPsBSSeq, c(args, list(covariates = column_name))),
            expected,
            info = column_name
        )
    }
})

test_that("plotting preserves phenotype column names from files", {
    fixture <- makeSyntheticPlotFixture()
    pheno <- fixture$pheno
    colnames(pheno)[colnames(pheno) == "Sample_Group"] <- "Sample group"
    pheno$Sample.group <- "different group"
    pheno_file <- tempfile(fileext = ".tsv")
    plot_file <- tempfile(fileext = ".pdf")
    circos_file <- tempfile(fileext = ".pdf")
    withr::defer(unlink(c(pheno_file, plot_file, circos_file)))
    write.table(pheno, pheno_file, sep = "\t", quote = FALSE, col.names = NA)

    plot <- plotDMR(
        fixture$dmrs, dmr_index = 1, beta = fixture$beta_handler,
        pheno = pheno_file, sample_group_col = "Sample group",
        genome = "hg38", plot_motif = FALSE, output_file = plot_file
    )
    expect_s3_class(plot, "gtable")

    local_mocked_bindings(
        .getCytobandData = function(...) {
            data.frame(
                V1 = c("chr1", "chr2"), V2 = 0L, V3 = 1e6L,
                V4 = "p", V5 = "gneg"
            )
        },
        .prepareCircosLinkData = function(...) NULL,
        .package = "CMEnt"
    )
    dmrs <- fixture$dmrs
    S4Vectors::mcols(dmrs)$pwm <- replicate(length(dmrs), matrix(0.25, 4, 10), simplify = FALSE)
    state <- CMEnt:::.prepareCircosPlotState(
        dmrs, beta = fixture$beta_handler, pheno = pheno_file,
        sample_group_col = "Sample group", max_num_samples_per_group = 0,
        query_components_with_jaspar = FALSE
    )
    expect_identical(colnames(state$reduced_pheno), colnames(pheno))
    expect_identical(state$reduced_pheno[["Sample group"]], pheno[["Sample group"]])
    expect_identical(state$reduced_pheno$Sample.group, pheno$Sample.group)

    expect_no_error(plotDMRsCircos(
        dmrs, beta = fixture$beta_handler, pheno = pheno_file,
        sample_group_col = "Sample group", max_num_samples_per_group = 0,
        query_components_with_jaspar = FALSE, output_file = circos_file
    ))
    expect_gt(file.info(circos_file)$size, 0)
})
