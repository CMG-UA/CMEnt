test_that("getSortedGenomicLocs reports missing annotation packages clearly", {
    local_mocked_bindings(
        .isPackageInstalled = function(pkg_name) FALSE,
        .package = "CMEnt"
    )
    withr::local_options(list(
        CMEnt.annotation_cache_dir = tempfile("cment-anno-cache-"),
        CMEnt.use_annotation_cache = FALSE
    ))

    expect_error(
        getSortedGenomicLocs(array = "450K", genome = "hg19"),
        "IlluminaHumanMethylation450kanno\\.ilmn12\\.hg19"
    )
    expect_error(
        getSortedGenomicLocs(array = "450K", genome = "hg19"),
        "BiocManager::install"
    )
})

test_that("getBetaHandler surfaces annotation dependency hints for matrix input", {
    local_mocked_bindings(
        .isPackageInstalled = function(pkg_name) FALSE,
        .package = "CMEnt"
    )

    beta <- matrix(runif(6), nrow = 3, dimnames = list(paste0("cg", 1:3), paste0("S", 1:2)))

    expect_error(
        getBetaHandler(beta = beta, array = "450K", genome = "hg19"),
        "getBetaHandler\\(\\) requires the following package"
    )
})

test_that("getDMRSequences reports missing BSgenome packages clearly", {
    local_mocked_bindings(
        .isPackageInstalled = function(pkg_name) FALSE,
        .package = "CMEnt"
    )

    dmrs <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 10000, end = 10050)
    )

    expect_error(
        getDMRSequences(dmrs, genome = "hg19"),
        "BSgenome\\.Hsapiens\\.UCSC\\.hg19"
    )
    expect_error(
        getDMRSequences(dmrs, genome = "hg19"),
        "Install with:"
    )
})

test_that("extractDMRMotifs aggregates sequence dependencies early", {
    local_mocked_bindings(
        .isPackageInstalled = function(pkg_name) FALSE,
        .package = "CMEnt"
    )

    dmrs <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 10000, end = 10100),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg19")
    )
    S4Vectors::mcols(dmrs)$seeds <- "cg00000029"

    err <- expect_error(
        extractDMRMotifs(dmrs, genome = "hg19", array = "450K"),
        class = "error"
    )
    expect_match(conditionMessage(err), "BSgenome\\.Hsapiens\\.UCSC\\.hg19")
    expect_match(conditionMessage(err), "IlluminaHumanMethylation450kanno\\.ilmn12\\.hg19")
})

test_that("annotateDMRsWithGenes reports TxDb and OrgDb requirements clearly", {
    local_mocked_bindings(
        .isPackageInstalled = function(pkg_name) FALSE,
        .package = "CMEnt"
    )

    dmrs <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 10000, end = 10100),
        seqinfo = GenomeInfoDb::Seqinfo(genome = "hg38")
    )

    expect_error(
        annotateDMRsWithGenes(dmrs, genome = "hg38"),
        "TxDb\\.Hsapiens\\.UCSC\\.hg38\\.knownGene"
    )
    expect_error(
        annotateDMRsWithGenes(dmrs, genome = "hg38"),
        "org\\.Hs\\.eg\\.db"
    )
})

test_that("buildDMRs aggregates optional dependencies before running", {
    local_mocked_bindings(
        .isPackageInstalled = function(pkg_name) FALSE,
        .package = "CMEnt"
    )

    beta <- matrix(
        runif(6),
        nrow = 3,
        dimnames = list(c("cg00000029", "cg00000108", "cg00000109"), c("Sample1", "Sample2"))
    )
    seeds <- data.frame(
        site_id = c("cg00000029", "cg00000108"),
        pval = c(0.01, 0.02),
        stringsAsFactors = FALSE
    )
    pheno <- data.frame(
        Sample_Group = c("control", "case"),
        row.names = c("Sample1", "Sample2"),
        stringsAsFactors = FALSE
    )

    err <- expect_error(
        buildDMRs(
            beta = beta,
            seeds = seeds,
            pheno = pheno,
            seeds_id_col = "site_id",
            sample_group_col = "Sample_Group",
            array = "450K",
            genome = "hg19",
            njobs = 1,
            verbose = 0
        ),
        class = "error"
    )
    expect_match(conditionMessage(err), "IlluminaHumanMethylation450kanno\\.ilmn12\\.hg19")
    expect_match(conditionMessage(err), "TxDb\\.Hsapiens\\.UCSC\\.hg19\\.knownGene")
    expect_match(conditionMessage(err), "org\\.Hs\\.eg\\.db")
    expect_match(conditionMessage(err), "BSgenome\\.Hsapiens\\.UCSC\\.hg19")
})

test_that("dependency assertions can auto-install missing packages", {
    installed <- new.env(parent = emptyenv())
    installed$pkgs <- "AlreadyPkg"
    installed$called <- character()

    local_mocked_bindings(
        .isPackageInstalled = function(pkg_name) pkg_name %in% installed$pkgs,
        .installDependencyPackages = function(pkg_names) {
            installed$called <- pkg_names
            installed$pkgs <- union(installed$pkgs, pkg_names)
            invisible(pkg_names)
        },
        .package = "CMEnt"
    )
    withr::local_options(list(CMEnt.auto_install_dep_if_missing = TRUE))

    reqs <- CMEnt:::.makeDependencyRequirements(c("MissingPkg", "AlreadyPkg"), "reason")
    pkgs <- CMEnt:::.assertDependencyRequirements(reqs, "test context")

    expect_equal(installed$called, "MissingPkg")
    expect_equal(pkgs, reqs$pkg_name)
})

test_that("findDMPsBSSeq reports DSS as a required dependency", {
    local_mocked_bindings(
        .isPackageInstalled = function(pkg_name) {
            !identical(pkg_name, "DSS")
        },
        .package = "CMEnt"
    )

    expect_error(
        findDMPsBSSeq(
            bsseq = "dummy.rds",
            samplesheet = data.frame(),
            sample_group_col = "Sample_Group"
        ),
        "DSS"
    )
})
