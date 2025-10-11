# Helper functions for tests

#' Create test data for DMR finding tests
#' @param n_dmps Number of DMPs to generate
#' @param n_cpgs Number of CpGs (rows) in beta value matrix
#' @param n_samples Number of samples to generate
#' @param seed Random seed for reproducibility
#' @param platform Array platform to simulate ("450K", "EPIC")
#' @return List containing test data components
create_test_data <- function(n_dmps = 100, n_cpgs = 10000, n_connections = 20, n_samples = 10, n_cases = 5, seed = 42, platform = "450K") {
    set.seed(seed)
    platform <- toupper(platform)
    # Annotation requirements
    if (platform == "450K") {
        if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
            testthat::skip("Illumina 450K annotation not installed")
        }
        locs <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
    } else if (platform == "EPIC") {
        if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
            testthat::skip("Illumina EPIC annotation not installed")
        }
        locs <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
    } else {
        stop("Unsupported platform. Use '450K' or 'EPIC'.")
    }
    all_ids <- rownames(locs)
    if (length(all_ids) < n_cpgs) n_cpgs <- length(all_ids)
    if (n_cpgs < n_dmps) stop("n_cpgs must be >= n_dmps")
    sel_ids <- sample(all_ids, n_cpgs)
    # Sort selected IDs by genomic position to satisfy function precondition
    sel_locs <- locs[sel_ids, , drop = FALSE]
    # Use the same ordering logic as findDMRsFromSeeds (stringr::str_order with numeric=TRUE
    # on a combined chr:pos key) so that the beta file passes the internal
    # sorted-by-position precondition check.
    ord <- stringr::str_order(paste0(sel_locs$chr, ":", sel_locs$pos), numeric = TRUE)
    cpg_ids <- sel_ids[ord]
    sel_locs <- sel_locs[ord, , drop = FALSE]

    beta_values <- matrix(runif(n_cpgs * n_samples, 0, 1),
        nrow = n_cpgs,
        dimnames = list(cpg_ids, paste0("sample", seq_len(n_samples)))
    )
    beta_file <- tempfile(fileext = ".txt")
    write.table(cbind(ID = rownames(beta_values), beta_values),
        file = beta_file,
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    pheno <- data.frame(
        row.names = colnames(beta_values),
        Sample_Group = factor(rep(c("Case", "Control"), c(n_cases, n_samples - n_cases))),
        casecontrol = as.integer(rep(c(1, 0), c(n_cases, n_samples - n_cases)))
    )
    dmp_ids <- sample(cpg_ids, n_dmps)
    dmp_locs <- sel_locs[dmp_ids, , drop = FALSE]
    ord <- stringr::str_order(paste0(dmp_locs$chr, ":", dmp_locs$pos), numeric = TRUE)
    dmp_ids <- dmp_ids[ord]
    dmp_locs <- dmp_locs[ord, , drop = FALSE]
    dmps <- data.frame(
        dmp = dmp_ids,
        chr = dmp_locs$chr,
        pos = dmp_locs$pos,
        pval_adj = runif(n_dmps, 0, 0.01),
        Sample_Group = rep("Case", n_dmps)
    )
    # Add statistical correlation on the beta values of neighboring DMPs in the same chromosome , limited by n_connections

    dmps_file <- tempfile(fileext = ".txt")
    write.table(dmps, file = dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
    list(beta_file = beta_file, beta_values = beta_values, dmps_file = dmps_file, dmps = dmps, pheno = pheno, cpg_ids = cpg_ids)
}
