# Helper functions for tests

#' Create test data for DMR finding tests
#' @param n_cpgs Number of CpGs to generate
#' @param n_samples Number of samples to generate
#' @param seed Random seed for reproducibility
#' @param platform Array platform to simulate ("450K", "EPIC", or "mixed")
#' @return List containing test data components
create_test_data <- function(n_cpgs = 10, n_samples = 10, seed = 42, platform = "450K") {
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
    # Treat any other value as 450K for now
    if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
      testthat::skip("Illumina 450K annotation not installed")
    }
    locs <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
  }
  all_ids <- rownames(locs)
  if (length(all_ids) < n_cpgs) n_cpgs <- length(all_ids)
  sel_ids <- sample(all_ids, n_cpgs)
  # Sort selected IDs by genomic position to satisfy function precondition
  sel_locs <- locs[sel_ids, , drop = FALSE]
  # Use the same ordering logic as findDMRsFromDMPs (stringr::str_order with numeric=TRUE
  # on a combined chr:pos key) so that the beta file passes the internal
  # sorted-by-position precondition check.
  ord <- stringr::str_order(paste0(sel_locs$chr, ":", sel_locs$pos), numeric = TRUE)
  cpg_ids <- sel_ids[ord]
  sel_locs <- sel_locs[ord, , drop = FALSE]

  n_half <- n_samples %/% 2
  beta_values <- matrix(runif(n_cpgs * n_samples, 0, 1), nrow = n_cpgs,
                        dimnames = list(cpg_ids, paste0("sample", seq_len(n_samples))))
  beta_file <- tempfile(fileext = ".txt")
  write.table(cbind(ID = rownames(beta_values), beta_values), file = beta_file,
              sep = "\t", quote = FALSE, row.names = FALSE)
  pheno <- data.frame(row.names = colnames(beta_values),
                      Sample_Group = factor(rep(c("Case", "Control"), c(n_half, n_samples - n_half))),
                      casecontrol = as.integer(rep(c(1,0), c(n_half, n_samples - n_half))))
  k <- max(2, floor(n_cpgs/2))
  dmps <- data.frame(dmp = cpg_ids[seq_len(k)],
                     chr = sel_locs$chr[seq_len(k)],
                     pos = sel_locs$pos[seq_len(k)],
                     pval = runif(k, 0, 0.01),
                     pval_adj = runif(k, 0, 0.01),
                     qval = runif(k, 0, 0.01),
                     delta_beta = runif(k, 0.2, 0.5),
                     cases_beta = runif(k, 0.6, 0.9),
                     controls_beta = runif(k, 0.2, 0.5),
                     cases_beta_sd = 0.1,
                     controls_beta_sd = 0.1,
                     cases_num = n_half,
                     controls_num = n_samples - n_half,
                     Sample_Group = rep("Case", k))
  dmps_file <- tempfile(fileext = ".txt")
  write.table(dmps, file = dmps_file, sep = "\t", quote = FALSE, row.names = FALSE)
  list(beta_file = beta_file, beta_values = beta_values, dmps_file = dmps_file, dmps = dmps, pheno = pheno, cpg_ids = cpg_ids)
}