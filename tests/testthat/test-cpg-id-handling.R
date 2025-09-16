context("CpG ID handling")

test_that("Initial regions preserve CpG IDs", {
  # Create test data with mixed 450k and EPIC probes
  test_dmps <- data.frame(
    chr = c("chr1", "chr1", "chr1", "chr1"),
    pos = c(1000, 1200, 1400, 1600),
    dmp = c("cg01", "cg00000002", "cg03", "cg00000004"),  # Mix of formats
    pval = c(0.01, 0.02, 0.01, 0.03),
    delta_beta = c(0.2, 0.3, 0.25, 0.15),
    stringsAsFactors = FALSE
  )
  
  # Test initial region finding
  regions <- .findInitialRegions(test_dmps, max.lookup.dist = 300, min.dmps = 2)
  
  # Verify CpG IDs are preserved
  expect_equal(regions$start_dmp, "cg01")
  expect_equal(regions$end_dmp, "cg00000004")  # EPIC format
  expect_equal(regions$n_dmps, 4)
  
  # Verify position-based ordering is used, not string-based
  expect_true(all(diff(regions$pos) > 0))
})

test_that("DMR expansion maintains correct CpG IDs", {
  # Create test data for initial regions
  test_regions <- data.frame(
    chr = "chr1",
    start = 1000,
    end = 1600,
    start_dmp = "cg01",
    end_dmp = "cg04",
    n_dmps = 4,
    mean_pval = 0.02,
    stringsAsFactors = FALSE
  )
  
  # Create test methylation data with mixed probe types
  test_meth_data <- data.frame(
    chr = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1"),
    pos = c(800, 1000, 1200, 1400, 1600, 1800),
    cpg = c("cg00000000", "cg01", "cg00000002", "cg03", "cg00000004", "cg00000005"),  # Mix of formats
    stringsAsFactors = FALSE
  )
  
  # Test DMR expansion
  expanded <- .expandDMRs(test_regions, test_meth_data, lookup.dist = 300)
  
  # Verify CpG IDs are maintained after expansion
  expect_true(all(c("start_dmp", "end_dmp") %in% names(expanded)))
  expect_equal(expanded$start_dmp, "cg00")  # Should include one CpG upstream
  expect_equal(expanded$end_dmp, "cg05")    # Should include one CpG downstream
})

test_that("findDMRsFromDMPs maintains CpG ID integrity", {
  # Create minimal test data
  test_dmps <- data.frame(
    chr = c("chr1", "chr1", "chr1"),
    pos = c(1000, 1200, 1400),
    dmp = c("cg01", "cg02", "cg03"),
    pval = c(0.01, 0.02, 0.01),
    delta_beta = c(0.2, 0.3, 0.25),
    stringsAsFactors = FALSE
  )
  
  # Create test methylation data
  test_meth_data <- data.frame(
    chr = c("chr1", "chr1", "chr1", "chr1", "chr1"),
    pos = c(800, 1000, 1200, 1400, 1600),
    cpg = c("cg00", "cg01", "cg02", "cg03", "cg04"),
    stringsAsFactors = FALSE
  )
  
  # Run full DMR finding process
  dmrs <- findDMRsFromDMPs(test_dmps, test_meth_data)
  
  # Verify CpG ID fields exist and contain correct values
  expect_true(all(c("start_dmp", "end_dmp") %in% names(dmrs)))
  expect_type(dmrs$start_dmp, "character")
  expect_type(dmrs$end_dmp, "character")
  expect_true(all(dmrs$start_dmp %in% test_meth_data$cpg))
  expect_true(all(dmrs$end_dmp %in% test_meth_data$cpg))
})