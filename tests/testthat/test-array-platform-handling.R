context("Array platform handling")

test_that("Package handles 450k probes", {
    # Create 450k-style test data
    test_dmps <- data.frame(
        chr = rep("chr1", 5),
        pos = seq(1000, 5000, by = 1000),
        dmp = paste0("cg", sprintf("%07d", 1:5)),
        pval = rep(0.01, 5),
        delta_beta = rep(0.3, 5),
        stringsAsFactors = FALSE
    )
    
    # Should work with 450k annotations
    expect_error(findDMRsFromDMPs(test_dmps, verbose = TRUE), NA)
})

test_that("Package handles EPIC probes", {
    # Create EPIC-style test data
    test_dmps <- data.frame(
        chr = rep("chr1", 5),
        pos = seq(1000, 5000, by = 1000),
        dmp = paste0("cg", sprintf("%08d", 1:5)),
        pval = rep(0.01, 5),
        delta_beta = rep(0.3, 5),
        stringsAsFactors = FALSE
    )
    
    # Should work with EPIC annotations
    expect_error(findDMRsFromDMPs(test_dmps, verbose = TRUE), NA)
})

test_that("Package handles mixed probe types", {
    # Create mixed 450k and EPIC test data
    test_dmps <- data.frame(
        chr = rep("chr1", 6),
        pos = seq(1000, 6000, by = 1000),
        dmp = c(
            paste0("cg", sprintf("%07d", 1:3)),    # 450k style
            paste0("cg", sprintf("%08d", 4:6))     # EPIC style
        ),
        pval = rep(0.01, 6),
        delta_beta = rep(0.3, 6),
        stringsAsFactors = FALSE
    )
    
    # Should work with both annotation types
    expect_error(findDMRsFromDMPs(test_dmps, verbose = TRUE), NA)
    
    # Verify position-based ordering
    result <- findDMRsFromDMPs(test_dmps, verbose = TRUE)
    expect_true(all(diff(result$pos) > 0))
})

test_that("Package handles custom CpG positions", {
    # Create test data with custom positions
    test_dmps <- data.frame(
        chr = rep("chr1", 3),
        pos = seq(1000, 3000, by = 1000),
        dmp = c("custom_1", "custom_2", "custom_3"),
        pval = rep(0.01, 3),
        delta_beta = rep(0.3, 3),
        stringsAsFactors = FALSE
    )
    
    # Should use positions from DMP data
    expect_error(findDMRsFromDMPs(test_dmps, verbose = TRUE), NA)
})