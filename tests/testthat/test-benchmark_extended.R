library(testthat)

test_that("Extended benchmark vignette components work", {
    skip_on_cran()
    skip_if_not_installed("DMRcate")
    skip_if_not_installed("bumphunter")
    
    beta <- loadExampleInputData("beta")
    pheno <- loadExampleInputData("pheno")
    array_type <- loadExampleInputData("array_type")
    
    beta_handler <- getBetaHandler(beta, array = array_type, genome = "hg19")
    beta_mat <- as.matrix(beta_handler$getBeta())
    locs <- beta_handler$getBetaLocs()
    mvalues <- log2(beta_mat / (1 - beta_mat + 1e-6) + 1e-6)
    
    expect_true(nrow(mvalues) > 0)
    expect_true(ncol(mvalues) == nrow(pheno))
})

test_that("DMR statistics calculation works", {
    example_output_path <- system.file("extdata", "example_output.rds", package = "DMRsegal")
    dmrs <- readRDS(example_output_path)
    
    get_dmr_stats <- function(dmrs, method) {
        if (length(dmrs) == 0) {
            return(data.frame(
                Method = method,
                Number_of_DMRs = 0,
                Mean_DMR_Width = NA,
                Median_DMR_Width = NA,
                Total_Coverage_bp = 0,
                Min_Width = NA,
                Max_Width = NA
            ))
        }
        
        data.frame(
            Method = method,
            Number_of_DMRs = length(dmrs),
            Mean_DMR_Width = mean(width(dmrs)),
            Median_DMR_Width = median(width(dmrs)),
            Total_Coverage_bp = sum(width(dmrs)),
            Min_Width = min(width(dmrs)),
            Max_Width = max(width(dmrs))
        )
    }
    
    stats <- get_dmr_stats(dmrs, "TestMethod")
    
    expect_equal(nrow(stats), 1)
    expect_true(stats$Number_of_DMRs > 0)
    expect_true(stats$Mean_DMR_Width > 0)
    expect_true(stats$Total_Coverage_bp > 0)
})

test_that("Overlap calculations work correctly", {
    example_output_path <- system.file("extdata", "example_output.rds", package = "DMRsegal")
    dmrs1 <- readRDS(example_output_path)
    dmrs2 <- dmrs1[1:min(5, length(dmrs1))]
    
    calculate_overlap <- function(gr1, gr2) {
        if (length(gr1) == 0 || length(gr2) == 0) {
            return(list(count = 0, percentage = 0))
        }
        overlap_count <- length(subsetByOverlaps(gr1, gr2))
        pct <- (overlap_count / length(gr1)) * 100
        list(count = overlap_count, percentage = pct)
    }
    
    overlap_result <- calculate_overlap(dmrs1, dmrs2)
    
    expect_true(overlap_result$count >= 0)
    expect_true(overlap_result$percentage >= 0)
    expect_true(overlap_result$percentage <= 100)
    expect_equal(overlap_result$count, length(dmrs2))
    expect_equal(overlap_result$percentage, (length(dmrs2) / length(dmrs1)) * 100)
})

test_that("Jaccard index calculation works", {
    example_output_path <- system.file("extdata", "example_output.rds", package = "DMRsegal")
    dmrs1 <- readRDS(example_output_path)
    dmrs2 <- dmrs1[1:min(5, length(dmrs1))]
    
    jaccard_index <- function(gr1, gr2) {
        if (length(gr1) == 0 || length(gr2) == 0) return(0)
        overlap_count <- length(subsetByOverlaps(gr1, gr2))
        total_unique <- length(unique(c(gr1, gr2)))
        if (total_unique == 0) return(0)
        overlap_count / total_unique
    }
    
    jaccard <- jaccard_index(dmrs1, dmrs2)
    
    expect_true(jaccard >= 0)
    expect_true(jaccard <= 1)
})

test_that("IoU calculation works", {
    example_output_path <- system.file("extdata", "example_output.rds", package = "DMRsegal")
    dmrs1 <- readRDS(example_output_path)
    dmrs2 <- dmrs1[1:min(5, length(dmrs1))]
    
    iou <- function(gr1, gr2) {
        if (length(gr1) == 0 || length(gr2) == 0) return(0)
        intersection <- sum(width(GenomicRanges::intersect(gr1, gr2)))
        union <- sum(width(GenomicRanges::union(gr1, gr2)))
        if (union == 0) return(0)
        intersection / union
    }
    
    iou_val <- iou(dmrs1, dmrs2)
    
    expect_true(iou_val >= 0)
    expect_true(iou_val <= 1)
})

test_that("Consensus DMR finding works", {
    example_output_path <- system.file("extdata", "example_output.rds", package = "DMRsegal")
    dmrs <- readRDS(example_output_path)
    
    dmrs1 <- dmrs[1:min(10, length(dmrs))]
    dmrs2 <- dmrs[5:min(15, length(dmrs))]
    dmrs3 <- dmrs[8:min(18, length(dmrs))]
    
    methods_list <- list(
        Method1 = dmrs1,
        Method2 = dmrs2,
        Method3 = dmrs3
    )
    
    find_consensus <- function(methods_list, min_methods = 2) {
        all_dmrs <- do.call(c, unname(methods_list))
        all_dmrs <- GenomicRanges::reduce(all_dmrs)
        
        consensus_counts <- vapply(seq_along(all_dmrs), function(i) {
            region <- all_dmrs[i]
            sum(sapply(methods_list, function(method_dmrs) {
                length(subsetByOverlaps(region, method_dmrs)) > 0
            }))
        }, integer(1))
        
        consensus_dmrs <- all_dmrs[consensus_counts >= min_methods]
        mcols(consensus_dmrs)$n_methods <- consensus_counts[consensus_counts >= min_methods]
        
        return(consensus_dmrs)
    }
    
    consensus_2plus <- find_consensus(methods_list, min_methods = 2)
    consensus_3plus <- find_consensus(methods_list, min_methods = 3)
    
    expect_true(length(consensus_2plus) >= length(consensus_3plus))
    expect_true(all(mcols(consensus_2plus)$n_methods >= 2))
    if (length(consensus_3plus) > 0) {
        expect_true(all(mcols(consensus_3plus)$n_methods >= 3))
    }
})

test_that("Method-specific DMR identification works", {
    example_output_path <- system.file("extdata", "example_output.rds", package = "DMRsegal")
    dmrs <- readRDS(example_output_path)
    
    dmrs1 <- dmrs[1:min(10, length(dmrs))]
    dmrs2 <- dmrs[6:min(15, length(dmrs))]
    
    methods_list <- list(
        Method1 = dmrs1,
        Method2 = dmrs2
    )
    
    find_specific_dmrs <- function(target_method, methods_list) {
        target_dmrs <- methods_list[[target_method]]
        if (length(target_dmrs) == 0) return(GRanges())
        
        other_methods <- methods_list[names(methods_list) != target_method]
        all_others <- do.call(c, unname(other_methods))
        
        if (length(all_others) == 0) {
            return(target_dmrs)
        }
        
        specific <- target_dmrs[countOverlaps(target_dmrs, all_others) == 0]
        return(specific)
    }
    
    specific1 <- find_specific_dmrs("Method1", methods_list)
    specific2 <- find_specific_dmrs("Method2", methods_list)
    
    expect_true(length(specific1) <= length(dmrs1))
    expect_true(length(specific2) <= length(dmrs2))
    
    if (length(specific1) > 0) {
        expect_equal(as.vector(countOverlaps(specific1, dmrs2)), rep(0, length(specific1)))
    }
})

test_that("Empty DMR sets are handled gracefully", {
    empty_dmrs <- GRanges()
    example_output_path <- system.file("extdata", "example_output.rds", package = "DMRsegal")
    dmrs <- readRDS(example_output_path)
    
    get_dmr_stats <- function(dmrs, method) {
        if (length(dmrs) == 0) {
            return(data.frame(
                Method = method,
                Number_of_DMRs = 0,
                Mean_DMR_Width = NA,
                Median_DMR_Width = NA,
                Total_Coverage_bp = 0,
                Min_Width = NA,
                Max_Width = NA
            ))
        }
        
        data.frame(
            Method = method,
            Number_of_DMRs = length(dmrs),
            Mean_DMR_Width = mean(width(dmrs)),
            Median_DMR_Width = median(width(dmrs)),
            Total_Coverage_bp = sum(width(dmrs)),
            Min_Width = min(width(dmrs)),
            Max_Width = max(width(dmrs))
        )
    }
    
    stats <- get_dmr_stats(empty_dmrs, "EmptyMethod")
    expect_equal(stats$Number_of_DMRs, 0)
    expect_equal(stats$Total_Coverage_bp, 0)
    expect_true(is.na(stats$Mean_DMR_Width))
    
    calculate_overlap <- function(gr1, gr2) {
        if (length(gr1) == 0 || length(gr2) == 0) {
            return(list(count = 0, percentage = 0))
        }
        overlap_count <- length(subsetByOverlaps(gr1, gr2))
        pct <- (overlap_count / length(gr1)) * 100
        list(count = overlap_count, percentage = pct)
    }
    
    overlap <- calculate_overlap(empty_dmrs, dmrs)
    expect_equal(overlap$count, 0)
    expect_equal(overlap$percentage, 0)
    
    overlap2 <- calculate_overlap(dmrs, empty_dmrs)
    expect_equal(overlap2$count, 0)
    expect_equal(overlap2$percentage, 0)
})
