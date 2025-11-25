test_that("getBeta returns big.matrix when size exceeds threshold with file input", {
    beta <- DMRsegal::loadExampleInputData("beta")[1:1000, ]
    
    temp_file <- tempfile(fileext = ".tsv")
    beta_df <- cbind(CpG = rownames(beta), as.data.frame(beta))
    write.table(beta_df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    beta_handler <- getBetaHandler(
        beta = temp_file,
        array = "450K",
        genome = "hg19",
        memory_threshold_mb = 0.1,
    )
    
    row_names_subset <- rownames(beta)[1:100]
    result <- beta_handler$getBeta(row_names = row_names_subset)
    
    expect_true(bigmemory::is.big.matrix(result))
    expect_equal(nrow(result), 100)
    expect_equal(ncol(result), ncol(beta))
    
    unlink(temp_file)
})

test_that("getBeta returns regular matrix when size below threshold", {
    beta <- DMRsegal::loadExampleInputData("beta")[1:1000, ]
    
    temp_file <- tempfile(fileext = ".tsv")
    beta_df <- cbind(CpG = rownames(beta), as.data.frame(beta))
    write.table(beta_df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    beta_handler <- getBetaHandler(
        beta = temp_file,
        array = "450K",
        genome = "hg19",
        memory_threshold_mb = 100,
    )
    
    row_names_subset <- rownames(beta)[1:100]
    result <- beta_handler$getBeta(row_names = row_names_subset)
    
    expect_false(bigmemory::is.big.matrix(result))
    expect_equal(nrow(result), 100)
    expect_equal(ncol(result), ncol(beta))
    
    unlink(temp_file)
})

test_that("getBeta big.matrix has correct values", {
    beta <- DMRsegal::loadExampleInputData("beta")[1:1000, ]
    
    temp_file <- tempfile(fileext = ".tsv")
    beta_df <- cbind(CpG = rownames(beta), as.data.frame(beta))
    write.table(beta_df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    beta_handler <- getBetaHandler(
        beta = temp_file,
        array = "450K",
        genome = "hg19",
        memory_threshold_mb = 0.1,
    )
    
    row_names_subset <- rownames(beta)[1:10]
    result_big <- beta_handler$getBeta(row_names = row_names_subset)
    
    beta_handler_normal <- getBetaHandler(
        beta = temp_file,
        array = "450K",
        genome = "hg19",
        memory_threshold_mb = 0.1
    )
    result_normal <- beta_handler_normal$getBeta(row_names = row_names_subset)
    
    expect_equal(as.matrix(result_big[1:5, 1:5]), as.matrix(result_normal[1:5, 1:5]))
    
    unlink(temp_file)
})

test_that("getBeta big.matrix with column subset", {
    beta <- DMRsegal::loadExampleInputData("beta")[1:1000, ]
    
    temp_file <- tempfile(fileext = ".tsv")
    beta_df <- cbind(CpG = rownames(beta), as.data.frame(beta))
    write.table(beta_df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    beta_handler <- getBetaHandler(
        beta = temp_file,
        array = "450K",
        genome = "hg19",
        memory_threshold_mb = 0.1
    )
    
    row_names_subset <- rownames(beta)[1:50]
    col_names_subset <- colnames(beta)[1:10]
    result <- beta_handler$getBeta(row_names = row_names_subset, col_names = col_names_subset)
    
    expect_true(bigmemory::is.big.matrix(result))
    expect_equal(nrow(result), 50)
    expect_equal(ncol(result), 10)
    expect_equal(colnames(result), col_names_subset)
    
    unlink(temp_file)
})

test_that("getBeta big.matrix handles all rows", {
    beta <- DMRsegal::loadExampleInputData("beta")[1:1000, ]
    
    temp_file <- tempfile(fileext = ".tsv")
    beta_df <- cbind(CpG = rownames(beta), as.data.frame(beta))
    write.table(beta_df, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    beta_handler <- getBetaHandler(
        beta = temp_file,
        array = "450K",
        genome = "hg19",
        memory_threshold_mb = 0,
        bigmatrix_threshold_mb = 0.1
    )
    
    result <- beta_handler$getBeta()
    
    expect_true(bigmemory::is.big.matrix(result))
    expect_equal(nrow(result), 1000)
    expect_equal(ncol(result), ncol(beta))
    
    unlink(temp_file)
})
