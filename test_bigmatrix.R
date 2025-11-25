library(DMRsegal)

beta_matrix <- loadExampleInputData("beta")
cat("Loaded beta matrix with", nrow(beta_matrix), "rows and", ncol(beta_matrix), "cols\n")

beta_handler <- getBetaHandler(
    beta = beta_matrix,
    array = "450K",
    genome = "hg19",
    memory_threshold_mb  = 0.001
)

row_names_to_test <- rownames(beta_matrix)[1:100]
result <- beta_handler$getBeta(row_names = row_names_to_test)

cat("Result class:", class(result), "\n")
cat("Result dimensions:", nrow(result), "x", ncol(result), "\n")

if (bigmemory::is.big.matrix(result)) {
    cat("SUCCESS: Result is a big.matrix as expected\n")
    cat("First few values:\n")
    print(result[1:3, 1:3])
} else {
    cat("Result is a regular matrix\n")
}

beta_handler2 <- getBetaHandler(
    beta = beta_matrix,
    array = "450K",
    genome = "hg19",
    memory_threshold_mb  = 1000
)

result2 <- beta_handler2$getBeta(row_names = row_names_to_test)
cat("\nWith high threshold - Result class:", class(result2), "\n")

if (!bigmemory::is.big.matrix(result2)) {
    cat("SUCCESS: Result is NOT a big.matrix when threshold is high\n")
} else {
    cat("FAIL: Result should not be a big.matrix\n")
}

