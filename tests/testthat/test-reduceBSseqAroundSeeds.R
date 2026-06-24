test_that("reduceBSSeqAroundSeeds returns a BSseq object with the expected dimensions", {
    # Create a mock BSseq object
    mock_bsseq <- bsseq::BSseq(
        chr = c("chr1", "chr1", "chr2"),
        pos = c(100, 200, 300),
        M = matrix(c(10, 20, 30, 40, 50, 60), nrow = 3),
        Cov = matrix(c(100, 200, 300, 400, 500, 600), nrow = 3),
        sampleNames = c("Sample1", "Sample2")
    )

    # Define seed positions
    seed_positions <- data.frame(
        chr = c("chr1", "chr2"),
        start = c(150, 250),
        end = c(250, 350)
    )

    # Call the function
    reduced_bsseq <- reduceBSseqAroundSeeds(mock_bsseq, seed_positions, expansion_window = 50)

    # Check that the output is a BSseq object
    expect_s4_class(reduced_bsseq, "BSseq")

    # Check that the dimensions of the output are as expected
    expect_equal(dim(reduced_bsseq), c(2, 2)) # Expecting 2 rows (one for each seed) and 2 columns (samples)
})