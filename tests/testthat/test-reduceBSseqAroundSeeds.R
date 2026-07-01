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

test_that("sparsifySeeds removes interior seeds from dense runs per chromosome", {
    seeds <- data.frame(
        chr = c("chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr2"),
        start = c(105, 100, 109, 130, 100, 109, 140),
        score = seq_len(7),
        row.names = paste0("seed", seq_len(7))
    )

    sparse <- sparsifySeeds(seeds, min_distance = 10)

    expect_identical(rownames(sparse), c("seed2", "seed3", "seed4", "seed5", "seed6", "seed7"))
    expect_identical(sparse$score, c(2L, 3L, 4L, 5L, 6L, 7L))
})

test_that("sparsifySeeds keeps close pairs and seeds at the exact minimum distance", {
    seeds <- data.frame(
        chr = c("chr1", "chr1", "chr1", "chr2", "chr2"),
        start = c(100, 110, 119, 100, 109),
        row.names = c("a", "b", "c", "d", "e")
    )

    sparse <- sparsifySeeds(seeds, min_distance = 10)

    expect_identical(rownames(sparse), c("a", "b", "c", "d", "e"))
})

test_that("sparsifySeeds leaves seeds unchanged when minimum distance is zero", {
    seeds <- data.frame(chr = c("chr1", "chr1"), start = c(100, 100))

    sparse <- sparsifySeeds(seeds, min_distance = 0)

    expect_identical(sparse, seeds)
})

test_that("sparsifySeeds validates required seed columns and distance", {
    seeds <- data.frame(chr = "chr1", start = 100)

    expect_error(
        sparsifySeeds(data.frame(chr = "chr1"), min_distance = 10),
        "Seeds data frame must contain 'chr' and 'start' columns",
        fixed = TRUE
    )
    expect_error(
        sparsifySeeds(seeds, min_distance = -1),
        "min_distance must be a non-negative numeric value",
        fixed = TRUE
    )
})
