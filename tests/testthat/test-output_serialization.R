test_that("non-tabular DMR output columns round-trip through TSV serialization", {
    dmrs_df <- data.frame(
        chr = c("chr1", "chr1"),
        start = c(100L, 200L),
        end = c(150L, 250L),
        stringsAsFactors = FALSE
    )
    dmrs_df$pwm <- list(
        matrix(c(0.2, 0.3, 0.3, 0.2), nrow = 2),
        matrix(c(0.1, 0.4, 0.4, 0.1), nrow = 2)
    )

    encoded <- CMEnt:::.encodeNonTabularColumns(dmrs_df)

    expect_equal(encoded$encoded_columns, "pwm")
    expect_true("pwm" %in% colnames(encoded$data))
    expect_true(is.character(encoded$data$pwm))
    expect_true(all(startsWith(encoded$data$pwm, CMEnt:::.serializedOutputPrefix)))

    tsv_file <- tempfile(fileext = ".tsv.gz")
    gz <- gzfile(tsv_file, "w")
    write.table(
        encoded$data,
        gz,
        sep = "\t",
        quote = FALSE,
        qmethod = "double",
        col.names = TRUE,
        row.names = FALSE
    )
    close(gz)

    reloaded_df <- read.delim(gzfile(tsv_file), check.names = FALSE, stringsAsFactors = FALSE)
    roundtrip <- CMEnt:::convertToGRanges(reloaded_df, genome = "hg38")

    expect_true("pwm" %in% names(S4Vectors::mcols(roundtrip)))
    expect_true(is.matrix(S4Vectors::mcols(roundtrip)$pwm[[1]]))
    expect_equal(
        S4Vectors::mcols(roundtrip)$pwm[[1]],
        dmrs_df$pwm[[1]]
    )
    expect_equal(
        S4Vectors::mcols(roundtrip)$pwm[[2]],
        dmrs_df$pwm[[2]]
    )

    unlink(tsv_file)
})
