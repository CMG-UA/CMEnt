test_that("getRegistry supports data.frame and TSV inputs", {
    df <- data.frame(
        id = sprintf("id%04d", 1:20),
        grp = rep(c("g1", "g2"), each = 10),
        value = seq_len(20),
        stringsAsFactors = FALSE
    )

    reg_mem <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        memory_threshold_mb = Inf
    )
    expect_equal(reg_mem$backend(), "memory")
    expect_equal(reg_mem$nrow(), nrow(df))
    expect_equal(reg_mem$ncol(), ncol(df))

    tsv <- tempfile(fileext = ".tsv")
    on.exit(unlink(tsv, force = TRUE), add = TRUE)
    write.table(df, file = tsv, sep = "\t", quote = FALSE, row.names = FALSE)

    reg_tsv <- getRegistry(
        data = tsv,
        indices = list(by_id = "id"),
        memory_threshold_mb = Inf
    )
    expect_equal(reg_tsv$backend(), "memory")
    expected_tsv <- df
    rownames(expected_tsv) <- expected_tsv$id
    expect_equal(as.data.frame(reg_tsv), expected_tsv)
})


test_that("indices accepts a single string column name", {
    df <- data.frame(
        id = sprintf("id%03d", 1:8),
        grp = rep(c("g1", "g2"), each = 4),
        value = seq_len(8),
        stringsAsFactors = FALSE
    )

    reg <- getRegistry(
        data = df,
        indices = "id",
        memory_threshold_mb = Inf
    )

    # Single-string indices are normalized to index name == column name.
    out <- reg$queryOnIndex("id", c("id006", "id001"))
    expect_equal(out$id, c("id006", "id001"))

    # Auto-selected master index should also make character row subset work.
    reg_view <- reg[c("id004", "id002"), c("id", "value")]
    bounded <- as.data.frame(reg_view)
    expect_equal(bounded$id, c("id004", "id002"))
    expect_equal(bounded$value, c(4, 2))
    expect_equal(nrow(reg), 8)
})


test_that("queryOnIndex preserves key input order and duplicates for single index", {
    df <- data.frame(
        id = c("a", "b", "c", "d"),
        grp = c("g1", "g2", "g1", "g2"),
        value = c(10, 20, 30, 40),
        stringsAsFactors = FALSE
    )

    reg <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        memory_threshold_mb = Inf
    )

    out <- reg$queryOnIndex("by_id", c("c", "a", "c", "x"))
    expect_equal(out$id, c("c", "a", "c"))
    expect_equal(out$value, c(30, 10, 30))
})


test_that("queryOnIndex works for composite indices in memory and sqlite", {
    df <- data.frame(
        k1 = c("A", "A", "B", "A", "B"),
        k2 = c("x", "x", "y", "x", "y"),
        payload = c(1, 2, 3, 4, 5),
        stringsAsFactors = FALSE
    )
    keys <- data.frame(
        k1 = c("A", "B", "A"),
        k2 = c("x", "y", "x"),
        stringsAsFactors = FALSE
    )

    reg_mem <- getRegistry(
        data = df,
        indices = list(by_pair = c("k1", "k2")),
        memory_threshold_mb = Inf
    )
    out_mem <- reg_mem$queryOnIndex("by_pair", keys)
    expect_equal(out_mem$payload, c(1, 2, 4, 3, 5, 1, 2, 4))

    reg_sql <- getRegistry(
        data = df,
        indices = list(by_pair = c("k1", "k2")),
        memory_threshold_mb = 0
    )
    out_sql <- reg_sql$queryOnIndex("by_pair", keys)
    expect_equal(out_sql$payload, c(1, 2, 4, 3, 5, 1, 2, 4))
    reg_sql$close()
})


test_that("sqlite queryOnIndex handles lookup cardinalities above SQLite variable limit", {
    n <- 2200L
    df <- data.frame(
        id = sprintf("id%05d", seq_len(n)),
        grp = rep(c("g1", "g2"), length.out = n),
        value = seq_len(n),
        stringsAsFactors = FALSE
    )
    keys <- sprintf("id%05d", sample.int(n, 1200))

    reg <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        memory_threshold_mb = 0
    )
    out <- reg$queryOnIndex("by_id", keys)
    expect_equal(nrow(out), length(keys))
    expect_equal(out$id, keys)
    reg$close()
})


test_that("subsetting returns a new bounded view and leaves original unchanged", {
    df <- data.frame(
        id = sprintf("id%03d", 1:30),
        grp = rep(letters[1:3], each = 10),
        value = seq_len(30),
        score = seq_len(30) * 2,
        stringsAsFactors = FALSE
    )

    reg <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        memory_threshold_mb = Inf
    )

    ret <- reg[5:15, c("score", "id")]
    expect_false(identical(ret, reg))

    orig_bounds <- reg$getActiveBounds()
    expect_equal(orig_bounds$row_idx, seq_len(nrow(df)))
    expect_equal(orig_bounds$col_names, colnames(df))

    bounds <- ret$getActiveBounds()
    expect_equal(bounds$row_idx, 5:15)
    expect_equal(bounds$col_names, c("score", "id"))

    out_df <- as.data.frame(ret)
    expected_subset <- df[5:15, c("score", "id"), drop = FALSE]
    rownames(expected_subset) <- expected_subset$id
    expect_equal(out_df, expected_subset)

    # Query is constrained by active rows and active columns.
    out_q <- ret$queryOnIndex("by_id", c("id001", "id007", "id015"))
    expect_equal(out_q$id, c("id007", "id015"))
    expect_equal(colnames(out_q), c("score", "id"))

    ret$resetBounds()
    reset_df <- as.data.frame(ret)
    expected_reset <- df
    rownames(expected_reset) <- expected_reset$id
    expect_equal(reset_df, expected_reset)
})


test_that("sqlite subset view reuses backend metadata and backing file path", {
    df <- data.frame(
        id = sprintf("id%03d", 1:25),
        grp = rep(c("a", "b", "c", "d", "e"), each = 5),
        value = seq_len(25),
        stringsAsFactors = FALSE
    )

    reg <- getRegistry(
        data = df,
        indices = list(by_id = "id", by_grp = "grp"),
        master_index = "by_id",
        memory_threshold_mb = 0
    )
    on.exit(reg$close(), add = TRUE)

    view <- reg[6:10, c("id", "value")]
    on.exit(view$close(), add = TRUE)

    expect_false(identical(view, reg))
    expect_equal(view$backend(), reg$backend())
    expect_equal(view$getMasterIndex(), reg$getMasterIndex())
    expect_equal(names(view$.__enclos_env__$private$.indices), names(reg$.__enclos_env__$private$.indices))
    expect_identical(
        view$.__enclos_env__$private$.sqlite_path,
        reg$.__enclos_env__$private$.sqlite_path
    )
    expect_equal(nrow(reg), nrow(df))
    expect_equal(nrow(view), 5)
})


test_that("data.frame-like accessors work on active bounds", {
    df <- data.frame(
        id = sprintf("id%03d", 1:12),
        grp = rep(c("a", "b", "c"), length.out = 12),
        value = seq_len(12),
        stringsAsFactors = FALSE
    )

    reg <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        memory_threshold_mb = Inf
    )

    reg_view <- reg[3:8, c("id", "value")]
    expect_equal(nrow(reg), 12)
    expect_equal(ncol(reg), 3)
    expect_equal(nrow(reg_view), 6)
    expect_equal(ncol(reg_view), 2)
    expect_equal(names(reg_view), c("id", "value"))
    expect_equal(colnames(reg_view), c("id", "value"))
    expect_equal(reg_view$id, df$id[3:8])
    expect_equal(reg_view[["value"]], df$value[3:8])
    expect_equal(reg_view[[2]], df$value[3:8])
    expect_equal(
        paste(reg_view$id, reg_view$value, sep = ":"),
        paste(df$id[3:8], df$value[3:8], sep = ":")
    )
    expect_true(is.function(reg_view$queryOnIndex))
    expect_null(reg_view$not_a_column)
})


test_that("character row subsetting uses master index", {
    df <- data.frame(
        id = c("id1", "id2", "id3", "id3", "id4"),
        grp = c("a", "a", "b", "b", "c"),
        value = c(1, 2, 3, 4, 5),
        stringsAsFactors = FALSE
    )

    reg <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = Inf
    )

    reg_view <- reg[c("id3", "id1", "id3"), c("id", "value")]
    bounds <- reg_view$getActiveBounds()
    expect_equal(bounds$row_idx, c(3L, 4L, 1L, 3L, 4L))
    expect_equal(as.data.frame(reg_view)$id, c("id3", "id3", "id1", "id3", "id3"))
    expect_equal(as.data.frame(reg_view)$value, c(3, 4, 1, 3, 4))
    reg_view_factor <- reg[factor(c("id3", "id1", "id3")), c("id", "value")]
    expect_equal(as.data.frame(reg_view_factor)$id, c("id3", "id3", "id1", "id3", "id3"))
    expect_equal(as.data.frame(reg_view_factor)$value, c(3, 4, 1, 3, 4))
    expect_equal(nrow(reg), 5)
})


test_that("master index range subsetting supports string1:string2 in memory backend", {
    df <- data.frame(
        id = sprintf("id%03d", 1:10),
        value = seq_len(10),
        stringsAsFactors = FALSE
    )
    reg <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = Inf
    )

    string1 <- "id003"
    string2 <- "id007"
    out <- as.data.frame(reg[string1:string2, c("id", "value")])
    expect_equal(out$id, sprintf("id%03d", 3:7))
    expect_equal(out$value, 3:7)

    out_rev <- as.data.frame(reg[string2:string1, c("id", "value")])
    expect_equal(out_rev$id, sprintf("id%03d", 7:3))
    expect_equal(out_rev$value, 7:3)
})


test_that("row names use master index values even if master column is not active", {
    df <- data.frame(
        id = c("id1", "id2", "id3", "id4"),
        grp = c("a", "a", "b", "b"),
        value = c(11, 12, 13, 14),
        stringsAsFactors = FALSE
    )

    reg <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = Inf
    )

    reg_view <- reg[c("id3", "id1"), "value"]
    expect_equal(row.names(reg_view), c("id3", "id1"))
    expect_equal(rownames(reg_view), c("id3", "id1"))
    expect_equal(dimnames(reg_view)[[1]], c("id3", "id1"))
    out_df <- as.data.frame(reg_view)
    expect_equal(rownames(out_df), c("id3", "id1"))
    expect_equal(out_df$value, c(13, 11))
})


test_that("returned data.frames use master index as row names for memory and sqlite", {
    df <- data.frame(
        id = sprintf("id%03d", 1:20),
        grp = rep(c("g1", "g2"), each = 10),
        value = seq_len(20),
        stringsAsFactors = FALSE
    )

    reg_mem <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = Inf
    )

    q_mem <- reg_mem$queryOnIndex("by_id", c("id004", "id002"), columns = "value")
    expect_equal(rownames(q_mem), c("id004", "id002"))
    expect_equal(q_mem$value, c(4, 2))

    reg_sql <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = 0
    )
    on.exit(reg_sql$close(), add = TRUE)

    q_sql <- reg_sql$queryOnIndex("by_id", c("id011", "id003"), columns = "value")
    expect_equal(rownames(q_sql), c("id011", "id003"))
    expect_equal(q_sql$value, c(11, 3))
})


test_that("character row subsetting auto-detects a single master index and errors when ambiguous", {
    df <- data.frame(
        id = c("id1", "id2", "id3"),
        grp = c("g1", "g1", "g2"),
        value = c(10, 20, 30),
        stringsAsFactors = FALSE
    )

    reg_auto <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        memory_threshold_mb = Inf
    )
    reg_auto_view <- reg_auto[c("id3", "id1"), ]
    expect_equal(as.data.frame(reg_auto_view)$id, c("id3", "id1"))
    expect_equal(as.data.frame(reg_auto)$id, df$id)

    reg_amb <- getRegistry(
        data = df,
        indices = list(by_id = "id", by_grp = "grp"),
        memory_threshold_mb = Inf
    )
    expect_error(
        reg_amb[c("id1"), ],
        "Character row selectors require a configured `master_index`"
    )
})


test_that("character row subsetting by master index works in sqlite backend", {
    df <- data.frame(
        id = sprintf("id%03d", 1:100),
        grp = rep(c("g1", "g2"), each = 50),
        value = seq_len(100),
        stringsAsFactors = FALSE
    )

    reg <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = 0
    )
    on.exit(reg$close(), add = TRUE)

    reg_view <- reg[c("id090", "id001", "id090"), c("id", "value")]
    out <- as.data.frame(reg_view)
    expect_equal(out$id, c("id090", "id001", "id090"))
    expect_equal(out$value, c(90, 1, 90))

    reg_view_factor <- reg[factor(c("id090", "id001", "id090")), c("id", "value")]
    out_factor <- as.data.frame(reg_view_factor)
    expect_equal(out_factor$id, c("id090", "id001", "id090"))
    expect_equal(out_factor$value, c(90, 1, 90))

    string1 <- "id010"
    string2 <- "id015"
    out_range <- as.data.frame(reg[string1:string2, c("id", "value")])
    expect_equal(out_range$id, sprintf("id%03d", 10:15))
    expect_equal(out_range$value, 10:15)
})


test_that("logical mask subsetting works for rows and columns in memory and sqlite", {
    df <- data.frame(
        id = sprintf("id%03d", 1:12),
        grp = rep(c("g1", "g2", "g3"), each = 4),
        value = seq_len(12),
        stringsAsFactors = FALSE
    )

    row_mask <- df$grp != "g2"
    col_mask <- c(TRUE, FALSE, TRUE)
    expected <- df[row_mask, col_mask, drop = FALSE]
    rownames(expected) <- expected$id

    reg_mem <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = Inf
    )
    out_mem <- as.data.frame(reg_mem[row_mask, col_mask])
    expect_equal(out_mem, expected)

    reg_sql <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = 0
    )
    on.exit(reg_sql$close(), add = TRUE)
    out_sql <- as.data.frame(reg_sql[row_mask, col_mask])
    expect_equal(out_sql, expected)
})


test_that("saveRDS/readRDS preserves sqlite-backed registry even if original db file is gone", {
    df <- data.frame(
        id = sprintf("id%03d", 1:40),
        grp = rep(c("g1", "g2"), each = 20),
        value = seq_len(40),
        stringsAsFactors = FALSE
    )

    reg <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = 0
    )

    reg_view <- reg[c("id020", "id001"), c("id", "value")]
    expected <- as.data.frame(reg_view)

    db_path <- reg_view$.__enclos_env__$private$.sqlite_path
    rds_file <- tempfile(fileext = ".rds")
    on.exit(unlink(rds_file, force = TRUE), add = TRUE)

    saveRDS(reg_view, rds_file)
    reg_view$close()
    reg$close()

    # Simulate missing original backing file; restored object should use embedded DB raw bytes.
    if (!is.null(db_path) && file.exists(db_path)) {
        unlink(db_path, force = TRUE)
    }

    reg2 <- readRDS(rds_file)
    on.exit(reg2$close(), add = TRUE)

    out <- as.data.frame(reg2)
    expect_equal(out, expected)

    q <- reg2$queryOnIndex("by_id", c("id005", "id020"))
    expect_equal(q$id, c("id020"))
})


test_that("select restricts loaded columns for data.frame and sqlite backend", {
    df <- data.frame(
        id = sprintf("id%03d", 1:20),
        grp = rep(c("a", "b"), each = 10),
        value = seq_len(20),
        score = seq_len(20) * 10,
        stringsAsFactors = FALSE
    )

    reg_mem <- getRegistry(
        data = df,
        select = c("id", "value"),
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = Inf
    )
    expect_equal(reg_mem$colnames(), c("id", "value"))
    expected_select <- df[, c("id", "value"), drop = FALSE]
    rownames(expected_select) <- expected_select$id
    expect_equal(as.data.frame(reg_mem), expected_select)

    reg_sql <- getRegistry(
        data = df,
        select = c("id", "value"),
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = 0
    )
    on.exit(reg_sql$close(), add = TRUE)
    expect_equal(reg_sql$colnames(), c("id", "value"))
    expect_equal(as.data.frame(reg_sql), expected_select)
    q <- reg_sql$queryOnIndex("by_id", c("id003", "id001"))
    expect_equal(q$id, c("id003", "id001"))
})


test_that("select with TSV input works and index must exist in selected columns", {
    df <- data.frame(
        id = sprintf("id%03d", 1:30),
        grp = rep(c("g1", "g2", "g3"), each = 10),
        value = seq_len(30),
        extra = seq_len(30) * 2,
        stringsAsFactors = FALSE
    )
    tsv <- tempfile(fileext = ".tsv")
    on.exit(unlink(tsv, force = TRUE), add = TRUE)
    write.table(df, file = tsv, sep = "\t", quote = FALSE, row.names = FALSE)

    reg <- getRegistry(
        data = tsv,
        select = c("id", "value"),
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = 0
    )
    on.exit(reg$close(), add = TRUE)
    expect_equal(reg$colnames(), c("id", "value"))
    expect_equal(reg[c("id010", "id001"), ]$colnames(), c("id", "value"))

    expect_error(
        getRegistry(
            data = tsv,
            select = c("value", "extra"),
            indices = list(by_id = "id"),
            memory_threshold_mb = 0
        ),
        "references missing columns"
    )
})


test_that("gzipped TSV input is supported in both memory and sqlite backends", {
    df <- data.frame(
        id = sprintf("id%03d", 1:40),
        grp = rep(c("g1", "g2"), each = 20),
        value = seq_len(40),
        extra = seq_len(40) * 2,
        stringsAsFactors = FALSE
    )

    tsv <- tempfile(fileext = ".tsv")
    gz_tsv <- paste0(tsv, ".gz")
    on.exit({
        unlink(tsv, force = TRUE)
        unlink(gz_tsv, force = TRUE)
    }, add = TRUE)
    write.table(df, file = tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    con <- gzfile(gz_tsv, open = "wt")
    writeLines(readLines(tsv), con)
    close(con)

    reg_mem <- getRegistry(
        data = gz_tsv,
        select = c("id", "value"),
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = Inf
    )
    expected_gz <- df[, c("id", "value"), drop = FALSE]
    rownames(expected_gz) <- expected_gz$id
    expect_equal(as.data.frame(reg_mem), expected_gz)
    q_mem <- reg_mem$queryOnIndex("by_id", c("id020", "id001"))
    expect_equal(q_mem$id, c("id020", "id001"))

    reg_sql <- getRegistry(
        data = gz_tsv,
        select = c("id", "value"),
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = 0
    )
    on.exit(reg_sql$close(), add = TRUE)
    reg_sql_view <- reg_sql[c("id040", "id001"), ]
    out_sql <- as.data.frame(reg_sql_view)
    expect_equal(out_sql$id, c("id040", "id001"))
    expect_equal(out_sql$value, c(40, 1))
})


test_that("rename and derive are applied at ingestion", {
    df <- data.frame(
        id = sprintf("id%03d", 1:12),
        grp = rep(c("g1", "g2"), each = 6),
        a = seq_len(12),
        b = seq_len(12) * 2,
        stringsAsFactors = FALSE
    )

    reg_mem <- getRegistry(
        data = df,
        select = c("id", "grp", "a", "b"),
        rename = c(id = "site_id"),
        derive = list(
            prod_ab = list(cols = c("a", "b"), fun = function(a, b) a * b),
            a_plus_b = function(d) d$a + d$b
        ),
        indices = list(by_site = "site_id"),
        master_index = "by_site",
        memory_threshold_mb = Inf
    )

    out_mem <- as.data.frame(reg_mem)
    expect_equal(colnames(out_mem), c("site_id", "grp", "a", "b", "prod_ab", "a_plus_b"))
    expect_equal(out_mem$site_id, df$id)
    expect_equal(out_mem$prod_ab, df$a * df$b)
    expect_equal(out_mem$a_plus_b, df$a + df$b)

    tsv <- tempfile(fileext = ".tsv")
    on.exit(unlink(tsv, force = TRUE), add = TRUE)
    write.table(df, file = tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    reg_sql <- getRegistry(
        data = tsv,
        select = c("id", "grp", "a", "b"),
        rename = c(id = "site_id"),
        derive = list(prod_ab = list(cols = c("a", "b"), fun = function(a, b) a * b)),
        indices = list(by_site = "site_id"),
        master_index = "by_site",
        memory_threshold_mb = 0
    )
    on.exit(reg_sql$close(), add = TRUE)
    out_sql <- as.data.frame(reg_sql)
    expect_equal(colnames(out_sql), c("site_id", "grp", "a", "b", "prod_ab"))
    expect_equal(out_sql$prod_ab, df$a * df$b)
})


test_that("addIndex adds queryable index to existing object", {
    df <- data.frame(
        id = sprintf("id%03d", 1:20),
        grp = rep(c("g1", "g2"), each = 10),
        value = seq_len(20),
        stringsAsFactors = FALSE
    )

    reg_sql <- getRegistry(
        data = df,
        indices = list(),
        memory_threshold_mb = 0
    )
    on.exit(reg_sql$close(), add = TRUE)

    expect_error(reg_sql$queryOnIndex("by_grp", c("g2")), "Unknown index")
    reg_sql$addIndex("by_grp", "grp")
    q <- reg_sql$queryOnIndex("by_grp", c("g2", "g1"), columns = c("id", "grp"))
    expect_equal(q$grp, c(rep("g2", 10), rep("g1", 10)))
})


test_that("paste works on single-column Registry views via as.character coercion", {
    df <- data.frame(
        id = sprintf("id%03d", 1:12),
        grp = rep(c("g1", "g2"), each = 6),
        value = seq_len(12),
        stringsAsFactors = FALSE
    )
    reg <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = Inf
    )

    combined <- paste(reg[, "id"], ":", reg[, "grp"])
    expect_equal(combined, paste(df$id, ":", df$grp))
})


test_that("arithmetic operators coerce Registry views to data.frame first", {
    df <- data.frame(
        id = sprintf("id%03d", 1:6),
        value = c(2, 4, 6, 8, 10, 12),
        stringsAsFactors = FALSE
    )

    reg_mem <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = Inf
    )
    v_mem <- reg_mem[, "value"]
    plus_mem <- v_mem + 1
    minus_mem <- v_mem - 1
    times_mem <- v_mem * 3
    left_plus_mem <- 1 + v_mem
    unary_minus_mem <- -v_mem

    expect_equal(plus_mem$value, df$value + 1)
    expect_equal(minus_mem$value, df$value - 1)
    expect_equal(times_mem$value, df$value * 3)
    expect_equal(left_plus_mem$value, df$value + 1)
    expect_equal(unary_minus_mem$value, -df$value)
    expect_equal(rownames(plus_mem), df$id)

    reg_sql <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = 0
    )
    on.exit(reg_sql$close(), add = TRUE)
    v_sql <- reg_sql[, "value"]
    plus_sql <- v_sql + 1
    minus_sql <- v_sql - 1
    times_sql <- v_sql * 3
    left_plus_sql <- 1 + v_sql
    unary_minus_sql <- -v_sql

    expect_equal(plus_sql$value, df$value + 1)
    expect_equal(minus_sql$value, df$value - 1)
    expect_equal(times_sql$value, df$value * 3)
    expect_equal(left_plus_sql$value, df$value + 1)
    expect_equal(unary_minus_sql$value, -df$value)
    expect_equal(rownames(plus_sql), df$id)
})


test_that("[<-.Registry supports adding/replacing columns in memory and sqlite backends", {
    df <- data.frame(
        id = sprintf("id%03d", 1:20),
        grp = rep(c("g1", "g2"), each = 10),
        value = seq_len(20),
        stringsAsFactors = FALSE
    )

    reg_mem <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = Inf
    )
    reg_mem[, "joined"] <- paste(reg_mem[, "id"], ":", reg_mem[, "grp"])
    out_mem <- as.data.frame(reg_mem)
    expect_true("joined" %in% colnames(out_mem))
    expect_equal(out_mem$joined, paste(df$id, ":", df$grp))

    reg_sql <- getRegistry(
        data = df,
        indices = list(by_id = "id"),
        master_index = "by_id",
        memory_threshold_mb = 0
    )
    on.exit(reg_sql$close(), add = TRUE)
    reg_sql[, "joined"] <- paste(reg_sql[, "id"], ":", reg_sql[, "grp"])
    out_sql <- as.data.frame(reg_sql)
    expect_true("joined" %in% colnames(out_sql))
    expect_equal(out_sql$joined, paste(df$id, ":", df$grp))

    reg_sql[c("id003", "id001"), "value"] <- c(300, 100)
    updated <- as.data.frame(reg_sql[c("id001", "id003"), c("id", "value")])
    expect_equal(updated$value, c(100, 300))
})
