#' Create a registry for indexed tabular retrieval
#'
#' @description
#' Builds an index-aware registry from a TSV file or a `data.frame`.
#' Depending on the estimated input size, data is kept in memory or persisted to
#' SQLite. The registry keeps active row/column bounds. Subsetting returns a
#' new bounded `Registry` view applied to all subsequent queries.
#'
#' @param data Path to a TSV file (plain or gzipped, e.g. `.tsv.gz`) or a `data.frame`.
#' @param indices Named list of predefined index definitions. Each list element
#'   must be a character vector with one or more column names.
#' @param select Optional column selector applied at ingestion time. Can be
#'   `NULL`, a character vector of column names, or a positive integer vector of
#'   column positions.
#' @param rename Optional rename map applied at ingestion time. Accepts a named
#'   character vector or named list. Preferred form is `old = "new"` (old
#'   column names as names, new names as values). A `new = "old"` style is also
#'   accepted when unambiguous.
#' @param derive Optional named list defining derived columns to add at ingestion
#'   time. Each element name is the new column name. Each element value can be:
#'   - a function taking the current data.frame chunk and returning a vector of
#'     length `nrow(chunk)`, or
#'   - a list with `cols` (character vector of input columns) and `fun`
#'     (function), called as `do.call(fun, chunk[cols])`.
#' @param master_index Optional name of the index used for character-based row
#'   subsetting (first position in `[`).
#'   If `NULL` and exactly one single-column index is available, that index is
#'   used automatically.
#' @param memory_threshold_mb Numeric. Inputs with estimated size greater than
#'   this threshold are stored in SQLite.
#' @param chunk_size Integer. Number of rows per chunk for SQLite ingestion.
#' @param sqlite_path Optional SQLite file path. If `NULL`, a temporary file is used.
#' @param query_chunk_size Integer. Number of rows per chunk for temporary query
#'   tables (lookup keys and active rows) in SQLite mode.
#'
#' @return A `Registry` object.
#' @examples
#' df <- data.frame(
#'     id = c("a", "b", "c"),
#'     group = c("g1", "g1", "g2"),
#'     value = c(1, 2, 3),
#'     stringsAsFactors = FALSE
#' )
#'
#' reg <- getRegistry(
#'     data = df,
#'     indices = list(id = "id", by_group = "group")
#' )
#'
#' reg$queryOnIndex("id", c("c", "a"))
#' reg_view <- reg[1:2, c("id", "value")]
#' as.data.frame(reg_view)
#' @export
getRegistry <- function(data,
                        indices = list(),
                        select = NULL,
                        rename = NULL,
                        derive = NULL,
                        master_index = NULL,
                        memory_threshold_mb = getOption("DMRsegal.registry_in_mem_threshold_mb", 500),
                        chunk_size = getOption("DMRsegal.registry_chunk_size", 50000),
                        sqlite_path = NULL,
                        query_chunk_size = getOption("DMRsegal.registry_query_chunk_size", 50000)) {
    if (inherits(data, "Registry")) {
        return(invisible(data))
    }
    invisible(Registry$new(
        data = data,
        indices = indices,
        select = select,
        rename = rename,
        derive = derive,
        master_index = master_index,
        memory_threshold_mb = memory_threshold_mb,
        chunk_size = chunk_size,
        sqlite_path = sqlite_path,
        query_chunk_size = query_chunk_size
    ))
}


.registry_resolve_row_selector <- function(selector_expr, eval_env) {
    if (is.call(selector_expr) &&
        length(selector_expr) == 3L &&
        identical(selector_expr[[1L]], as.name(":"))) {
        lhs <- tryCatch(eval(selector_expr[[2L]], envir = eval_env), error = function(e) NULL)
        rhs <- tryCatch(eval(selector_expr[[3L]], envir = eval_env), error = function(e) NULL)

        lhs_chr <- if (is.factor(lhs)) as.character(lhs) else lhs
        rhs_chr <- if (is.factor(rhs)) as.character(rhs) else rhs

        lhs_is_chr <- is.character(lhs_chr)
        rhs_is_chr <- is.character(rhs_chr)
        if (lhs_is_chr || rhs_is_chr) {
            if (!(lhs_is_chr && rhs_is_chr)) {
                stop("Master-index range selectors cannot mix character and numeric boundaries.")
            }
            if (length(lhs_chr) != 1L || length(rhs_chr) != 1L) {
                stop("Master-index range selectors require scalar character boundaries.")
            }
            return(structure(
                list(start = lhs_chr[[1L]], end = rhs_chr[[1L]]),
                class = "registry_master_range"
            ))
        }

        if (is.numeric(lhs) && is.numeric(rhs) && length(lhs) == 1L && length(rhs) == 1L) {
            if (anyNA(c(lhs, rhs)) || !all(is.finite(c(lhs, rhs)))) {
                stop("Numeric range selector cannot contain NA/NaN/Inf values.")
            }
            return(seq.int(lhs[[1L]], rhs[[1L]]))
        }
    }
    eval(selector_expr, envir = eval_env)
}


#' @rdname getRegistry
#' @param x A `Registry` object.
#' @param i Row selector relative to the current active row bounds. Character
#'   selectors use `master_index`.
#' @param j Column selector relative to the current active column bounds.
#' @param drop Included for compatibility. Ignored.
#' @param ... Passed to `Registry$as.data.frame()`.
#' @export
`[.Registry` <- function(x, i, j, drop = FALSE) {
    i_sel <- if (missing(i)) NULL else .registry_resolve_row_selector(substitute(i), parent.frame())
    j_sel <- if (missing(j)) NULL else j
    .registry_subset_view(x = x, i = i_sel, j = j_sel, drop = drop)
}


#' @rdname getRegistry
#' @param value Replacement values for `[<-`. Must be a scalar or a vector of
#'   length matching the selected rows.
#' @export
`[<-.Registry` <- function(x, i, j, value) {
    i_sel <- if (missing(i)) NULL else .registry_resolve_row_selector(substitute(i), parent.frame())
    j_sel <- if (missing(j)) NULL else j
    if (is.null(j_sel)) {
        stop("Column selector `j` must be provided for `[<-.Registry`.")
    }

    private <- x$.__enclos_env__$private
    target_rows <- private$.active_row_idx
    if (!is.null(i_sel)) {
        if (inherits(i_sel, "registry_master_range")) {
            target_rows <- private$.subset_rows_by_master_index_range(i_sel$start, i_sel$end)
        } else if (is.character(i_sel) || is.factor(i_sel)) {
            target_rows <- private$.subset_rows_by_master_index(as.character(i_sel))
        } else {
            target_rows <- private$.subset_absolute_indices(
                current_idx = private$.active_row_idx,
                selector = i_sel,
                axis_names = NULL,
                axis = "rows"
            )
        }
    }

    active_col_names <- private$.base_col_names[private$.active_col_idx]
    is_new_col <- FALSE
    if (is.numeric(j_sel)) {
        if (length(j_sel) != 1L || is.na(j_sel)) {
            stop("Numeric column selector for `[<-` must be a single non-NA value.")
        }
        mapped <- seq_along(active_col_names)[j_sel]
        if (length(mapped) != 1L || is.na(mapped)) {
            stop("Numeric column selector for `[<-` is out of active column bounds.")
        }
        target_col <- active_col_names[[mapped]]
    } else if (is.character(j_sel)) {
        if (length(j_sel) != 1L || is.na(j_sel) || j_sel == "") {
            stop("Character column selector for `[<-` must be a single non-empty string.")
        }
        target_col <- j_sel
        is_new_col <- !(target_col %in% private$.base_col_names)
    } else {
        stop("Unsupported column selector type for `[<-.Registry`.")
    }

    if (inherits(value, "Registry")) {
        value_df <- as.data.frame(value)
        if (ncol(value_df) != 1L) {
            stop("`value` Registry for `[<-` must have exactly one active column.")
        }
        value <- value_df[[1]]
    } else if (is.data.frame(value)) {
        if (ncol(value) != 1L) {
            stop("`value` data.frame for `[<-` must have exactly one column.")
        }
        value <- value[[1]]
    } else if (is.matrix(value)) {
        if (ncol(value) != 1L) {
            stop("`value` matrix for `[<-` must have exactly one column.")
        }
        value <- as.vector(value[, 1])
    }

    n_target <- length(target_rows)
    if (n_target > 0L) {
        if (length(value) == 1L) {
            value <- rep.int(value, n_target)
        }
        if (length(value) != n_target) {
            stop(
                "Replacement length mismatch in `[<-.Registry`: got ",
                length(value),
                ", expected ",
                n_target,
                "."
            )
        }
    }

    if (n_target > 0L) {
        keep <- !duplicated(target_rows, fromLast = TRUE)
        target_rows <- as.integer(target_rows[keep])
        value <- value[keep]
    } else {
        target_rows <- integer(0)
    }

    if (private$.backend == "memory") {
        base_df <- private$.data_in_memory
        if (is_new_col) {
            base_df[[target_col]] <- rep.int(NA, nrow(base_df))
            private$.base_col_names <- c(private$.base_col_names, target_col)
            private$.active_col_idx <- c(private$.active_col_idx, length(private$.base_col_names))
        }
        if (length(target_rows) > 0L) {
            base_df[target_rows, target_col] <- value
        }
        private$.data_in_memory <- base_df
        private$.memory_index_cache <- list()
        return(x)
    }

    private$.ensure_sqlite_connection()
    conn <- private$.conn
    table_q <- .registry_quote_ident(conn, private$.table_name)
    col_q <- .registry_quote_ident(conn, target_col)

    if (is_new_col) {
        DBI::dbExecute(
            conn,
            paste0("ALTER TABLE ", table_q, " ADD COLUMN ", col_q)
        )
        private$.base_col_names <- c(private$.base_col_names, target_col)
        private$.active_col_idx <- c(private$.active_col_idx, length(private$.base_col_names))
    }

        if (length(target_rows) > 0L) {
            tmp_assign <- .registry_new_temp_table("tmp_registry_assign")
            tmp_assign_q <- .registry_quote_ident(conn, tmp_assign)
            on.exit({
                if (DBI::dbExistsTable(conn, tmp_assign)) {
                DBI::dbExecute(conn, paste0("DROP TABLE ", tmp_assign_q))
            }
        }, add = TRUE)

            template <- data.frame(
                .registry_rowid = integer(0),
                .value = value[0],
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
        DBI::dbWriteTable(conn, tmp_assign, template, temporary = TRUE, row.names = FALSE)
        for (start in seq.int(1L, length(target_rows), by = private$.query_chunk_size)) {
            end <- min(start + private$.query_chunk_size - 1L, length(target_rows))
            chunk <- data.frame(
                .registry_rowid = target_rows[start:end],
                .value = value[start:end],
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
            DBI::dbWriteTable(conn, tmp_assign, chunk, append = TRUE, row.names = FALSE)
        }

        sql <- paste0(
            "UPDATE ", table_q,
            " SET ", col_q, " = (",
            "SELECT t.\".value\" FROM ", tmp_assign_q, " t ",
            "WHERE t.\".registry_rowid\" = ", table_q, ".\".registry_rowid\"",
            ") WHERE \".registry_rowid\" IN (SELECT \".registry_rowid\" FROM ", tmp_assign_q, ")"
        )
        DBI::dbExecute(conn, sql)
    }
    private$.snapshot_sqlite_file()
    x
}


.registry_subset_view <- function(x, i = NULL, j = NULL, drop = FALSE) {
    private <- x$.__enclos_env__$private
    new_row_idx <- private$.active_row_idx
    new_col_idx <- private$.active_col_idx

    if (!is.null(i)) {
        if (inherits(i, "registry_master_range")) {
            new_row_idx <- private$.subset_rows_by_master_index_range(i$start, i$end)
        } else if (is.character(i) || is.factor(i)) {
            new_row_idx <- private$.subset_rows_by_master_index(as.character(i))
        } else {
            new_row_idx <- private$.subset_absolute_indices(
                current_idx = new_row_idx,
                selector = i,
                axis_names = NULL,
                axis = "rows"
            )
        }
    }
    if (!is.null(j)) {
        active_col_names <- private$.base_col_names[new_col_idx]
        new_col_idx <- private$.subset_absolute_indices(
            current_idx = new_col_idx,
            selector = j,
            axis_names = active_col_names,
            axis = "columns"
        )
    }
    private$.spawn_view(new_row_idx = new_row_idx, new_col_idx = new_col_idx)
}


#' @rdname getRegistry
#' @method as.data.frame Registry
#' @rawNamespace S3method(as.data.frame,Registry)
as.data.frame.Registry <- function(x, ...) {
    x$as.data.frame(...)
}


#' @rdname getRegistry
#' @method as.character Registry
#' @rawNamespace S3method(as.character,Registry)
as.character.Registry <- function(x, ...) {
    df <- x$as.data.frame(...)
    if (ncol(df) == 0L) {
        return(rep.int("", nrow(df)))
    }
    if (ncol(df) == 1L) {
        return(as.character(df[[1L]]))
    }
    apply(df, 1L, function(row) paste(as.character(row), collapse = "\t"))
}

#' @rdname getRegistry
#' @method as.integer Registry
#' @rawNamespace S3method(as.integer,Registry)
as.integer.Registry <- function(x, ...) {
    df <- x$as.data.frame(...)
    if (ncol(df) == 1L) {
        return(as.integer(df[[1L]]))
    }
    as.integer(df)
}


#' @rdname getRegistry
#' @method Ops Registry
#' @rawNamespace S3method(Ops,Registry)
Ops.Registry <- function(e1, e2 = NULL) {
    op_fun <- get(.Generic, envir = baseenv(), mode = "function")

    if (nargs() == 1L) {
        lhs <- if (inherits(e1, "Registry")) as.data.frame(e1) else e1
        return(op_fun(lhs))
    }

    lhs <- if (inherits(e1, "Registry")) as.data.frame(e1) else e1
    rhs <- if (inherits(e2, "Registry")) as.data.frame(e2) else e2
    op_fun(lhs, rhs)
}


#' @rdname getRegistry
#' @method dim Registry
#' @rawNamespace S3method(dim,Registry)
dim.Registry <- function(x) {
    bounds <- x$getActiveBounds()
    c(length(bounds$row_idx), length(bounds$col_idx))
}


#' @rdname getRegistry
#' @method names Registry
#' @rawNamespace S3method(names,Registry)
names.Registry <- function(x) {
    x$getActiveBounds()$col_names
}

#' @rdname getRegistry
#' @method head Registry
#' @rawNamespace S3method(head,Registry)
head.Registry <- function(x, n = 6L, ...) {
    bounds <- x$getActiveBounds()
    head_idx <- head(bounds$row_idx, n)
    x$subsetRowsCols(i = head_idx)
}


#' @rdname getRegistry
#' @method dimnames Registry
#' @rawNamespace S3method(dimnames,Registry)
dimnames.Registry <- function(x) {
    list(row.names(x), x$getActiveBounds()$col_names)
}


#' @rdname getRegistry
#' @method colnames Registry
#' @rawNamespace S3method(colnames,Registry)
colnames.Registry <- function(x, do.NULL = TRUE, prefix = "col") { # nolint
    x$getActiveBounds()$col_names
}

#' @rdname getRegistry
#' @method row.names Registry
#' @description The row names of a `Registry` correspond to the master index
#'   values for the currently active row bounds. If no master index is defined,
#'   row names are generated as `row1`, `row2`, etc.
#' @rawNamespace S3method(row.names,Registry)
row.names.Registry <- function(x, do.NULL = TRUE, prefix = "row") { # nolint
    bounds <- x$getActiveBounds()
    if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix)) {
        stop("`prefix` must be a single non-missing character string.")
    }
    if (length(bounds$row_idx) == 0L) {
        return(character(0))
    }

    master_idx <- x$getMasterIndex()
    if (is.null(master_idx)) {
        return(paste0(prefix, bounds$row_idx))
    }

    private <- NULL
    if (!is.null(x$.__enclos_env__) && !is.null(x$.__enclos_env__$private)) {
        private <- x$.__enclos_env__$private
    }
    if (is.null(private) || is.null(private$.indices)) {
        return(paste0(prefix, bounds$row_idx))
    }

    idx_cols <- private$.indices[[master_idx]]
    if (!is.character(idx_cols) || length(idx_cols) != 1L) {
        return(paste0(prefix, bounds$row_idx))
    }
    idx_col <- idx_cols[[1]]

    # Prefer private bounded fetch; this works even when the master-index
    # column is not part of the active column bounds.
    fetch_fun <- private$.fetch_rows_with_columns
    if (is.function(fetch_fun)) {
        vals <- fetch_fun(bounds$row_idx, idx_col)[[idx_col]]
        return(as.character(vals))
    }

    if (idx_col %in% bounds$col_names) {
        return(as.character(x[[idx_col]]))
    }

    paste0(prefix, bounds$row_idx)
}


#' @rdname getRegistry
#' @method [[ Registry
#' @rawNamespace S3method("[[",Registry)
`[[.Registry` <- function(x, i, ...) {
    if (length(i) != 1) {
        stop("`[[` expects a single column selector.")
    }

    if (is.character(i)) {
        if (exists(i, envir = x, inherits = FALSE)) {
            return(.subset2(x, i))
        }
        col_names <- x$getActiveBounds()$col_names
        if (!(i %in% col_names)) {
            stop("Unknown column `", i, "` in active bounds.")
        }
        private <- x$.__enclos_env__$private
        vals <- private$.fetch_rows_with_columns(private$.active_row_idx, i)
        return(vals[[i]])
    }

    if (is.numeric(i)) {
        col_names <- x$getActiveBounds()$col_names
        if (is.na(i) || i < 1 || i > length(col_names) || as.integer(i) != i) {
            stop("Numeric `[[` selector must be a valid single column position.")
        }
        col <- col_names[[as.integer(i)]]
        private <- x$.__enclos_env__$private
        vals <- private$.fetch_rows_with_columns(private$.active_row_idx, col)
        return(vals[[col]])
    }

    stop("Unsupported selector type for `[[`.")
}


#' @rdname getRegistry
#' @method $ Registry
#' @rawNamespace S3method("$",Registry)
`$.Registry` <- function(x, name) {
    if (exists(name, envir = x, inherits = FALSE)) {
        return(.subset2(x, name))
    }
    col_names <- x$getActiveBounds()$col_names
    if (name %in% col_names) {
        return(x[[name]])
    }
    NULL
}


.registry_escape_string <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("\\\\", "\\\\\\\\", x)
    x <- gsub("\"", "\\\"", x, fixed = TRUE)
    x <- gsub("\n", "\\n", x, fixed = TRUE)
    x <- gsub("\r", "\\r", x, fixed = TRUE)
    x
}


.registry_encode_single_key <- function(x) {
    vals <- as.character(x)
    is_na <- is.na(vals)
    vals[is_na] <- ""
    paste0(ifelse(is_na, "N", "V"), vals)
}


.registry_make_composite_key <- function(df, cols) {
    if (length(cols) == 0) {
        return(rep.int("", nrow(df)))
    }
    parts <- lapply(cols, function(col) {
        vals <- as.character(df[[col]])
        is_na <- is.na(vals)
        vals[is_na] <- ""
        vals <- .registry_escape_string(vals)
        paste0(ifelse(is_na, "N", "V"), nchar(vals), ":", vals, ";")
    })
    Reduce(paste0, parts)
}


.registry_make_index_keys <- function(df, cols) {
    if (length(cols) == 1L) {
        return(.registry_encode_single_key(df[[cols[[1L]]]]))
    }
    .registry_make_composite_key(df, cols)
}


.registry_quote_ident <- function(conn, x) {
    as.character(DBI::dbQuoteIdentifier(conn, x))
}


.registry_validate_index_spec <- function(indices, col_names) {
    if (is.character(indices)) {
        if (length(indices) != 1L || is.na(indices) || indices == "") {
            stop("Character `indices` input must be a single non-empty column name.")
        }
        idx_name <- names(indices)
        if (is.null(idx_name) || length(idx_name) != 1L || is.na(idx_name) || idx_name == "") {
            idx_name <- indices
        }
        indices <- setNames(list(indices), idx_name)
    }

    if (!is.list(indices)) {
        stop("`indices` must be a named list of character vectors or a single character column name.")
    }
    if (length(indices) == 0) {
        return(list())
    }
    idx_names <- names(indices)
    if (is.null(idx_names) || anyNA(idx_names) || any(idx_names == "")) {
        stop("`indices` must be a named list with non-empty names.")
    }
    for (idx_name in idx_names) {
        cols <- indices[[idx_name]]
        if (!is.character(cols) || length(cols) == 0) {
            stop("Index `", idx_name, "` must be a non-empty character vector.")
        }
        if (anyNA(cols) || any(cols == "")) {
            stop("Index `", idx_name, "` contains invalid column names.")
        }
        missing_cols <- setdiff(cols, col_names)
        if (length(missing_cols) > 0) {
            stop(
                "Index `", idx_name, "` references missing columns: ",
                paste(missing_cols, collapse = ", ")
            )
        }
    }
    indices
}


.registry_new_temp_table <- function(prefix) {
    ts <- gsub("[^0-9]", "", format(Sys.time(), "%Y%m%d%H%M%OS6"))
    paste0(prefix, "_", ts, "_", sample.int(1000000, 1))
}


.registry_open_read_connection <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        return(gzfile(path, open = "rt"))
    }
    file(path, open = "rt")
}


.registry_seq_bounds <- function(x) {
    n <- length(x)
    if (n == 0L || anyNA(x)) {
        return(NULL)
    }
    x_int <- as.integer(x)
    if (n == 1L) {
        return(c(x_int[[1L]], x_int[[1L]]))
    }
    if (all(diff(x_int) == 1L)) {
        return(c(x_int[[1L]], x_int[[n]]))
    }
    NULL
}


.registry_safe_row_names <- function(x, prefix = "row") {
    out <- as.character(x)
    bad <- is.na(out) | out == ""
    if (any(bad)) {
        out[bad] <- prefix
    }
    make.unique(out, sep = ".")
}


#' Registry Class
#'
#' @description Internal R6 class implementing indexed tabular retrieval with
#' in-memory and SQLite backends.
#' @keywords internal
Registry <- R6::R6Class("Registry", # nolint
    public = list(
        initialize = function(data,
                              indices = list(),
                              select = NULL,
                              rename = NULL,
                              derive = NULL,
                              master_index = NULL,
                              memory_threshold_mb = getOption("DMRsegal.registry_in_mem_threshold_mb", 500),
                              chunk_size = getOption("DMRsegal.registry_chunk_size", 50000),
                              sqlite_path = NULL,
                              query_chunk_size = getOption("DMRsegal.registry_query_chunk_size", 50000)) {
            init_ok <- FALSE
            on.exit({
                if (!init_ok) {
                    suppressWarnings(try(self$close(), silent = TRUE))
                }
            }, add = TRUE)

            if (!is.numeric(memory_threshold_mb) || length(memory_threshold_mb) != 1 || is.na(memory_threshold_mb) || memory_threshold_mb < 0) {
                stop("`memory_threshold_mb` must be a single non-negative number.")
            }
            if (!is.numeric(chunk_size) || length(chunk_size) != 1 || chunk_size < 1) {
                stop("`chunk_size` must be a positive integer.")
            }
            if (!is.numeric(query_chunk_size) || length(query_chunk_size) != 1 || query_chunk_size < 1) {
                stop("`query_chunk_size` must be a positive integer.")
            }

            private$.chunk_size <- as.integer(chunk_size)
            private$.query_chunk_size <- as.integer(query_chunk_size)
            private$.load_data(
                data = data,
                memory_threshold_mb = as.numeric(memory_threshold_mb),
                sqlite_path = sqlite_path,
                select = select,
                rename = rename,
                derive = derive
            )
            private$.indices <- .registry_validate_index_spec(indices, private$.base_col_names)
            for (idx_name in names(private$.indices)) {
                private$.create_index(idx_name, private$.indices[[idx_name]])
            }
            private$.master_index <- private$.resolve_master_index(master_index)
            private$.snapshot_sqlite_file()
            init_ok <- TRUE
        },

        addIndex = function(index_name, columns) {
            if (!is.character(index_name) || length(index_name) != 1 || is.na(index_name) || index_name == "") {
                stop("`index_name` must be a non-empty character scalar.")
            }
            validated <- .registry_validate_index_spec(
                setNames(list(columns), index_name),
                private$.base_col_names
            )
            private$.indices[[index_name]] <- validated[[index_name]]
            private$.create_index(index_name, private$.indices[[index_name]])
            private$.snapshot_sqlite_file()
            invisible(self)
        },

        setMasterIndex = function(master_index) {
            private$.master_index <- private$.resolve_master_index(master_index)
            private$.snapshot_sqlite_file()
            invisible(self)
        },

        getMasterIndex = function() {
            private$.master_index
        },

        queryOnIndex = function(index_name, values, columns = NULL) {
            if (!is.character(index_name) || length(index_name) != 1 || is.na(index_name) || index_name == "") {
                stop("`index_name` must be a non-empty character scalar.")
            }
            if (is.null(private$.indices[[index_name]])) {
                stop("Unknown index `", index_name, "`. Available indices: ", paste(names(private$.indices), collapse = ", "))
            }

            index_cols <- private$.indices[[index_name]]
            lookup_df <- private$.normalize_lookup_values(values, index_cols)
            if (nrow(lookup_df) == 0) {
                return(private$.empty_result_df(private$.resolve_output_columns(columns)))
            }

            out_cols <- private$.resolve_output_columns(columns)
            if (private$.backend == "memory") {
                return(private$.query_memory(index_name, index_cols, lookup_df, out_cols))
            }
            private$.query_sqlite(index_cols, lookup_df, out_cols)
        },

        subsetRowsCols = function(i = NULL, j = NULL, drop = FALSE) {
            .registry_subset_view(x = self, i = i, j = j, drop = drop)
        },

        as.data.frame = function(...) {
            selected_cols <- private$.base_col_names[private$.active_col_idx]
            out <- private$.fetch_rows_with_columns(private$.active_row_idx, selected_cols)
            private$.apply_master_row_names(out, row_ids = private$.active_row_idx)
        },

        resetBounds = function() {
            private$.active_row_idx <- seq_len(private$.n_rows)
            private$.active_col_idx <- seq_along(private$.base_col_names)
            invisible(self)
        },

        getActiveBounds = function() {
            list(
                row_idx = private$.active_row_idx,
                col_idx = private$.active_col_idx,
                col_names = private$.base_col_names[private$.active_col_idx]
            )
        },

        getRowNames = function(prefix = "row") {
            private$.active_row_names(prefix = prefix)
        },

        backend = function() {
            private$.backend
        },

        nrow = function() {
            length(private$.active_row_idx)
        },

        ncol = function() {
            length(private$.active_col_idx)
        },

        colnames = function() {
            private$.base_col_names[private$.active_col_idx]
        },

        close = function() {
            if (!is.null(private$.conn)) {
                is_valid <- isTRUE(tryCatch(DBI::dbIsValid(private$.conn), error = function(e) FALSE))
                if (is_valid && !isTRUE(private$.shared_conn)) {
                    suppressWarnings(try(DBI::dbDisconnect(private$.conn), silent = TRUE))
                }
                private$.conn <- NULL
                private$.shared_conn <- FALSE
            }
            if (isTRUE(private$.owns_db_file) && !is.null(private$.sqlite_path) && file.exists(private$.sqlite_path)) {
                suppressWarnings(unlink(private$.sqlite_path, force = TRUE))
            }
            invisible(self)
        }
    ),
    private = list(
        .backend = NULL,
        .table_name = "registry_data",
        .data_in_memory = NULL,
        .conn = NULL,
        .shared_conn = FALSE,
        .owns_db_file = FALSE,
        .sqlite_path = NULL,
        .n_rows = 0L,
        .base_col_names = character(0),
        .indices = list(),
        .memory_index_cache = list(),
        .master_index = NULL,
        .sqlite_file_raw = NULL,
        .active_row_idx = integer(0),
        .active_col_idx = integer(0),
        .chunk_size = 50000L,
        .query_chunk_size = 50000L,
        finalize = function() {
            try(self$close(), silent = TRUE)
            invisible()
        },

        .estimate_input_size_mb = function(data) {
            if (is.character(data) && length(data) == 1 && file.exists(data)) {
                return(as.numeric(file.info(data)$size) / (1024^2))
            }
            if (is.data.frame(data)) {
                return(as.numeric(object.size(data)) / (1024^2))
            }
            stop("`data` must be either an existing TSV path or a data.frame.")
        },

        .load_data = function(data, memory_threshold_mb, sqlite_path, select = NULL, rename = NULL, derive = NULL) {
            size_mb <- private$.estimate_input_size_mb(data)
            use_memory <- size_mb <= memory_threshold_mb

            if (use_memory) {
                private$.backend <- "memory"
                private$.load_to_memory(data, select = select, rename = rename, derive = derive)
            } else {
                private$.backend <- "sqlite"
                private$.load_to_sqlite(data, sqlite_path, select = select, rename = rename, derive = derive)
            }
            private$.active_row_idx <- seq_len(private$.n_rows)
            private$.active_col_idx <- seq_along(private$.base_col_names)
        },

        .load_to_memory = function(data, select = NULL, rename = NULL, derive = NULL) {
            df <- private$.coerce_input_to_df(data, select = select, rename = rename, derive = derive)
            if (ncol(df) == 0) {
                stop("Input must contain at least one column.")
            }
            df$.registry_rowid <- seq_len(nrow(df))
            private$.data_in_memory <- df
            private$.n_rows <- nrow(df)
            private$.base_col_names <- setdiff(colnames(df), ".registry_rowid")
        },

        .load_to_sqlite = function(data, sqlite_path, select = NULL, rename = NULL, derive = NULL) {
            db_file <- sqlite_path
            owns_file <- FALSE
            if (is.null(db_file)) {
                db_file <- tempfile(pattern = "dmrsegal_registry_", fileext = ".sqlite")
                owns_file <- TRUE
            }
            private$.conn <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_file)
            private$.shared_conn <- FALSE
            private$.configure_sqlite_connection()
            private$.sqlite_path <- db_file
            private$.owns_db_file <- owns_file

            if (is.data.frame(data)) {
                private$.ingest_df_to_sqlite(data, select = select, rename = rename, derive = derive)
            } else if (is.character(data) && length(data) == 1 && file.exists(data)) {
                private$.ingest_tsv_to_sqlite(data, select = select, rename = rename, derive = derive)
            } else {
                stop("`data` must be either an existing TSV path or a data.frame.")
            }

            private$.create_base_rowid_index()
        },

        .coerce_input_to_df = function(data, select = NULL, rename = NULL, derive = NULL) {
            if (is.data.frame(data)) {
                df <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
                spec <- private$.build_ingest_spec(
                    available_col_names = colnames(df),
                    select = select,
                    rename = rename,
                    derive = derive
                )
                df <- df[, spec$selected_source_cols, drop = FALSE]
                return(private$.apply_ingest_spec(df, spec))
            }
            if (is.character(data) && length(data) == 1 && file.exists(data)) {
                spec <- private$.build_ingest_spec_from_path(
                    path = data,
                    select = select,
                    rename = rename,
                    derive = derive
                )
                df <- data.table::fread(
                    file = data,
                    select = spec$selected_source_cols,
                    data.table = FALSE,
                    showProgress = getOption("DMRsegal.verbose", 0) > 1
                )
                df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
                if (!identical(colnames(df), spec$selected_source_cols)) {
                    colnames(df) <- spec$selected_source_cols
                }
                return(private$.apply_ingest_spec(df, spec))
            }
            stop("`data` must be either an existing TSV path or a data.frame.")
        },

        .ingest_df_to_sqlite = function(df, select = NULL, rename = NULL, derive = NULL) {
            df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
            spec <- private$.build_ingest_spec(
                available_col_names = colnames(df),
                select = select,
                rename = rename,
                derive = derive
            )
            if (length(spec$final_cols) == 0) {
                stop("Input must contain at least one column.")
            }
            private$.base_col_names <- spec$final_cols
            n <- nrow(df)
            next_rowid <- 1L
            table_exists <- FALSE
            if (n == 0) {
                empty <- df[0, spec$selected_source_cols, drop = FALSE]
                empty <- private$.apply_ingest_spec(empty, spec)
                empty$.registry_rowid <- integer(0)
                empty <- empty[, c(".registry_rowid", private$.base_col_names), drop = FALSE]
                DBI::dbWriteTable(private$.conn, private$.table_name, empty, row.names = FALSE)
                private$.n_rows <- 0L
                return(invisible())
            }
            for (start in seq.int(1L, n, by = private$.chunk_size)) {
                end <- min(start + private$.chunk_size - 1L, n)
                chunk <- df[start:end, spec$selected_source_cols, drop = FALSE]
                chunk <- private$.apply_ingest_spec(chunk, spec)
                chunk$.registry_rowid <- seq.int(next_rowid, length.out = nrow(chunk))
                next_rowid <- next_rowid + nrow(chunk)
                chunk <- chunk[, c(".registry_rowid", private$.base_col_names), drop = FALSE]
                DBI::dbWriteTable(
                    private$.conn,
                    private$.table_name,
                    chunk,
                    append = table_exists,
                    row.names = FALSE
                )
                table_exists <- TRUE
            }
            private$.n_rows <- n
        },

        .ingest_tsv_to_sqlite = function(path, select = NULL, rename = NULL, derive = NULL) {
            con <- .registry_open_read_connection(path)
            on.exit(close(con), add = TRUE)

            header_line <- readLines(con, n = 1, warn = FALSE)
            if (length(header_line) == 0) {
                stop("Input TSV file is empty: ", path)
            }
            file_col_names <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
            spec <- private$.build_ingest_spec(
                available_col_names = file_col_names,
                select = select,
                rename = rename,
                derive = derive
            )
            private$.base_col_names <- spec$final_cols
            if (length(private$.base_col_names) == 0) {
                stop("Input TSV file does not contain valid column names.")
            }

            rows_written <- 0L
            table_exists <- FALSE
            repeat {
                lines <- readLines(con, n = private$.chunk_size, warn = FALSE)
                if (length(lines) == 0) {
                    break
                }
                chunk <- data.table::fread(
                    text = paste(c(header_line, lines), collapse = "\n"),
                    sep = "\t",
                    header = TRUE,
                    select = spec$selected_source_cols,
                    data.table = FALSE,
                    showProgress = FALSE
                )
                chunk <- as.data.frame(chunk, stringsAsFactors = FALSE, check.names = FALSE)
                if (!identical(colnames(chunk), spec$selected_source_cols)) {
                    colnames(chunk) <- spec$selected_source_cols
                }
                chunk <- private$.apply_ingest_spec(chunk, spec)
                chunk$.registry_rowid <- seq.int(rows_written + 1L, length.out = nrow(chunk))
                chunk <- chunk[, c(".registry_rowid", private$.base_col_names), drop = FALSE]
                DBI::dbWriteTable(
                    private$.conn,
                    private$.table_name,
                    chunk,
                    append = table_exists,
                    row.names = FALSE
                )
                table_exists <- TRUE
                rows_written <- rows_written + nrow(chunk)
            }

            if (!table_exists) {
                empty <- as.data.frame(setNames(vector("list", length(private$.base_col_names)), private$.base_col_names), stringsAsFactors = FALSE)
                empty$.registry_rowid <- integer(0)
                empty <- empty[, c(".registry_rowid", private$.base_col_names), drop = FALSE]
                DBI::dbWriteTable(private$.conn, private$.table_name, empty, row.names = FALSE)
            }
            private$.n_rows <- as.integer(rows_written)
        },

        .build_ingest_spec_from_path = function(path, select = NULL, rename = NULL, derive = NULL) {
            con <- .registry_open_read_connection(path)
            on.exit(close(con), add = TRUE)
            header_line <- readLines(con, n = 1, warn = FALSE)
            if (length(header_line) == 0) {
                stop("Input TSV file is empty: ", path)
            }
            col_names <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
            private$.build_ingest_spec(
                available_col_names = col_names,
                select = select,
                rename = rename,
                derive = derive
            )
        },

        .build_ingest_spec = function(available_col_names, select = NULL, rename = NULL, derive = NULL) {
            selected_source_cols <- unname(private$.resolve_selected_columns(select, available_col_names))
            rename_map <- private$.normalize_rename_map(rename, selected_source_cols)
            renamed_cols <- unname(vapply(selected_source_cols, function(col) {
                if (col %in% names(rename_map)) rename_map[[col]] else col
            }, character(1)))
            if (anyDuplicated(renamed_cols)) {
                stop("Renaming produced duplicate output column names.")
            }
            derive_spec <- private$.normalize_derive_spec(derive, renamed_cols)
            final_cols <- unname(c(renamed_cols, names(derive_spec)))
            if (anyDuplicated(final_cols)) {
                stop("Final output columns contain duplicates.")
            }
            list(
                selected_source_cols = selected_source_cols,
                rename_map = rename_map,
                renamed_cols = renamed_cols,
                derive = derive_spec,
                final_cols = final_cols
            )
        },

        .apply_ingest_spec = function(df, spec) {
            if (!identical(colnames(df), spec$selected_source_cols)) {
                colnames(df) <- spec$selected_source_cols
            }

            colnames(df) <- spec$renamed_cols

            if (length(spec$derive) > 0) {
                for (new_col in names(spec$derive)) {
                    ds <- spec$derive[[new_col]]
                    if (identical(ds$mode, "df")) {
                        val <- ds$fun(df)
                    } else {
                        args <- lapply(ds$cols, function(cc) df[[cc]])
                        val <- do.call(ds$fun, args)
                    }
                    if (is.data.frame(val)) {
                        if (ncol(val) != 1L) {
                            stop("Derived column function for `", new_col, "` must return a vector or single-column data.frame.")
                        }
                        val <- val[[1]]
                    }
                    if (length(val) != nrow(df)) {
                        stop(
                            "Derived column `", new_col, "` has invalid length: ",
                            length(val), " (expected ", nrow(df), ")."
                        )
                    }
                    df[[new_col]] <- val
                }
            }

            df[, spec$final_cols, drop = FALSE]
        },

        .normalize_rename_map = function(rename = NULL, selected_source_cols) {
            if (is.null(rename)) {
                return(setNames(character(0), character(0)))
            }

            if (is.list(rename) && !is.character(rename)) {
                rename <- unlist(rename, use.names = TRUE)
            }
            if (!is.character(rename) || is.null(names(rename))) {
                stop("`rename` must be a named character vector or named list.")
            }
            if (anyNA(names(rename)) || any(names(rename) == "") || anyNA(rename) || any(rename == "")) {
                stop("`rename` contains invalid empty/NA names or values.")
            }

            names_match <- all(names(rename) %in% selected_source_cols)
            values_match <- all(as.character(rename) %in% selected_source_cols)

            map_from_names <- setNames(as.character(rename), names(rename))
            map_from_values <- setNames(names(rename), as.character(rename))

            rename_map <- NULL
            if (names_match && !values_match) {
                rename_map <- map_from_names
            } else if (!names_match && values_match) {
                rename_map <- map_from_values
            } else if (names_match && values_match) {
                if (identical(map_from_names, map_from_values)) {
                    rename_map <- map_from_names
                } else {
                    stop(
                        "`rename` is ambiguous (both names and values match source columns). ",
                        "Use unambiguous mapping."
                    )
                }
            } else {
                stop("`rename` does not reference selected source columns.")
            }

            if (anyDuplicated(unname(rename_map))) {
                stop("`rename` maps multiple source columns to the same target name.")
            }
            unchanged <- setdiff(selected_source_cols, names(rename_map))
            if (any(unname(rename_map) %in% unchanged)) {
                stop("`rename` target names collide with unchanged selected columns.")
            }
            rename_map
        },

        .normalize_derive_spec = function(derive = NULL, available_after_rename) {
            if (is.null(derive)) {
                return(list())
            }
            if (!is.list(derive) || is.null(names(derive)) || anyNA(names(derive)) || any(names(derive) == "")) {
                stop("`derive` must be a named list.")
            }
            if (anyDuplicated(names(derive))) {
                stop("`derive` contains duplicated output column names.")
            }
            if (any(names(derive) %in% available_after_rename)) {
                stop("`derive` output names must be new columns and cannot overwrite existing columns.")
            }

            normalized <- list()
            for (new_col in names(derive)) {
                spec <- derive[[new_col]]
                if (is.function(spec)) {
                    normalized[[new_col]] <- list(mode = "df", fun = spec, cols = character(0))
                    next
                }
                if (is.list(spec) && is.function(spec$fun)) {
                    cols <- spec$cols
                    if (is.null(cols)) {
                        cols <- character(0)
                    }
                    if (!is.character(cols) || anyNA(cols) || any(cols == "")) {
                        stop("`derive[[\"", new_col, "\"]]$cols` must be a character vector of column names.")
                    }
                    if (!all(cols %in% available_after_rename)) {
                        missing_cols <- setdiff(cols, available_after_rename)
                        stop("`derive[[\"", new_col, "\"]]` references missing columns: ", paste(missing_cols, collapse = ", "))
                    }
                    normalized[[new_col]] <- list(mode = "args", fun = spec$fun, cols = cols)
                    next
                }
                stop(
                    "`derive[[\"", new_col, "\"]]` must be a function or list(cols=..., fun=...)."
                )
            }
            normalized
        },

        .resolve_selected_columns = function(select = NULL, available_col_names) {
            if (is.null(select)) {
                return(unname(available_col_names))
            }

            if (is.character(select)) {
                if (anyNA(select) || any(select == "")) {
                    stop("`select` contains invalid empty/NA column names.")
                }
                if (!all(select %in% available_col_names)) {
                    missing_cols <- setdiff(select, available_col_names)
                    stop("`select` references missing columns: ", paste(missing_cols, collapse = ", "))
                }
                return(unname(unique(select)))
            }

            if (is.numeric(select)) {
                if (length(select) == 0L) {
                    return(character(0))
                }
                if (anyNA(select) || any(select < 1) || any(as.integer(select) != select) || any(select > length(available_col_names))) {
                    stop("Numeric `select` must contain valid positive integer column positions.")
                }
                return(unname(available_col_names[unique(as.integer(select))]))
            }

            stop("`select` must be NULL, a character vector, or a numeric vector.")
        },

        .create_index = function(index_name, columns) {
            if (private$.backend != "sqlite") {
                return(invisible())
            }
            private$.ensure_sqlite_connection()
            idx_sql_name <- gsub("[^A-Za-z0-9_]", "_", paste0("idx_registry_", index_name))
            table_quoted <- .registry_quote_ident(private$.conn, private$.table_name)
            idx_quoted <- .registry_quote_ident(private$.conn, idx_sql_name)
            cols_sql <- paste(vapply(columns, function(col) .registry_quote_ident(private$.conn, col), character(1)), collapse = ", ")
            sql <- paste0(
                "CREATE INDEX IF NOT EXISTS ",
                idx_quoted,
                " ON ",
                table_quoted,
                " (",
                cols_sql,
                ")"
            )
            DBI::dbExecute(private$.conn, sql)
            invisible()
        },

        .create_base_rowid_index = function() {
            if (private$.backend != "sqlite") {
                return(invisible())
            }
            private$.ensure_sqlite_connection()
            idx_sql_name <- "idx_registry_rowid"
            table_quoted <- .registry_quote_ident(private$.conn, private$.table_name)
            idx_quoted <- .registry_quote_ident(private$.conn, idx_sql_name)
            col_quoted <- .registry_quote_ident(private$.conn, ".registry_rowid")
            sql <- paste0(
                "CREATE INDEX IF NOT EXISTS ",
                idx_quoted,
                " ON ",
                table_quoted,
                " (",
                col_quoted,
                ")"
            )
            DBI::dbExecute(private$.conn, sql)
            invisible()
        },

        .ensure_sqlite_connection = function() {
            if (private$.backend != "sqlite") {
                return(invisible())
            }

            is_valid <- FALSE
            if (!is.null(private$.conn)) {
                is_valid <- isTRUE(tryCatch(DBI::dbIsValid(private$.conn), error = function(e) FALSE))
                if (!is_valid) {
                    private$.conn <- NULL
                }
            }
            if (!is.null(private$.conn)) {
                return(invisible())
            }

            if (!is.null(private$.sqlite_path) && file.exists(private$.sqlite_path)) {
                private$.conn <- DBI::dbConnect(RSQLite::SQLite(), dbname = private$.sqlite_path)
                private$.shared_conn <- FALSE
                private$.configure_sqlite_connection()
                return(invisible())
            }

            if (is.null(private$.sqlite_file_raw) || length(private$.sqlite_file_raw) == 0L) {
                stop("SQLite backend has no available database file to restore from.")
            }

            restored_path <- tempfile(pattern = "dmrsegal_registry_restore_", fileext = ".sqlite")
            private$.write_raw_file(restored_path, private$.sqlite_file_raw)
            private$.sqlite_path <- restored_path
            private$.owns_db_file <- TRUE
            private$.conn <- DBI::dbConnect(RSQLite::SQLite(), dbname = private$.sqlite_path)
            private$.shared_conn <- FALSE
            private$.configure_sqlite_connection()
            invisible()
        },

        .configure_sqlite_connection = function() {
            if (private$.backend != "sqlite" || is.null(private$.conn)) {
                return(invisible())
            }
            suppressWarnings(try(DBI::dbExecute(private$.conn, "PRAGMA temp_store = MEMORY"), silent = TRUE))
            suppressWarnings(try(DBI::dbExecute(private$.conn, "PRAGMA synchronous = NORMAL"), silent = TRUE))
            suppressWarnings(try(DBI::dbExecute(private$.conn, "PRAGMA foreign_keys = OFF"), silent = TRUE))
            suppressWarnings(try(DBI::dbExecute(private$.conn, "PRAGMA busy_timeout = 5000"), silent = TRUE))
            invisible()
        },

        .snapshot_sqlite_file = function() {
            if (private$.backend != "sqlite") {
                return(invisible())
            }
            private$.ensure_sqlite_connection()
            if (is.null(private$.sqlite_path) || !file.exists(private$.sqlite_path)) {
                return(invisible())
            }
            private$.sqlite_file_raw <- private$.read_file_raw(private$.sqlite_path)
            invisible()
        },

        .read_file_raw = function(path) {
            size <- file.info(path)$size
            if (is.na(size) || size <= 0) {
                return(raw(0))
            }
            con <- file(path, open = "rb")
            on.exit(close(con), add = TRUE)
            remaining <- as.double(size)
            out <- raw(0)
            chunk <- 1024 * 1024
            while (remaining > 0) {
                n <- as.integer(min(chunk, remaining))
                buf <- readBin(con, what = "raw", n = n)
                if (length(buf) == 0) {
                    break
                }
                out <- c(out, buf)
                remaining <- remaining - length(buf)
            }
            out
        },

        .write_raw_file = function(path, raw_data) {
            con <- file(path, open = "wb")
            on.exit(close(con), add = TRUE)
            writeBin(raw_data, con)
            invisible()
        },

        .resolve_master_index = function(master_index = NULL) {
            if (is.null(master_index)) {
                single_col_indices <- names(private$.indices)[
                    vapply(private$.indices, function(cols) length(cols) == 1L, logical(1))
                ]
                if (length(single_col_indices) == 1L) {
                    return(single_col_indices[[1]])
                }
                return(NULL)
            }

            if (!is.character(master_index) || length(master_index) != 1L || is.na(master_index) || master_index == "") {
                stop("`master_index` must be NULL or a non-empty character scalar.")
            }

            if (!is.null(private$.indices[[master_index]])) {
                if (length(private$.indices[[master_index]]) != 1L) {
                    stop("`master_index` must resolve to a single-column index.")
                }
                return(master_index)
            }

            if (master_index %in% private$.base_col_names) {
                private$.indices[[master_index]] <- master_index
                private$.create_index(master_index, master_index)
                return(master_index)
            }

            stop("`master_index` must match an existing index name or column name.")
        },

        .subset_rows_by_master_index = function(selector) {
            if (is.null(private$.master_index)) {
                stop(
                    "Character row selectors require a configured `master_index`. ",
                    "Set it in `getRegistry(master_index = ...)` or call `$setMasterIndex()`."
                )
            }
            if (anyNA(selector)) {
                stop("Character row selectors cannot contain NA.")
            }

            if (length(selector) == 0L || length(private$.active_row_idx) == 0L) {
                return(integer(0))
            }

            idx_cols <- private$.indices[[private$.master_index]]
            if (length(idx_cols) != 1L) {
                stop("Internal error: `master_index` does not resolve to a single column.")
            }
            idx_col <- idx_cols[[1]]
            selector_chr <- as.character(selector)

            if (private$.backend == "memory") {
                idx_cache <- private$.get_memory_index_map(private$.master_index, idx_cols)
                selector_keys <- .registry_encode_single_key(selector_chr)
                seq_bounds <- .registry_seq_bounds(private$.active_row_idx)

                if (!is.null(seq_bounds)) {
                    lo <- as.integer(seq_bounds[[1L]])
                    hi <- as.integer(seq_bounds[[2L]])
                    if (identical(idx_cache$mode, "unique_single")) {
                        selected_row_ids <- unname(idx_cache$map[selector_keys])
                        selected_row_ids <- selected_row_ids[!is.na(selected_row_ids)]
                        selected_row_ids <- selected_row_ids[selected_row_ids >= lo & selected_row_ids <= hi]
                    } else {
                        selected_row_ids <- unlist(lapply(selector_keys, function(key) {
                            rows <- idx_cache$map[[key]]
                            if (is.null(rows) || length(rows) == 0L) {
                                return(integer(0))
                            }
                            rows[rows >= lo & rows <= hi]
                        }), use.names = FALSE)
                    }
                    if (length(selected_row_ids) == 0L) {
                        return(integer(0))
                    }
                    return(as.integer(selected_row_ids))
                }

                active_pos_by_rowid <- split(seq_along(private$.active_row_idx), private$.active_row_idx)
                if (identical(idx_cache$mode, "unique_single")) {
                    selected_row_ids <- unname(idx_cache$map[selector_keys])
                    selected_row_ids <- selected_row_ids[!is.na(selected_row_ids)]
                    if (length(selected_row_ids) == 0L) {
                        return(integer(0))
                    }
                    selected_rel_pos <- unlist(lapply(selected_row_ids, function(row_id) {
                        active_pos_by_rowid[[as.character(row_id)]]
                    }), use.names = FALSE)
                } else {
                    selected_rel_pos <- unlist(lapply(selector_keys, function(key) {
                        rows <- idx_cache$map[[key]]
                        if (is.null(rows) || length(rows) == 0L) {
                            return(integer(0))
                        }
                        rel <- unlist(lapply(rows, function(row_id) {
                            active_pos_by_rowid[[as.character(row_id)]]
                        }), use.names = FALSE)
                        if (length(rel) > 1L) {
                            rel <- sort.int(rel, method = "quick")
                        }
                        rel
                    }), use.names = FALSE)
                }
                if (length(selected_rel_pos) == 0L) {
                    return(integer(0))
                }
                return(private$.active_row_idx[selected_rel_pos])
            }

            private$.match_active_rowids_sqlite(index_col = idx_col, selector = selector_chr)
        },

        .subset_rows_by_master_index_range = function(start_value, end_value) {
            if (is.null(private$.master_index)) {
                stop(
                    "Character row selectors require a configured `master_index`. ",
                    "Set it in `getRegistry(master_index = ...)` or call `$setMasterIndex()`."
                )
            }

            if (!is.character(start_value) || length(start_value) != 1L || is.na(start_value) ||
                !is.character(end_value) || length(end_value) != 1L || is.na(end_value)) {
                stop("Master index range selectors must be two non-missing character scalars.")
            }

            if (length(private$.active_row_idx) == 0L) {
                return(integer(0))
            }

            active_names <- private$.active_row_names(prefix = "row")
            start_pos <- match(start_value, active_names)
            end_pos <- match(end_value, active_names)

            if (is.na(start_pos) || is.na(end_pos)) {
                missing_vals <- c()
                if (is.na(start_pos)) missing_vals <- c(missing_vals, start_value)
                if (is.na(end_pos)) missing_vals <- c(missing_vals, end_value)
                stop(
                    "Master index range selector value(s) not found in active rows: ",
                    paste(missing_vals, collapse = ", ")
                )
            }

            if (start_pos <= end_pos) {
                rel <- seq.int(start_pos, end_pos)
            } else {
                rel <- seq.int(start_pos, end_pos, by = -1L)
            }
            private$.active_row_idx[rel]
        },

        .active_row_names = function(prefix = "row") {
            if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix)) {
                stop("`prefix` must be a single non-missing character string.")
            }
            if (length(private$.active_row_idx) == 0L) {
                return(character(0))
            }
            if (is.null(private$.master_index)) {
                return(paste0(prefix, private$.active_row_idx))
            }

            idx_cols <- private$.indices[[private$.master_index]]
            if (length(idx_cols) != 1L) {
                stop("Internal error: `master_index` does not resolve to a single column.")
            }
            idx_col <- idx_cols[[1]]
            if (identical(private$.backend, "memory")) {
                out <- private$.data_in_memory[[idx_col]][private$.active_row_idx]
            } else {
                out <- private$.fetch_rows_with_columns(private$.active_row_idx, idx_col)[[idx_col]]
            }
            as.character(out)
        },

        .apply_master_row_names = function(df, row_ids = NULL, prefix = "row") {
            if (!is.data.frame(df) || nrow(df) == 0L) {
                return(df)
            }
            if (is.null(private$.master_index)) {
                return(df)
            }

            idx_cols <- private$.indices[[private$.master_index]]
            if (!is.character(idx_cols) || length(idx_cols) != 1L) {
                return(df)
            }
            idx_col <- idx_cols[[1L]]

            vals <- NULL
            if (idx_col %in% colnames(df)) {
                vals <- df[[idx_col]]
            } else if (!is.null(row_ids)) {
                if (length(row_ids) != nrow(df)) {
                    return(df)
                }
                if (identical(private$.backend, "memory")) {
                    vals <- private$.data_in_memory[[idx_col]][as.integer(row_ids)]
                } else {
                    vals <- private$.fetch_rows_with_columns(as.integer(row_ids), idx_col)[[idx_col]]
                }
            }

            if (is.null(vals) || length(vals) != nrow(df)) {
                return(df)
            }

            rownames(df) <- .registry_safe_row_names(vals, prefix = prefix)
            df
        },

        .spawn_view = function(new_row_idx, new_col_idx) {
            saved_conn <- private$.conn
            saved_shared <- private$.shared_conn
            private$.conn <- NULL
            private$.shared_conn <- FALSE
            on.exit({
                private$.conn <- saved_conn
                private$.shared_conn <- saved_shared
            }, add = TRUE)

            view <- self$clone(deep = FALSE)
            view_private <- view$.__enclos_env__$private
            view_private$.active_row_idx <- new_row_idx
            view_private$.active_col_idx <- new_col_idx

            # SQLite views share the same backing DB metadata. They must never
            # delete the shared source file.
            if (identical(view_private$.backend, "sqlite")) {
                view_private$.owns_db_file <- FALSE
                view_private$.conn <- saved_conn
                view_private$.shared_conn <- !is.null(saved_conn)
            }
            view
        },

        .match_active_rowids_sqlite = function(index_col, selector) {
            private$.ensure_sqlite_connection()
            conn <- private$.conn
            data_table_q <- .registry_quote_ident(conn, private$.table_name)
            tmp_keys <- .registry_new_temp_table("tmp_registry_master_keys")
            tmp_keys_q <- .registry_quote_ident(conn, tmp_keys)
            tmp_active <- NULL
            tmp_active_q <- NULL

            on.exit({
                if (DBI::dbExistsTable(conn, tmp_keys)) {
                    DBI::dbExecute(conn, paste0("DROP TABLE ", tmp_keys_q))
                }
                if (!is.null(tmp_active) && DBI::dbExistsTable(conn, tmp_active)) {
                    DBI::dbExecute(conn, paste0("DROP TABLE ", tmp_active_q))
                }
            }, add = TRUE)

            key_template <- data.frame(
                .query_ord = integer(0),
                .lookup_key = character(0),
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
            DBI::dbWriteTable(conn, tmp_keys, key_template, temporary = TRUE, row.names = FALSE)
            ord <- seq_along(selector)
            for (start in seq.int(1L, length(selector), by = private$.query_chunk_size)) {
                end <- min(start + private$.query_chunk_size - 1L, length(selector))
                chunk <- data.frame(
                    .query_ord = ord[start:end],
                    .lookup_key = selector[start:end],
                    stringsAsFactors = FALSE,
                    check.names = FALSE
                )
                DBI::dbWriteTable(conn, tmp_keys, chunk, append = TRUE, row.names = FALSE)
            }

            idx_col_q <- .registry_quote_ident(conn, index_col)
            active_bounds <- .registry_seq_bounds(private$.active_row_idx)
            if (!is.null(active_bounds)) {
                sql <- paste0(
                    "SELECT d.\".registry_rowid\"",
                    " FROM ", tmp_keys_q, " q",
                    " JOIN ", data_table_q, " d ON d.", idx_col_q, " = q.\".lookup_key\"",
                    " WHERE d.\".registry_rowid\" BETWEEN ",
                    as.integer(active_bounds[[1L]]),
                    " AND ",
                    as.integer(active_bounds[[2L]]),
                    " ORDER BY q.\".query_ord\", d.\".registry_rowid\""
                )
            } else {
                tmp_active <- .registry_new_temp_table("tmp_registry_master_active")
                tmp_active_q <- .registry_quote_ident(conn, tmp_active)
                active_template <- data.frame(
                    .active_ord = integer(0),
                    .registry_rowid = integer(0),
                    stringsAsFactors = FALSE,
                    check.names = FALSE
                )
                DBI::dbWriteTable(conn, tmp_active, active_template, temporary = TRUE, row.names = FALSE)
                active_ord <- seq_along(private$.active_row_idx)
                for (start in seq.int(1L, length(private$.active_row_idx), by = private$.query_chunk_size)) {
                    end <- min(start + private$.query_chunk_size - 1L, length(private$.active_row_idx))
                    chunk <- data.frame(
                        .active_ord = active_ord[start:end],
                        .registry_rowid = private$.active_row_idx[start:end],
                        stringsAsFactors = FALSE,
                        check.names = FALSE
                    )
                    DBI::dbWriteTable(conn, tmp_active, chunk, append = TRUE, row.names = FALSE)
                }
                sql <- paste0(
                    "SELECT a.\".registry_rowid\"",
                    " FROM ", tmp_keys_q, " q",
                    " JOIN ", data_table_q, " d ON d.", idx_col_q, " = q.\".lookup_key\"",
                    " JOIN ", tmp_active_q, " a ON a.\".registry_rowid\" = d.\".registry_rowid\"",
                    " ORDER BY q.\".query_ord\", a.\".active_ord\""
                )
            }
            res <- DBI::dbGetQuery(conn, sql)
            if (nrow(res) == 0L) {
                return(integer(0))
            }
            as.integer(res[[1]])
        },

        .subset_absolute_indices = function(current_idx, selector, axis_names = NULL, axis = c("rows", "columns")) {
            axis <- match.arg(axis)
            if (is.null(selector)) {
                return(current_idx)
            }
            if (is.character(selector)) {
                if (axis != "columns") {
                    stop("Character selectors are only supported for columns.")
                }
                if (is.null(axis_names)) {
                    stop("Internal error: missing axis names for column selection.")
                }
                mapped <- match(selector, axis_names)
                if (anyNA(mapped)) {
                    missing_names <- unique(selector[is.na(mapped)])
                    stop("Unknown ", axis, " names: ", paste(missing_names, collapse = ", "))
                }
                return(current_idx[mapped])
            }
            if (is.logical(selector)) {
                if (length(selector) == 0L || length(current_idx) == 0L) {
                    return(integer(0))
                }
                if (anyNA(selector)) {
                    stop("Logical ", axis, " selector cannot contain NA.")
                }
                n <- length(current_idx)
                if (length(selector) != n) {
                    if (length(selector) == 1L || (n %% length(selector) == 0L)) {
                        selector <- rep_len(selector, n)
                    } else {
                        stop(
                            "Logical ",
                            axis,
                            " selector length (",
                            length(selector),
                            ") is not compatible with active ",
                            axis,
                            " length (",
                            n,
                            ")."
                        )
                    }
                }
                return(current_idx[selector])
            }
            pos_idx <- seq_along(current_idx)[selector]
            if (anyNA(pos_idx)) {
                stop("Invalid ", axis, " selector: out-of-range or NA indices are not allowed.")
            }
            current_idx[pos_idx]
        },

        .resolve_output_columns = function(columns) {
            active_cols <- private$.base_col_names[private$.active_col_idx]
            if (is.null(columns)) {
                return(active_cols)
            }
            if (is.character(columns)) {
                if (!all(columns %in% active_cols)) {
                    missing_cols <- setdiff(columns, active_cols)
                    stop("Requested output columns are outside active bounds: ", paste(missing_cols, collapse = ", "))
                }
                return(columns)
            }
            pos <- seq_along(active_cols)[columns]
            if (anyNA(pos)) {
                stop("Invalid output column selector.")
            }
            active_cols[pos]
        },

        .get_memory_index_map = function(index_name, index_cols) {
            if (private$.backend != "memory") {
                stop("Internal error: memory index map requested for non-memory backend.")
            }
            cache <- private$.memory_index_cache[[index_name]]
            if (!is.null(cache)) {
                return(cache)
            }
            row_ids <- private$.data_in_memory$.registry_rowid
            source_df <- private$.data_in_memory[, index_cols, drop = FALSE]
            keys <- .registry_make_index_keys(source_df, index_cols)

            if (length(index_cols) == 1L && anyDuplicated(keys) == 0L) {
                cache <- list(
                    mode = "unique_single",
                    map = setNames(row_ids, keys)
                )
            } else {
                cache <- list(
                    mode = "generic",
                    map = split(row_ids, keys)
                )
            }

            private$.memory_index_cache[[index_name]] <- cache
            cache
        },

        .normalize_lookup_values = function(values, index_cols) {
            if (length(index_cols) == 1) {
                if (is.data.frame(values)) {
                    if (ncol(values) != 1) {
                        stop("For single-column index `", index_cols, "`, `values` must have one column.")
                    }
                    lookup_df <- values[, 1, drop = FALSE]
                    colnames(lookup_df) <- index_cols
                } else {
                    lookup_df <- data.frame(setNames(list(values), index_cols), stringsAsFactors = FALSE, check.names = FALSE)
                }
            } else {
                if (!is.data.frame(values)) {
                    stop("For multi-column indices, `values` must be a data.frame.")
                }
                if (!setequal(colnames(values), index_cols)) {
                    stop(
                        "Lookup data.frame columns must match index columns exactly.\n",
                        "Expected: ", paste(index_cols, collapse = ", "), "\n",
                        "Got: ", paste(colnames(values), collapse = ", ")
                    )
                }
                lookup_df <- values[, index_cols, drop = FALSE]
            }
            lookup_df <- as.data.frame(lookup_df, stringsAsFactors = FALSE, check.names = FALSE)
            for (col in index_cols) {
                lookup_df[[col]] <- as.character(lookup_df[[col]])
            }
            if (anyNA(lookup_df)) {
                stop("Lookup values cannot contain NA.")
            }
            lookup_df
        },

        .empty_result_df = function(cols) {
            if (length(cols) == 0) {
                return(data.frame(row.names = integer(0)))
            }
            out <- as.data.frame(setNames(vector("list", length(cols)), cols), stringsAsFactors = FALSE, check.names = FALSE)
            out[0, , drop = FALSE]
        },

        .fetch_rows_with_columns = function(row_ids, cols) {
            if (length(row_ids) == 0) {
                return(private$.empty_result_df(cols))
            }
            if (private$.backend == "memory") {
                base_df <- private$.data_in_memory
                if (anyNA(row_ids) || any(row_ids < 1L) || any(row_ids > nrow(base_df))) {
                    stop("Invalid active row bounds detected.")
                }
                if (length(cols) == 0) {
                    return(data.frame(row.names = seq_along(row_ids)))
                }
                return(base_df[row_ids, cols, drop = FALSE])
            }
            private$.fetch_rows_sqlite(row_ids, cols)
        },

        .fetch_rows_sqlite = function(row_ids, cols) {
            private$.ensure_sqlite_connection()
            conn <- private$.conn
            table_quoted <- .registry_quote_ident(conn, private$.table_name)
            unique_cols <- unique(cols)
            if (length(unique_cols) == 0) {
                return(data.frame(row.names = seq_along(row_ids)))
            }

            seq_bounds <- .registry_seq_bounds(row_ids)
            if (!is.null(seq_bounds)) {
                select_parts <- vapply(unique_cols, function(col) paste0("d.", .registry_quote_ident(conn, col)), character(1))
                sql <- paste0(
                    "SELECT ",
                    paste(select_parts, collapse = ", "),
                    " FROM ",
                    table_quoted,
                    " d WHERE d.\".registry_rowid\" BETWEEN ",
                    as.integer(seq_bounds[[1L]]),
                    " AND ",
                    as.integer(seq_bounds[[2L]]),
                    " ORDER BY d.\".registry_rowid\""
                )
                res <- DBI::dbGetQuery(conn, sql)
                return(res[, cols, drop = FALSE])
            }

            tmp_active <- .registry_new_temp_table("tmp_registry_active_rows")
            tmp_active_q <- .registry_quote_ident(conn, tmp_active)
            on.exit({
                if (DBI::dbExistsTable(conn, tmp_active)) {
                    DBI::dbExecute(conn, paste0("DROP TABLE ", tmp_active_q))
                }
            }, add = TRUE)

            active_template <- data.frame(
                .active_ord = integer(0),
                .registry_rowid = integer(0),
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
            DBI::dbWriteTable(conn, tmp_active, active_template, temporary = TRUE, row.names = FALSE)
            ord <- seq_along(row_ids)
            for (start in seq.int(1L, length(row_ids), by = private$.query_chunk_size)) {
                end <- min(start + private$.query_chunk_size - 1L, length(row_ids))
                chunk <- data.frame(
                    .active_ord = ord[start:end],
                    .registry_rowid = row_ids[start:end],
                    stringsAsFactors = FALSE,
                    check.names = FALSE
                )
                DBI::dbWriteTable(conn, tmp_active, chunk, append = TRUE, row.names = FALSE)
            }

            select_parts <- c(
                "a.\".active_ord\"",
                vapply(unique_cols, function(col) paste0("d.", .registry_quote_ident(conn, col)), character(1))
            )
            sql <- paste0(
                "SELECT ",
                paste(select_parts, collapse = ", "),
                " FROM ",
                tmp_active_q,
                " a JOIN ",
                table_quoted,
                " d ON d.\".registry_rowid\" = a.\".registry_rowid\"",
                " ORDER BY a.\".active_ord\""
            )
            res <- DBI::dbGetQuery(conn, sql)
            res <- res[, unique_cols, drop = FALSE]
            res[, cols, drop = FALSE]
        },

        .query_memory = function(index_name, index_cols, lookup_df, out_cols) {
            if (length(private$.active_row_idx) == 0L) {
                return(private$.empty_result_df(out_cols))
            }

            query_keys <- .registry_make_index_keys(lookup_df, index_cols)
            idx_cache <- private$.get_memory_index_map(index_name, index_cols)
            seq_bounds <- .registry_seq_bounds(private$.active_row_idx)
            if (!is.null(seq_bounds)) {
                lo <- as.integer(seq_bounds[[1L]])
                hi <- as.integer(seq_bounds[[2L]])
                if (identical(idx_cache$mode, "unique_single")) {
                    matched_row_ids <- unname(idx_cache$map[query_keys])
                    matched_row_ids <- matched_row_ids[!is.na(matched_row_ids)]
                    matched_row_ids <- matched_row_ids[matched_row_ids >= lo & matched_row_ids <= hi]
                } else {
                    matched_row_ids <- unlist(lapply(query_keys, function(k) {
                        rows <- idx_cache$map[[k]]
                        if (is.null(rows) || length(rows) == 0L) {
                            return(integer(0))
                        }
                        rows[rows >= lo & rows <= hi]
                    }), use.names = FALSE)
                }
                if (length(matched_row_ids) == 0L) {
                    return(private$.empty_result_df(out_cols))
                }
                out <- private$.fetch_rows_with_columns(matched_row_ids, out_cols)
                return(private$.apply_master_row_names(out, row_ids = matched_row_ids))
            }

            if (identical(idx_cache$mode, "unique_single")) {
                matched_row_ids <- unname(idx_cache$map[query_keys])
                matched_row_ids <- matched_row_ids[!is.na(matched_row_ids)]
                if (length(matched_row_ids) == 0L) {
                    return(private$.empty_result_df(out_cols))
                }
                active_pos_by_rowid <- split(seq_along(private$.active_row_idx), private$.active_row_idx)
                matched_positions <- unlist(lapply(matched_row_ids, function(row_id) {
                    active_pos_by_rowid[[as.character(row_id)]]
                }), use.names = FALSE)
                if (length(matched_positions) == 0L) {
                    return(private$.empty_result_df(out_cols))
                }
                matched_row_ids <- private$.active_row_idx[matched_positions]
                if (length(out_cols) == 0L) {
                    out <- data.frame(row.names = seq_along(matched_positions))
                    return(private$.apply_master_row_names(out, row_ids = matched_row_ids))
                }
                out <- private$.fetch_rows_with_columns(matched_row_ids, out_cols)
                return(private$.apply_master_row_names(out, row_ids = matched_row_ids))
            }

            needed_cols <- unique(c(index_cols, out_cols))
            active_data <- private$.fetch_rows_with_columns(private$.active_row_idx, needed_cols)
            if (nrow(active_data) == 0) {
                return(private$.empty_result_df(out_cols))
            }
            active_keys <- .registry_make_index_keys(active_data, index_cols)
            pos_by_key <- split(seq_along(active_keys), active_keys)
            matched_positions <- unlist(lapply(query_keys, function(k) pos_by_key[[k]]), use.names = FALSE)
            if (length(matched_positions) == 0L) {
                return(private$.empty_result_df(out_cols))
            }
            matched_row_ids <- private$.active_row_idx[matched_positions]
            if (length(out_cols) == 0L) {
                out <- data.frame(row.names = seq_along(matched_positions))
                return(private$.apply_master_row_names(out, row_ids = matched_row_ids))
            }
            out <- active_data[matched_positions, out_cols, drop = FALSE]
            private$.apply_master_row_names(out, row_ids = matched_row_ids)
        },

        .query_sqlite = function(index_cols, lookup_df, out_cols) {
            private$.ensure_sqlite_connection()
            if (length(private$.active_row_idx) == 0L) {
                return(private$.empty_result_df(out_cols))
            }
            conn <- private$.conn
            data_table_q <- .registry_quote_ident(conn, private$.table_name)
            tmp_keys <- .registry_new_temp_table("tmp_registry_keys")
            tmp_keys_q <- .registry_quote_ident(conn, tmp_keys)
            tmp_active <- NULL
            tmp_active_q <- NULL

            on.exit({
                if (DBI::dbExistsTable(conn, tmp_keys)) {
                    DBI::dbExecute(conn, paste0("DROP TABLE ", tmp_keys_q))
                }
                if (!is.null(tmp_active) && DBI::dbExistsTable(conn, tmp_active)) {
                    DBI::dbExecute(conn, paste0("DROP TABLE ", tmp_active_q))
                }
            }, add = TRUE)

            key_template <- lookup_df[0, , drop = FALSE]
            key_template$.query_ord <- integer(0)
            key_template <- key_template[, c(".query_ord", index_cols), drop = FALSE]
            for (col in index_cols) {
                key_template[[col]] <- as.character(key_template[[col]])
            }
            DBI::dbWriteTable(conn, tmp_keys, key_template, temporary = TRUE, row.names = FALSE)

            if (nrow(lookup_df) > 0) {
                ord <- seq_len(nrow(lookup_df))
                for (start in seq.int(1L, nrow(lookup_df), by = private$.query_chunk_size)) {
                    end <- min(start + private$.query_chunk_size - 1L, nrow(lookup_df))
                    chunk <- lookup_df[start:end, , drop = FALSE]
                    for (col in index_cols) {
                        chunk[[col]] <- as.character(chunk[[col]])
                    }
                    chunk$.query_ord <- ord[start:end]
                    chunk <- chunk[, c(".query_ord", index_cols), drop = FALSE]
                    DBI::dbWriteTable(conn, tmp_keys, chunk, append = TRUE, row.names = FALSE)
                }
            }

            join_pred <- paste(vapply(index_cols, function(col) {
                col_q <- .registry_quote_ident(conn, col)
                paste0("d.", col_q, " = q.", col_q)
            }, character(1)), collapse = " AND ")

            unique_out_cols <- unique(out_cols)
            active_bounds <- .registry_seq_bounds(private$.active_row_idx)
            if (!is.null(active_bounds)) {
                select_parts <- c("q.\".query_ord\"", "d.\".registry_rowid\"")
                if (length(unique_out_cols) > 0) {
                    select_parts <- c(
                        select_parts,
                        vapply(unique_out_cols, function(col) paste0("d.", .registry_quote_ident(conn, col)), character(1))
                    )
                }
                sql <- paste0(
                    "SELECT ",
                    paste(select_parts, collapse = ", "),
                    " FROM ",
                    tmp_keys_q,
                    " q JOIN ",
                    data_table_q,
                    " d ON ",
                    join_pred,
                    " WHERE d.\".registry_rowid\" BETWEEN ",
                    as.integer(active_bounds[[1L]]),
                    " AND ",
                    as.integer(active_bounds[[2L]]),
                    " ORDER BY q.\".query_ord\", d.\".registry_rowid\""
                )
            } else {
                tmp_active <- .registry_new_temp_table("tmp_registry_active")
                tmp_active_q <- .registry_quote_ident(conn, tmp_active)
                active_template <- data.frame(
                    .active_ord = integer(0),
                    .registry_rowid = integer(0),
                    stringsAsFactors = FALSE,
                    check.names = FALSE
                )
                DBI::dbWriteTable(conn, tmp_active, active_template, temporary = TRUE, row.names = FALSE)
                ord <- seq_along(private$.active_row_idx)
                for (start in seq.int(1L, length(private$.active_row_idx), by = private$.query_chunk_size)) {
                    end <- min(start + private$.query_chunk_size - 1L, length(private$.active_row_idx))
                    chunk <- data.frame(
                        .active_ord = ord[start:end],
                        .registry_rowid = private$.active_row_idx[start:end],
                        stringsAsFactors = FALSE,
                        check.names = FALSE
                    )
                    DBI::dbWriteTable(conn, tmp_active, chunk, append = TRUE, row.names = FALSE)
                }

                select_parts <- c("q.\".query_ord\"", "a.\".active_ord\"", "d.\".registry_rowid\"")
                if (length(unique_out_cols) > 0) {
                    select_parts <- c(
                        select_parts,
                        vapply(unique_out_cols, function(col) paste0("d.", .registry_quote_ident(conn, col)), character(1))
                    )
                }
                sql <- paste0(
                    "SELECT ",
                    paste(select_parts, collapse = ", "),
                    " FROM ",
                    tmp_keys_q,
                    " q JOIN ",
                    data_table_q,
                    " d ON ",
                    join_pred,
                    " JOIN ",
                    tmp_active_q,
                    " a ON a.\".registry_rowid\" = d.\".registry_rowid\"",
                    " ORDER BY q.\".query_ord\", a.\".active_ord\""
                )
            }

            res <- DBI::dbGetQuery(conn, sql)
            if (nrow(res) == 0) {
                return(private$.empty_result_df(out_cols))
            }
            row_ids <- as.integer(res[[".registry_rowid"]])
            if (length(unique_out_cols) == 0) {
                out <- data.frame(row.names = seq_len(nrow(res)))
                return(private$.apply_master_row_names(out, row_ids = row_ids))
            }
            res <- res[, unique_out_cols, drop = FALSE]
            out <- res[, out_cols, drop = FALSE]
            private$.apply_master_row_names(out, row_ids = row_ids)
        }
    )
)
