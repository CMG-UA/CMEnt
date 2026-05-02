#' Get File Hash for Caching
#'
#' @description Internal helper to compute MD5 hash of a file for cache validation
#'
#' @param file_path Path to the file
#' @return MD5 hash string
#' @keywords internal
#' @noRd
.getFileHash <- function(file_path) {
    tryCatch(
        {
            tools::md5sum(file_path)
        },
        error = function(e) {
            # Fallback: use file size and modification time
            info <- file.info(file_path)
            paste0(info$size, "_", as.numeric(info$mtime))
        }
    )
}


.getBetaColNamesAndInds <- function(beta_file, beta_col_names = NULL, is_tabix = FALSE) {
    if (endsWith(beta_file, "gz")) {
        conn <- gzfile(beta_file, "r")
    } else {
        conn <- file(beta_file, "r")
    }
    file_beta_col_names <- scan(conn,
        sep = "\t",
        what = character(),
        nlines = 1,
        quiet = TRUE
    )
    close(conn)
    if (is_tabix) {
        file_beta_col_names <- file_beta_col_names[7:length(file_beta_col_names)]
    }
    if (is.null(beta_col_names)) {
        beta_col_names <- file_beta_col_names[-1]
        cols_inds <- seq(2, length(beta_col_names))
    } else {
        cols_inds <- match(beta_col_names, file_beta_col_names)
        cols_inds <- cols_inds[!is.na(cols_inds), drop = FALSE]
        if (length(cols_inds) == 0) {
            stop(
                "Beta file does not contain any phenotype rownames ",
                "as column name. First 5 supplied phenotype rownames: ",
                paste(beta_col_names[seq_len(min(5, length(beta_col_names)))], sep = ",")
            )
        }
        beta_col_names <- file_beta_col_names[cols_inds]
    }
    ret <- list(
        beta_col_names = beta_col_names, beta_col_inds = cols_inds,
        file_beta_col_names = file_beta_col_names[-1]
    )
    invisible(ret)
}


.subsetBetaFile <- function(beta_file,
                            sites,
                            beta_row_names = NULL,
                            beta_col_names = NULL) {
    if (is.null(beta_row_names)) {
        beta_row_names <- unlist(data.table::fread(
            file = beta_file,
            select = 1,
            header = TRUE
        ))
    }

    ret <- .getBetaColNamesAndInds(beta_file, beta_col_names)
    beta_col_names <- ret$beta_col_names
    cols_inds <- ret$beta_col_inds
    if (endsWith(beta_file, "gz")) {
        conn <- gzfile(beta_file, "r")
    } else {
        conn <- file(beta_file, "r")
    }
    if (!is.null(sites)) {
        sites_inds <- which(beta_row_names %in% sites)
        sites_inds_steps <- diff(c(-1, sites_inds)) - 1
        beta_sites <- lapply(sites_inds_steps, function(step) {
            l <- scan(
                conn,
                what = character(),
                skip = step,
                nlines = 1,
                quiet = TRUE
            )
            as.numeric(l[cols_inds])
        })
    } else {
        beta_sites <- readLines(conn)[2:(length(beta_row_names) + 1)]
    }
    close(conn)
    beta_sites <- do.call(rbind, beta_sites)
    rownames(beta_sites) <- sites
    colnames(beta_sites) <- beta_col_names
    beta_sites
}


.readSamplesheet <- function(samplesheet_file,
                             samplesheet_file_sep,
                             sample_group_col,
                             sample_group_case,
                             sample_group_control,
                             target_col,
                             subset = NULL) {
    ref <- read.table(
        samplesheet_file,
        header = TRUE,
        sep = samplesheet_file_sep,
        quote = "",
        comment.char = "",
        check.names = FALSE
    )
    if (is.null(sample_group_col)) {
        return(ref)
    }
    if (is.null(sample_group_case) && (is.null(sample_group_control))) {
        stop("Both sample_group_case and sample_group_control provided are NULL, cannot continue.")
    }
    if (!(sample_group_col %in% colnames(ref))) {
        stop(
            "Sample group column ",
            sample_group_col,
            " not found in samplesheet columns: ",
            paste0(colnames(ref), collapse = ",")
        )
    }
    if (!(target_col %in% colnames(ref))) {
        stop(
            "Target column ",
            target_col,
            " not found in samplesheet columns: ",
            paste0(colnames(ref), collapse = ",")
        )
    }
    if (!is.null(subset)) {
        stopifnot(any(ref[, target_col] %in% subset))
        ref <- ref[ref[, target_col] %in% subset, ]
    }
    if (!is.null(sample_group_control)) {
        control_samplesheet <- ref[ref[, sample_group_col] %in% sample_group_control, c(target_col, sample_group_col),
            drop =
                FALSE
        ]
        if (nrow(control_samplesheet) == 0) {
            stop(
                paste0(
                    "The control label ",
                    sample_group_control,
                    " was not found in the provided samplesheet sample_group_col ",
                    sample_group_col,
                    " which has the following values: ",
                    paste(unique(ref[, sample_group_col]), collapse = ",")
                )
            )
        }
    }

    if (!is.null(sample_group_case)) {
        case_samplesheet <- ref[ref[, sample_group_col] %in% sample_group_case, c(target_col, sample_group_col),
            drop =
                FALSE
        ]

        if (nrow(case_samplesheet) == 0) {
            stop(
                paste0(
                    "The case label ",
                    sample_group_case,
                    " was not found in the provided samplesheet sample_group_col ",
                    sample_group_col,
                    " which has the following values: ",
                    paste(unique(ref[, sample_group_col]), collapse = ",")
                )
            )
        }
    }

    if (!is.null(sample_group_case) || !is.null(sample_group_control)) {
        if (!is.null(sample_group_case) && !is.null(sample_group_control)) {
            subset_samplesheet <- rbind(control_samplesheet, case_samplesheet)
        } else {
            if (is.null(sample_group_case)) {
                subset_samplesheet <- control_samplesheet
            } else {
                subset_samplesheet <- case_samplesheet
            }
        }
    } else {
        subset_samplesheet <- ref
    }

    rownames(subset_samplesheet) <- subset_samplesheet[, target_col]
    subset_samplesheet <- subset_samplesheet[, sample_group_col,
        drop =
            FALSE
    ]
    subset_samplesheet <- subset_samplesheet[
        stringr::str_order(rownames(subset_samplesheet),
            numeric = TRUE
        ), ,
        drop = FALSE
    ]
    stopifnot(ncol(subset_samplesheet) != 0)
    stopifnot(nrow(subset_samplesheet) != 0)
    .log_info(
        "Read samplesheet head:\n\t",
        paste(capture.output(print(
            head(subset_samplesheet)
        )), collapse = "\n\t")
    )
    if (!is.null(sample_group_col)) {
        if (!is.null(sample_group_control)) {
            subset_samplesheet[, "casecontrol"] <- 1
            if (length(sample_group_control) == 1) {
                sample_group_control <- base::strsplit(sample_group_control, ",")[[1]]
            }
            subset_samplesheet[subset_samplesheet[, sample_group_col] %in% sample_group_control, "casecontrol"] <- 0
        } else {
            subset_samplesheet[, "casecontrol"] <- 0
            if (length(sample_group_case) == 1) {
                sample_group_case <- base::strsplit(sample_group_case, ",")[[1]]
            }
            subset_samplesheet[subset_samplesheet[, sample_group_col] %in% sample_group_case, "casecontrol"] <- 1
        }
    }
    subset_samplesheet
}

# Lightweight styled logging helpers -----------------------------------------

# Internal state for timing steps
.CMEnt_log_env <- local({ # nolint
    e <- new.env(parent = emptyenv())
    e$last_step_time <- list()
    e
})

.fmt_dur <- function(start_time) {
    if (is.null(start_time)) {
        return("")
    }
    secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (is.na(secs) || secs < 0) {
        return("")
    }
    if (secs < 0.001) {
        x <- " (took %.2f\u03bcs)"
        Encoding(x) <- "UTF-8"
        sprintf(x, secs * 1000000)
    } else if (secs < 1) {
        sprintf(" (took %.2fms)", secs * 1000)
    } else if (secs < 60) {
        sprintf(" (took %.2fs)", secs)
    } else {
        sprintf(" (took %dm %02ds)", floor(secs / 60), round(secs %% 60))
    }
}

.has_ansi <- function() {
    isFALSE(nzchar(Sys.getenv("NO_COLOR"))) && cli::num_ansi_colors() > 1L
}

.col <- function(x, col = c("cyan", "green", "yellow", "blue", "red")) {
    col <- strex::match_arg(col, ignore_case = TRUE)
    if (!.has_ansi()) {
        return(x)
    }
    switch(col,
        cyan   = cli::col_cyan(x),
        green  = cli::col_green(x),
        yellow = cli::col_yellow(x),
        blue   = cli::col_blue(x),
        red    = cli::col_red(x)
    )
}



#' @keywords internal
#' @noRd
.format_log_output <- function(msg, lead, level, width = getOption("width")) {
    prefix <- paste(rep("\t", max(0, level - 1)), collapse = "")
    paste(strwrap(prefix = paste(" ", prefix), initial =  paste(prefix, lead, " ", sep = ""), width = (0.9 - 0.05 * max(0, level - 1)) * getOption("width"), msg), collapse = "\n")
}

#' @keywords internal
#' @noRd
.log_error <- function(..., .envir = parent.frame()) {
    msg <- paste0(..., collapse = "")
    lead <- .col(cli::symbol$cross, "red")
    stop(.format_log_output(msg, lead = lead, level = 0), call. = FALSE)
}

#' @keywords internal
#' @noRd
.log_warn <- function(..., .envir = parent.frame()) {
    msg <- paste0(..., collapse = "")
    lead <- .col(cli::symbol$warning, "yellow")
    warning(.format_log_output(msg, lead = lead, level = 0))
    invisible()
}

.node_size <- function() {
    bit <- 8L * .Machine$sizeof.pointer
    if (!(bit == 32L || bit == 64L)) {
        stop("Unknown architecture", call. = FALSE)
    }

    if (bit == 32L) 28L else 56L
}

.mem_used <- function() {
    sum(gc()[, 1] * c(.node_size(), 8))
}


.format_mem_used <- function(digits = 3, ...) {
    x <- .mem_used()
    power <- min(floor(log(abs(x), 1000)), 4)
    if (power < 1) {
        unit <- "B"
    } else {
        unit <- c("kB", "MB", "GB", "TB")[[power]]
        x <- x / (1000^power)
    }

    formatted <- format(signif(x, digits = digits),
        big.mark = ",",
        scientific = FALSE
    )

    paste(formatted, unit)
}

#' @keywords internal
#' @noRd
.log_success <- function(..., .envir = parent.frame(), level = 1) {
    if (getOption("CMEnt.verbose", 0) < level) {
        return(invisible())
    }
    if (level <= length(.CMEnt_log_env$last_step_time)) {
        dur <- .fmt_dur(.CMEnt_log_env$last_step_time[[level]])
    } else {
        .log_warn("No previous step time recorded for level ", level, " to calculate duration.")
        dur <- ""
    }
    msg <- paste0(paste0(..., collapse = ""), dur)
    # if level is equal or greater than 2, report memory usage in MBs as well
    if (level >= 2) {
        msg <- paste0(msg, " [mem: ", .format_mem_used(), "]")
    }
    lead <- .col(cli::symbol$tick, "green")
    message(.format_log_output(msg, lead = lead, level = level))
    invisible()
}


#' @keywords internal
#' @noRd
.log_info <- function(..., .envir = parent.frame(), level = 1) {
    if (getOption("CMEnt.verbose", 0) < level) {
        return(invisible())
    }
    # Suppress output from parallel workers
    if (!is.null(getOption("future.fork.enable")) && getOption("future.fork.enable", TRUE)) {
        if (exists(".Random.seed", envir = .GlobalEnv) && !identical(Sys.getpid(), getOption("future.main.pid", Sys.getpid()))) {
            return(invisible())
        }
    }
    msg <- paste0(..., collapse = "")
    lead <- .col(cli::symbol$info, "blue")
    message(.format_log_output(msg, lead = lead, level = level))
    invisible()
}


#' @keywords internal
#' @noRd
.log_step <- function(..., .envir = parent.frame(), level = 1) {
    if (getOption("CMEnt.verbose", 0) < level) {
        return(invisible())
    }
    .CMEnt_log_env$last_step_time[[level]] <- Sys.time() # nolint
    msg <- paste0(..., collapse = "")
    lead <- .col(cli::symbol$arrow_right, "cyan")
    message(.format_log_output(msg, lead = lead, level = level))
    invisible()
}

#' @keywords internal
#' @noRd
.makeBiocParallelParam <- function(njobs, n_tasks = NULL, progressbar = FALSE) {
    workers <- suppressWarnings(as.integer(njobs))
    if (length(workers) == 0L || is.na(workers) || workers < 1L) {
        stop("njobs must be a positive integer.", call. = FALSE)
    }

    if (!is.null(n_tasks)) {
        n_tasks <- suppressWarnings(as.integer(n_tasks))
        if (length(n_tasks) > 0L && !is.na(n_tasks) && n_tasks > 0L) {
            workers <- min(workers, n_tasks)
        }
    }

    if (workers <= 1L) {
        return(BiocParallel::SerialParam(progressbar = progressbar))
    }

    if (identical(.Platform$OS.type, "unix")) {
        BiocParallel::MulticoreParam(workers = workers, progressbar = progressbar)
    } else {
        BiocParallel::SnowParam(workers = workers, type = "SOCK", progressbar = progressbar)
    }
}

#' @keywords internal
#' @noRd
.processSamplesheet <- function(args,
                                subset = NULL) {
    samplesheet_file <- args$samplesheet
    target_col <- args$target_col
    sample_group_col <- args$sample_group_col
    sample_group_control <- args$sample_group_control
    sample_group_case <- args$sample_group_case
    if (!is.null(sample_group_case)) {
        sample_group_case <- base::strsplit(sample_group_case, ",")[[1]]
    }
    if (!is.null(sample_group_control)) {
        sample_group_control <- base::strsplit(sample_group_control, ",")[[1]]
    }

    if (endsWith(samplesheet_file, ".csv")) {
        samplesheet_file_sep <- ","
    } else if (endsWith(samplesheet_file, ".tsv")) {
        samplesheet_file_sep <- "\t"
    }
    subset_samplesheet <- .readSamplesheet(
        samplesheet_file = samplesheet_file,
        samplesheet_file_sep = samplesheet_file_sep,
        sample_group_col = sample_group_col,
        sample_group_control = sample_group_control,
        sample_group_case = sample_group_case,
        target_col = target_col,
        subset = subset
    )
    ret <- list(samplesheet = subset_samplesheet[, c(sample_group_col, "casecontrol")])
    ret
}


.getTabixCacheDir <- function(output_dir) {
    cache_dir <- if (is.null(output_dir)) tempdir() else output_dir
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    cache_dir
}

createH5file <- function(input_file, output_h5file = tempfile(fileext = ".h5"), dataset_name = "data", select = NULL,
                         chunk_size = 100000, sep = "\t") {
    stopifnot(is.character(input_file), length(input_file) == 1, file.exists(input_file))
    stopifnot(is.character(output_h5file), length(output_h5file) == 1)

    if (file.exists(output_h5file)) file.remove(output_h5file)
    dir.create(dirname(output_h5file), recursive = TRUE, showWarnings = FALSE)
    rhdf5::h5createFile(output_h5file)

    row_offset <- 0L
    initialized <- FALSE
    p <- NULL

    cb <- readr::DataFrameCallback$new(function(df, pos) {
        if (!is.null(select)) {
            df <- df[, select, drop = FALSE]
        }
        if (!initialized) {
            p <<- ncol(df)
            rhdf5::h5createDataset(
                file = output_h5file,
                dataset = dataset_name,
                dims = c(0, p),
                maxdims = c(rhdf5::H5Sunlimited(), p),
                chunk = c(min(chunk_size, max(1L, nrow(df))), p),
                storage.mode = "character",
                level = 7
            )
            rhdf5::h5write(names(df), output_h5file, paste0(dataset_name, "_colnames"))
            initialized <<- TRUE
        }

        n_new <- nrow(df)
        if (n_new == 0L) {
            return(invisible(NULL))
        }

        # Convert whole chunk to character matrix
        # (keeps NA as NA_character_)
        mat <- as.matrix(data.frame(lapply(df, as.character), check.names = FALSE))

        new_total <- row_offset + n_new

        # Extend then write
        rhdf5::h5set_extent(output_h5file, dataset_name, c(new_total, p))

        idx_rows <- (row_offset + 1L):new_total
        rhdf5::h5write(
            mat,
            file = output_h5file,
            name = dataset_name,
            index = list(idx_rows, 1:p)
        )

        row_offset <<- new_total
        invisible(NULL)
    })

    readr::read_tsv_chunked(
        file = input_file,
        callback = cb,
        chunk_size = chunk_size,
        show_col_types = FALSE,
        progress = FALSE
    )

    rhdf5::h5write(row_offset, output_h5file, paste0(dataset_name, "_nrows"))
    invisible(output_h5file)
}

.postProcessRegistry <- function(df, select = NULL, rename = NULL, derive = NULL, indices = NULL) {
    if (!is.null(select)) {
        df <- df[, select, drop = FALSE]
    }
    if (!is.null(rename)) {
        for (name in names(rename)) {
            colnames(df)[colnames(df) == name] <- rename[[name]]
        }
    }
    if (!is.null(derive)) {
        for (new_col in names(derive)) {
            col_info <- derive[[new_col]]
            if (!all(col_info$cols %in% colnames(df))) {
                stop(
                    "Cannot derive column ", new_col, " because not all required columns are present. ",
                    "Required columns: ", paste(col_info$cols, collapse = ", "), ". ",
                    "Available columns: ", paste(colnames(df), collapse = ", ")
                )
            }
            df[[new_col]] <- do.call(col_info$fun, as.data.frame(df[, col_info$cols]))
        }
    }
    if (!is.null(indices)) {
        missing_indices <- setdiff(indices, colnames(df))
        if (length(missing_indices) > 0) {
            stop(
                "Cannot set indices because the following specified index columns are missing: ",
                paste(missing_indices, collapse = ", "), ". ",
                "Available columns: ", paste(colnames(df), collapse = ", ")
            )
        }
        if (length(indices) == 1) {
            rownames(df) <- df[[indices]]
        } else {
            rownames(df) <- do.call(paste, c(df[, indices], sep = ":"))
        }
    }
    df
}

getRegistry <- function(obj, indices = NULL, select = NULL, rename = NULL, derive = NULL,
                        chunk_size = 100000, output_h5file = NULL) {
    if (is.data.frame(obj)) {
        return(.postProcessRegistry(obj, select = select, rename = rename, derive = derive, indices = indices))
    }
    if (is.null(output_h5file)) {
        output_h5file <- tempfile(fileext = ".h5")
    }
    createH5file(
        input_file = obj,
        output_h5file = output_h5file,
        dataset_name = "data",
        select = select,
        chunk_size = chunk_size
    )
    da <- HDF5Array::HDF5Array(output_h5file, "data")
    x <- DelayedDataFrame::DelayedDataFrame(da)
    colnames(x) <- rhdf5::h5read(output_h5file, "data_colnames")
    .postProcessRegistry(x, select = NULL, rename = rename, derive = derive, indices = indices)
}


#' Create Genomic Location Registry from Tabix BED File
#'
#' @description This function creates a Registry from a Tabix-indexed BED file.
#' @param input_tabix Character. Path to the Tabix-indexed BED file.
#' @param output_dir Character. Directory used for temporary or explicit derived files.
#' @param num_rows Integer. Number of rows in the BED file. If NULL, the function will compute it automatically (default: NULL)
#' @param hash Character. Hash string used for deterministic temporary file names.
#' @param chunk_size Integer. Number of rows to process in each chunk for memory efficiency (default: 50000)
#' @return Returns a DelayedDataFrame object
#' @keywords internal
#' @noRd
genomicLocsFromTabix <- function(input_tabix, output_dir = NULL, num_rows = NULL, hash = NULL,
                                 chunk_size = 50000, use_id_as_rownames = FALSE,
                                 chrom_col = "#chrom", start_col = "start",
                                 output_h5file = NULL) { # nolint
    renaming <- c("chr", "start")
    names(renaming) <- c(chrom_col, start_col)
    if (is.null(output_h5file)) {
        output_dir <- .getTabixCacheDir(output_dir)
        if (is.null(hash)) {
            hash <- .getFileHash(input_tabix)
        }
        output_h5file <- file.path(output_dir, paste0("bed_locations_", hash, ".h5"))
    }
    if (!use_id_as_rownames) {
        sorted_locs <- getRegistry(
            input_tabix,
            select = c(chrom_col, start_col),
            rename = renaming,
            derive = list(
                index = list(
                    cols = c("chr", "start"),
                    fun = function(chr, start) paste0(chr, ":", start)
                )
            ),
            indices = "index",
            chunk_size = chunk_size,
            output_h5file = output_h5file
        )
    } else {
        sorted_locs <- getRegistry(
            input_tabix,
            select = c(chrom_col, start_col, "end", "id"),
            rename = renaming,
            indices = "id",
            chunk_size = chunk_size,
            output_h5file = output_h5file
        )
    }
    sorted_locs
}


#' Read and Process Custom Methylation BED Data
#'
#' @description Reads methylation data from a custom BED file format, converts it to
#' a tabix-indexed format for efficient random access, and creates genomic location
#' indices. This function is designed to handle custom methylation array data or
#' sequencing-based methylation data in BED format, making it compatible with the
#' CMEnt workflow.
#'
#' @param bed_file Character. Path to the input BED file containing methylation data.
#'   The file should have chromosome and position columns, plus sample columns with
#'   methylation values. Can be gzipped (default: NULL)
#' @param pheno Data frame. Phenotype data with sample IDs as rownames. Only samples
#'   present in both the pheno rownames and BED file header will be processed
#' @param genome Character. Genome version to use (e.g., "hg38", "hg19", "hs1") (default: "hg38")
#' @param chrom_col Character. Name of the chromosome column in the BED file
#'   (default: "#chrom")
#' @param start_col Character. Name of the start position column in the BED file
#'   (default: "start")
#' @param output_dir Character. Directory for caching processed files. If NULL, uses
#'   a temporary working directory unless `output_prefix` is provided (default: NULL)
#' @param chunk_size Integer. Number of rows to process in each chunk for memory
#'   efficiency (default: 50000)
#' @param output_prefix Character. Optional prefix used to persist derived BED/tabix
#'   artifacts next to analysis outputs.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item tabix_file: Character path to the created tabix-indexed BED file
#'   \item locations: Disk-backed genomic location registry
#' }
#'
#' @details
#' The function performs the following workflow:
#' \enumerate{
#'   \item Validates that tabix and bgzip are available in the system PATH
#'   \item Checks the BED file header for required columns and sample IDs
#'   \item Processes the BED file in chunks to minimize memory usage
#'   \item Normalizes the BED format with standard BED6 columns (#chrom, start, end, id, score, strand)
#'   \item Converts chromosomes to integer factors for efficient sorting
#'   \item Creates a tabix-indexed compressed file for fast random access
#'   \item Persists derived artifacts under `output_prefix` when provided
#' }
#'
#' @section Requirements:
#' This function requires tabix and bgzip command-line tools to be installed and
#' available in the system PATH. These tools are part of the HTSlib/samtools suite.
#'
#' @section Memory Management:
#' The function uses chunk-based processing to handle large BED files without
#' loading the entire dataset into memory. The genomic locations are stored in
#' a Registry object that can exceed available RAM by using disk-backed
#' storage.
#'
#' @examples
#' # Create a simple phenotype data frame
#' pheno <- data.frame(
#'     sample_group = c("case", "case", "control", "control"),
#'     row.names = c("Sample1", "Sample2", "Sample3", "Sample4")
#' )
#'
#' # Process a custom BED file
#' result <- readCustomMethylationBedData(
#'     bed_file = "custom_methylation.bed.gz",
#'     genome = "hg38",
#'     pheno = pheno
#' )
#'
#' # Use custom column names
#' result <- readCustomMethylationBedData(
#'     bed_file = "custom_methylation.bed",
#'     pheno = pheno,
#'     genome = "hg38",
#'     chrom_col = "chromosome",
#'     start_col = "position"
#' )
#'
#' # Specify custom output directory
#' result <- readCustomMethylationBedData(
#'     bed_file = "custom_methylation.bed.gz",
#'     pheno = pheno,
#'     genome = "hg38",
#'     output_dir = "/path/to/output"
#' )
#'
#' # Access the processed files
#' tabix_file <- result$tabix_file
#' locations <- result$locations
#'
#' @seealso
#' \code{\link{convertBetaToTabix}} for converting standard beta files to tabix format
#' \code{\link{getBetaHandler}} for creating a BetaHandler object from processed files
#'
#' @export
readCustomMethylationBedData <- function(bed_file, pheno, genome = "hg38", chrom_col = "#chrom",
                                         start_col = "start", output_dir = NULL, chunk_size = 50000,
                                         output_prefix = NULL) {
    tabix_available <- tryCatch(
        {
            system2("which", "tabix", stdout = FALSE, stderr = FALSE)
            system2("which", "bgzip", stdout = FALSE, stderr = FALSE)
            TRUE
        },
        error = function(e) FALSE,
        warning = function(w) FALSE
    )

    if (!tabix_available) {
        stop("tabix/bgzip not found in PATH. Cannot process BED file.")
    }

    cache_dir <- .getTabixCacheDir(output_dir)
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    hash <- .getFileHash(bed_file)
    normalized_bed_file <- file.path(tempdir(), paste0("bed_", hash, ".tsv"))

    # Read BED file header
    bed_header <- base::strsplit(readLines(bed_file, n = 1), "\t")[[1]]
    # Ensure required columns are present
    required_cols <- c(chrom_col, start_col)
    missing_cols <- setdiff(required_cols, bed_header)
    if (length(missing_cols) > 0) {
        stop("Missing required columns in BED file: ", paste(missing_cols, collapse = ", "), ". Available columns: ", paste(bed_header, collapse = ", "))
    }
    sample_ids <- rownames(pheno)
    existing_ids <- intersect(sample_ids, bed_header)

    # Map existing sample IDs to BED file rows
    id_mapping <- match(existing_ids, bed_header)
    id_mapping <- id_mapping[!is.na(id_mapping)]
    if (length(id_mapping) == 0) {
        stop("None of the provided sample IDs were found in the BED file header.")
    }
    if (length(id_mapping) < length(sample_ids)) {
        missing_ids <- setdiff(sample_ids, bed_header)
        .log_warn(length(missing_ids), " out of ", length(sample_ids), " sample IDs were not found in the BED file header and will be ignored. The IDs are: ", paste(missing_ids, collapse = ", "))
    }

    # Quickly read number of rows in BED file
    tmp_con <- if (endsWith(bed_file, ".gz")) gzfile(bed_file, "r") else file(bed_file, "r")
    num_rows <- sum(sapply(readLines(tmp_con), function(x) nchar(x) > 0)) - 1
    close(tmp_con)
    .log_info("Processing BED file with ", num_rows, " rows and ", length(existing_ids), " matching sample IDs.", level = 2)


    # Read chunks of the BED file to minimize memory usage
    con <- if (endsWith(bed_file, ".gz")) gzfile(bed_file, "r") else file(bed_file, "r")
    # Write normalized header to new BED file
    norm_bed_header <- c("#chrom", "start", "end", "id", "score", "strand", existing_ids)
    writeLines(paste(norm_bed_header, collapse = "\t"), normalized_bed_file)
    # Skip header line
    readLines(con, n = 1)
    count <- 0
    while (length(chunk <- readLines(con, n = chunk_size)) > 0) {
        bed_data <- data.table::fread(paste(chunk, collapse = "\n"), sep = "\t", header = FALSE, data.table = FALSE)
        colnames(bed_data) <- bed_header
        bed_data$chr <- bed_data[[chrom_col]]
        bed_data$start <- as.integer(bed_data[[start_col]])
        bed_data$end <- bed_data$start + 1
        bed_data$score <- "."
        bed_data$id <- seq(count + 1, count + nrow(bed_data))
        bed_data$strand <- "*"

        # Write normalized BED data
        bed_subset <- bed_data[, c("chr", "start", "end", "id", "score", "strand", existing_ids), drop = FALSE]
        data.table::fwrite(
            bed_subset,
            file = normalized_bed_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE
        )
        count <- count + nrow(bed_data)
    }
    close(con)

    hash <- .getFileHash(bed_file)
    tabix_file_path <- .getDerivedOutputPath(output_prefix, ".input_beta.tabix.bed.gz")
    if (is.null(tabix_file_path)) {
        tabix_file_path <- file.path(cache_dir, paste0("bed_beta_", hash, ".bed.gz"))
    }
    # Convert to tabix
    convertBetaToTabix(
        .bed_file = normalized_bed_file,
        output_file = tabix_file_path,
        chunk_size = chunk_size,
        njobs = 1,
        output_prefix = output_prefix
    )
    locations_h5_file <- .getDerivedOutputPath(output_prefix, ".input_beta.locations.h5")
    locations <- genomicLocsFromTabix(
        tabix_file_path,
        num_rows = num_rows,
        hash = hash,
        output_dir = cache_dir,
        chunk_size = chunk_size,
        output_h5file = locations_h5_file
    )


    list(
        tabix_file = tabix_file_path,
        locations = locations
    )
}


#' Convert Beta File to Tabix-Indexed Format
#'
#' @description Converts a methylation beta values file to a tabix-indexed BED format
#' for faster random access during DMR analysis. The function uses a memory-efficient
#' chunk-based approach to handle large files and can persist the derived tabix file
#' next to analysis outputs when `output_prefix` is supplied.
#'
#' @param beta_file Character. Path to the input beta values file
#' @param sorted_locs Data frame with genomic locations containing 'chr' and 'start' columns.
#'   If NULL, will be retrieved automatically using getSortedGenomicLocs() (default: NULL)
#' @param array Character. Array platform type. Only used if sorted_locs is NULL (default: "450K")
#' @param genome Character. Genome version. Only used if  sorted_locs is NULL (default: "hg38")
#' @param locations_file Character. Optional path to an explicit genomic locations file passed through to `getSortedGenomicLocs()`.
#' @param output_file Character. Path for the output tabix file. If NULL, a temporary
#'   file is used unless `output_prefix` is supplied.
#' @param chunk_size Integer. Number of rows to process in each chunk (default: 50000)
#' @param njobs Integer. Number of parallel jobs for sorting (default: 1)
#' @param .bed_file Character. Internal precomputed BED path used to skip beta-to-BED conversion.
#' @param output_prefix Character. Optional prefix used to persist derived tabix
#'   artifacts next to analysis outputs.
#'
#' @return Character. Path to the created tabix file, or NULL if conversion failed
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks if tabix and bgzip tools are available in the system PATH
#'   \item Processes the beta file in chunks (50,000 rows at a time) to minimize memory usage
#'   \item Converts beta values to BED format with genomic coordinates
#'   \item Sorts, compresses (bgzip), and indexes (tabix) the file
#'   \item Persists the derived file if an explicit output path or `output_prefix` is provided
#' }
#'
#' @examples
#' # Convert a beta file to tabix format
#' tabix_file <- convertBetaToTabix(
#'     beta_file = "methylation_beta.txt",
#'     array = "450K"
#' )
#'
#' # Use custom output location
#' tabix_file <- convertBetaToTabix(
#'     beta_file = "methylation_beta.txt",
#'     output_file = "my_custom_location.bed.gz",
#'     array = "EPIC"
#' )
#'
#' @export
convertBetaToTabix <- function(beta_file,
                               sorted_locs = NULL,
                               array = c("450K", "27K", "EPIC", "EPICv2"),
                               genome = "hg38",
                               locations_file = NULL,
                               output_file = NULL,
                               chunk_size = 50000,
                               njobs = 1,
                               .bed_file = NULL,
                               output_prefix = NULL) {
    # Check if tabix/bgzip are available
    tabix_available <- tryCatch(
        {
            system2("which", "tabix", stdout = FALSE, stderr = FALSE)
            system2("which", "bgzip", stdout = FALSE, stderr = FALSE)
            TRUE
        },
        error = function(e) FALSE,
        warning = function(w) FALSE
    )

    if (!tabix_available) {
        .log_warn("tabix/bgzip not found in PATH. Skipping tabix conversion.")
        return(NULL)
    }

    if (is.null(output_file)) {
        output_file <- .getDerivedOutputPath(output_prefix, ".input_beta.tabix.bed.gz")
        if (is.null(output_file)) {
            beta_hash <- .getFileHash(beta_file)
            output_file <- file.path(tempdir(), paste0("beta_", beta_hash, ".bed.gz"))
        }
    }
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

    .log_step("Converting beta file to tabix format...", level = 1)

    tryCatch(
        {
            if (is.null(.bed_file)) {
                array <- strex::match_arg(array, ignore_case = TRUE)
                # Get sorted locations if not provided
                if (is.null(sorted_locs)) {
                    sorted_locs <- getSortedGenomicLocs(array = array, genome = genome)
                }
                # Read header to get column names
                .log_step("Reading beta file header...", level = 2)

                header_conn <- if (endsWith(beta_file, ".gz")) gzfile(beta_file, "r") else file(beta_file, "r")
                header_line <- readLines(header_conn, n = 1)
                close(header_conn)
                col_names <- base::strsplit(header_line, "\t")[[1]]

                # Get total number of rows for progress tracking
                .log_step("Counting rows in beta file...", level = 2)

                # Count lines efficiently (cross-platform)
                if (endsWith(beta_file, ".gz")) {
                    conn <- gzfile(beta_file, "r")
                } else {
                    conn <- file(beta_file, "r")
                }
                n_lines <- length(readLines(conn))
                close(conn)
                n_rows <- n_lines - 1 # Exclude header

                .log_info("Processing ", n_rows, " site sites...", level = 2)

                # Create temporary BED file for writing chunks
                temp_bed <- tempfile(fileext = ".bed")
                withr::defer(unlink(temp_bed))

                # Write header to temp BED file with 6 mandatory BED columns
                bed_header <- c("#chrom", "start", "end", "id", "score", "strand", col_names[-1])
                writeLines(paste(bed_header, collapse = "\t"), temp_bed)

                # Process file in chunks to avoid memory issues

                skip_rows <- 1 # Start after header
                rows_processed <- 0

                while (rows_processed < n_rows) {
                    .log_info("Processing rows ", rows_processed + 1, " to ",
                        min(rows_processed + chunk_size, n_rows), "...",
                        level = 3
                    )

                    # Read chunk
                    chunk_data <- data.table::fread(
                        beta_file,
                        header = FALSE,
                        skip = skip_rows,
                        nrows = chunk_size,
                        data.table = FALSE,
                        showProgress = FALSE
                    )

                    if (nrow(chunk_data) == 0) break

                    # Set column names
                    colnames(chunk_data) <- col_names

                    site_ids <- chunk_data[[1]]

                    # Match with genomic locations
                    common_sites <- intersect(site_ids, rownames(sorted_locs))

                    if (length(common_sites) > 0) {
                        # Create BED format for this chunk with 6 mandatory columns
                        bed_chunk <- as.data.frame(sorted_locs[common_sites, c("chr", "start"), drop = FALSE])
                        # Tabix requires plain integer coordinates; fwrite may emit
                        # scientific notation for doubles like 45000000 ("4.5e+07"),
                        # which tabix then parses as 4 and rejects as an invalid BED interval.
                        bed_chunk$start <- as.integer(round(bed_chunk$start))
                        bed_chunk$end <- bed_chunk$start + 1L
                        bed_chunk$id <- rownames(bed_chunk)
                        bed_chunk$score <- 0
                        bed_chunk$strand <- "*"
                        bed_chunk <- bed_chunk[, c("chr", "start", "end", "id", "score", "strand")]

                        # Add beta values as additional columns
                        beta_subset <- chunk_data[match(common_sites, site_ids), -1, drop = FALSE]
                        bed_chunk <- cbind(bed_chunk, beta_subset)

                        # Append to temp BED file
                        data.table::fwrite(
                            bed_chunk,
                            file = temp_bed,
                            sep = "\t",
                            quote = FALSE,
                            row.names = FALSE,
                            col.names = FALSE,
                            append = TRUE
                        )
                    }

                    rows_processed <- rows_processed + nrow(chunk_data)
                    skip_rows <- skip_rows + nrow(chunk_data)
                }

                .log_success("Processed ", rows_processed, " rows", level = 2)

                # Check if any data was written
                if (file.info(temp_bed)$size <= length(paste(bed_header, collapse = "\t")) + 1) {
                    .log_warn("No common sites found between beta file and genomic locations")
                    return(NULL)
                }
            } else {
                temp_bed <- .bed_file
                .log_info("Using provided BED file for tabix conversion: ", temp_bed, level = 2)
            }

            # Sort, compress with bgzip, and index with tabix
            .log_step("Sorting BED file...", level = 2)
            temp_sorted <- tempfile(fileext = ".bed")

            # Platform-specific sorting
            is_windows <- .Platform$OS.type == "windows"

            if (is_windows) {
                # Windows: use external sorting for large files
                # Read and write header first
                header_line <- readLines(temp_bed, n = 1)
                writeLines(header_line, temp_sorted)

                # Process file in chunks, sort each chunk, write to temp files
                chunk_files <- character()
                skip <- 1
                chunk_num <- 0
                chunk_sort_size <- 100000 # rows per chunk

                repeat {
                    chunk <- data.table::fread(temp_bed,
                        skip = skip, nrows = chunk_sort_size,
                        header = FALSE, data.table = TRUE
                    )
                    if (nrow(chunk) == 0) break

                    # Sort chunk by chr (col 1) then position (col 2)
                    data.table::setorderv(chunk, cols = c(1, 2), order = c(1, 1))

                    # Write sorted chunk to temp file
                    chunk_num <- chunk_num + 1
                    chunk_file <- tempfile(fileext = paste0("_chunk", chunk_num, ".txt"))
                    data.table::fwrite(chunk, chunk_file,
                        sep = "\t", quote = FALSE,
                        row.names = FALSE, col.names = FALSE
                    )
                    chunk_files <- c(chunk_files, chunk_file)

                    skip <- skip + nrow(chunk)
                }

                # K-way merge of sorted chunks
                if (length(chunk_files) > 0) {
                    # Helper function to compare two BED lines (chr:start)
                    compare_bed_lines <- function(line1, line2) {
                        parts1 <- base::strsplit(line1, "\t", fixed = TRUE)[[1]]
                        parts2 <- base::strsplit(line2, "\t", fixed = TRUE)[[1]]

                        chr1 <- parts1[1]
                        chr2 <- parts2[1]

                        # Compare chromosomes
                        if (chr1 != chr2) {
                            return(chr1 < chr2)
                        }

                        # Same chromosome, compare positions
                        pos1 <- as.numeric(parts1[2])
                        pos2 <- as.numeric(parts2[2])
                        pos1 < pos2
                    }

                    # Open all chunk files and read first line from each
                    connections <- lapply(chunk_files, function(f) file(f, "r"))

                    heap <- list()

                    for (i in seq_along(connections)) {
                        line <- readLines(connections[[i]], n = 1)
                        if (length(line) > 0) {
                            heap[[length(heap) + 1]] <- list(line = line, chunk_idx = i)
                        }
                    }

                    out_conn <- file(temp_sorted, "a")

                    # Merge: repeatedly extract minimum, write it, and refill from same chunk
                    while (length(heap) > 0) {
                        # Find minimum element in heap
                        min_idx <- 1
                        for (i in seq_along(heap)) {
                            if (i > 1 && compare_bed_lines(heap[[i]]$line, heap[[min_idx]]$line)) {
                                min_idx <- i
                            }
                        }

                        # Write minimum line
                        writeLines(heap[[min_idx]]$line, out_conn)

                        # Read next line from the same chunk
                        chunk_idx <- heap[[min_idx]]$chunk_idx
                        next_line <- readLines(connections[[chunk_idx]], n = 1)

                        if (length(next_line) > 0) {
                            # Replace with new line from same chunk
                            heap[[min_idx]]$line <- next_line
                        } else {
                            # This chunk is exhausted, remove from heap
                            heap[[min_idx]] <- NULL
                        }
                    }

                    close(out_conn)
                    # Close all chunk file connections
                    lapply(connections, close)
                    # Clean up chunk files
                    lapply(chunk_files, unlink)
                }
            } else {
                # Unix/Linux/Mac: use efficient system sort
                sort_cmd <- sprintf(
                    "(head -n 1 %s && tail -n +2 %s | sort --parallel=%d -V -k1,1 -k2,2n) > %s",
                    shQuote(temp_bed), shQuote(temp_bed), njobs, shQuote(temp_sorted)
                )
                system(sort_cmd)
            }
            unlink(temp_bed)
            # Compress with bgzip
            .log_step("Compressing with bgzip...", level = 2)
            .log_info("Expected output compressed file: ", output_file, level = 3)
            error_file <- tempfile(fileext = ".log")
            status_code <- system2("bgzip", args = c("-c", shQuote(temp_sorted)), stdout = output_file, stderr = error_file)
            unlink(temp_sorted)
            if (status_code != 0) {
                con <- file(error_file, "r")
                error <- readLines(con)
                close(con)
                stop("bgzip compression failed with exit code ", status_code, ": ", error)
            }
            # Index with tabix
            .log_step("Creating tabix index...", level = 2)
            error_file <- tempfile(fileext = ".log")

            status_code <- system2("tabix", args = c("-f", "-p", "bed", shQuote(output_file)), stderr = error_file, stdout = NULL)
            if (status_code != 0) {
                con <- file(error_file, "r")
                error <- readLines(con)
                close(con)
                stop("tabix indexing failed with exit code ", status_code, ": ", error)
            }

            # Clean up temp files

            if (file.exists(output_file) && file.exists(paste0(output_file, ".tbi"))) {
                .log_success("Tabix file created: ", output_file, level = 1)
                return(output_file)
            } else {
                .log_warn("Failed to create tabix index")
                NULL
            }
        },
        error = function(e) {
            .log_warn("Error converting to tabix: ", e$message)
            NULL
        }
    )
}


#' Sort Beta File by Genomic Coordinates
#'
#' @description This helper function sorts a methylation beta values file by genomic coordinates
#' (chromosome and position) as required by the findDMRsFromSeeds function. The function reads
#' the beta file, sorts the site sites according to their genomic positions using array annotation,
#' and writes the sorted data to a new file.
#'
#' @param beta_file Character. Path to the input beta values file to be sorted
#' @param output_file Character. Path for the output sorted beta file (default: adds "_sorted" suffix)
#' @param array Character. Array platform type (default: "450K")
#' @param genome Character. Genome version (default: "hg38")
#' @param genomic_locs Data frame. Optional pre-computed genomic locations. If NULL, locations will be retrieved automatically (default: NULL)
#' @param overwrite Logical. Whether to overwrite existing output file (default: FALSE)
#'
#' @return Character. Path to the sorted output file
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Reads the beta values file
#'   \item Loads the appropriate array annotation (450K or EPIC)
#'   \item Sorts site sites by genomic coordinates (chr:start)
#'   \item Writes the sorted data to a new file
#'   \item Validates that the output is properly sorted
#' }
#' @note If you want to convert to tabix, consider using the convertBetaToTabix function instead directly, sorting is done internally.
#'
#' @examples
#' # Sort a beta file for 450K array
#' sorted_file <- sortBetaFileByCoordinates(
#'     beta_file = "unsorted_beta.txt",
#'     output_file = "sorted_beta.txt",
#'     array = "450K"
#' )
#'
#' # Sort an EPIC array beta file with default output name
#' sorted_file <- sortBetaFileByCoordinates(
#'     beta_file = "epic_beta.txt",
#'     array = "EPIC"
#' )
#'
#' @export
sortBetaFileByCoordinates <- function(beta_file,
                                      output_file = NULL,
                                      array = c("450K", "27K", "EPIC", "EPICv2"),
                                      genome = "hg38",
                                      genomic_locs = NULL,
                                      overwrite = FALSE) {
    # Validate inputs
    if (!file.exists(beta_file)) {
        stop("Beta file does not exist: ", beta_file)
    }

    # Set default output file name
    if (is.null(output_file)) {
        file_ext <- tools::file_ext(beta_file)
        file_base <- tools::file_path_sans_ext(beta_file)
        output_file <- paste0(file_base, "_sorted.", file_ext)
    }

    # Check if output file exists
    if (file.exists(output_file) && !overwrite) {
        stop(
            "Output file already exists: ", output_file,
            ". Set overwrite=TRUE to overwrite or choose a different output_file name."
        )
    }

    .log_step("Reading beta file", beta_file, level = 2)
    # Read the beta file
    beta_data <- data.table::fread(beta_file, header = TRUE, data.table = FALSE, showProgress = getOption("CMEnt.verbose", 0) > 1)

    # Get row names (site IDs) from first column
    site_ids <- beta_data[[1]]
    beta_values <- beta_data[, -1, drop = FALSE]
    rownames(beta_values) <- site_ids
    .log_success("Beta loaded: ", nrow(beta_values), " sites across ", ncol(beta_values), " samples", level = 2)

    sorted_locs <- genomic_locs
    if (is.null(sorted_locs)) {
        array <- strex::match_arg(array, ignore_case = TRUE)
        sorted_locs <- getSortedGenomicLocs(array = array, genome = genome)
    }


    # Find sites that are present in both the beta file and array annotation
    common_sites <- intersect(site_ids, rownames(sorted_locs))
    missing_from_annotation <- setdiff(site_ids, rownames(sorted_locs))
    if (length(missing_from_annotation) > 0) {
        stop(
            "Found ", length(missing_from_annotation), " site sites in beta file that are not in ",
            array, " annotation. First 5 missing: ", paste(head(missing_from_annotation, 5), collapse = ", ")
        )
    }

    missing_from_beta <- setdiff(rownames(sorted_locs), site_ids)
    if (length(missing_from_beta) > 0) {
        .log_info("Note: ", length(missing_from_beta), " sites in ", array, " annotation are missing from beta file", level = 2)
    }

    final_order <- rownames(sorted_locs)[rownames(sorted_locs) %in% common_sites]

    # Reorder beta values
    sorted_beta_values <- beta_values[final_order, , drop = FALSE]

    # Prepare output data frame
    output_data <- data.frame(
        ID = rownames(sorted_beta_values),
        sorted_beta_values,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    .log_step("Writing sorted beta file", output_file, level = 2)

    # Write sorted file
    data.table::fwrite(
        output_data,
        file = output_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
    )

    output_file
}

.prepareCovariateModel <- function(pheno, covariates = NULL) {
    if (is.null(covariates) || length(covariates) == 0L || is.null(pheno)) {
        return(NULL)
    }
    if (!all(covariates %in% colnames(pheno))) {
        stop("Not all covariates are present in pheno.")
    }
    xc <- as.data.frame(pheno[, covariates, drop = FALSE])
    xc <- data.frame(`(Intercept)` = 1, xc, check.names = FALSE)
    for (col in colnames(xc)) {
        if (is.character(xc[[col]]) || is.factor(xc[[col]])) {
            xc[[col]] <- as.numeric(as.factor(xc[[col]]))
        }
    }
    xc <- as.matrix(xc)
    storage.mode(xc) <- "double"
    xtx_inv <- tryCatch(solve(crossprod(xc)), error = function(e) NULL)
    pseudo_solution <- if (is.null(xtx_inv)) NULL else xtx_inv %*% t(xc)
    list(
        covariate_matrix = xc,
        t_covariate_matrix = t(xc),
        pseudo_solution = pseudo_solution,
        is_singular = is.null(xtx_inv)
    )
}

.remove_confounder_effect <- function(signal, covariate_matrix, pseudo_solution = NULL, t_covariate_matrix = NULL) {
    if (is.null(covariate_matrix) || ncol(covariate_matrix) == 0L) {
        return(signal)
    }
    if (is.null(pseudo_solution)) {
        xtx <- crossprod(covariate_matrix)
        xtx_inv <- tryCatch(solve(xtx), error = function(e) NULL)
        if (is.null(xtx_inv)) {
            return(signal)
        }
        pseudo_solution <- xtx_inv %*% t(covariate_matrix)
    }
    if (is.null(t_covariate_matrix)) {
        t_covariate_matrix <- t(covariate_matrix)
    }
    effect <- t(pseudo_solution %*% t(signal))
    fitted <- effect %*% t_covariate_matrix
    if (inherits(signal, "DelayedArray")) {
        fitted <- DelayedArray::DelayedArray(fitted)
    }
    signal - fitted
}

.transformBeta <- function(beta, pheno, covariates = NULL, covariate_model = NULL) {
    if (inherits(beta, "DelayedDataFrame")) {
        beta <- DelayedArray::DelayedArray(beta)
    }
    m_values <- log2(beta / (1 - beta + 1e-6) + 1e-6)
    if (is.null(covariate_model)) {
        covariate_model <- .prepareCovariateModel(pheno = pheno, covariates = covariates)
    }
    if (!is.null(covariate_model)) {
        if (isTRUE(covariate_model$is_singular)) {
            return(m_values)
        }
        m_values <- .remove_confounder_effect(
            m_values,
            covariate_matrix = covariate_model$covariate_matrix,
            pseudo_solution = covariate_model$pseudo_solution,
            t_covariate_matrix = covariate_model$t_covariate_matrix
        )
    }
    m_values[is.na(m_values)] <- 0
    m_values
}


.liftOverFromGenomeToGenome <- function(granges, from_genome, to_genome) {
    if (from_genome == to_genome) {
        return(granges)
    }
    cache_dir <- getOption(
        "CMEnt.annotation_cache_dir",
        .getOSCacheDir(file.path("R", "CMEnt", "annotation_cache"))
    )
    bfc <- .getBiocFileCache(cache_dir)
    chain_name <- paste0(from_genome, "To", stringr::str_to_title(to_genome), ".over.chain")
    chain_file <- .getBiocFileCachePath(
        bfc,
        rname = paste0("liftOver_", chain_name),
        ext = ".over.chain"
    )
    if (!file.exists(chain_file)) {
        temp_gz <- tempfile(fileext = ".over.chain.gz")
        on.exit(unlink(temp_gz), add = TRUE)
        .downloadFirstAvailable(
            urls = paste0(
                "https://hgdownload.soe.ucsc.edu/goldenPath/",
                from_genome, "/liftOver/", chain_name, ".gz"
            ),
            destfile = temp_gz,
            mode = "wb"
        )
        R.utils::gunzip(
            filename = temp_gz,
            destname = chain_file,
            overwrite = TRUE,
            remove = FALSE
        )
    }
    chain <- rtracklayer::import.chain(chain_file)
    lifted <- rtracklayer::liftOver(granges, chain)
    lifted_unlisted <- unlist(lifted)
    if (length(lifted_unlisted) > 0) {
        GenomeInfoDb::genome(lifted_unlisted) <- to_genome
    }
    lifted_unlisted
}


.loadAnnotationLocations <- function(pkg_name, source_genome) {
    pkg_dir <- find.package(pkg_name)
    data_env <- new.env(parent = emptyenv())
    lazyLoad(file.path(pkg_dir, "data", "Rdata"), envir = data_env)

    if (!exists("Locations", envir = data_env, inherits = FALSE)) {
        stop("Annotation package '", pkg_name, "' does not contain a 'Locations' dataset.")
    }

    locs <- get("Locations", envir = data_env, inherits = FALSE)
    locs <- as.data.frame(locs, stringsAsFactors = FALSE)

    required_cols <- c("chr", "pos")
    missing_cols <- setdiff(required_cols, colnames(locs))
    if (length(missing_cols) > 0) {
        stop(
            "Annotation package '", pkg_name, "' is missing required columns in 'Locations': ",
            paste(missing_cols, collapse = ", ")
        )
    }

    gr <- GenomicRanges::GRanges(
        seqnames = locs$chr,
        ranges = IRanges::IRanges(start = locs$pos, width = 1)
    )
    names(gr) <- rownames(locs)
    GenomeInfoDb::genome(gr) <- source_genome
    gr
}


.isPackageInstalled <- function(pkg_name) {
    nzchar(system.file(package = pkg_name))
}


.packageInstallSpec <- function(pkg_name) {
    if (identical(pkg_name, "IlluminaMouseMethylationanno.12.v1.mm10")) {
        return(list(type = "github", value = "chiaraherzog/IlluminaMouseMethylationanno.12.v1.mm10"))
    }
    if (pkg_name %in% c("optparse", "devtools")) {
        return(list(type = "cran", value = pkg_name))
    }
    list(type = "bioc", value = pkg_name)
}


.formatInstallInstructions <- function(pkg_names) {
    if (length(pkg_names) == 0L) {
        return(character(0))
    }

    specs <- lapply(pkg_names, .packageInstallSpec)
    spec_types <- vapply(specs, `[[`, character(1), "type")
    spec_values <- vapply(specs, `[[`, character(1), "value")
    lines <- character(0)

    cran_pkgs <- unique(spec_values[spec_types == "cran"])
    bioc_pkgs <- unique(spec_values[spec_types == "bioc"])
    github_repos <- unique(spec_values[spec_types == "github"])

    if (length(cran_pkgs) > 0L) {
        lines <- c(
            lines,
            paste0(
                "install.packages(c(",
                paste(sprintf("\"%s\"", cran_pkgs), collapse = ", "),
                "))"
            )
        )
    }

    if (length(bioc_pkgs) > 0L) {
        lines <- c(
            lines,
            "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")",
            paste0(
                "BiocManager::install(c(",
                paste(sprintf("\"%s\"", bioc_pkgs), collapse = ", "),
                "))"
            )
        )
    }

    if (length(github_repos) > 0L) {
        if (!"devtools" %in% cran_pkgs) {
            lines <- c(lines, "install.packages(\"devtools\")")
        }
        lines <- c(
            lines,
            vapply(github_repos, function(repo) {
                paste0("devtools::install_github(\"", repo, "\")")
            }, character(1))
        )
    }

    lines
}


.assertPackagesInstalled <- function(pkg_names, context, reason = NULL) {
    pkg_names <- unique(pkg_names[!is.na(pkg_names) & nzchar(pkg_names)])
    if (length(pkg_names) == 0L) {
        return(invisible(character(0)))
    }

    missing_pkgs <- pkg_names[!vapply(pkg_names, .isPackageInstalled, logical(1))]
    if (length(missing_pkgs) == 0L) {
        return(invisible(pkg_names))
    }

    msg <- c(
        paste0(
            context,
            " requires the following package",
            if (length(missing_pkgs) > 1L) "s" else "",
            ": ",
            paste(sprintf("'%s'", missing_pkgs), collapse = ", "),
            "."
        )
    )
    if (!is.null(reason) && nzchar(reason)) {
        msg <- c(msg, reason)
    }
    msg <- c(msg, "Install with:", .formatInstallInstructions(missing_pkgs))

    stop(paste(msg, collapse = "\n"), call. = FALSE)
}


.makeDependencyRequirements <- function(pkg_names, reason) {
    pkg_names <- unique(pkg_names[!is.na(pkg_names) & nzchar(pkg_names)])
    if (length(pkg_names) == 0L) {
        return(data.frame(pkg_name = character(0), reason = character(0), stringsAsFactors = FALSE))
    }
    data.frame(
        pkg_name = pkg_names,
        reason = rep(reason, length(pkg_names)),
        stringsAsFactors = FALSE
    )
}


.combineDependencyRequirements <- function(...) {
    parts <- list(...)
    parts <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, parts)
    if (length(parts) == 0L) {
        return(data.frame(pkg_name = character(0), reason = character(0), stringsAsFactors = FALSE))
    }
    combined <- do.call(rbind, parts)
    combined[!duplicated(combined), , drop = FALSE]
}


.assertDependencyRequirements <- function(requirements, context) {
    if (!is.data.frame(requirements) || nrow(requirements) == 0L) {
        return(invisible(character(0)))
    }

    missing_mask <- !vapply(requirements$pkg_name, .isPackageInstalled, logical(1))
    if (!any(missing_mask)) {
        return(invisible(requirements$pkg_name))
    }

    missing_reqs <- requirements[missing_mask, , drop = FALSE]
    missing_pkgs <- unique(missing_reqs$pkg_name)
    missing_reasons <- unique(missing_reqs$reason[nzchar(missing_reqs$reason)])
    msg <- c(
        paste0(
            context,
            " requires the following package",
            if (length(missing_pkgs) > 1L) "s" else "",
            ": ",
            paste(sprintf("'%s'", missing_pkgs), collapse = ", "),
            "."
        )
    )
    if (length(missing_reasons) > 0L) {
        msg <- c(msg, "Required by this run:", paste0("- ", missing_reasons))
    }
    msg <- c(msg, "Install with:", .formatInstallInstructions(missing_pkgs))

    stop(paste(msg, collapse = "\n"), call. = FALSE)
}


.getArrayAnnotationPackage <- function(array, genome) {
    array <- strex::match_arg(array, choices = c("450K", "27K", "EPIC", "EPICv2", "Mouse"), ignore_case = TRUE)
    genome <- strex::match_arg(genome, choices = c("hg38", "hg19", "hs1", "mm10", "mm39"), ignore_case = TRUE)

    if (genome %in% c("hg19", "hg38", "hs1")) {
        return(switch(tolower(array),
            "450k" = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
            "epic" = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
            "epicv2" = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
            "27k" = "IlluminaHumanMethylation27kanno.ilmn12.hg19",
            stop(
                paste0(
                    "Incorrect array and genome combination was provided. ",
                    "For hg19, hg38, and hs1, ('450K','EPIC','EPICv2','27K') arrays are supported."
                ),
                call. = FALSE
            )
        ))
    }

    if (genome %in% c("mm10", "mm39")) {
        if (!identical(array, "Mouse")) {
            stop(
                "Incorrect array and genome combination was provided. For mm10 and mm39 only 'Mouse' array is supported.",
                call. = FALSE
            )
        }
        return("IlluminaMouseMethylationanno.12.v1.mm10")
    }

    stop("Unsupported array/genome combination: ", array, "/", genome, call. = FALSE)
}


.assertArrayAnnotationPackageInstalled <- function(array, genome, context) {
    pkg_name <- .getArrayAnnotationPackage(array = array, genome = genome)
    .assertDependencyRequirements(
        requirements = .makeDependencyRequirements(
            pkg_names = pkg_name,
            reason = paste0(
                "The requested array/genome combination ('", array, "', '", genome,
                "') needs this annotation package to resolve probe coordinates."
            )
        ),
        context = context
    )
    pkg_name
}


.arrayAnnotationDependencyRequirements <- function(array, genome, reason = NULL) {
    pkg_name <- .getArrayAnnotationPackage(array = array, genome = genome)
    .makeDependencyRequirements(
        pkg_names = pkg_name,
        reason = if (is.null(reason)) {
            paste0(
                "The requested array/genome combination ('", array, "', '", genome,
                "') needs this annotation package to resolve probe coordinates."
            )
        } else {
            reason
        }
    )
}

supportedOrganisms <- function() {
    filename <- system.file(
        "extdata",
        "supportedOrganisms.csv",
        package = "CMEnt",
        mustWork = FALSE
    )
    if (!nzchar(filename)) {
        filename <- file.path("inst", "extdata", "supportedOrganisms.csv")
    }
    if (!file.exists(filename)) {
        stop(
            "Could not locate 'supportedOrganisms.csv' in the installed package or source tree.",
            call. = FALSE
        )
    }
    read.csv(filename, header = TRUE, stringsAsFactors = FALSE)
}

.getGeneAnnotationPackages <- function(genome) {
    target_genome <- tolower(genome)
    annotation_source_genome <- if (target_genome == "hs1") "hg38" else target_genome
    supported_organisms <- supportedOrganisms()
    matched_row <- supported_organisms[grepl(annotation_source_genome, tolower(supported_organisms$TxDb)), , drop = FALSE]
    if (nrow(matched_row) == 0L) {
        stop("Unsupported genome: ", genome, call. = FALSE)
    }
    c(txdb = matched_row$TxDb[1], orgdb = matched_row$OrgDb[1])
}


.assertGeneAnnotationPackagesInstalled <- function(genome, context) {
    pkgs <- .getGeneAnnotationPackages(genome)
    .assertDependencyRequirements(
        requirements = .makeDependencyRequirements(
            pkg_names = unname(pkgs),
            reason = paste0(
                "Gene annotation for genome '", genome,
                "' requires both the TxDb and OrgDb annotation packages."
            )
        ),
        context = context
    )
    pkgs
}


.geneAnnotationDependencyRequirements <- function(genome, reason = NULL) {
    pkgs <- .getGeneAnnotationPackages(genome)
    .makeDependencyRequirements(
        pkg_names = unname(pkgs),
        reason = if (is.null(reason)) {
            paste0(
                "Gene annotation for genome '", genome,
                "' requires both the TxDb and OrgDb annotation packages."
            )
        } else {
            reason
        }
    )
}


.resolveBSGenomePackage <- function(genome) {
    if (genome == "hg19") {
        return("BSgenome.Hsapiens.UCSC.hg19")
    }
    if (genome == "hg38") {
        return("BSgenome.Hsapiens.UCSC.hg38")
    }
    if (genome == "hs1") {
        return("BSgenome.Hsapiens.UCSC.hs1")
    }
    if (genome == "mm10") {
        return("BSgenome.Mmusculus.UCSC.mm10")
    }
    if (genome == "mm39") {
        return("BSgenome.Mmusculus.UCSC.mm39")
    }
    NULL
}


.assertBSGenomePackageInstalled <- function(genome, context) {
    pkg_name <- .resolveBSGenomePackage(genome)
    if (is.null(pkg_name)) {
        stop(
            "Unsupported genome '", genome,
            "' for local sequence extraction. Set use_online = TRUE or provide one of: hg19, hg38, hs1, mm10, mm39.",
            call. = FALSE
        )
    }
    .assertDependencyRequirements(
        requirements = .makeDependencyRequirements(
            pkg_names = pkg_name,
            reason = paste0(
                "Local sequence extraction for genome '", genome,
                "' uses the matching BSgenome package."
            )
        ),
        context = context
    )
    pkg_name
}


.bsgenomeDependencyRequirements <- function(genome, reason = NULL) {
    pkg_name <- .resolveBSGenomePackage(genome)
    if (is.null(pkg_name)) {
        stop(
            "Unsupported genome '", genome,
            "' for local sequence extraction. Set use_online = TRUE or provide one of: hg19, hg38, hs1, mm10, mm39.",
            call. = FALSE
        )
    }
    .makeDependencyRequirements(
        pkg_names = pkg_name,
        reason = if (is.null(reason)) {
            paste0(
                "Local sequence extraction for genome '", genome,
                "' uses the matching BSgenome package."
            )
        } else {
            reason
        }
    )
}


.jasparDependencyRequirements <- function(reason = NULL) {
    tax_group <- getOption("CMEnt.jaspar_tax_group", "vertebrates")
    jaspar_version <- getOption("CMEnt.jaspar_version", 2024)
    jaspar_pkg <- paste0("JASPAR", jaspar_version)
    cache_dir <- getOption(
        "CMEnt.jaspar_cache_dir",
        .getOSCacheDir(file.path("R", "CMEnt", "jaspar_cache"))
    )
    cache_key <- paste0("jaspar", jaspar_version, "_", tax_group, "_pwms")
    if (.hasBiocFileCacheRDS(cache_dir, cache_key)) {
        return(data.frame(pkg_name = character(0), reason = character(0), stringsAsFactors = FALSE))
    }
    .makeDependencyRequirements(
        pkg_names = jaspar_pkg,
        reason = if (is.null(reason)) {
            paste0(
                "JASPAR motif matching was requested and no cached PWM set was available for tax_group='",
                tax_group, "' and version ", jaspar_version, "."
            )
        } else {
            reason
        }
    )
}


.experimentHubDependencyRequirements <- function(resource, reason = NULL) {
    .makeDependencyRequirements(
        pkg_names = c("ExperimentHub", "AnnotationHub"),
        reason = if (is.null(reason)) {
            paste0(
                "The bundled example resource '", resource,
                "' was not found locally, so CMEnt needs ExperimentHub support to fetch it."
            )
        } else {
            reason
        }
    )
}


.findDMRsDependencyRequirements <- function(beta, array, genome,
                                           annotate_with_genes = TRUE,
                                           extract_motifs = TRUE,
                                           bed_provided = FALSE) {
    is_tabix_input <- is.character(beta) &&
        length(beta) == 1L &&
        file.exists(beta) &&
        file_is_tabix(beta)
    is_bed_input <- isTRUE(bed_provided) ||
        (is.character(beta) && length(beta) == 1L && file.exists(beta) && identical(tolower(tools::file_ext(beta)), "bed"))
    needs_array_annotations <- !inherits(beta, "BetaHandler") &&
        !is_bsseq(beta) &&
        !is_tabix_input &&
        !is_bed_input &&
        !is.null(array)

    .combineDependencyRequirements(
        if (needs_array_annotations) {
            .arrayAnnotationDependencyRequirements(
                array = array,
                genome = genome,
                reason = paste0(
                    "The supplied beta input needs array annotations for array='", array,
                    "' and genome='", genome, "' before DMR finding can start."
                )
            )
        },
        if (isTRUE(annotate_with_genes)) {
            .geneAnnotationDependencyRequirements(
                genome = genome,
                reason = paste0(
                    "'annotate_with_genes = TRUE' requires gene annotation packages for genome '",
                    genome, "'."
                )
            )
        },
        if (isTRUE(extract_motifs)) {
            .bsgenomeDependencyRequirements(
                genome = genome,
                reason = paste0(
                    "'extract_motifs = TRUE' requires the BSgenome package for genome '",
                    genome, "'."
                )
            )
        }
    )
}


.motifDependencyRequirements <- function(genome, array = NULL, beta_locs = NULL, context = "motif extraction") {
    if (is.null(beta_locs) && is.null(array)) {
        stop(
            context,
            " requires either 'beta_locs' or a non-NULL 'array' value to resolve site coordinates.",
            call. = FALSE
        )
    }
    .combineDependencyRequirements(
        .bsgenomeDependencyRequirements(
            genome = genome,
            reason = paste0(context, " requires the BSgenome package for genome '", genome, "'.")
        ),
        if (!is.null(array) && (is.null(beta_locs) || (is.character(beta_locs) && length(beta_locs) == 1L && file.exists(beta_locs)))) {
            .arrayAnnotationDependencyRequirements(
                array = array,
                genome = genome,
                reason = paste0(
                    context, " needs array annotations for array='", array,
                    "' and genome='", genome, "'."
                )
            )
        }
    )
}


#' Get Sorted Array Locations
#'
#' @description Retrieves and sorts genomic location annotations for the specified
#' methylation array platform and genome version. Performs liftOver if necessary.
#' The function caches the results.
#'
#' @param array Character. Array platform type (supported: "450K", "EPIC", "EPICv2", "27K", "Mouse", 'NULL'), ignored when locations_file is provided. Must be 'NULL' when the experiment is not array-based.
#' @param genome Character. Genome version (supported: "hg38", "hg19", "hs1", "mm10", "mm39"), ignored if locations_file is provided
#' @param locations_file Character. Optional path to a precomputed locations file (RDS format). If provided, this file will be used directly (default: NULL)
#'
#' @return A data frame containing sorted genomic locations with rownames as site IDs and columns:
#' \itemize{
#'   \item chr: Chromosome
#'   \item start: Genomic position
#'   \item start: Start position (same as start)
#'   \item end: End position (start + 1)
#' }
#'
#' @examples
#' # Get sorted locations for 450K array (hg38)
#' locs_450k <- getSortedGenomicLocs("450K")
#'
#' # Get sorted locations for EPIC array with hg38
#' locs_epic <- getSortedGenomicLocs("EPIC", "hg38")
#'
#' # Get sorted locations for EPICv2 array
#' locs_epicv2 <- getSortedGenomicLocs("EPICv2", "hg38")
#'
#' @export
getSortedGenomicLocs <- function(array = c("450K", "27K", "EPIC", "EPICv2", "Mouse"), genome = c("hg38", "hg19", "hs1", "mm10", "mm39"), locations_file = NULL) {
    if (!is.null(locations_file) && file.exists(locations_file)) {
        locs <- readRDS(locations_file)
        return(locs)
    }
    genome <- strex::match_arg(genome, ignore_case = TRUE)
    array_based <- !is.null(array)
    if (!array_based) {
        stop("Provided array is NULL but locations file was not provided.")
    }
    array <- strex::match_arg(array, ignore_case = TRUE)
    cache_dir <- getOption(
        "CMEnt.annotation_cache_dir",
        .getOSCacheDir(file.path("R", "CMEnt", "annotation_cache"))
    )

    array <- tolower(array)
    genome <- tolower(genome)
    cache_key <- paste0(
        array, "_", genome,
        "_locations"
    )
    locs <- if (getOption("CMEnt.use_annotation_cache", TRUE)) {
        .readBiocFileCacheRDS(cache_dir, cache_key)
    } else {
        NULL
    }
    if (!is.null(locs)) {
        .log_info("Using cached annotation file: ", cache_key, level = 3)
        return(locs)
    }
    pkg_name <- .assertArrayAnnotationPackageInstalled(
        array = array,
        genome = genome,
        context = "getSortedGenomicLocs()"
    )
    source_genome <- NULL
    if (tolower(array) %in% c("450k", "epic", "27k")) {
        source_genome <- "hg19"
    } else if (tolower(array) == "epicv2") {
        source_genome <- "hg38"
    } else if (tolower(array) == "mouse") {
        source_genome <- "mm10"
    }
    locs <- .loadAnnotationLocations(pkg_name, source_genome = source_genome)
    from_genome <- NULL
    if (!is.null(source_genome) && !identical(genome, source_genome)) {
        from_genome <- source_genome
    }
    if (!is.null(from_genome)) {
        locs <- .liftOverFromGenomeToGenome(locs, from_genome, genome)
    }
    locs <- convertToDataFrame(locs)
    ord <- stringr::str_order(paste0(locs[, "chr"], ":", locs[, "start"]), numeric = TRUE)
    locs <- locs[ord, , drop = FALSE]
    locs <- locs[!duplicated(rownames(locs)), ]
    if (!"start" %in% colnames(locs)) {
        locs[, "start"] <- locs[, "start"]
    }
    if (!"start" %in% colnames(locs)) {
        locs[, "start"] <- locs[, "start"]
    }
    if (!"end" %in% colnames(locs)) {
        locs[, "end"] <- locs[, "start"] + 1
    }
    locs[locs[, "end"] == locs[, "start"], "end"] <- locs[locs[
        ,
        "end"
    ] == locs[, "start"], "start"] + 1
    locs$name <- rownames(locs)
    locs <- getRegistry(locs, "name")
    tryCatch(
        {
            .saveBiocFileCacheRDS(locs, cache_dir, cache_key)
        },
        warning = function(w) {
            .log_warn(
                "Could not write annotation cache entry '", cache_key,
                "' (warning: ", conditionMessage(w), "). Continuing without cache persistence."
            )
        },
        error = function(e) {
            .log_warn(
                "Could not write annotation cache entry '", cache_key,
                "' (error: ", conditionMessage(e), "). Continuing without cache persistence."
            )
        }
    )
    locs
}

#' Orders a vector of indices according to their corresponding genomic
#' locations (chromosome and position). This function is useful for sorting site
#' sites or other genomic features by their physical positions.
#'
#' @param x Character or integer vector. Indices or identifiers to be ordered
#' @param array Character. Array platform type, either "450K" or "EPIC" (default: "450K")
#' @param genome Character. Genome version, either "hg38", "hg19", or "hs1" (default: "hg38")
#' @param genomic_locs Data frame. Optional pre-computed genomic locations. If NULL,
#' locations will be retrieved using getSortedGenomicLocs (default: NULL)
#'
#' @return Integer vector of ordered indices
#'
#' @examples
#' # Order site indices by genomic location
#' site_ids <- c("cg00000029", "cg00000108", "cg00000109")
#' ordered_indices <- orderByLoc(site_ids, array = "450K")
#'
#' # Order using pre-computed genomic locations
#' locs <- getSortedGenomicLocs("EPIC", "hg38")
#' ordered_indices <- orderByLoc(site_ids, genomic_locs = locs)
#'
#' @export
orderByLoc <- function(x,
                       array = c("450K", "27K", "EPIC", "EPICv2"),
                       genome = c("hg38", "hg19", "hs1", "mm10", "mm39"),
                       genomic_locs = NULL) {
    if (is.null(genomic_locs)) {
        genomic_locs <- getSortedGenomicLocs(array, genome)
    }
    order(match(x, rownames(genomic_locs)), method = "radix")
}

#' @keywords internal
#' @noRd
.splitCsvValues <- function(x) {
    if (length(x) == 0 || is.null(x)) {
        return(character(0))
    }
    x <- x[[1]]
    if (is.na(x) || !nzchar(as.character(x))) {
        return(character(0))
    }
    vals <- unlist(base::strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE)
    vals <- trimws(vals)
    vals <- vals[nzchar(vals)]
    vals[!is.na(vals)]
}


#' @keywords internal
#' @noRd
.splitCsvIndices <- function(x) {
    vals <- .splitCsvValues(x)
    if (length(vals) == 0) {
        return(integer(0))
    }
    suppressWarnings({
        inds <- as.integer(vals)
    })
    inds[!is.na(inds)]
}


#' @keywords internal
#' @noRd
.downsampleFlankIndices <- function(indices, max_sup_sites_per_dmr_side) {
    if (is.null(max_sup_sites_per_dmr_side) || max_sup_sites_per_dmr_side <= 0) {
        return(indices)
    }
    if (length(indices) <= max_sup_sites_per_dmr_side) {
        return(indices)
    }
    step <- ceiling(length(indices) / max_sup_sites_per_dmr_side)
    indices[seq.int(1L, length(indices), by = step)]
}


.extendGRangesWithFlanks <- function(granges, uflank_size = 0, dflank_size = 0) {
    if (uflank_size > 0 || dflank_size > 0) {
        flanked_start <- GenomicRanges::start(granges) - uflank_size
        flanked_end <- GenomicRanges::end(granges) + dflank_size
        start_off_limit <- (flanked_start < 1) * (1 - flanked_start)
        seqls <- GenomeInfoDb::seqlengths(granges)[
            as.character(GenomeInfoDb::seqnames(granges))
        ]
        end_off_limit <- (flanked_end > seqls) * (flanked_end - seqls)
        BiocGenerics::start(granges) <- pmax(flanked_start, 1)
        BiocGenerics::end(granges) <- pmin(flanked_end, seqls)
        ret <- list(granges = granges, start_off_limit = start_off_limit, end_off_limit = end_off_limit)
    } else {
        ret <- list(granges = granges, start_off_limit = rep(0, length(granges)), end_off_limit = rep(0, length(granges)))
    }
    ret
}

.getBSGenomePackage <- function(genome) {
    pkg_name <- .resolveBSGenomePackage(genome)
    if (is.null(pkg_name) || !.isPackageInstalled(pkg_name)) {
        return(NULL)
    }
    pkg_name
}


.getOSCacheDir <- function(prefix) {
    x <- normalizePath(tools::R_user_dir(prefix, which = "cache"), mustWork = FALSE)
    try(dir.create(x, recursive = TRUE, showWarnings = FALSE))
    x
}

.getBiocFileCache <- function(cache_dir) {
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    BiocFileCache::BiocFileCache(cache_dir, ask = FALSE)
}

.getBiocFileCachePath <- function(bfc, rname, ext = "", create = TRUE) {
    cache_records <- BiocFileCache::bfcquery(
        bfc,
        query = rname,
        field = "rname",
        exact = TRUE
    )
    if (nrow(cache_records) > 0) {
        existing_paths <- cache_records$rpath[file.exists(cache_records$rpath)]
        if (length(existing_paths) > 0) {
            return(existing_paths[[1]])
        }
        if (!is.na(cache_records$rpath[[1]]) && nzchar(cache_records$rpath[[1]])) {
            return(cache_records$rpath[[1]])
        }
    }

    if (!isTRUE(create)) {
        return(NULL)
    }

    unname(BiocFileCache::bfcnew(bfc, rname = rname, ext = ext))
}

.readBiocFileCacheRDS <- function(cache_dir, rname) {
    bfc <- .getBiocFileCache(cache_dir)
    cache_file <- .getBiocFileCachePath(bfc, rname = rname, ext = ".rds", create = FALSE)
    if (is.null(cache_file) || !file.exists(cache_file)) {
        return(NULL)
    }
    readRDS(cache_file)
}

.hasBiocFileCacheRDS <- function(cache_dir, rname) {
    bfc <- .getBiocFileCache(cache_dir)
    cache_file <- .getBiocFileCachePath(bfc, rname = rname, ext = ".rds", create = FALSE)
    !is.null(cache_file) && file.exists(cache_file)
}

.saveBiocFileCacheRDS <- function(object, cache_dir, rname) {
    bfc <- .getBiocFileCache(cache_dir)
    cache_file <- .getBiocFileCachePath(bfc, rname = rname, ext = ".rds")
    saveRDS(object, cache_file)
    invisible(cache_file)
}

.getDerivedOutputPath <- function(output_prefix, suffix) {
    if (is.null(output_prefix) || !nzchar(output_prefix)) {
        return(NULL)
    }
    paste0(output_prefix, suffix)
}

.downloadFirstAvailable <- function(urls, destfile, mode = "wb", quiet = TRUE) {
    attempts <- character(length(urls))

    for (i in seq_along(urls)) {
        error_message <- tryCatch(
            {
                utils::download.file(urls[[i]], destfile, quiet = quiet, mode = mode)
                NULL
            },
            error = function(e) {
                conditionMessage(e)
            }
        )

        if (is.null(error_message)) {
            return(urls[[i]])
        }

        attempts[[i]] <- paste0(urls[[i]], " (", error_message, ")")
    }

    stop(
        "Failed to download resource from all candidate URLs: ",
        paste(attempts, collapse = "; ")
    )
}

#' Extract DNA Sequences for DMRs
#'
#' @description Retrieves the DNA sequences corresponding to genomic regions
#' specified in a GRanges object. This function is useful for extracting the
#' actual DNA sequence of identified DMRs for downstream analyses such as
#' motif finding or sequence composition analysis.
#'
#' @param dmrs GRanges object containing genomic coordinates of DMRs
#' @param genome Character. Genome version to use for sequence extraction, .e.g. "hg38" or "hs1".
#' @param use_online Logical. If TRUE, forces use of online UCSC API instead of
#'   BSgenome packages. If FALSE (default), uses BSgenome packages with online
#'   fallback when packages are unavailable (default: FALSE)
#' @param uflank_size Integer. Number of base pairs to add as flanking regions
#'   upstream of each DMR (default: 0)
#' @param dflank_size Integer. Number of base pairs to add as flanking regions
#'   downstream of each DMR (default: 0)
#' @param batch_size Integer. For online API, number of regions to process per batch (default: 100)
#' @param njobs Integer. For online API, number of cores for parallel processing (default: 1)
#' @return A Character vector containing DNA sequences for each DMR
#'
#' @details
#' The function first attempts to use genome-appropriate BSgenome packages:
#' \itemize{
#'   \item hg19: BSgenome.Hsapiens.UCSC.hg19
#'   \item hg38: BSgenome.Hsapiens.UCSC.hg38
#'   \item hs1: BSgenome.Hsapiens.UCSC.hs1
#'   \item mm10: BSgenome.Mmusculus.UCSC.mm10
#'   \item mm39: BSgenome.Mmusculus.UCSC.mm39
#' }
#'
#' If the required BSgenome package is not installed, the function raises an
#' error with installation instructions. Set `use_online = TRUE` to query
#' sequences from the UCSC Genome Browser REST API instead. The online method
#' processes sequences in batches with optional parallel processing for
#' improved performance with large datasets.
#'
#' For large numbers of DMRs (>10k), consider using parallel processing by setting
#' njobs > 1 when using the online API, or install the appropriate BSgenome package
#' for much faster local sequence retrieval.
#'
#' @examples
#' # Extract sequences for DMRs using BSgenome packages
#' sequences <- getDMRSequences(dmrs, "hg19")
#'
#' # Force use of online UCSC API with parallel processing
#' sequences <- getDMRSequences(dmrs, "hg19", use_online = TRUE, njobs = 4)
#'
#' # Calculate GC content
#' gc_content <- sapply(sequences, function(s) {
#'     (stringr::str_count(s, "G") + stringr::str_count(s, "C")) / nchar(s)
#' })
#'
#' @importFrom BSgenome getSeq
#' @importFrom rtracklayer import.chain liftOver
#' @export
getDMRSequences <- function(dmrs, genome, use_online = FALSE, uflank_size = 0, dflank_size = 0,
                            batch_size = 100, njobs = 1) {
    if (!use_online) {
        pkg_name <- .assertBSGenomePackageInstalled(genome, context = "getDMRSequences()")
        use_bsgenome <- TRUE
    } else {
        pkg_name <- NULL
        use_bsgenome <- FALSE
    }
    extended_ret <- .extendGRangesWithFlanks(dmrs, uflank_size, dflank_size)
    dmrs <- extended_ret$granges
    add_na_to_the_start <- extended_ret$start_off_limit
    add_na_to_the_end <- extended_ret$end_off_limit

    if (use_bsgenome) {
        .log_info("Querying sequences using BSgenome package...", level = 2)
        if (!isNamespaceLoaded(pkg_name)) {
            loadNamespace(pkg_name)
        }
        seq_db <- getExportedValue(pkg_name, pkg_name)
        sequences <- Biostrings::getSeq(seq_db, dmrs, as.character = TRUE)
        if (is.list(sequences)) {
            sequences <- sapply(sequences, function(x) paste(x, collapse = ""))
        }
    } else {
        .log_info("Querying sequences from UCSC Genome Browser API...", level = 2)
        sequences <- .getSequencesFromUCSC(dmrs, genome, batch_size = batch_size, njobs = njobs)
    }
    off_bound_mask <- add_na_to_the_start > 0 | add_na_to_the_end > 0
    if (any(off_bound_mask)) {
        .log_info("  Found ", sum(off_bound_mask), " regions with out-of-bound flanking extensions", level = 2)
        .log_info("Adding 'N' padding for out-of-bound flanking regions...", level = 2)
        sequences[off_bound_mask] <- mapply(function(seq, na_start, na_end) {
            if (is.na(seq)) {
                return(NA_character_)
            }
            seq <- paste0(
                paste(rep("N", na_start), collapse = ""),
                seq,
                paste(rep("N", na_end), collapse = "")
            )
            seq
        }, sequences[off_bound_mask], add_na_to_the_start[off_bound_mask], add_na_to_the_end[off_bound_mask], SIMPLIFY = TRUE)
    }
    sequences
}


#' @export
getSiteBackgroundCounts <- function(regions, genome, njobs = 1, canonical_chr = TRUE) {
    pkg_name <- .getBSGenomePackage(genome)
    if (is.null(pkg_name)) {
        sequences <- getDMRSequences(regions, genome, use_online = TRUE, njobs = njobs)
        site_counts <- sapply(sequences, function(seq) {
            if (is.na(seq)) {
                return(NA_integer_)
            }
            stringr::str_count(seq, "CG")
        })
        return(unlist(site_counts))
    }
    cache_dir <- getOption(
        "CMEnt.annotation_cache_dir",
        .getOSCacheDir(file.path("R", "CMEnt", "annotation_cache"))
    )
    cache_key <- paste0(pkg_name, "_site_positions")
    if (canonical_chr) {
        cache_key <- paste0(pkg_name, "_canonical_chr_site_positions")
    }
    sites <- .readBiocFileCacheRDS(cache_dir, cache_key)
    if (!is.null(sites)) {
        sites <- data.table::as.data.table(sites)
        data.table::setkey(sites, chr, start, end)
    } else {
        .log_info("Missing site positions cache. Generating site positions from BSgenome package...", level = 2)
        if (!isNamespaceLoaded(pkg_name)) {
            loadNamespace(pkg_name)
        }
        seq_db <- getExportedValue(pkg_name, pkg_name)
        chrs <- names(seq_db)
        if (canonical_chr) {
            chrs <- chrs[grepl("^chr[0-9XYM]+$", chrs)]
        }
        cgs <- lapply(chrs, function(x) start(Biostrings::matchPattern("CG", seq_db[[x]])))
        names(cgs) <- chrs
        suppressWarnings(
            sites <- do.call(
                c, lapply(
                    seq_along(chrs),
                    function(x) {
                        GenomicRanges::GRanges(
                            seqnames = chrs[x],
                            ranges = IRanges::IRanges(cgs[[x]], width = 2)
                        )
                    }
                )
            )
        )
        sites <- data.table::as.data.table(as.data.frame(sites, stringsAsFactors = FALSE))[, 1:2]
        colnames(sites) <- c("chr", "start")
        sites[, "end"] <- sites[, "start"]
        data.table::setkey(sites, chr, start, end)
        .saveBiocFileCacheRDS(sites, cache_dir, cache_key)
    }
    regions <- as.data.frame(regions, stringsAsFactors = FALSE)[, 1:3]
    regions <- data.table::as.data.table(regions)
    colnames(regions) <- c("rchr", "rstart", "rend")
    regions[, "id"] <- seq_len(nrow(regions))
    data.table::setkey(regions, rchr, rstart, rend)
    site_counts <- data.table::foverlaps(regions,
        sites,
        by.x = c("rchr", "rstart", "rend"),
        by.y = c("chr", "start", "end"),
    )[, .N, by = id]
    site_counts <- site_counts[order(site_counts$id), "N"]
    unlist(site_counts)
    # }
}

#' Query DNA Sequences from UCSC Genome Browser API
#'
#' @description Internal helper function to retrieve DNA sequences from the
#' UCSC Genome Browser REST API when BSgenome packages are not available.
#' Uses parallel processing and batched requests for improved performance.
#'
#' @param dmrs GRanges object containing genomic coordinates
#' @param genome Character. Genome version (e.g., "hg38", "hs1", "mm39")
#' @param batch_size Integer. Number of regions to process per batch (default: 100)
#' @param njobs Integer. Number of cores for parallel processing (default: 1)
#'
#' @return Character vector of DNA sequences
#'
#' @keywords internal
#' @noRd
.getSequencesFromUCSC <- function(dmrs, genome, batch_size = 100, njobs = 1) {
    base_url <- "https://api.genome.ucsc.edu/getData/sequence"
    n_dmrs <- length(dmrs)

    if (n_dmrs == 0) {
        return(character(0))
    }

    n_batches <- ceiling(n_dmrs / batch_size)

    if (n_dmrs > 100) {
        .log_info(sprintf(
            "Retrieving sequences for %d regions in %d batches using %d core(s)...",
            n_dmrs, n_batches, njobs
        ), level = 2)
    }

    fetch_batch <- function(batch_idx) {
        start_idx <- (batch_idx - 1) * batch_size + 1
        end_idx <- min(batch_idx * batch_size, n_dmrs)
        batch_dmrs <- dmrs[start_idx:end_idx]
        batch_sequences <- character(length(batch_dmrs))

        for (j in seq_along(batch_dmrs)) {
            chr <- as.character(GenomeInfoDb::seqnames(batch_dmrs[j]))
            start <- GenomicRanges::start(batch_dmrs[j])
            end <- GenomicRanges::end(batch_dmrs[j])

            url <- sprintf(
                "%s?genome=%s;chrom=%s;start=%d;end=%d",
                base_url, genome, chr, start - 1, end
            )

            tryCatch(
                {
                    response <- readLines(url, warn = FALSE)
                    json_data <- jsonlite::fromJSON(paste(response, collapse = ""))

                    if (!is.null(json_data$dna)) {
                        batch_sequences[j] <- toupper(json_data$dna)
                    } else {
                        batch_sequences[j] <- NA_character_
                    }
                },
                error = function(e) {
                    batch_sequences[j] <- NA_character_
                }
            )

            if (j %% 10 == 0) {
                Sys.sleep(0.05)
            }
        }

        batch_sequences
    }

    if (njobs > 1 && n_batches > 1) {
        bp_param <- .makeBiocParallelParam(njobs, n_tasks = n_batches)
        batch_results <- BiocParallel::bplapply(
            seq_len(n_batches),
            fetch_batch,
            BPPARAM = bp_param
        )
        sequences <- unlist(batch_results)
    } else {
        sequences <- character(n_dmrs)
        for (i in seq_len(n_batches)) {
            batch_seqs <- fetch_batch(i)
            start_idx <- (i - 1) * batch_size + 1
            end_idx <- min(i * batch_size, n_dmrs)
            sequences[start_idx:end_idx] <- batch_seqs

            if (n_dmrs > 100 && i %% 10 == 0) {
                .log_info(sprintf("  Progress: %d/%d batches completed", i, n_batches), level = 2)
            }
        }
    }

    sequences
}

#' @keywords internal
#' @noRd
.annotateDMRsWithGeneFeature <- function(dmrs, features, orgdb_pkg,
                                         feature_type = c("promoter", "gene_body")) {
    feature_type <- match.arg(feature_type)
    annotations <- rep(NA_character_, length(dmrs))
    if (length(dmrs) == 0L || length(features) == 0L) {
        return(annotations)
    }

    overlaps <- GenomicRanges::findOverlaps(dmrs, features)
    if (length(overlaps) == 0L) {
        return(annotations)
    }

    if (!isNamespaceLoaded(orgdb_pkg)) {
        loadNamespace(orgdb_pkg)
    }
    orgdb <- getExportedValue(orgdb_pkg, orgdb_pkg)
    feature_hits <- S4Vectors::subjectHits(overlaps)
    entrez_ids <- switch(feature_type,
        promoter = as.character(S4Vectors::mcols(features[feature_hits])$name),
        gene_body = as.character(names(features)[feature_hits])
    )
    valid_entrez_ids <- unique(entrez_ids[!is.na(entrez_ids) & nzchar(entrez_ids)])
    if (length(valid_entrez_ids) == 0L) {
        return(annotations)
    }

    symbols <- suppressMessages(AnnotationDbi::mapIds(orgdb,
        keys = valid_entrez_ids,
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
    ))
    names(symbols) <- valid_entrez_ids
    symbols_by_dmr <- split(
        symbols[entrez_ids],
        S4Vectors::queryHits(overlaps)
    )
    collapsed_symbols <- vapply(symbols_by_dmr, function(x) {
        genes_vec <- unique(stats::na.omit(as.character(x)))
        if (length(genes_vec) == 0L) {
            return(NA_character_)
        }
        paste(genes_vec, collapse = ",")
    }, character(1))
    annotations[as.integer(names(collapsed_symbols))] <- unname(collapsed_symbols)
    annotations
}

#' Annotate DMRs with Gene Information
#'
#' @description Annotates DMRs with overlapping gene promoters and gene bodies
#' using TxDb annotations. For each DMR, identifies genes whose promoters or
#' gene bodies overlap with the DMR coordinates.
#'
#' @param dmrs Dataframe or GRanges object containing DMR coordinates
#' @param genome Character. Genome version to use for gene annotation. (default: "hg38")
#' @param promoter_upstream Integer. Number of base pairs upstream of TSS to
#'   define promoter region (default: 2000)
#' @param promoter_downstream Integer. Number of base pairs downstream of TSS
#'   to define promoter region (default: 200)
#' @param njobs Integer. Number of parallel jobs used to annotate promoter and
#'   gene-body overlaps (default: `getOption("CMEnt.njobs")`)
#'
#' @return The input Dataframe/GRanges object with additional metadata columns:
#' \itemize{
#'   \item in_promoter_of: Character vector of gene symbols with promoters overlapping the DMR (comma-separated)
#'   \item in_gene_body_of: Character vector of gene symbols with gene bodies overlapping the DMR (comma-separated)
#' }
#'
#' @details
#' The function uses genome-appropriate TxDb packages. For `hs1`, CMEnt
#' uses hg38 gene models and lifts them to `hs1` before computing overlaps.
#' Gene symbols are retrieved from the appropriate org.*.eg.db package.
#' Multiple overlapping genes are concatenated with commas.
#'
#' @examples
#' # Annotate DMRs with gene information
#' dmrs <- data.frame(
#'     chr = c("chr1", "chr2"),
#'     start = c(1000000, 2000000),
#'     end = c(1001000, 2001000)
#' )
#' dmrs_annotated <- annotateDMRsWithGenes(dmrs, genome = "hg38")
#'
#' # Use custom promoter definition
#' dmrs_annotated <- annotateDMRsWithGenes(
#'     dmrs,
#'     genome = "hg38",
#'     promoter_upstream = 5000,
#'     promoter_downstream = 1000,
#'     njobs = 2
#' )
#'
#' @export
annotateDMRsWithGenes <- function(dmrs, genome = "hg38",
                                  promoter_upstream = 2000,
                                  promoter_downstream = 200,
                                  njobs = getOption("CMEnt.njobs", min(8, future::availableCores() - 1))) {
    cache_dir <- getOption(
        "CMEnt.annotation_cache_dir",
        .getOSCacheDir(file.path("R", "CMEnt", "annotation_cache"))
    )
    dmrs_df_provided <- is.data.frame(dmrs)
    dmrs <- convertToGRanges(dmrs, genome)

    target_genome <- tolower(genome)
    annotation_source_genome <- if (target_genome == "hs1") "hg38" else target_genome
    annotation_pkgs <- .assertGeneAnnotationPackagesInstalled(
        genome = genome,
        context = "annotateDMRsWithGenes()"
    )
    txdb_pkg <- unname(annotation_pkgs[["txdb"]])
    orgdb_pkg <- unname(annotation_pkgs[["orgdb"]])
    if (annotation_source_genome != target_genome) {
        .log_info(
            "Using ", annotation_source_genome,
            " gene models lifted to ", target_genome, " for annotation.",
            level = 2
        )
    }
    .log_step("Loading gene annotations for ", genome, "...", level = 2)

    # Load TxDb namespace
    if (!isNamespaceLoaded(txdb_pkg)) {
        loadNamespace(txdb_pkg)
    }

    # Load TxDb - the main object has the same name as the package
    txdb <- getExportedValue(txdb_pkg, txdb_pkg)

    # Get genes and promoters
    # Load them from cache if available
    suppressMessages({
        genes_key <- paste0("genes_", genome)
        genes <- if (getOption("CMEnt.use_annotation_cache", TRUE)) {
            .readBiocFileCacheRDS(cache_dir, genes_key)
        } else {
            NULL
        }
        if (!is.null(genes)) {
            .log_info("Loading cached genes from ", genes_key, level = 2)
        } else {
            genes <- GenomicFeatures::genes(txdb)
            # get genes only within standard chromosomes
            std_chroms <- GenomeInfoDb::standardChromosomes(GenomeInfoDb::seqinfo(txdb))
            genes <- genes[as.character(GenomeInfoDb::seqnames(genes)) %in% std_chroms]
            if (annotation_source_genome != target_genome) {
                genes <- .liftOverFromGenomeToGenome(genes, annotation_source_genome, target_genome)
                target_std_chroms <- GenomeInfoDb::standardChromosomes(GenomeInfoDb::Seqinfo(genome = target_genome))
                genes <- genes[as.character(GenomeInfoDb::seqnames(genes)) %in% target_std_chroms]
            }
            tryCatch(
                {
                    .saveBiocFileCacheRDS(genes, cache_dir, genes_key)
                },
                warning = function(w) {
                    .log_warn(
                        "Could not write annotation cache entry '", genes_key,
                        "' (warning: ", conditionMessage(w), "). Continuing without cache persistence."
                    )
                },
                error = function(e) {
                    .log_warn(
                        "Could not write annotation cache entry '", genes_key,
                        "' (error: ", conditionMessage(e), "). Continuing without cache persistence."
                    )
                }
            )
        }
        promoters_key <- paste0("promoters_", genome)
        promoters <- if (getOption("CMEnt.use_annotation_cache", TRUE)) {
            .readBiocFileCacheRDS(cache_dir, promoters_key)
        } else {
            NULL
        }
        if (!is.null(promoters)) {
            .log_info("Loading cached promoters from ", promoters_key, level = 2)
        } else {
            transcripts_by_gene <- GenomicFeatures::transcriptsBy(txdb, by = "gene")
            transcripts_by_gene <- transcripts_by_gene[names(transcripts_by_gene) %in% names(genes)]
            promoters <- GenomicFeatures::promoters(transcripts_by_gene, upstream = promoter_upstream, downstream = promoter_downstream)
            promoters <- stack(promoters)
            if (annotation_source_genome != target_genome) {
                promoters <- .liftOverFromGenomeToGenome(promoters, annotation_source_genome, target_genome)
                target_std_chroms <- GenomeInfoDb::standardChromosomes(GenomeInfoDb::Seqinfo(genome = target_genome))
                promoters <- promoters[as.character(GenomeInfoDb::seqnames(promoters)) %in% target_std_chroms]
            }
            tryCatch(
                {
                    .saveBiocFileCacheRDS(promoters, cache_dir, promoters_key)
                },
                warning = function(w) {
                    .log_warn(
                        "Could not write annotation cache entry '", promoters_key,
                        "' (warning: ", conditionMessage(w), "). Continuing without cache persistence."
                    )
                },
                error = function(e) {
                    .log_warn(
                        "Could not write annotation cache entry '", promoters_key,
                        "' (error: ", conditionMessage(e), "). Continuing without cache persistence."
                    )
                }
            )
        }
    })
    .log_success("Gene annotations loaded: ", length(genes), " genes", level = 2)
    .log_step("Finding overlaps with promoters and gene bodies...", level = 2)
    .log_step("Mapping overlapping Entrez IDs to gene symbols...", level = 2)
    annotation_specs <- list(
        list(
            column = "in_promoter_of",
            features = promoters,
            feature_type = "promoter"
        ),
        list(
            column = "in_gene_body_of",
            features = genes,
            feature_type = "gene_body"
        )
    )
    if (!is.null(njobs) && is.finite(njobs) && as.integer(njobs) > 1L) {
        withr::local_options(list(CMEnt.njobs = as.integer(njobs)))
        .setupParallel()
        on.exit(.finalizeParallel(), add = TRUE)
        annotation_results <- future.apply::future_lapply(
            annotation_specs,
            function(spec) {
                list(
                    column = spec$column,
                    values = .annotateDMRsWithGeneFeature(
                        dmrs = dmrs,
                        features = spec$features,
                        orgdb_pkg = orgdb_pkg,
                        feature_type = spec$feature_type
                    )
                )
            },
            future.seed = TRUE,
            future.stdout = NA,
            future.globals = c(".annotateDMRsWithGeneFeature", "dmrs", "orgdb_pkg")
        )
    } else {
        annotation_results <- lapply(annotation_specs, function(spec) {
            list(
                column = spec$column,
                values = .annotateDMRsWithGeneFeature(
                    dmrs = dmrs,
                    features = spec$features,
                    orgdb_pkg = orgdb_pkg,
                    feature_type = spec$feature_type
                )
            )
        })
    }
    for (annotation_result in annotation_results) {
        S4Vectors::mcols(dmrs)[[annotation_result$column]] <- annotation_result$values
    }
    .log_success("Gene symbols mapped", level = 2)
    if (dmrs_df_provided) {
        dmrs <- convertToDataFrame(dmrs)
    }
    dmrs
}

convertToGRanges <- function(obj, genome) {
    input_is_df <- !inherits(obj, "GRanges")
    if (input_is_df) {
        obj <- .decodeSerializedOutputColumns(obj)
        # If the chr prefix is missing, add it (e.g. "1" -> "chr1")
        if (!any(grepl("^chr", obj[[1]]))) {
            obj[[1]] <- paste0("chr", obj[[1]])
            obj$chr_prefix_added <- TRUE
        }
        # if the chromosome appears in the form of chr1:1230, save the original location in a separate column and parse the location into chr, start, end
        if (any(grepl(":", obj[[1]]))) {
            obj$original_location <- obj[[1]]
            loc_split <- base::strsplit(as.character(obj[[1]]), ":", fixed = TRUE)
            obj[[1]] <- sapply(loc_split, function(x) x[1])
        }
        obj <- GenomicRanges::makeGRangesFromDataFrame(obj,
            keep.extra.columns = TRUE,
            seqinfo = GenomeInfoDb::Seqinfo(genome = genome),
            na.rm = TRUE
        )
    } else {
        if (!is(obj, "GRanges")) {
            stop("dmrs must be a data.frame or GRanges object")
        }
        # if the genome info in the gr is different from the specified genome, update the locations with liftOver
        grs_genome <- GenomeInfoDb::genome(GenomeInfoDb::seqinfo(obj))[[1]]
        if (is.na(grs_genome)) {
            .log_warn("Input GRanges object has no genome information. Assuming genome: ", genome)
            grs_genome <- genome
            GenomeInfoDb::seqinfo(obj) <- GenomeInfoDb::Seqinfo(genome = genome)
        }
        if (grs_genome != genome) {
            obj <- .liftOverFromGenomeToGenome(obj, grs_genome, genome)
        }
    }
    obj
}

.serializedOutputPrefix <- "CMEnt:serialized_base64:"

#' @keywords internal
#' @noRd
.serializeOutputValue <- function(x) {
    encoded <- jsonlite::base64_enc(serialize(x, connection = NULL, version = 2))
    encoded <- gsub("[\r\n]", "", encoded)
    paste0(
        .serializedOutputPrefix,
        encoded
    )
}

#' @keywords internal
#' @noRd
.isSerializedOutputColumn <- function(x) {
    is.character(x) &&
        length(x) > 0L &&
        any(!is.na(x) & startsWith(x, .serializedOutputPrefix)) &&
        all(is.na(x) | startsWith(x, .serializedOutputPrefix))
}

#' @keywords internal
#' @noRd
.encodeNonTabularColumns <- function(df) {
    stopifnot(is.data.frame(df))

    encoded_columns <- names(df)[vapply(df, is.list, logical(1))]
    if (length(encoded_columns) == 0L) {
        return(list(data = df, encoded_columns = character(0)))
    }

    encoded_df <- df
    for (col in encoded_columns) {
        encoded_df[[col]] <- vapply(encoded_df[[col]], .serializeOutputValue, character(1))
    }

    list(data = encoded_df, encoded_columns = encoded_columns)
}

#' @keywords internal
#' @noRd
.decodeSerializedOutputColumns <- function(df) {
    stopifnot(is.data.frame(df))

    serialized_columns <- names(df)[vapply(df, .isSerializedOutputColumn, logical(1))]
    if (length(serialized_columns) == 0L) {
        return(df)
    }

    decoded_df <- df
    for (col in serialized_columns) {
        decoded_df[[col]] <- lapply(decoded_df[[col]], function(x) {
            if (is.na(x) || !startsWith(x, .serializedOutputPrefix)) {
                return(x)
            }

            payload <- sub(
                paste0("^", .serializedOutputPrefix),
                "",
                x
            )
            unserialize(jsonlite::base64_dec(payload))
        })
    }

    decoded_df
}

convertToDataFrame <- function(gr) {
    if (is.data.frame(gr)) {
        return(gr)
    }
    chr_prefix_added <- FALSE
    if ("chr_prefix_added" %in% names(mcols(gr))) {
        chr_prefix_added <- TRUE
    }
    df <- as.data.frame(gr, stringsAsFactors = FALSE)
    colnames(df)[colnames(df) == "seqnames"] <- "chr"
    if ("original_location" %in% colnames(df)) {
        df <- df[, c("original_location", setdiff(colnames(df), c("chr", "original_location")))]
        colnames(df)[colnames(df) == "original_location"] <- "chr"
    }
    if (chr_prefix_added) {
        df$chr <- sub("^chr", "", df$chr)
    }
    df
}

.already_logged_dir <- tempdir()
.already_logged_file <- file.path(.already_logged_dir, "cment_already_logged_parallel.txt")
#' @keywords internal
#' @noRd
.cleanupParallelState <- function() {
    # Reset plan first so current backend gets torn down by future.
    tryCatch(
        future::plan(future::sequential),
        error = function(e) invisible(NULL)
    )

    # Optional deep cleanup for stale multisession clusters.
    # Disabled by default because stopping dead clusters may block.
    if (isTRUE(getOption("CMEnt.force_cluster_cleanup", FALSE))) {
        reg <- tryCatch(
            getFromNamespace("clusterRegistry", "future"),
            error = function(e) NULL
        )
        if (is.list(reg) && is.function(reg$stopCluster)) {
            timeout_sec <- getOption("CMEnt.cluster_cleanup_timeout_sec", 2)
            # Keep cleanup bounded so we don't stall before any progress is shown.
            setTimeLimit(elapsed = timeout_sec, transient = TRUE)
            on.exit(setTimeLimit(elapsed = Inf, transient = FALSE), add = TRUE)
            tryCatch(
                reg$stopCluster(),
                error = function(e) invisible(NULL)
            )
        }
    }

    gc(verbose = FALSE)
    invisible()
}

#' @keywords internal
#' @noRd
.setupParallel <- function() {
    if (!dir.exists(.already_logged_dir)) {
        dir.create(.already_logged_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (nzchar(Sys.getenv("CI"))) {
        # CI-safe future backend
        .log_info("Running in CI environment, using sequential processing to avoid parallel issues in CI", level = 2)
        future::plan(future::sequential)
        return()
    }
    njobs <- getOption("CMEnt.njobs")
    if (njobs < 0) {
        njobs <- future::availableCores() + njobs
    }
    if (njobs > 1) {
        parallel_backend <- getOption("CMEnt.parallel_backend", "auto")
        parallel_backend <- as.character(parallel_backend)[1]
        if (is.na(parallel_backend) || !nzchar(parallel_backend)) {
            parallel_backend <- "auto"
        }
        parallel_backend <- tolower(parallel_backend)
        if (!parallel_backend %in% c("auto", "multicore", "multisession")) {
            warning("Unsupported CMEnt.parallel_backend='", parallel_backend, "'. Falling back to 'auto'.")
            parallel_backend <- "auto"
        }
        use_multicore <- parallel_backend == "multicore" ||
            (parallel_backend == "auto" && future::availableCores("multicore") > 1L)
        if (use_multicore) {
            if (!file.exists(.already_logged_file)) {
                .log_info("Using multicore parallelization with ", njobs, " workers", level = 2)
                writeLines("TRUE", con = .already_logged_file)
            }
            future::plan(future::multicore, workers = njobs)
        } else {
            if (!file.exists(.already_logged_file)) {
                .log_info("Using multisession parallelization with ", njobs, " workers", level = 2)
                writeLines("TRUE", con = .already_logged_file)
            }
            future::plan(future::multisession, workers = njobs)
        }
        withr::defer(future::plan(future::sequential), envir = parent.frame())
    } else {
        if (!file.exists(.already_logged_file)) {
            .log_info("Using sequential processing (njobs=1)", level = 2)
            writeLines("TRUE", con = .already_logged_file)
        }
        future::plan(future::sequential)
    }
}

#' @keywords internal
#' @noRd
.finalizeParallel <- function() {
    future::plan(future::sequential)
}

#' @keywords internal
#' @noRd
.calculateBetaStats <- function(beta_values, pheno, aggfun) {
    is_case <- pheno[, "__casecontrol__"] == 1
    cases <- beta_values[, is_case, drop = FALSE]
    cases <- as.matrix(cases)

    if (identical(aggfun, mean)) {
        cases_beta <- matrixStats::rowMeans2(cases, na.rm = TRUE)
    } else if (identical(aggfun, stats::median)) {
        cases_beta <- matrixStats::rowMedians(cases, na.rm = TRUE)
    } else {
        cases_beta <- apply(cases, 1, aggfun, na.rm = TRUE)
    }
    cases_sd <- matrixStats::rowSds(cases, na.rm = TRUE)
    cases_num <- matrixStats::rowCounts(!is.na(cases))
    rm(cases)

    is_ctrl <- !is_case
    ctrl <- beta_values[, is_ctrl, drop = FALSE]
    ctrl <- as.matrix(ctrl)
    if (identical(aggfun, mean)) {
        controls_beta <- matrixStats::rowMeans2(ctrl, na.rm = TRUE)
    } else if (identical(aggfun, stats::median)) {
        controls_beta <- matrixStats::rowMedians(ctrl, na.rm = TRUE)
    } else {
        controls_beta <- apply(ctrl, 1, aggfun, na.rm = TRUE)
    }
    controls_sd <- matrixStats::rowSds(ctrl, na.rm = TRUE)
    controls_num <- matrixStats::rowCounts(!is.na(ctrl))
    rm(ctrl)

    list(
        cases_beta = cases_beta,
        controls_beta = controls_beta,
        cases_beta_sd = cases_sd,
        controls_beta_sd = controls_sd,
        cases_num = cases_num,
        controls_num = controls_num
    )
}


#' @keywords internal
#' @noRd
.winsorizeVector <- function(x, probs = c(0.05, 0.95)) {
    if (!is.numeric(x) || length(x) == 0L) {
        return(x)
    }
    q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8)
    x[x < q[1]] <- q[1]
    x[x > q[2]] <- q[2]
    x
}

#' @keywords internal
#' @noRd
.momentStats <- function(x) {
    if (!is.numeric(x) || length(x) < 3L) {
        return(c(skew = NA_real_, excess_kurtosis = NA_real_))
    }
    s <- stats::sd(x)
    if (!is.finite(s) || s <= 1e-12) {
        return(c(skew = 0, excess_kurtosis = 0))
    }
    z <- (x - mean(x)) / s
    c(
        skew = mean(z^3),
        excess_kurtosis = mean(z^4) - 3
    )
}

#' @keywords internal
#' @noRd
.summarizeCorrelationAssumptions <- function(x_mat, y_mat, n_valid) {
    min_nvalid_q10 <- getOption("CMEnt.auto_pval_min_nvalid_q10", 10)
    corr_delta_threshold <- getOption("CMEnt.auto_pval_corr_delta_threshold", 0.10)
    skew_threshold <- getOption("CMEnt.auto_pval_skew_threshold", 2)
    excess_kurtosis_threshold <- getOption("CMEnt.auto_pval_excess_kurtosis_threshold", 7)
    pilot_max_pairs <- as.integer(getOption("CMEnt.auto_pval_pilot_pairs", 2000L))
    pilot_max_pairs <- max(1L, pilot_max_pairs)

    q10_n_valid <- if (length(n_valid) > 0L) {
        as.numeric(stats::quantile(n_valid, probs = 0.10, na.rm = TRUE, names = FALSE, type = 7))
    } else {
        NA_real_
    }
    if (!is.finite(q10_n_valid) || q10_n_valid < min_nvalid_q10) {
        return(list(
            use_parametric = FALSE,
            q10_n_valid = q10_n_valid,
            median_abs_delta_spearman = NA_real_,
            median_abs_delta_winsorized = NA_real_,
            median_abs_skew = NA_real_,
            median_excess_kurtosis = NA_real_,
            n_pairs_used = 0L,
            reason = "low effective sample size"
        ))
    }

    valid_idx <- which(n_valid >= 3L)
    if (length(valid_idx) == 0L) {
        return(list(
            use_parametric = FALSE,
            q10_n_valid = q10_n_valid,
            median_abs_delta_spearman = NA_real_,
            median_abs_delta_winsorized = NA_real_,
            median_abs_skew = NA_real_,
            median_excess_kurtosis = NA_real_,
            n_pairs_used = 0L,
            reason = "no valid pairs for diagnostics"
        ))
    }

    if (length(valid_idx) > pilot_max_pairs) {
        sel_pos <- unique(as.integer(round(seq(1, length(valid_idx), length.out = pilot_max_pairs))))
        valid_idx <- valid_idx[sel_pos]
    }

    delta_spearman <- rep(NA_real_, length(valid_idx))
    delta_winsorized <- rep(NA_real_, length(valid_idx))
    abs_skew <- rep(NA_real_, 2L * length(valid_idx))
    excess_kurtosis <- rep(NA_real_, 2L * length(valid_idx))
    k <- 0L
    for (i in seq_along(valid_idx)) {
        r <- valid_idx[i]
        xv <- x_mat[r, ]
        yv <- y_mat[r, ]
        ok <- !is.na(xv) & !is.na(yv)
        if (sum(ok) < 3L) {
            next
        }
        xv <- xv[ok]
        yv <- yv[ok]
        r_pearson <- suppressWarnings(stats::cor(xv, yv, method = "pearson"))
        r_spearman <- suppressWarnings(stats::cor(xv, yv, method = "spearman"))
        xw <- .winsorizeVector(xv)
        yw <- .winsorizeVector(yv)
        r_winsor <- suppressWarnings(stats::cor(xw, yw, method = "pearson"))
        if (is.finite(r_pearson) && is.finite(r_spearman)) {
            delta_spearman[i] <- abs(r_pearson - r_spearman)
        }
        if (is.finite(r_pearson) && is.finite(r_winsor)) {
            delta_winsorized[i] <- abs(r_pearson - r_winsor)
        }
        mx <- .momentStats(xv)
        my <- .momentStats(yv)
        k <- k + 1L
        abs_skew[(2L * k) - 1L] <- abs(mx["skew"])
        abs_skew[2L * k] <- abs(my["skew"])
        excess_kurtosis[(2L * k) - 1L] <- mx["excess_kurtosis"]
        excess_kurtosis[2L * k] <- my["excess_kurtosis"]
    }

    median_abs_delta_spearman <- median(delta_spearman, na.rm = TRUE)
    median_abs_delta_winsorized <- median(delta_winsorized, na.rm = TRUE)
    median_abs_skew <- median(abs_skew, na.rm = TRUE)
    median_excess_kurtosis <- median(excess_kurtosis, na.rm = TRUE)
    if (!is.finite(median_abs_delta_spearman)) {
        median_abs_delta_spearman <- Inf
    }
    if (!is.finite(median_abs_delta_winsorized)) {
        median_abs_delta_winsorized <- Inf
    }
    if (!is.finite(median_abs_skew)) {
        median_abs_skew <- Inf
    }
    if (!is.finite(median_excess_kurtosis)) {
        median_excess_kurtosis <- Inf
    }

    use_parametric <- (median_abs_delta_spearman <= corr_delta_threshold) &&
        (median_abs_delta_winsorized <= corr_delta_threshold) &&
        (median_abs_skew <= skew_threshold) &&
        (median_excess_kurtosis <= excess_kurtosis_threshold)
    reason <- if (use_parametric) {
        "assumptions look acceptable"
    } else {
        "assumptions not met"
    }
    list(
        use_parametric = use_parametric,
        q10_n_valid = q10_n_valid,
        median_abs_delta_spearman = median_abs_delta_spearman,
        median_abs_delta_winsorized = median_abs_delta_winsorized,
        median_abs_skew = median_abs_skew,
        median_excess_kurtosis = median_excess_kurtosis,
        n_pairs_used = k,
        reason = reason
    )
}

#' Load CMEnt Example Resources
#'
#' @description Helper function to load lightweight example resources bundled as
#' package data. The shipped local example bundle is a compact chr5/chr11 subset
#' intended for package examples, tests, and vignettes. If a local resource is not
#' found, the function can fall back to ExperimentHub.
#'
#' @param resource Character. Name of the resource to load. Available resources:
#' \itemize{
#'   \item "beta": Example beta values matrix
#'   \item "pheno": Example phenotype data
#'   \item "dmps": Example differentially methylated positions
#'   \item "array_type": Example array type annotation
#' }
#' @param use_experiment_hub Logical. Whether to attempt loading from ExperimentHub
#'   if the resource is not found locally (default: TRUE)
#'
#' @return The requested example object
#'
#' @details
#' The function follows this priority order:
#' \enumerate{
#'   \item Load the packaged example resource from the installed package `data/`
#'   \item If use_experiment_hub is TRUE, query ExperimentHub
#'   \item If all fail, raise an informative error
#' }
#'
#' Available local resources correspond to `.rda` files in `data/`:
#' \itemize{
#'   \item `beta.rda` - Example methylation beta values
#'   \item `pheno.rda` - Example phenotype/sample information
#'   \item `dmps.rda` - Example differential methylation results
#'   \item `array_type.rda` - Array platform type
#' }
#'
#' @examples
#' # Load beta values
#' beta <- loadExampleInputData("beta")
#'
#' # Load phenotype data
#' pheno <- loadExampleInputData("pheno")
#'
#' @export
loadExampleInputData <- function(resource, use_experiment_hub = TRUE) {
    available_resources <- c("beta", "pheno", "dmps", "array_type")
    if (!resource %in% available_resources) {
        stop(
            "Unknown resource: ", resource, "\n",
            "Available resources: ", paste(available_resources, collapse = ", ")
        )
    }

    verbose_setting <- getOption("CMEnt.verbose", 1)
    resource_files <- c(
        beta = "beta.rda",
        pheno = "pheno.rda",
        dmps = "dmps.rda",
        array_type = "array_type.rda"
    )
    resource_file <- unname(resource_files[[resource]])
    local_candidates <- c(
        system.file("data", resource_file, package = "CMEnt", mustWork = FALSE),
        file.path("data", resource_file)
    )
    local_path <- local_candidates[file.exists(local_candidates)][1]
    if (!is.na(local_path)) {
        if (verbose_setting >= 2) {
            .log_info("Loading ", resource, " from existing local path", level = 2)
        }
        load(local_path)
        return(get(resource))
    }

    # Fall back to ExperimentHub when local example files are unavailable.
    if (use_experiment_hub) {
        if (verbose_setting >= 2) {
            .log_info("Resource not found locally, attempting to load from ExperimentHub...", level = 2)
        }

        .assertDependencyRequirements(
            requirements = .experimentHubDependencyRequirements(resource = resource),
            context = "loadExampleInputData()"
        )

        tryCatch(
            {
                cache <- ExperimentHub::getExperimentHubOption("CACHE")
                dir.create(cache, showWarnings = FALSE, recursive = TRUE)
                eh <- ExperimentHub::ExperimentHub()
                cment_resources <- AnnotationHub::query(eh, "DMRsegaldata")
                resource_match <- cment_resources[grepl(resource, cment_resources$title, ignore.case = TRUE)]
                if (length(resource_match) > 0) {
                    if (verbose_setting >= 2) {
                        .log_success("Found resource in ExperimentHub", level = 2)
                    }
                    return(resource_match[[1]])
                }
                if (verbose_setting >= 2) {
                    .log_warn("Resource '", resource, "' not found in ExperimentHub")
                }
            },
            error = function(e) {
                if (verbose_setting >= 2) {
                    .log_warn("Failed to query ExperimentHub: ", e$message)
                }
            }
        )
    }

    stop(
        "Resource '", resource, "' was not found locally",
        if (isTRUE(use_experiment_hub)) " or in ExperimentHub" else "",
        ".\nAvailable resources: ", paste(available_resources, collapse = ", "),
        call. = FALSE
    )
}

#' @rdname loadExampleInputData
#' @description Compatibility wrapper returning the chr5/chr11 example subset.
#' The packaged local example resources already correspond to chr5/chr11, so this
#' wrapper is mainly kept for backwards compatibility with earlier releases and
#' still subsets ExperimentHub-backed resources when needed.
#' @export
loadExampleInputDataChr5And11 <- function(resource, use_experiment_hub = TRUE) {
    ret <- loadExampleInputData(resource, use_experiment_hub = use_experiment_hub)
    if (resource == "beta") {
        getBetaHandler(ret, array = "450k")$getBeta(chr = c("chr5", "chr11"))
    } else if (resource == "pheno") {
        ret
    } else if (resource == "dmps") {
        locs <- getSortedGenomicLocs(array = "450k")
        locs <- locs[locs$chr %in% c("chr5", "chr11"), ]
        locs <- locs[rownames(locs) %in% rownames(ret), , drop = FALSE]
        ret[rownames(locs), , drop = FALSE]
    } else if (resource == "array_type") {
        ret
    }
}
