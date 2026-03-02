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
                sample_group_control <- strsplit(sample_group_control, ",")[[1]]
            }
            subset_samplesheet[subset_samplesheet[, sample_group_col] %in% sample_group_control, "casecontrol"] <- 0
        } else {
            subset_samplesheet[, "casecontrol"] <- 0
            if (length(sample_group_case) == 1) {
                sample_group_case <- strsplit(sample_group_case, ",")[[1]]
            }
            subset_samplesheet[subset_samplesheet[, sample_group_col] %in% sample_group_case, "casecontrol"] <- 1
        }
    }
    subset_samplesheet
}

# Lightweight styled logging helpers -----------------------------------------

# Internal state for timing steps
.DMRsegal_log_env <- local({ # nolint
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

.col <- function(x, col = c("cyan", "green", "yellow", "blue")) {
    col <- strex::match_arg(col, ignore_case = TRUE)
    if (!.has_ansi()) {
        return(x)
    }
    switch(col,
        cyan   = cli::col_cyan(x),
        green  = cli::col_green(x),
        yellow = cli::col_yellow(x),
        blue   = cli::col_blue(x)
    )
}

#' Internal logging helpers using cli

#' @keywords internal
#' @noRd
.log_warn <- function(..., .envir = parent.frame()) {
    msg <- paste0(..., collapse = "")
    lead <- .col(cli::symbol$warning, "yellow")
    warning(paste(lead, msg))
    invisible()
}

#' @keywords internal
#' @noRd
.log_success <- function(..., .envir = parent.frame(), level = 1) {
    if (getOption("DMRsegal.verbose", 0) < level) {
        return(invisible())
    }
    if (level <= length(.DMRsegal_log_env$last_step_time)) {
        dur <- .fmt_dur(.DMRsegal_log_env$last_step_time[[level]])
    } else {
        .log_warn("No previous step time recorded for level ", level, " to calculate duration.")
        dur <- ""
    }
    msg <- paste0(paste0(..., collapse = ""), dur)
    lead <- paste(rep(" ", level - 1), .col(cli::symbol$tick, "green"), sep = "")
    message(paste(lead, msg))
    invisible()
}

#' @keywords internal
#' @noRd
.log_info <- function(..., .envir = parent.frame(), level = 1) {
    if (getOption("DMRsegal.verbose", 0) < level) {
        return(invisible())
    }
    # Suppress output from parallel workers
    if (!is.null(getOption("future.fork.enable")) && getOption("future.fork.enable", TRUE)) {
        if (exists(".Random.seed", envir = .GlobalEnv) && !identical(Sys.getpid(), getOption("future.main.pid", Sys.getpid()))) {
            return(invisible())
        }
    }
    msg <- paste0(..., collapse = "")
    lead <- paste(rep(" ", level - 1), .col(cli::symbol$info, "blue"), sep = "")
    message(paste(lead, msg))
    invisible()
}


#' @keywords internal
#' @noRd
.log_step <- function(..., .envir = parent.frame(), level = 1) {
    if (getOption("DMRsegal.verbose", 0) < level) {
        return(invisible())
    }
    .DMRsegal_log_env$last_step_time[[level]] <- Sys.time() # nolint
    msg <- paste0(..., collapse = "")
    lead <- paste(rep(" ", level - 1), .col(cli::symbol$arrow_right, "cyan"), sep = "")
    message(paste(lead, msg))
    invisible()
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
        sample_group_case <- strsplit(sample_group_case, ",")[[1]]
    }
    if (!is.null(sample_group_control)) {
        sample_group_control <- strsplit(sample_group_control, ",")[[1]]
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
    if (is.null(output_dir)) {
        use_cache <- getOption("DMRsegal.bed_cache_dir", NULL)
        if (!is.null(use_cache) && !isFALSE(use_cache)) {
            if (is.character(use_cache)) {
                cache_dir <- use_cache
            } else {
                cache_dir <- file.path(path.expand("~"), ".cache", "DMRsegal", "tabix_cache")
            }
        } else {
            cache_dir <- tempdir()
        }
    } else {
        cache_dir <- output_dir
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

getRegistry <- function(obj, indices = NULL, select = NULL, rename = NULL, derive = NULL, chunk_size = 100000) {
    if (is.data.frame(obj)) {
        return(.postProcessRegistry(obj, select = select, rename = rename, derive = derive, indices = indices))
    }
    cache_dir <- getOption("DMRsegal.h5_cache_dir")
    h5_file <- file.path(cache_dir, paste0(.getFileHash(obj), ".h5"))
    if (!file.exists(h5_file)) {
        createH5file(
            input_file = obj,
            output_h5file = h5_file,
            dataset_name = "data",
            select = select,
            chunk_size = chunk_size
        )
    }
    da <- HDF5Array::HDF5Array(h5_file, "data")
    x <- DelayedDataFrame::DelayedDataFrame(da)
    colnames(x) <- rhdf5::h5read(h5_file, "data_colnames")
    .postProcessRegistry(x, select = NULL, rename = rename, derive = derive, indices = indices)
}


#' Create Genomic Location Registry from Tabix BED File
#'
#' @description This function creates a Registry from a Tabix-indexed BED file.
#' @param input_tabix Character. Path to the Tabix-indexed BED file.
#' @param output_dir Character. Directory for caching processed files. If NULL, uses a default cache directory at `~/.cache/R/DMRsegal/tabix_cache/` (default: NULL)
#' @param num_rows Integer. Number of rows in the BED file. If NULL, the function will compute it automatically (default: NULL)
#' @param hash Character. Hash string for caching. If NULL, the function will compute it from the input file (default: NULL)
#' @param chunk_size Integer. Number of rows to process in each chunk for memory efficiency (default: 50000)
#' @return Returns a DelayedDataFrame object
#' @keywords internal
#' @noRd
genomicLocsFromTabix <- function(input_tabix, output_dir = NULL, num_rows = NULL, hash = NULL, chunk_size = 50000, use_id_as_rownames = FALSE, chrom_col = "#chrom", start_col = "start") { # nolint
    cache_dir <- .getTabixCacheDir(output_dir)
    if (is.null(hash)) {
        hash <- .getFileHash(input_tabix)
    }
    cache_file <- file.path(cache_dir, paste0("bed_locations_", hash, ".rds"))

    if (file.exists(cache_file) && getOption("DMRsegal.use_tabix_cache", FALSE)) {
        return(readRDS(cache_file))
    }
    renaming <- c("chr", "start")
    names(renaming) <- c(chrom_col, start_col)
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
            chunk_size = chunk_size
        )
    } else {
        sorted_locs <- getRegistry(
            input_tabix,
            select = c(chrom_col, start_col, "end", "id"),
            rename = renaming,
            indices = "id",
            chunk_size = chunk_size
        )
    }
    if (getOption("DMRsegal.use_tabix_cache", FALSE)) {
        saveRDS(sorted_locs, file = cache_file)
    }
    sorted_locs
}


#' Read and Process Custom Methylation BED Data
#'
#' @description Reads methylation data from a custom BED file format, converts it to
#' a tabix-indexed format for efficient random access, and creates genomic location
#' indices. This function is designed to handle custom methylation array data or
#' sequencing-based methylation data in BED format, making it compatible with the
#' DMRsegal workflow.
#'
#' @param bed_file Character. Path to the input BED file containing methylation data.
#'   The file should have chromosome and position columns, plus sample columns with
#'   methylation values. Can be gzipped (default: NULL)
#' @param pheno Data frame. Phenotype data with sample IDs as rownames. Only samples
#'   present in both the pheno rownames and BED file header will be processed
#' @param genome Character. Genome version to use (e.g., "hg19", "hg38") (default: "hg19")
#' @param chrom_col Character. Name of the chromosome column in the BED file
#'   (default: "#chrom")
#' @param start_col Character. Name of the start position column in the BED file
#'   (default: "start")
#' @param output_dir Character. Directory for caching processed files. If NULL, uses
#'   a default cache directory at `~/.cache/R/DMRsegal/tabix_cache/` (default: NULL)
#' @param chunk_size Integer. Number of rows to process in each chunk for memory
#'   efficiency (default: 50000)
#'
#' @return A list with two elements:
#' \itemize{
#'   \item tabix_file: Character path to the created tabix-indexed BED file
#'   \item locations_file: Character path to the RDS file containing genomic location indices
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
#'   \item Caches all processed files based on input file hash for reuse
#' }
#'
#' The cache directory structure includes:
#' \itemize{
#'   \item `bed_beta_<hash>.bed.gz`: Tabix-indexed BED file with methylation values
#'   \item `bed_beta_<hash>.bed.gz.tbi`: Tabix index file
#'   \item `bed_locations_<hash>.rds`: Serialized descriptor for loading locations
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
#'     output_dir = "/path/to/cache"
#' )
#'
#' # Access the processed files
#' tabix_file <- result$tabix_file
#' locations_file <- result$locations_file
#'
#' @seealso
#' \code{\link{convertBetaToTabix}} for converting standard beta files to tabix format
#' \code{\link{getBetaHandler}} for creating a BetaHandler object from processed files
#'
#' @export
readCustomMethylationBedData <- function(bed_file, pheno, genome = "hg19", chrom_col = "#chrom",
                                         start_col = "start", output_dir = NULL, chunk_size = 50000) {
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
    bed_header <- strsplit(readLines(bed_file, n = 1), "\t")[[1]]
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
    tabix_file_path <- file.path(cache_dir, paste0("bed_beta_", hash, ".bed.gz"))
    # Convert to tabix
    convertBetaToTabix(
        .bed_file = normalized_bed_file,
        output_file = tabix_file_path,
        chunk_size = chunk_size,
        njobs = 1
    )
    locations <- genomicLocsFromTabix(
        tabix_file_path,
        num_rows = num_rows,
        hash = hash,
        output_dir = cache_dir,
        chunk_size = chunk_size
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
#' chunk-based approach to handle large files and automatically uses a cache directory
#' to avoid redundant conversions of the same beta file.
#'
#' @param beta_file Character. Path to the input beta values file
#' @param sorted_locs Data frame with genomic locations containing 'chr' and 'start' columns.
#'   If NULL, will be retrieved automatically using getSortedGenomicLocs() (default: NULL)
#' @param array Character. Array platform type. Only used if sorted_locs is NULL (default: "450K")
#' @param genome Character. Genome version. Only used if  sorted_locs is NULL (default: "hg19")
#' @param output_file Character. Path for the output tabix file. If NULL, uses a cache
#'   directory in tempdir() with hash-based naming (default: NULL)
#' @param chunk_size Integer. Number of rows to process in each chunk (default: 50000)
#' @param njobs Integer. Number of parallel jobs for sorting (default: 1)
#'
#' @return Character. Path to the created tabix file, or NULL if conversion failed
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks if tabix and bgzip tools are available in the system PATH
#'   \item Computes a hash of the input beta file for cache lookup
#'   \item Checks if a cached tabix file exists; if so, returns it immediately
#'   \item Processes the beta file in chunks (50,000 rows at a time) to minimize memory usage
#'   \item Converts beta values to BED format with genomic coordinates
#'   \item Sorts, compresses (bgzip), and indexes (tabix) the file
#'   \item Saves to cache directory for future reuse
#' }
#'
#' The chunk-based processing ensures that even very large beta files (millions of CpGs)
#' can be converted without running out of memory. The cache directory is located at
#' \code{tempdir()/DMRsegal_tabix_cache/} and persists for the duration of the c session.
#' Files are named based on the MD5 hash of the input beta file, ensuring that identical
#' files reuse the same cached version.
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
                               genome = "hg19",
                               locations_file = NULL,
                               output_file = NULL,
                               chunk_size = 50000,
                               njobs = 1,
                               .bed_file = NULL) {
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

    # Set default output file name - use cache directory
    if (is.null(output_file)) {
        # Create cache directory in temp folder
        cache_dir <- getOption("DMRsegal.tabix_cache_dir", file.path(path.expand("~"), ".cache", "R", "DMRsegal", "tabix_cache"))
        if (!dir.exists(cache_dir)) {
            dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        }
        # Compute hash of beta file
        beta_hash <- .getFileHash(beta_file)

        # Create cache filename based on hash
        output_file <- file.path(cache_dir, paste0("beta_", beta_hash, ".bed.gz"))

        # Check if tabix file already exists in cache
        if (getOption("DMRsegal.use_tabix_cache", TRUE) && file.exists(output_file) && file.exists(paste0(output_file, ".tbi"))) {
            .log_info("Using cached tabix file: ", basename(output_file), level = 2)
            return(output_file)
        }
        if (!getOption("DMRsegal.use_tabix_cache", TRUE)) {
            output_file <- file.path(tempdir(), paste0("beta_", beta_hash, ".bed.gz"))
            .log_info("Tabix caching disabled; will create new tabix file at a temporary location: ", basename(output_file), level = 2)
        }
    }

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
                col_names <- strsplit(header_line, "\t")[[1]]

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

                .log_info("Processing ", n_rows, " CpG sites...", level = 2)

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

                    cpg_ids <- chunk_data[[1]]

                    # Match with genomic locations
                    common_cpgs <- intersect(cpg_ids, rownames(sorted_locs))

                    if (length(common_cpgs) > 0) {
                        # Create BED format for this chunk with 6 mandatory columns
                        bed_chunk <- as.data.frame(sorted_locs[common_cpgs, c("chr", "start"), drop = FALSE])
                        bed_chunk$end <- bed_chunk$start + 1
                        bed_chunk$id <- rownames(bed_chunk)
                        bed_chunk$score <- 0
                        bed_chunk$strand <- "*"
                        bed_chunk <- bed_chunk[, c("chr", "start", "end", "id", "score", "strand")]

                        # Add beta values as additional columns
                        beta_subset <- chunk_data[match(common_cpgs, cpg_ids), -1, drop = FALSE]
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
                    .log_warn("No common CpGs found between beta file and genomic locations")
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
                        parts1 <- strsplit(line1, "\t", fixed = TRUE)[[1]]
                        parts2 <- strsplit(line2, "\t", fixed = TRUE)[[1]]

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
#' the beta file, sorts the CpG sites according to their genomic positions using array annotation,
#' and writes the sorted data to a new file.
#'
#' @param beta_file Character. Path to the input beta values file to be sorted
#' @param output_file Character. Path for the output sorted beta file (default: adds "_sorted" suffix)
#' @param array Character. Array platform type (default: "450K")
#' @param genome Character. Genome version (default: "hg19")
#' @param genomic_locs Data frame. Optional pre-computed genomic locations. If NULL, locations will be retrieved automatically (default: NULL)
#' @param overwrite Logical. Whether to overwrite existing output file (default: FALSE)
#'
#' @return Character. Path to the sorted output file
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Reads the beta values file
#'   \item Loads the appropriate array annotation (450K or EPIC)
#'   \item Sorts CpG sites by genomic coordinates (chr:start)
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
                                      genome = "hg19",
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
    beta_data <- data.table::fread(beta_file, header = TRUE, data.table = FALSE, showProgress = getOption("DMRsegal.verbose", 0) > 1)

    # Get row names (CpG IDs) from first column
    cpg_ids <- beta_data[[1]]
    beta_values <- beta_data[, -1, drop = FALSE]
    rownames(beta_values) <- cpg_ids
    .log_success("Beta loaded: ", nrow(beta_values), " CpGs across ", ncol(beta_values), " samples", level = 2)

    sorted_locs <- genomic_locs
    if (is.null(sorted_locs)) {
        array <- strex::match_arg(array, ignore_case = TRUE)
        sorted_locs <- getSortedGenomicLocs(array = array, genome = genome)
    }


    # Find CpGs that are present in both the beta file and array annotation
    common_cpgs <- intersect(cpg_ids, rownames(sorted_locs))
    missing_from_annotation <- setdiff(cpg_ids, rownames(sorted_locs))
    if (length(missing_from_annotation) > 0) {
        stop(
            "Found ", length(missing_from_annotation), " CpG sites in beta file that are not in ",
            array, " annotation. First 5 missing: ", paste(head(missing_from_annotation, 5), collapse = ", ")
        )
    }

    missing_from_beta <- setdiff(rownames(sorted_locs), cpg_ids)
    if (length(missing_from_beta) > 0) {
        .log_info("Note: ", length(missing_from_beta), " CpGs in ", array, " annotation are missing from beta file", level = 2)
    }

    final_order <- rownames(sorted_locs)[rownames(sorted_locs) %in% common_cpgs]

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
    signal - effect %*% t_covariate_matrix
}

.transformBeta <- function(beta, pheno, covariates = NULL, covariate_model = NULL) {
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
    cache_dir <- getOption("DMRsegal.annotation_cache_dir", file.path(
        path.expand("~"),
        ".cache", "R", "DMRsegal", "annotations"
    ))
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    chain_name <- paste0(from_genome, "To", stringr::str_to_title(to_genome), ".over.chain")
    chain_file <- file.path(cache_dir, chain_name)
    if (!file.exists(chain_file)) {
        utils::download.file(
            url = paste0(
                "http://hgdownload.soe.ucsc.edu/goldenPath/",
                from_genome, "/liftOver/", chain_name, ".gz"
            ),
            destfile = paste0(chain_file, ".gz"), mode = "wb"
        )
        R.utils::gunzip(paste0(chain_file, ".gz"), overwrite = TRUE)
    }
    chain <- rtracklayer::import.chain(chain_file)
    lifted <- rtracklayer::liftOver(granges, chain)
    lifted_unlisted <- unlist(lifted)
    lifted_unlisted
}


#' Get Sorted Array Locations
#'
#' @description Retrieves and sorts genomic location annotations for the specified
#' methylation array platform and genome version. Performs liftOver if necessary.
#' The function caches the results.
#'
#' @param array Character. Array platform type (supported: "450K", "EPIC", "EPICv2", "27K", "Mouse", 'NULL'), ignored when locations_file is provided. Must be 'NULL' when the experiment is not array-based.
#' @param genome Character. Genome version (supported: "hg19", "hg38", "mm10", "mm39"), ignored if locations_file is provided
#' @param locations_file Character. Optional path to a precomputed locations file (RDS format). If provided, this file will be used directly (default: NULL)
#'
#' @return A data frame containing sorted genomic locations with rownames as CpG IDs and columns:
#' \itemize{
#'   \item chr: Chromosome
#'   \item start: Genomic position
#'   \item start: Start position (same as start)
#'   \item end: End position (start + 1)
#' }
#'
#' @examples
#' # Get sorted locations for 450K array (hg19)
#' locs_450k <- getSortedGenomicLocs("450K")
#'
#' # Get sorted locations for EPIC array with hg38
#' locs_epic <- getSortedGenomicLocs("EPIC", "hg38")
#'
#' # Get sorted locations for EPICv2 array
#' locs_epicv2 <- getSortedGenomicLocs("EPICv2", "hg38")
#'
#' @export
getSortedGenomicLocs <- function(array = c("450K", "27K", "EPIC", "EPICv2", "Mouse"), genome = c("hg19", "hg38", "mm10", "mm39"), locations_file = NULL) {
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
    cache_dir <- getOption("DMRsegal.annotation_cache_dir", file.path(
        path.expand("~"),
        ".cache", "R", "DMRsegal", "annotations"
    ))
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }

    array <- tolower(array)
    genome <- tolower(genome)
    cache_file <- file.path(cache_dir, paste0(
        array, "_", genome,
        "_locations.rds"
    ))
    if (getOption("DMRsegal.use_annotation_cache", TRUE) && file.exists(cache_file)) {
        .log_info("Using cached annotation file: ", basename(cache_file), level = 3)
        locs <- readRDS(cache_file)
        return(locs)
    }
    pkg_name <- NULL
    if (genome %in% c("hg19", "hg38")) {
        if (array == "450k") {
            pkg_name <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
        } else if (array == "epic") {
            pkg_name <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
        } else if (array == "epicv2") {
            pkg_name <- "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
        } else if (array == "27k") {
            pkg_name <- "IlluminaHumanMethylation27kanno.ilmn12.hg19"
        } else {
            stop(
                paste0(
                    "Incorrect array and genome combination was provided. ",
                    "For hg19 and hg38, ('450K','EPIC','EPICv2','27K') arrays are supported."
                )
            )
        }
    } else if (genome %in% c("mm10", "mm39")) {
        if (array != "Mouse") {
            stop("Incorrect array and genome combination was provided. For mm10 and mm39 only 'Mouse' array is supported.")
        }
        pkg_name <- "IlluminaMouseMethylationanno.12.v1.mm10"
        if (!requireNamespace(pkg_name, quietly = TRUE)) {
            if (!requireNamespace("devtools", quietly = TRUE)) {
                install.packages("devtools")
            }
            devtools::install_github(paste0("chiaraherzog/", pkg_name))
        }
    }
    if (is.null(pkg_name)) {
        stop(
            "Unsupported array/genome combination: ", array,
            "/", genome
        )
    }
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
        .log_info(
            "Installing required annotation package: ",
            pkg_name
        )
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        BiocManager::install(pkg_name)
    }
    locs <- minfi::getLocations(pkg_name)
    from_genome <- NULL
    if (genome == "mm39") {
        from_genome <- "mm10"
    }
    if (genome == "hg38") {
        if (tolower(array) != "epicv2") {
            from_genome <- "hg19"
        }
    } else {
        if (tolower(array) == "epicv2") {
            from_genome <- "hg38"
        }
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
            saveRDS(locs, cache_file)
        },
        warning = function(w) {
            .log_warn(
                "Could not write annotation cache file '", cache_file,
                "' (warning: ", conditionMessage(w), "). Continuing without cache persistence."
            )
        },
        error = function(e) {
            .log_warn(
                "Could not write annotation cache file '", cache_file,
                "' (error: ", conditionMessage(e), "). Continuing without cache persistence."
            )
        }
    )
    locs
}


#' Order Indices by Genomic Location
#'
#' @description Orders a vector of indices according to their corresponding genomic
#' locations (chromosome and position). This function is useful for sorting CpG
#' sites or other genomic features by their physical positions.
#'
#' @param x Character or integer vector. Indices or identifiers to be ordered
#' @param array Character. Array platform type, either "450K" or "EPIC" (default: "450K")
#' @param genome Character. Genome version, either "hg19" or "hg38" (default: "hg19")
#' @param genomic_locs Data frame. Optional pre-computed genomic locations. If NULL,
#' locations will be retrieved using getSortedGenomicLocs (default: NULL)
#'
#' @return Integer vector of ordered indices
#'
#' @examples
#' # Order CpG indices by genomic location
#' cpg_ids <- c("cg00000029", "cg00000108", "cg00000109")
#' ordered_indices <- orderByLoc(cpg_ids, array = "450K")
#'
#' # Order using pre-computed genomic locations
#' locs <- getSortedGenomicLocs("EPIC", "hg38")
#' ordered_indices <- orderByLoc(cpg_ids, genomic_locs = locs)
#'
#' @export
orderByLoc <- function(x,
                       array = c("450K", "27K", "EPIC", "EPICv2"),
                       genome = c("hg19", "hg38", "mm10", "mm39"),
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
    vals <- unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE)
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
.downsampleFlankIndices <- function(indices, max_sup_cpgs_per_dmr_side) {
    if (is.null(max_sup_cpgs_per_dmr_side) || max_sup_cpgs_per_dmr_side <= 0) {
        return(indices)
    }
    if (length(indices) <= max_sup_cpgs_per_dmr_side) {
        return(indices)
    }
    step <- ceiling(length(indices) / max_sup_cpgs_per_dmr_side)
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
    if (genome == "hg19") {
        pkg_name <- "BSgenome.Hsapiens.UCSC.hg19"
    } else if (genome == "hg38") {
        pkg_name <- "BSgenome.Hsapiens.UCSC.hg38"
    } else if (genome == "mm10") {
        pkg_name <- "BSgenome.Mmusculus.UCSC.mm10"
    } else if (genome == "mm39") {
        pkg_name <- "BSgenome.Mmusculus.UCSC.mm39"
    } else {
        if (!requireNamespace("BSgenome", quietly = TRUE)) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager")
            }
            BiocManager::install("BSgenome")
        }
        available_genomes <- BSgenome::available.genomes()
        matched_genome <- grep(paste0("^BSgenome\\..*\\.UCSC\\.", genome, "$"), available_genomes, value = TRUE)
        if (length(matched_genome) > 0) {
            pkg_name <- matched_genome[1]
        } else {
            pkg_name <- NULL
        }
    }
    if (is.null(pkg_name)) {
        return(NULL)
    }
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
        if (!requireNamespace(pkg_name, quietly = TRUE)) {
            .log_warn("BSgenome package not available: ", pkg_name)
            .log_info("Attempting to install...")
            tryCatch(
                {
                    if (!requireNamespace("BiocManager", quietly = TRUE)) {
                        install.packages("BiocManager")
                    }
                    BiocManager::install(pkg_name, update = FALSE)
                },
                error = function(e) {
                    .log_warn("Installation failed")
                    pkg_name <<- NULL
                }
            )
        }
    }
    pkg_name
}

#' Extract DNA Sequences for DMRs
#'
#' @description Retrieves the DNA sequences corresponding to genomic regions
#' specified in a GRanges object. This function is useful for extracting the
#' actual DNA sequence of identified DMRs for downstream analyses such as
#' motif finding or sequence composition analysis.
#'
#' @param dmrs GRanges object containing genomic coordinates of DMRs
#' @param genome Character. Genome version to use for sequence extraction, .e.g. "hg19".
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
#'   \item mm10: BSgenome.Mmusculus.UCSC.mm10
#'   \item mm39: BSgenome.Mmusculus.UCSC.mm39
#' }
#'
#' If the required BSgenome package is not installed and installation fails,
#' the function will automatically fall back to querying sequences from the
#' UCSC Genome Browser REST API. The online method processes sequences in
#' batches with optional parallel processing for improved performance with large datasets.
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
        pkg_name <- .getBSGenomePackage(genome)
        use_bsgenome <- !is.null(pkg_name)
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
getCpGBackgroundCounts <- function(regions, genome, njobs = 1, canonical_chr = TRUE) {
    pkg_name <- .getBSGenomePackage(genome)
    if (is.null(pkg_name)) {
        sequences <- getDMRSequences(regions, genome, use_online = TRUE, njobs = njobs)
        cpg_counts <- sapply(sequences, function(seq) {
            if (is.na(seq)) {
                return(NA_integer_)
            }
            stringr::str_count(seq, "CG")
        })
        return(unlist(cpg_counts))
    }
    cache_dir <- getOption(
        "DMRsegal.annotation_cache_dir", file.path(
            path.expand("~"),
            ".cache", "R", "DMRsegal", "annotations"
        )
    )
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    cpg_positions_file <- file.path(cache_dir, paste0(
        pkg_name, "_cpg_positions.rds"
    ))
    if (canonical_chr) {
        cpg_positions_file <- file.path(cache_dir, paste0(
            pkg_name, "_canonical_chr_cpg_positions.rds"
        ))
    }
    if (file.exists(cpg_positions_file)) {
        cpgs <- readRDS(cpg_positions_file)
    } else {
        .log_info("Missing CpG positions cache. Generating CpG positions from BSgenome package...", level = 2)
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
            cpgs <- do.call(
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
        cpgs <- data.table::as.data.table(as.data.frame(cpgs, stringsAsFactors = FALSE))[, 1:2]
        colnames(cpgs) <- c("chr", "start")
        cpgs[, "end"] <- cpgs[, "start"]
        data.table::setkey(cpgs, chr, start, end)
        saveRDS(cpgs, cpg_positions_file)
    }
    regions <- as.data.frame(regions, stringsAsFactors = FALSE)[, 1:3]
    regions <- data.table::as.data.table(regions)
    colnames(regions) <- c("rchr", "rstart", "rend")
    regions[, "id"] <- seq_len(nrow(regions))
    data.table::setkey(regions, rchr, rstart, rend)
    cpg_counts <- data.table::foverlaps(regions,
        cpgs,
        by.x = c("rchr", "rstart", "rend"),
        by.y = c("chr", "start", "end"),
    )[, .N, by = id]
    cpg_counts <- cpg_counts[order(cpg_counts$id), "N"]
    unlist(cpg_counts)
    # }
}

#' Query DNA Sequences from UCSC Genome Browser API
#'
#' @description Internal helper function to retrieve DNA sequences from the
#' UCSC Genome Browser REST API when BSgenome packages are not available.
#' Uses parallel processing and batched requests for improved performance.
#'
#' @param dmrs GRanges object containing genomic coordinates
#' @param genome Character. Genome version (e.g., "hg19", "mm39")
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
        if (!requireNamespace("parallel", quietly = TRUE)) {
            warning("parallel package not available, using single core processing")
            njobs <- 1
        }
    }

    if (njobs > 1 && n_batches > 1) {
        cl <- parallel::makeCluster(min(njobs, n_batches))
        on.exit(parallel::stopCluster(cl), add = TRUE)

        parallel::clusterExport(cl, c("base_url", "genome", "batch_size", "n_dmrs", "dmrs"),
            envir = environment()
        )
        parallel::clusterEvalQ(cl, {
            library(GenomicRanges)
            library(GenomeInfoDb)
            library(jsonlite)
        })

        batch_results <- parallel::parLapply(cl, seq_len(n_batches), fetch_batch)
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

#' Annotate DMRs with Gene Information
#'
#' @description Annotates DMRs with overlapping gene promoters and gene bodies
#' using TxDb annotations. For each DMR, identifies genes whose promoters or
#' gene bodies overlap with the DMR coordinates.
#'
#' @param dmrs Dataframe or GRanges object containing DMR coordinates
#' @param genome Character. Genome version to use for gene annotation. (default: "hg19")
#' @param promoter_upstream Integer. Number of base pairs upstream of TSS to
#'   define promoter region (default: 2000)
#' @param promoter_downstream Integer. Number of base pairs downstream of TSS
#'   to define promoter region (default: 200)
#'
#' @return The input Dataframe/GRanges object with additional metadata columns:
#' \itemize{
#'   \item in_promoter_of: Character vector of gene symbols with promoters overlapping the DMR (comma-separated)
#'   \item in_gene_body_of: Character vector of gene symbols with gene bodies overlapping the DMR (comma-separated)
#' }
#'
#' @details
#' The function uses genome-appropriate TxDb packages.
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
#' dmrs_annotated <- annotateDMRsWithGenes(dmrs, genome = "hg19")
#'
#' # Use custom promoter definition
#' dmrs_annotated <- annotateDMRsWithGenes(
#'     dmrs,
#'     genome = "hg38",
#'     promoter_upstream = 5000,
#'     promoter_downstream = 1000
#' )
#'
#' @export
annotateDMRsWithGenes <- function(dmrs, genome = "hg19",
                                  promoter_upstream = 2000,
                                  promoter_downstream = 200) {
    cache_dir <- getOption("DMRsegal.annotation_cache_dir", file.path(
        path.expand("~"),
        ".cache", "R", "DMRsegal", "annotations"
    ))
    dmrs_df_provided <- is.data.frame(dmrs)
    dmrs <- convertToGRanges(dmrs, genome)

    supported_organisms <- Organism.dplyr::supportedOrganisms()
    # find row matching the genome
    matched_row <- supported_organisms[grepl(tolower(genome), tolower(supported_organisms$TxDb)), , drop = FALSE]
    if (length(matched_row) == 0) {
        stop("Unsupported genome: ", genome)
    }
    txdb_pkg <- matched_row$TxDb
    orgdb_pkg <- matched_row$OrgDb

    # Load required packages
    if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        BiocManager::install(txdb_pkg)
    }

    if (!requireNamespace(orgdb_pkg, quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        BiocManager::install(orgdb_pkg)
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
        genes_file <- file.path(cache_dir, paste0("genes_", genome, ".rds"))
        if (file.exists(genes_file) && getOption("DMRsegal.use_annotation_cache", TRUE)) {
            .log_info("Loading cached genes from ", genes_file, level = 2)
            genes <- readRDS(genes_file)
        } else {
            genes <- GenomicFeatures::genes(txdb)
            # get genes only within standard chromosomes
            std_chroms <- GenomeInfoDb::standardChromosomes(GenomeInfoDb::seqinfo(txdb))
            genes <- genes[as.character(GenomeInfoDb::seqnames(genes)) %in% std_chroms]
            tryCatch(
                {
                    saveRDS(genes, genes_file)
                },
                warning = function(w) {
                    .log_warn(
                        "Could not write annotation cache file '", genes_file,
                        "' (warning: ", conditionMessage(w), "). Continuing without cache persistence."
                    )
                },
                error = function(e) {
                    .log_warn(
                        "Could not write annotation cache file '", genes_file,
                        "' (error: ", conditionMessage(e), "). Continuing without cache persistence."
                    )
                }
            )
        }
        promoters_file <- file.path(cache_dir, paste0("promoters_", genome, ".rds"))
        if (file.exists(promoters_file) && getOption("DMRsegal.use_annotation_cache", TRUE)) {
            .log_info("Loading cached promoters from ", promoters_file, level = 2)
            promoters <- readRDS(promoters_file)
        } else {
            transcripts_by_gene <- GenomicFeatures::transcriptsBy(txdb, by = "gene")
            transcripts_by_gene <- transcripts_by_gene[names(transcripts_by_gene) %in% names(genes)]
            promoters <- GenomicFeatures::promoters(transcripts_by_gene, upstream = promoter_upstream, downstream = promoter_downstream)
            promoters <- stack(promoters)
            tryCatch(
                {
                    saveRDS(promoters, promoters_file)
                },
                warning = function(w) {
                    .log_warn(
                        "Could not write annotation cache file '", promoters_file,
                        "' (warning: ", conditionMessage(w), "). Continuing without cache persistence."
                    )
                },
                error = function(e) {
                    .log_warn(
                        "Could not write annotation cache file '", promoters_file,
                        "' (error: ", conditionMessage(e), "). Continuing without cache persistence."
                    )
                }
            )
        }
    })
    .log_success("Gene annotations loaded: ", length(genes), " genes", level = 2)
    .log_step("Finding overlaps with promoters and gene bodies...", level = 2)

    # Find overlaps
    promoter_overlaps <- GenomicRanges::findOverlaps(dmrs, promoters)
    gene_body_overlaps <- GenomicRanges::findOverlaps(dmrs, genes)

    # Load org.db namespace
    if (!isNamespaceLoaded(orgdb_pkg)) {
        loadNamespace(orgdb_pkg)
    }

    # Get gene symbols - the main object has the same name as the package
    orgdb <- getExportedValue(orgdb_pkg, orgdb_pkg)

    # Extract Entrez IDs
    promoter_regions <- promoters[S4Vectors::subjectHits(promoter_overlaps)]

    promoter_entrez <- as.character(mcols(promoter_regions)$name)
    gene_body_entrez <- names(genes)[S4Vectors::subjectHits(gene_body_overlaps)]

    # Initialize annotation columns
    n_dmrs <- length(dmrs)
    dmrs$in_promoter_of <- rep(NA_character_, n_dmrs)
    dmrs$in_gene_body_of <- rep(NA_character_, n_dmrs)
    .log_step("Mapping overlapping Entrez IDs to gene symbols...", level = 2)
    # Convert Entrez IDs to symbols only if there are overlaps
    promoter_symbols <- character(0)
    if (length(promoter_entrez) > 0 && any(!is.na(promoter_entrez))) {
        valid_promoter_entrez <- promoter_entrez[!is.na(promoter_entrez) & promoter_entrez != ""]
        if (length(valid_promoter_entrez) > 0) {
            promoter_symbols <- suppressMessages(AnnotationDbi::mapIds(orgdb,
                keys = valid_promoter_entrez,
                column = "SYMBOL",
                keytype = "ENTREZID",
                multiVals = "first"
            ))
            names(promoter_symbols) <- valid_promoter_entrez
        }
    }

    gene_body_symbols <- character(0)
    if (length(gene_body_entrez) > 0 && any(!is.na(gene_body_entrez))) {
        valid_gene_body_entrez <- gene_body_entrez[!is.na(gene_body_entrez) & gene_body_entrez != ""]
        if (length(valid_gene_body_entrez) > 0) {
            gene_body_symbols <- suppressMessages(AnnotationDbi::mapIds(orgdb,
                keys = valid_gene_body_entrez,
                column = "SYMBOL",
                keytype = "ENTREZID",
                multiVals = "first"
            ))
            names(gene_body_symbols) <- valid_gene_body_entrez
        }
    }
    .log_success("Gene symbols mapped", level = 2)
    # Aggregate gene symbols for each DMR
    if (length(promoter_overlaps) > 0 && length(promoter_symbols) > 0) {
        promoter_entrez_char <- as.character(promoter_entrez)
        promoter_symbols_mapped <- promoter_symbols[promoter_entrez_char]
        promoter_by_dmr <- split(
            promoter_symbols_mapped,
            S4Vectors::queryHits(promoter_overlaps)
        )

        for (i in names(promoter_by_dmr)) {
            idx <- as.integer(i)
            genes_vec <- unique(na.omit(promoter_by_dmr[[i]]))
            if (length(genes_vec) > 0) {
                dmrs$in_promoter_of[idx] <- paste(genes_vec, collapse = ",")
            }
        }
    }

    if (length(gene_body_overlaps) > 0 && length(gene_body_symbols) > 0) {
        gene_body_entrez_char <- as.character(gene_body_entrez)
        gene_body_symbols_mapped <- gene_body_symbols[gene_body_entrez_char]
        gene_body_by_dmr <- split(
            gene_body_symbols_mapped,
            S4Vectors::queryHits(gene_body_overlaps)
        )

        for (i in names(gene_body_by_dmr)) {
            idx <- as.integer(i)
            genes_vec <- unique(na.omit(gene_body_by_dmr[[i]]))
            if (length(genes_vec) > 0) {
                dmrs$in_gene_body_of[idx] <- paste(genes_vec, collapse = ",")
            }
        }
    }
    if (dmrs_df_provided) {
        dmrs <- convertToDataFrame(dmrs)
    }
    dmrs
}

convertToGRanges <- function(obj, genome) {
    input_is_df <- !inherits(obj, "GRanges")
    if (input_is_df) {
        # If the chr prefix is missing, add it (e.g. "1" -> "chr1")
        if (!any(grepl("^chr", obj[[1]]))) {
            obj[[1]] <- paste0("chr", obj[[1]])
            obj$chr_prefix_added <- TRUE
        }
        # if the chromosome appears in the form of chr1:1230, save the original location in a separate column and parse the location into chr, start, end
        if (any(grepl(":", obj[[1]]))) {
            obj$original_location <- obj[[1]]
            loc_split <- strsplit(as.character(obj[[1]]), ":", fixed = TRUE)
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
.already_logged_file <- file.path(.already_logged_dir, "dmrsegal_already_logged_parallel.txt")
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
    if (isTRUE(getOption("DMRsegal.force_cluster_cleanup", FALSE))) {
        reg <- tryCatch(
            getFromNamespace("clusterRegistry", "future"),
            error = function(e) NULL
        )
        if (is.list(reg) && is.function(reg$stopCluster)) {
            timeout_sec <- getOption("DMRsegal.cluster_cleanup_timeout_sec", 2)
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

    njobs <- getOption("DMRsegal.njobs")
    if (njobs < 0) {
        njobs <- future::availableCores() + njobs
    }
    if (njobs > 1) {
        if (future::availableCores("multicore") > 1L) {
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
    cases <- as.matrix(cases, ncol = ncol(cases))

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
    ctrl <- as.matrix(ctrl, ncol = ncol(ctrl))
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
.resolveAdaptiveMinCpgDeltaBeta <- function(seeds_beta, pheno, aggfun, base_threshold = 0.1) {
    if (!is.finite(base_threshold)) {
        return(base_threshold)
    }
    case_mask <- pheno[, "__casecontrol__"] == 1
    ctrl_mask <- pheno[, "__casecontrol__"] == 0
    if (sum(case_mask, na.rm = TRUE) == 0L || sum(ctrl_mask, na.rm = TRUE) == 0L) {
        return(base_threshold)
    }
    beta_mat <- as.matrix(seeds_beta)
    max_rows <- as.integer(getOption("DMRsegal.adaptive_min_cpg_delta_beta_max_rows", 50000L))
    max_rows <- max(1L, max_rows)
    if (nrow(beta_mat) > max_rows) {
        sel <- unique(as.integer(round(seq(1, nrow(beta_mat), length.out = max_rows))))
        beta_mat <- beta_mat[sel, , drop = FALSE]
    }

    if (identical(aggfun, mean)) {
        cases <- matrixStats::rowMeans2(beta_mat[, case_mask, drop = FALSE], na.rm = TRUE)
        ctrls <- matrixStats::rowMeans2(beta_mat[, ctrl_mask, drop = FALSE], na.rm = TRUE)
    } else if (identical(aggfun, stats::median)) {
        cases <- matrixStats::rowMedians(beta_mat[, case_mask, drop = FALSE], na.rm = TRUE)
        ctrls <- matrixStats::rowMedians(beta_mat[, ctrl_mask, drop = FALSE], na.rm = TRUE)
    } else {
        cases <- apply(beta_mat[, case_mask, drop = FALSE], 1, aggfun, na.rm = TRUE)
        ctrls <- apply(beta_mat[, ctrl_mask, drop = FALSE], 1, aggfun, na.rm = TRUE)
    }
    deltas <- abs(cases - ctrls)
    deltas <- deltas[is.finite(deltas)]
    if (length(deltas) == 0L) {
        return(base_threshold)
    }
    qprob <- getOption("DMRsegal.adaptive_min_cpg_delta_beta_quantile", 0.25)
    adaptive <- as.numeric(stats::quantile(deltas, probs = qprob, na.rm = TRUE, names = FALSE, type = 7))
    if (!is.finite(adaptive)) {
        return(base_threshold)
    }
    max(base_threshold, adaptive)
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
    min_nvalid_q10 <- getOption("DMRsegal.auto_pval_min_nvalid_q10", 10)
    corr_delta_threshold <- getOption("DMRsegal.auto_pval_corr_delta_threshold", 0.10)
    skew_threshold <- getOption("DMRsegal.auto_pval_skew_threshold", 2)
    excess_kurtosis_threshold <- getOption("DMRsegal.auto_pval_excess_kurtosis_threshold", 7)
    pilot_max_pairs <- as.integer(getOption("DMRsegal.auto_pval_pilot_pairs", 2000L))
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

#' Load DMRsegal Data Resources
#'
#' @description Helper function to load data resources from the package data folder.
#' Resources are lazy loaded when the package is loaded, so they exist in the package
#' namespace. If a resource is not found in the package, it will attempt to query
#' ExperimentHub. If that also fails, it provides instructions to generate the data.
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
#' @return The requested data object, or NULL if not found
#'
#' @details
#' The function follows this priority order:
#' \enumerate{
#'   \item Check if the resource exists in the package namespace (lazy loaded)
#'   \item If not found, try loading directly from the data file
#'   \item If use_experiment_hub is TRUE, query ExperimentHub
#'   \item If all fail, provide instructions for generating the data
#' }
#'
#' Available resources correspond to .rda files in the data/ folder:
#' \itemize{
#'   \item beta.rda - Example methylation beta values
#'   \item pheno.rda - Example phenotype/sample information
#'   \item dmps.rda - Example differential methylation results
#'   \item array_type.rda - Array platform type
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

    # First, try using data() to load the resource
    # This works even during covr::package_coverage()
    verbose_setting <- getOption("DMRsegal.verbose", 1)
    tryCatch(
        {
            # Suppress data() output
            invisible(utils::data(list = resource, package = "DMRsegal", envir = environment()))
            if (exists(resource, inherits = FALSE)) {
                if (verbose_setting >= 2) {
                    .log_info("Loaded ", resource, " using data()", level = 2)
                }
                return(get(resource, inherits = FALSE))
            }
        },
        error = function(e) {
            # Silently continue to next method
        }
    )

    # Second, check if the resource exists in the package namespace (lazy loaded)
    if (exists(resource, envir = asNamespace("DMRsegal"), inherits = FALSE)) {
        if (verbose_setting >= 2) {
            .log_info("Loading ", resource, " from package namespace (lazy loaded)", level = 2)
        }
        return(get(resource, envir = asNamespace("DMRsegal"), inherits = FALSE))
    }

    # Third, try to load from data file directly
    data_file <- paste0(resource, ".rda")
    data_path <- system.file("data", data_file, package = "DMRsegal", mustWork = FALSE)
    if (file.exists(data_path)) {
        if (verbose_setting >= 2) {
            .log_info("Loading ", resource, " from package data file", level = 2)
        }
        env <- new.env()
        load(data_path, envir = env)
        if (exists(resource, envir = env)) {
            return(get(resource, envir = env))
        }
    }

    # Fourth, try ExperimentHub if enabled
    if (use_experiment_hub) {
        if (verbose_setting >= 2) {
            .log_info("Resource not found locally, attempting to load from ExperimentHub...", level = 2)
        }

        if (!requireNamespace("ExperimentHub", quietly = TRUE)) {
            if (verbose_setting >= 2) {
                .log_warn("ExperimentHub package not available. Install with: BiocManager::install('ExperimentHub')")
            }
        } else {
            tryCatch(
                {
                    cache <- ExperimentHub::getExperimentHubOption("CACHE")
                    dir.create(cache, showWarnings = FALSE, recursive = TRUE)
                    eh <- ExperimentHub::ExperimentHub()
                    # Query for DMRsegaldata resources
                    dmrsegal_resources <- ExperimentHub::query(eh, "DMRsegaldata")
                    # Find the specific resource
                    resource_match <- dmrsegal_resources[grepl(resource, dmrsegal_resources$title, ignore.case = TRUE)]
                    if (length(resource_match) > 0) {
                        if (verbose_setting >= 2) {
                            .log_success("Found resource in ExperimentHub", level = 2)
                        }
                        return(resource_match[[1]])
                    } else {
                        if (verbose_setting >= 2) {
                            .log_warn("Resource '", resource, "' not found in ExperimentHub")
                        }
                    }
                },
                error = function(e) {
                    if (verbose_setting >= 2) {
                        .log_warn("Failed to query ExperimentHub: ", e$message)
                    }
                }
            )
        }
    }

    # If we get here, the resource was not found
    stop(
        "Resource '", resource, "' not found in package or ExperimentHub.\n",
        "Available resources: ", paste(available_resources, collapse = ", ")
    )
}

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
        ret[rownames(ret) %in% rownames(locs), , drop = FALSE]
    } else if (resource == "array_type") {
        ret
    }
}
