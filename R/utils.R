#' Get File Hash for Caching
#'
#' @description Internal helper to compute MD5 hash of a file for cache validation
#'
#' @param file_path Path to the file
#' @return MD5 hash string
#' @keywords internal
.getFileHash <- function(file_path) {
    tryCatch({
        tools::md5sum(file_path)
    }, error = function(e) {
        # Fallback: use file size and modification time
        info <- file.info(file_path)
        paste0(info$size, "_", as.numeric(info$mtime))
    })
}

#' Convert Beta File to Tabix-Indexed Format
#'
#' @description Converts a methylation beta values file to a tabix-indexed BED format
#' for faster random access during DMR analysis. The function uses a memory-efficient
#' chunk-based approach to handle large files and automatically uses a cache directory
#' to avoid redundant conversions of the same beta file.
#'
#' @param beta_file Character. Path to the input beta values file
#' @param sorted_locs Data frame with genomic locations containing 'chr' and 'pos' columns.
#'   If NULL, will be retrieved automatically using getSortedGenomicLocs() (default: NULL)
#' @param array Character. Array platform type, either "450K" or "EPIC". Only used if
#'   sorted_locs is NULL (default: "450K")
#' @param output_file Character. Path for the output tabix file. If NULL, uses a cache
#'   directory in tempdir() with hash-based naming (default: NULL)
#' @param chunk_size Integer. Number of rows to process in each chunk (default: 50000)
#' @param verbose Logical. Whether to print progress messages (default: TRUE)

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
#' \code{tempdir()/DMRSegal_tabix_cache/} and persists for the duration of the R session.
#' Files are named based on the MD5 hash of the input beta file, ensuring that identical
#' files reuse the same cached version.
#'
#' @examples
#' \dontrun{
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
#' }
#'
#' @export
convertBetaToTabix <- function(beta_file, sorted_locs = NULL, array = c("450K", "EPIC"), 
                               output_file = NULL,
                               chunk_size = 50000,
                               verbose = TRUE) {
    # Get sorted locations if not provided
    if (is.null(sorted_locs)) {
        array <- match.arg(array)
        sorted_locs <- getSortedGenomicLocs(array = array)
    }
    # Check if tabix/bgzip are available
    tabix_available <- tryCatch({
        system2("which", "tabix", stdout = FALSE, stderr = FALSE)
        system2("which", "bgzip", stdout = FALSE, stderr = FALSE)
        TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    
    if (!tabix_available) {
        if (verbose) {
            .log_warn("tabix/bgzip not found in PATH. Skipping tabix conversion.")
        }
        return(NULL)
    }
    
    # Set default output file name - use cache directory
    if (is.null(output_file)) {
        # Create cache directory in temp folder
        cache_dir <- file.path(tempdir(), "DMRSegal_tabix_cache")
        dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
        
        # Compute hash of beta file
        beta_hash <- .getFileHash(beta_file)
        
        # Create cache filename based on hash
        output_file <- file.path(cache_dir, paste0("beta_", beta_hash, ".bed.gz"))
        
        # Check if tabix file already exists in cache
        if (file.exists(output_file) && file.exists(paste0(output_file, ".tbi"))) {
            if (verbose) {
                .log_info("Using cached tabix file: ", basename(output_file))
            }
            return(output_file)
        }
    }
    
    if (verbose) {
        .log_step("Converting beta file to tabix format...")
    }
    
    tryCatch({
        # Read header to get column names
        if (verbose) {
            .log_step("Reading beta file header...", level = 2)
        }
        
        header_conn <- if (endsWith(beta_file, ".gz")) gzfile(beta_file, "r") else file(beta_file, "r")
        header_line <- readLines(header_conn, n = 1)
        close(header_conn)
        col_names <- strsplit(header_line, "\t")[[1]]
        
        # Get total number of rows for progress tracking
        if (verbose) {
            .log_step("Counting rows in beta file...", level = 2)
        }
        
        # Count lines efficiently
        if (endsWith(beta_file, ".gz")) {
            n_lines <- as.numeric(system(sprintf("zcat %s | wc -l", shQuote(beta_file)), intern = TRUE))
        } else {
            n_lines <- as.numeric(system(sprintf("wc -l < %s", shQuote(beta_file)), intern = TRUE))
        }
        n_rows <- n_lines - 1  # Exclude header
        
        if (verbose) {
            .log_info("Processing ", n_rows, " CpG sites...", level = 2)
        }
        
        # Create temporary BED file for writing chunks
        temp_bed <- tempfile(fileext = ".bed")
        
        # Write header to temp BED file
        bed_header <- c("chr", "start", "end", "name", col_names[-1])
        writeLines(paste(bed_header, collapse = "\t"), temp_bed)
        
        # Process file in chunks to avoid memory issues
        
        skip_rows <- 1  # Start after header
        rows_processed <- 0
        
        while (rows_processed < n_rows) {
            if (verbose && getOption("DMRSegal.verbose", 1) > 1) {
                .log_info("Processing rows ", rows_processed + 1, " to ", 
                         min(rows_processed + chunk_size, n_rows), "...", level = 3)
            }
            
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
                # Create BED format for this chunk
                bed_chunk <- sorted_locs[common_cpgs, c("chr", "pos"), drop = FALSE]
                bed_chunk$start <- bed_chunk$pos
                bed_chunk$end <- bed_chunk$pos + 1
                bed_chunk$name <- rownames(bed_chunk)
                bed_chunk <- bed_chunk[, c("chr", "start", "end", "name")]
                
                # Add beta values
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
        
        if (verbose) {
            .log_success("Processed ", rows_processed, " rows", level = 2)
        }
        
        # Check if any data was written
        if (file.info(temp_bed)$size <= length(paste(bed_header, collapse = "\t")) + 1) {
            if (verbose) {
                .log_warn("No common CpGs found between beta file and genomic locations")
            }
            unlink(temp_bed)
            return(NULL)
        }
        
        # Sort, compress with bgzip, and index with tabix
        if (verbose) {
            .log_step("Sorting BED file...", level = 2)
        }
        temp_sorted <- tempfile(fileext = ".bed")
        sort_cmd <- sprintf("(head -n 1 %s && tail -n +2 %s | sort -k1,1 -k2,2n) > %s", 
                           shQuote(temp_bed), shQuote(temp_bed), shQuote(temp_sorted))
        system(sort_cmd)
        
        # Compress with bgzip
        if (verbose) {
            .log_step("Compressing with bgzip...", level = 2)
        }
        bgzip_cmd <- sprintf("bgzip -c %s > %s", shQuote(temp_sorted), shQuote(output_file))
        system(bgzip_cmd)
        
        # Index with tabix
        if (verbose) {
            .log_step("Creating tabix index...", level = 2)
        }
        tabix_cmd <- sprintf("tabix -p bed %s", shQuote(output_file))
        system(tabix_cmd)
        
        # Clean up temp files
        unlink(temp_bed)
        unlink(temp_sorted)
        
        if (file.exists(output_file) && file.exists(paste0(output_file, ".tbi"))) {
            if (verbose) {
                .log_success("Tabix file created: ", output_file)
            }
            return(output_file)
        } else {
            if (verbose) {
                .log_warn("Failed to create tabix index")
            }
            return(NULL)
        }
    }, error = function(e) {
        if (verbose) {
            .log_warn("Error converting to tabix: ", e$message)
        }
        return(NULL)
    })
}

.readSamplesheet <- function(samplesheet_file,
                             samplesheet_file_sep,
                             sample_group_col,
                             sample_group_case,
                             sample_group_control,
                             target_col,
                             subset = NULL,
                             max_samples = NULL,
                             max_case_samples = NULL,
                             max_control_samples = NULL) {
    require(stringr)
    require(data.table)

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
        if (!is.null(max_control_samples) && max_control_samples < nrow(control_samplesheet)) {
            inds <- sort(sample(seq_len(nrow(control_samplesheet)), max_control_samples, replace = FALSE))
            control_samplesheet <- control_samplesheet[inds, ]
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
        if (!is.null(max_case_samples) && max_case_samples < nrow(case_samplesheet)) {
            inds <- sort(sample(seq_len(nrow(case_samplesheet)), max_case_samples, replace = FALSE))
            case_samplesheet <- case_samplesheet[inds, ]
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

    if (!is.null(max_samples) && max_samples < nrow(subset_samplesheet)) {
        inds <- sort(sample(seq_len(nrow(subset_samplesheet)), max_samples, replace = FALSE))
        subset_samplesheet <- subset_samplesheet[inds, ]
    }

    rownames(subset_samplesheet) <- subset_samplesheet[, target_col]
    subset_samplesheet <- subset_samplesheet[, sample_group_col,
        drop =
            FALSE
    ]
    subset_samplesheet <- subset_samplesheet[
        str_order(rownames(subset_samplesheet),
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
.dmrsegal_log_env <- local({
    e <- new.env(parent = emptyenv())
    e$last_step_time <- list()
    e
})

.fmt_dur <- function(start_time) {
    if (is.null(start_time)) {
        return("")
    }
    secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (is.na(secs)) {
        return("")
    }
    if (secs < 60) sprintf(" (took %.2fs)", secs) else sprintf(" (took %dm %02ds)", floor(secs / 60), round(secs %% 60))
}

.has_ansi <- function() {
    isFALSE(nzchar(Sys.getenv("NO_COLOR"))) && cli::num_ansi_colors() > 1L
}

.col <- function(x, col = c("cyan", "green", "yellow", "blue")) {
    col <- match.arg(col)
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
.log_info <- function(..., .envir = parent.frame(), level = 1) {
    if (getOption("DMRSegal.verbose", 1) < level) {
        return(invisible())
    }
    msg <- paste0(..., collapse = "")
    lead <- paste(rep("\t", level - 1), .col(cli::symbol$info, "blue"), sep = "")
    cat(paste(lead, msg, '\n'), file = stderr())
    invisible()
}

#' @keywords internal
.log_warn <- function(..., .envir = parent.frame()) {
    msg <- paste0(..., collapse = "")
    lead <- .col(cli::symbol$warning, "yellow")
    cat(paste(lead, msg, '\n'), file = stderr())
    invisible()
}

#' @keywords internal
.log_success <- function(..., .envir = parent.frame(), level = 1) {
    if (getOption("DMRSegal.verbose", 1) < level) {
        return(invisible())
    }
    dur <- .fmt_dur(.dmrsegal_log_env$last_step_time[[level]])
    msg <- paste0(paste0(..., collapse = ""), dur)
    lead <- paste(rep("\t", level - 1), .col(cli::symbol$tick, "green"), sep = "")
    cat(paste(lead, msg, '\n'), file = stderr())
    invisible()
}

#' @keywords internal
.log_step <- function(title, ..., .envir = parent.frame(), level = 1) {
    if (getOption("DMRSegal.verbose", 1) < level) {
        return(invisible())
    }
    .dmrsegal_log_env$last_step_time[level:max(1, length(.dmrsegal_log_env$last_step_time))] <- Sys.time()
    subtitle <- paste0(..., collapse = "")
    header <- title
    if (nzchar(subtitle)) header <- paste0(header, ": ", subtitle)
    lead <- paste(rep("\t", level - 1), .col(cli::symbol$arrow_right, "cyan"), sep = "")
    cat(paste(lead, header, '\n'), file = stderr())
    invisible()
}

.processSamplesheet <- function(args,
                                subset = NULL) {
    suppressWarnings(suppressMessages({
        require(stringr)
    }))
    samplesheet_file <- args$samplesheet
    target_col <- args$target_col
    sample_group_col <- args$sample_group_col
    sample_group_control <- args$sample_group_control
    sample_group_case <- args$sample_group_case
    max_missing_ratio <- args$max_missing_cov_ratio
    if (is.null(max_missing_ratio)) {
        max_missing_ratio <- 0.3
    }
    if (!is.null(sample_group_case)) {
        sample_group_case <- strsplit(sample_group_case, ",")[[1]]
    }
    if (!is.null(sample_group_control)) {
        sample_group_control <- strsplit(sample_group_control, ",")[[1]]
    }
    max_samples <- args$max_samples
    max_case_samples <- args$max_case_samples
    max_control_samples <- args$max_control_samples

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
        subset = subset,
        max_samples = max_samples,
        max_case_samples = max_case_samples,
        max_control_samples = max_control_samples
    )
    ret <- list(samplesheet = subset_samplesheet[, c(sample_group_col, "casecontrol")])
    ret
}

#' Get Sorted Genomic Locations for Array Platform
#'
#' @description Retrieves and sorts genomic location annotations for the specified
#' methylation array platform (450K or EPIC). The locations are sorted by
#' chromosome and position to ensure proper genomic ordering.
#'
#' @param array Character. Array platform type, either "450K" or "EPIC"
#'
#' @return A data frame containing sorted genomic locations with columns:
#' \itemize{
#'   \item chr: Chromosome
#'   \item pos: Genomic position
#'   \item start: Start position (same as pos)
#'   \item end: End position (pos + 1)
#' }
#'
#' @examples
#' \dontrun{
#' # Get sorted locations for 450K array
#' locs_450k <- getSortedGenomicLocs("450K")
#'
#' # Get sorted locations for EPIC array
#' locs_epic <- getSortedGenomicLocs("EPIC")
#' }
#'
#' @export
getSortedGenomicLocs <- function(array) {
    if (array == "450K") {
        sorted_locs <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
    } else if (array == "EPIC") {
        sorted_locs <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
    } else {
        stop("Unknown array type: ", array)
    }
    sorted_locs <- sorted_locs[str_order(paste0(sorted_locs[, "chr"], ":", sorted_locs[, "pos"]), numeric = TRUE), ]
    sorted_locs[, "start"] <- sorted_locs[, "pos"]
    sorted_locs[, "end"] <- sorted_locs[, "pos"] + 1
    sorted_locs
}

#' Order Indices by Genomic Location
#'
#' @description Orders a vector of indices according to their corresponding genomic
#' locations (chromosome and position). This function is useful for sorting CpG
#' sites or other genomic features by their physical positions.
#'
#' @param x Character or integer vector. Indices or identifiers to be ordered
#' @param array Character. Array platform type, either "450K" or "EPIC" (default: "450K")
#' @param genomic_locs Data frame. Optional pre-computed genomic locations. If NULL,
#' locations will be retrieved using getSortedGenomicLocs (default: NULL)
#'
#' @return Integer vector of ordered indices
#'
#' @examples
#' \dontrun{
#' # Order CpG indices by genomic location
#' cpg_ids <- c("cg00000029", "cg00000108", "cg00000109")
#' ordered_indices <- orderByLoc(cpg_ids, array = "450K")
#'
#' # Order using pre-computed genomic locations
#' locs <- getSortedGenomicLocs("EPIC")
#' ordered_indices <- orderByLoc(cpg_ids, genomic_locs = locs)
#' }
#'
#' @export
orderByLoc <- function(x, array = "450K", genomic_locs = NULL) {
    if (is.null(genomic_locs)) {
        genomic_locs <- getSortedGenomicLocs(array)
    }
    str_order(paste0(genomic_locs[x, "chr"], ":", genomic_locs[x, "pos"]), numeric = TRUE)
}
