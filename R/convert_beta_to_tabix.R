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
#' if (nzchar(Sys.which("tabix")) && nzchar(Sys.which("bgzip"))) {
#'     beta_file <- tempfile(fileext = ".tsv")
#'     writeLines(c("\tsample1", "cg1\t0.5"), beta_file)
#'     locs <- data.frame(chr = "chr1", start = 100L, row.names = "cg1")
#'     tabix_file <- convertBetaToTabix(beta_file, sorted_locs = locs)
#' }
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
    if (!all(nzchar(Sys.which(c("tabix", "bgzip"))))) {
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
            owns_temp_bed <- is.null(.bed_file)
            if (owns_temp_bed) {
                array <- strex::match_arg(array, ignore_case = TRUE)
                # Get sorted locations if not provided
                if (is.null(sorted_locs)) {
                    sorted_locs <- getSortedGenomicLocs(array = array, genome = genome)
                }
                # Read header to get column names
                .log_step("Reading beta file header...", level = 2)

                header_info <- .getBetaFileHeaderInfo(beta_file)
                col_names <- header_info$file_beta_col_names

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
                bed_header <- c("#chrom", "start", "end", "id", "score", "strand", header_info$sample_col_names)
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
                        beta_subset <- chunk_data[match(common_sites, site_ids), header_info$sample_col_names, drop = FALSE]
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
                system2("sh", c("-c", shQuote(sort_cmd)))
            }
            if (owns_temp_bed) {
                unlink(temp_bed)
            }
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
