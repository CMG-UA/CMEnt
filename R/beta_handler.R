is_file <- function(file_path) {
    is.character(file_path) && length(file_path) == 1 && file.exists(file_path)
}
file_is_tabix <- function(file_path) {
    endsWith(file_path, ".gz") && file.exists(paste0(file_path, ".tbi"))
}

#' Beta Handler Class
#'
#' @description An R6 class to handle methylation beta value files efficiently,
#' with support for in-memory loading, tabix indexing, and various file formats.
#'
#'
#' @importFrom R6 R6Class
#' @keywords internal
BetaHandler <- R6::R6Class("BetaHandler", # nolint
    public = list(
        #' @field beta Path to beta values file, or a tabix file, or in-memory beta matrix
        beta = NULL,
        #' @field genome Reference genome
        genome = "hg19",
        #' @field array Array platform, ignore for mouse genomes or when sorted_locs provided
        array = "450K",
        #' @field beta_row_names_file Path to row names file
        beta_row_names_file = NULL,
        #' @field sorted_locs Sorted genomic locations
        sorted_locs = NULL,
        #' @field njobs Number of parallel jobs
        njobs = 1,
        #' @description Create a new BetaHandler object
        #' @param beta Path to beta values file, or a tabix, or a beta matrix
        #' @param array Array platform type. Ignored if sorted_locs, or a tabix file have been provided.
        #' @param genome Reference genome version, eg. hg19. Only human and mouse genomes are supported. Ignored if sorted_locs, or a tabix file have been provided.
        #' @param beta_row_names_file Path to row names file. If NULL, row names will be read from input `beta`.
        #' @param sorted_locs Sorted genomic locations data frame. If NULL, will be retrieved automatically
        #' @param njobs Number of parallel jobs
        #' @return A new BetaHandler object
        initialize = function(beta = NULL,
                              array = c("450K", "27K", "EPIC", "EPICv2"),
                              genome = "hg19",
                              beta_row_names_file = NULL,
                              sorted_locs = NULL,
                              njobs = 1) {
            # Validate inputs
            if (is.null(beta)) {
                stop("Beta values must be provided")
            }
            if (!is.null(beta) && is.character(beta) && length(beta) == 1 && !file.exists(beta)) {
                stop("Provided beta file does not exist: ", beta)
            }
            if (is.null(sorted_locs) && !(is_file(beta) && file_is_tabix(beta))) {
                array <- strex::match_arg(array, ignore_case = TRUE)
                self$array <- array
                self$genome <- genome
            }


            # Set fields
            self$beta <- beta
            self$njobs <- njobs

            self$beta_row_names_file <- beta_row_names_file
            if (!is.null(sorted_locs)) {
                if (is.character(sorted_locs) && length(sorted_locs) == 1 && file.exists(sorted_locs)) {
                    sorted_locs <- readRDS(sorted_locs)
                    try(sorted_locs <- bigmemory::attach.big.matrix(sorted_locs), silent = TRUE)
                }
            } else if (is_file(beta) && file_is_tabix(beta)) {
                .log_info("Loading genomic locations from tabix beta file...", level = 2)
                sorted_locs <- bigmemory::attach(readRDS(genomicLocsFromTabixToDescriptor(
                    input_tabix = beta
                )))
            }
            self$sorted_locs <- sorted_locs
            private$.sorted_locs_is_bigmatrix <- bigmemory::is.big.matrix(sorted_locs)

            # Initialize private fields
            private$.beta_col_names <- NULL
            private$.beta_row_names <- NULL
            private$.loaded <- FALSE
            private$.validated <- FALSE

            # Validate and load

            invisible(self$validate())
        },

        #' @description Load beta file data into memory or prepare for file-based access
        #' @return Self (invisibly)
        load = function() {
            if (private$.loaded) {
                return(invisible(self))
            }
            if (!is.character(self$beta) && length(self$beta) > 0) {
                private$.beta_file_in_memory <- self$beta
                self$beta <- NULL
                private$.loaded <- TRUE
                return(invisible(self))
            }

            # Check if beta file is tabix
            if (is_file(self$beta)) {
                if (file_is_tabix(self$beta)) {
                    .log_step("Beta file appears to be tabix-indexed. Using as tabix file...", level = 1)
                    private$.tabix_file <- self$beta
                } else {
                    private$.beta_file <- self$beta
                }
            }
            if (!is.null(private$.beta_file)) {
                # Check file size
                file_size_mb <- file.info(private$.beta_file)$size / (1024^2)

                mem_thres <- getOption("DMRsegal.beta_in_mem_threshold_mb", 500)
                if (file_size_mb < mem_thres) {
                    .log_step(
                        "Beta file is small (", round(file_size_mb, 1),
                        " MB). Loading into memory for faster access..."
                    )

                    private$.beta_file_in_memory <- tryCatch(
                        {
                            temp_data <- data.table::fread(
                                private$.beta_file,
                                header = TRUE,
                                data.table = FALSE,
                                showProgress = getOption("DMRsegal.verbose", 0) > 1
                            )

                            # Set rownames from first column
                            rownames(temp_data) <- temp_data[[1]]
                            temp_data <- temp_data[, -1, drop = FALSE]

                            .log_success(
                                "Beta file loaded into memory (", nrow(temp_data),
                                " CpGs x ", ncol(temp_data), " samples)"
                            )

                            temp_data
                        },
                        error = function(e) {
                            .log_warn(
                                "Failed to load beta file into memory: ",
                                e$message, ". Will use alternative method."
                            )
                            NULL
                        }
                    )
                    private$.beta_file <- NULL
                } else {
                    .log_step(
                        "Beta file is large (", round(file_size_mb, 1),
                        " MB). Checking for tabix conversion opportunity..."
                    )
                }

                # If file wasn't loaded into memory, try tabix conversion
                if (is.null(private$.beta_file_in_memory)) {
                    sorted_locs <- self$getGenomicLocs()

                    # Attempt conversion
                    converted_tabix <- convertBetaToTabix(
                        beta_file = private$.beta_file,
                        sorted_locs = sorted_locs,
                        output_file = NULL
                    )

                    if (!is.null(converted_tabix)) {
                        .log_success("Beta file converted to tabix format for improved performance")
                        private$.tabix_file <- converted_tabix
                        private$.beta_file <- NULL
                    } else {
                        .log_info("Continuing with standard beta file (tabix conversion not available)")
                    }
                } else {
                    private$.beta_file <- NULL
                }
            }
            private$.loaded <- TRUE
            invisible(self)
        },

        #' @description Get **all** provided sorted genomic locations
        #' @return data.frame or matrix of the genomic locations
        getGenomicLocs = function() {
            if (is.null(self$sorted_locs)) {
                self$sorted_locs <- getSortedGenomicLocs(
                    array = self$array,
                    genome = self$genome
                )
            }
            self$sorted_locs
        },

        #' @description Get row names (CpG IDs) from the beta file
        #' @return Character vector of row names
        getBetaRowNames = function() {
            self$load()

            if (!is.null(private$.beta_row_names)) {
                return(private$.beta_row_names)
            }

            if (!is.null(self$beta_row_names_file) && file.exists(self$beta_row_names_file)) {
                .log_step("Reading row names from beta row names file: ",
                    self$beta_row_names_file, " ...",
                    level = 2
                )
                private$.beta_row_names <- unlist(read.table(
                    self$beta_row_names_file,
                    header = FALSE,
                    comment.char = "",
                    quote = ""
                ))
            } else {
                .log_step("Reading row names from input...", level = 2)
                if (!is.null(private$.beta_file_in_memory)) {
                    .log_info("Reading from beta matrix...", level = 3)
                    private$.beta_row_names <- rownames(private$.beta_file_in_memory)
                } else if (!is.null(private$.beta_file)) {
                    .log_info("Reading from beta file...", level = 3)
                    private$.beta_row_names <- unlist(data.table::fread(
                        file = private$.beta_file,
                        select = 1,
                        sep = "\t",
                        header = TRUE,
                        showProgress = TRUE,
                        nThread = self$njobs
                    ))
                } else {
                    .log_info("Reading from tabix file...", level = 3)
                    private$.beta_row_names <- unlist(data.table::fread(
                        file = private$.tabix_file,
                        select = 4,
                        sep = "\t",
                        header = TRUE,
                        showProgress = TRUE,
                        nThread = self$njobs
                    ))
                }
                if (!private$.sorted_locs_is_bigmatrix) {
                    sorted_locs <- self$getGenomicLocs()
                    private$.beta_row_names <- private$.beta_row_names[private$.beta_row_names %in% rownames(sorted_locs)]
                } else {
                    private$.beta_row_names <- as.integer(private$.beta_row_names)
                }
                if (!is.null(self$beta_row_names_file)) {
                    writeLines(
                        paste(private$.beta_row_names, collapse = "\n"),
                        self$beta_row_names_file
                    )
                }
            }
            .log_success("Row names read: ", length(private$.beta_row_names), level = 2)
            private$.beta_row_names
        },

        #' @description Get column names (sample IDs) from the beta file
        #' @return Character vector of column names
        getBetaColNames = function() {
            self$load()

            if (!is.null(private$.beta_col_names)) {
                return(private$.beta_col_names)
            }

            if (!is.null(private$.beta_file_in_memory)) {
                private$.beta_col_names <- colnames(private$.beta_file_in_memory)
            } else {
                if (!is.null(private$.beta_file)) {
                    if (endsWith(private$.beta_file, ".gz")) {
                        conn <- gzfile(private$.beta_file, "r")
                    } else {
                        conn <- file(private$.beta_file, "r")
                    }
                } else {
                    conn <- gzfile(private$.tabix_file, "r")
                }

                file_beta_col_names <- scan(conn,
                    sep = "\t",
                    what = character(),
                    nlines = 1,
                    quiet = TRUE
                )
                close(conn)

                if (!is.null(private$.beta_file)) {
                    # expected first column: cpg ID
                    file_beta_col_names <- file_beta_col_names[-1]
                } else {
                    # expected first columns: #chr start end id score strand
                    file_beta_col_names <- file_beta_col_names[7:length(file_beta_col_names)]
                }

                private$.beta_col_names <- file_beta_col_names
            }

            private$.beta_col_names
        },

        #' @description Validate that the beta file is properly sorted and formatted
        #' @return Self (invisibly)
        validate = function() {
            if (private$.validated) {
                return(invisible(self))
            }

            self$load()

            # Skip validation if no genomic locations available yet
            # (will be validated when actually used)
            if (is.null(self$sorted_locs)) {
                private$.validated <- TRUE
                return(invisible(self))
            }

            sorted_locs <- self$getGenomicLocs()
            beta_row_names <- self$getBetaRowNames()
            if (private$.sorted_locs_is_bigmatrix) {
                if (is.null(private$.tabix_file)) {
                    stop("When using big.matrix for sorted_locs, beta must be provided as tabix file.")
                }
            }
            if (is.null(private$.beta_file_in_memory)) {
                if (!private$.sorted_locs_is_bigmatrix) {
                    .log_step("Validating beta file sorting by position...", level = 2)

                    # Validate that file is sorted
                    if (!all(beta_row_names[
                        orderByLoc(beta_row_names,
                            genome = self$genome,
                            genomic_locs = sorted_locs
                        )
                    ] == beta_row_names)) {
                        stop("Provided beta file is not sorted by position!")
                    }
                    .log_success("Beta file sorting validated", level = 2)
                }
            } else {
                # Sort in-memory beta data
                private$.beta_file_in_memory <- private$.beta_file_in_memory[
                    orderByLoc(rownames(private$.beta_file_in_memory),
                        genome = self$genome,
                        genomic_locs = sorted_locs
                    ), ,
                    drop = FALSE
                ]
                rownames(private$.beta_file_in_memory) <- beta_row_names[
                    orderByLoc(beta_row_names,
                        genome = self$genome,
                        genomic_locs = sorted_locs
                    )
                ]
            }

            private$.validated <- TRUE
            invisible(self)
        },

        #' @description Get genomic locations for beta values
        #' @return Data frame of genomic locations
        getBetaLocs = function() {
            if (!is.null(private$.beta_locs)) {
                return(private$.beta_locs)
            }
            sorted_locs <- self$getGenomicLocs()
            if (private$.sorted_locs_is_bigmatrix) {
                private$.beta_locs <- sorted_locs
                return(sorted_locs)
            }
            beta_row_names <- self$getBetaRowNames()
            private$.beta_locs <- sorted_locs[beta_row_names, , drop = FALSE]
            private$.beta_locs
        },

        #' @description Extract beta values for specific CpG sites and samples
        #' @param row_names Character vector of CpG IDs to extract. If numeric, treated as row indices.
        #' @param col_names Character vector of sample IDs to extract (default: NULL for all)
        #' @param allow_missing Logical. If TRUE, missing CpG sites will be ignored instead of throwing an error (default: FALSE)
        #' @param check_mem Logical. If TRUE, checks memory usage and may return a big.matrix if size exceeds threshold (default: TRUE)
        #' @return Matrix of beta values, or big.matrix if estimated size exceeds mem_thres
        getBeta = function(row_names = NULL, col_names = NULL, allow_missing = FALSE, check_mem = FALSE) {
            self$validate()
            if (!is.null(private$.beta_file_in_memory)) {
                .log_step("Subsetting from in-memory beta data..", level = 3)
                if (!is.null(row_names)) {
                    if (is.numeric(row_names)) {
                        rcmp <- seq_len(nrow(private$.beta_file_in_memory))
                    } else {
                        rcmp <- rownames(private$.beta_file_in_memory)
                    }
                    if (allow_missing) {
                        row_names <- intersect(row_names, rcmp)
                    } else {
                        missing_rows <- setdiff(row_names, rcmp)
                        if (length(missing_rows) > 0) {
                            stop(
                                "Requested CpG sites not found in beta data: ",
                                paste(missing_rows, collapse = ", ")
                            )
                        }
                    }
                    if (is.null(col_names)) {
                        beta_subset <- private$.beta_file_in_memory[row_names, , drop = FALSE]
                    } else {
                        beta_subset <- private$.beta_file_in_memory[row_names, col_names, drop = FALSE]
                    }
                } else {
                    if (is.null(col_names)) {
                        beta_subset <- private$.beta_file_in_memory
                    } else {
                        beta_subset <- private$.beta_file_in_memory[, col_names, drop = FALSE]
                    }
                }
                .log_success("Beta values subsetted: ", nrow(beta_subset),
                    " CpGs x ", ncol(beta_subset), " samples",
                    level = 3
                )
                return(beta_subset)
            }
            if (check_mem){
            first_row <- if (!is.null(private$.beta_file)) {
                data.table::fread(private$.beta_file, nrows = 1, data.table = FALSE)
            } else {
                if (is.null(row_names)) {
                    locs <- self$getBetaLocs()
                } else {
                    locs <- self$getBetaLocs()[row_names, , drop = FALSE]
                }
                region_first <- data.frame(
                    chr = as.character(locs[1, "chr"]),
                    start = as.integer(locs[1, "start"]),
                    end = as.integer(locs[1, "end"])
                )
                bedr::tabix(region_first, private$.tabix_file,
                    check.valid = FALSE,
                    check.sort = FALSE, check.chr = FALSE, verbose = FALSE
                )
            }

            n_rows <- if (is.null(row_names)) length(self$getBetaRowNames()) else length(row_names)
            estimated_size_mb <- as.numeric(object.size(first_row)) * n_rows / (1024^2)
            mem_thres <- getOption("DMRsegal.subset_beta_as_bigmem_mb", 500)
            if (estimated_size_mb > mem_thres) {
                .log_info("Estimated size (", round(estimated_size_mb, 1),
                    " MB) exceeds threshold. Loading to big.matrix...",
                    level = 3
                )
                chunk_size <- max(2, ceiling(mem_thres * 1024^2 / as.numeric(object.size(first_row))))
                return(private$.loadToBigMatrix(row_names, col_names, allow_missing, chunk_size))
            }
            }

            if (!is.null(private$.beta_file)) {
                .log_step("Subsetting from beta file..", level = 3)
                beta_subset <- .subsetBetaFile(
                    private$.beta_file,
                    row_names,
                    beta_row_names = private$.beta_row_names,
                    beta_col_names = private$.beta_col_names
                )
                if (!is.null(row_names) && !allow_missing) {
                    missing_rows <- setdiff(row_names, rownames(beta_subset))
                    if (length(missing_rows) > 0) {
                        stop(
                            "Requested CpG sites not found in beta file: ",
                            paste(missing_rows, collapse = ", ")
                        )
                    }
                }
                if (!is.null(col_names)) {
                    beta_subset <- beta_subset[, col_names, drop = FALSE]
                }
            } else {
                .log_step("Subsetting from tabix file..", level = 3)
                if (is.null(row_names)) {
                    locs <- self$getBetaLocs()
                } else {
                    locs <- self$getBetaLocs()[row_names, , drop = FALSE]
                }
                regions <- data.frame(
                    chr = as.character(locs[, "chr"]),
                    start = as.integer(locs[, "start"]),
                    end = as.integer(locs[, "end"])
                )
                beta_subset <- bedr::tabix(regions, private$.tabix_file,
                    check.valid = FALSE,
                    check.sort = FALSE, check.chr = FALSE, verbose = FALSE
                )
                if (is.null(beta_subset) || (!allow_missing && !is.null(row_names) && nrow(beta_subset) < length(row_names))) {
                    stop("Requested CpG sites not found in beta tabix file")
                }
                beta_subset <- as.data.frame(sapply(beta_subset[, 7:ncol(beta_subset), drop = FALSE], as.numeric))
                if (!is.null(col_names)) {
                    beta_subset <- beta_subset[, col_names, drop = FALSE]
                }
            }

            .log_success("Beta values subsetted: ", nrow(beta_subset),
                " CpGs x ", ncol(beta_subset), " samples",
                level = 3
            )
            return(beta_subset)
        }
    ),
    private = list(
        .beta_col_names = NULL,
        .beta_row_names = NULL,
        .beta_locs = NULL,
        .loaded = FALSE,
        .validated = FALSE,
        .beta_file = NULL,
        .tabix_file = NULL,
        .beta_file_in_memory = NULL,
        .sorted_locs_is_bigmatrix = FALSE,
        .loadToBigMatrix = function(row_names = NULL, col_names = NULL, allow_missing, chunk_size) {
            beta_row_names <- self$getBetaRowNames()
            if (is.null(row_names)) {
                row_names <- beta_row_names
            }
            beta_col_names <- self$getBetaColNames()
            if (is.null(col_names)) {
                col_names <- beta_col_names
            }
            temp_file <- tempfile(fileext = ".bmat")
            beta_big <- bigmemory::big.matrix(
                nrow = length(row_names), ncol = length(col_names), type = "double",
                backingfile = basename(temp_file), backingpath = dirname(temp_file),
                descriptorfile = paste0(basename(temp_file), ".desc")
            )
            options(bigmemory.allow.dimnames=TRUE)
            if (!is.null(private$.beta_file)) {

                chunks <- split(seq_along(row_names), ceiling(seq_along(row_names) / chunk_size))
                for (chunk_idx in seq_along(chunks)) {
                    chunk_rows <- row_names[chunks[[chunk_idx]]]
                    chunk_data <- .subsetBetaFile(
                        private$.beta_file,
                        sites = chunk_rows,
                        beta_row_names = beta_row_names,
                        beta_col_names = beta_col_names
                    )
                    if (!is.null(col_names)) {
                        chunk_data <- chunk_data[, col_names, drop = FALSE]
                    }
                    beta_big[chunks[[chunk_idx]], ] <- as.matrix(chunk_data)
                }
            } else {
                if (is.null(row_names)) {
                    locs <- self$getBetaLocs()
                } else {
                    locs <- self$getBetaLocs()[row_names, , drop = FALSE]
                }
                chunks <- split(seq_len(nrow(locs)), ceiling(seq_len(nrow(locs)) / chunk_size))
                for (i in seq_along(chunks)) {
                    chunk_locs <- locs[chunks[[i]], , drop = FALSE]
                    regions_chunk <- data.frame(
                        chr = as.character(chunk_locs[, "chr"]),
                        start = as.integer(chunk_locs[, "start"]),
                        end = as.integer(chunk_locs[, "end"])
                    )
                    chunk <- bedr::tabix(regions_chunk, private$.tabix_file,
                        check.valid = FALSE,
                        check.sort = FALSE, check.chr = FALSE, verbose = FALSE
                    )
                    if (!is.null(chunk) && nrow(chunk) > 0) {
                        chunk <- as.data.frame(sapply(chunk[, 7:ncol(chunk), drop = FALSE], as.numeric))
                        if (!is.null(col_names)) {
                            chunk <- chunk[, col_names, drop = FALSE]
                        }
                        beta_big[chunks[[i]], ] <- as.matrix(chunk)
                    }
                }
            }
            rownames(beta_big) <- row_names
            colnames(beta_big) <- col_names
            .log_success("Loaded to big.matrix: ", nrow(beta_big), " CpGs x ", ncol(beta_big), " samples", level = 3)
            beta_big
        }
    )
)

#' Create a BetaHandler object for efficient beta value access
#'
#' @name getBetaHandler
#' @description Create a new BetaHandler object that manages methylation beta values
#' from various input formats (files, matrices, tabix) with memory-efficient access patterns.
#'
#' @param beta Path to beta values file, or a tabix, or a beta matrix
#' @param array Array platform type, **ignored** if sorted_locs have been provided
#' @param genome Reference genome version, eg. hg19. Only human and mouse genomes are supported. **Ignored** if sorted_locs have been provided.
#' @param beta_row_names_file Path to row names file
#' @param sorted_locs Data frame with genomic locations containing 'chr' and 'start' and 'end' columns, sorted by genomic position. If NULL, will be retrieved automatically using genome and array information.
#' @param njobs Number of parallel jobs
#' @return A new BetaHandler object
#'
#' @examples
#' beta_matrix <- beta <- loadExampleInputData("beta")
#'
#' beta_handler <- getBetaHandler(
#'     beta = beta_matrix,
#'     array = "450K",
#'     genome = "hg19"
#' )
#'
#' beta_locs <- beta_handler$getBetaLocs()
#' head(beta_locs)
#'
#' beta_values <- beta_handler$getBeta()
#' head(beta_values[, 1:5])
#'
#' @export
getBetaHandler <- function(beta, array = c("450K", "27K", "EPIC", "EPICv2"),
                           genome = c("hg19", "hg38", "mm10", "mm39"),
                           beta_row_names_file = NULL,
                           sorted_locs = NULL,
                           njobs = 1) {
    if (inherits(beta, "BetaHandler")) {
        .log_warn("Provided beta is already a BetaHandler instance. Returning it directly.")
        return(invisible(beta))
    }
    invisible(BetaHandler$new(
        beta = beta,
        array = array,
        genome = genome,
        beta_row_names_file = beta_row_names_file,
        njobs = njobs,
        sorted_locs = sorted_locs
    ))
}
