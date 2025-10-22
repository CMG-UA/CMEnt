#' Beta File Handler Class
#'
#' @description An R6 class to handle methylation beta value files efficiently,
#' with support for in-memory loading, tabix indexing, and various file formats.
#'
#'
#' @importFrom R6 R6Class
#' @export
BetaFileHandler <- R6::R6Class("BetaFileHandler",
    public = list(
        #' @field beta_file Path to beta values file
        beta_file = NULL,
        #' @field tabix_file Path to tabix-indexed file
        tabix_file = NULL,
        #' @field genome Reference genome
        genome = "hg19",
        #' @field array Array platform
        array = "450K",
        #' @field beta_file_in_memory In-memory beta data
        beta_file_in_memory = NULL,
        #' @field beta_row_names_file Path to row names file
        beta_row_names_file = NULL,
        #' @field sorted_locs Sorted genomic locations
        sorted_locs = NULL,
        #' @field verbose Verbosity level
        verbose = 0,
        #' @field memory_threshold_mb Memory threshold in MB
        memory_threshold_mb = 500,
        #' @field njobs Number of parallel jobs
        njobs = 1,

        #' @description Create a new BetaFileHandler object
        #' @param beta_file Path to beta values file
        #' @param tabix_file Path to tabix-indexed file
        #' @param array Array platform type
        #' @param genome Reference genome version
        #' @param beta_row_names_file Path to row names file
        #' @param verbose Verbosity level
        #' @param memory_threshold_mb Memory threshold in MB
        #' @param njobs Number of parallel jobs
        #' @return A new BetaFileHandler object
        initialize = function(beta_file = NULL,
                              tabix_file = NULL,
                              array = c("450K", "27K", "EPIC", "EPICv2"),
                              genome = c("hg19", "hg38", "mm10", "mm39"),
                              beta_row_names_file = NULL,
                              verbose = 0,
                              memory_threshold_mb = 500,
                              njobs = 1) {
            # Validate inputs
            if (is.null(beta_file) && is.null(tabix_file)) {
                stop("Either beta_file or tabix_file must be provided")
            }
            if (!is.null(beta_file) && !file.exists(beta_file)) {
                stop("Provided beta_file does not exist: ", beta_file)
            } else if (!is.null(tabix_file) && !file.exists(tabix_file)) {
                stop("Provided tabix_file does not exist: ", tabix_file)
            }

            array <- match.arg(array)
            genome <- match.arg(genome)

            # Set fields
            self$beta_file <- beta_file
            self$tabix_file <- tabix_file
            self$array <- array
            self$genome <- genome
            self$beta_row_names_file <- beta_row_names_file
            self$verbose <- verbose
            self$memory_threshold_mb <- memory_threshold_mb
            self$njobs <- njobs

            # Initialize private fields
            private$.beta_col_names <- NULL
            private$.beta_row_names <- NULL
            private$.loaded <- FALSE
            private$.validated <- FALSE

            # Validate and load
            self$validate()

            invisible(self)
        },

        #' @description Load beta file data into memory or prepare for file-based access
        #' @return Self (invisibly)
        load = function() {
            if (private$.loaded) {
                return(invisible(self))
            }

            if (!is.null(self$beta_file) && is.null(self$tabix_file)) {
                # Check file size
                file_size_mb <- file.info(self$beta_file)$size / (1024^2)

                if (file_size_mb < self$memory_threshold_mb) {
                    .log_step(
                        "Beta file is small (", round(file_size_mb, 1),
                        " MB). Loading into memory for faster access..."
                    )

                    self$beta_file_in_memory <- tryCatch(
                        {
                            temp_data <- data.table::fread(
                                self$beta_file,
                                header = TRUE,
                                data.table = FALSE,
                                showProgress = self$verbose > 1
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
                } else {
                    .log_step(
                        "Beta file is large (", round(file_size_mb, 1),
                        " MB). Checking for tabix conversion opportunity..."
                    )
                }

                # If file wasn't loaded into memory, try tabix conversion
                if (is.null(self$beta_file_in_memory)) {
                    sorted_locs <- private$get_sorted_locs()

                    # Attempt conversion
                    converted_tabix <- convertBetaToTabix(
                        beta_file = self$beta_file,
                        sorted_locs = sorted_locs,
                        output_file = NULL,
                        verbose = self$verbose > 0
                    )

                    if (!is.null(converted_tabix)) {
                        .log_success("Beta file converted to tabix format for improved performance")
                        self$tabix_file <- converted_tabix
                        self$beta_file <- NULL
                    } else {
                        .log_info("Continuing with standard beta file (tabix conversion not available)")
                    }
                }
            }

            private$.loaded <- TRUE
            invisible(self)
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
                .log_step("Reading row names from beta file...", level = 2)
                if (!is.null(self$beta_file)) {
                    private$.beta_row_names <- unlist(data.table::fread(
                        file = self$beta_file,
                        select = 1,
                        sep = "\t",
                        header = TRUE,
                        showProgress = TRUE,
                        nThread = self$njobs
                    ))
                } else {
                    private$.beta_row_names <- unlist(data.table::fread(
                        file = self$tabix_file,
                        select = 4,
                        sep = "\t",
                        header = TRUE,
                        showProgress = TRUE,
                        nThread = self$njobs
                    ))
                }

                if (!is.null(self$beta_row_names_file)) {
                    writeLines(
                        paste(private$.beta_row_names, collapse = "\n"),
                        self$beta_row_names_file
                    )
                }
            }

            .log_success("Row names read: ", length(private$.beta_row_names), level = 2)
            return(private$.beta_row_names)
        },

        #' @description Get column names (sample IDs) from the beta file
        #' @return Character vector of column names
        getBetaColNames = function() {
            self$load()

            if (!is.null(private$.beta_col_names)) {
                return(private$.beta_col_names)
            }

            if (!is.null(self$beta_file_in_memory)) {
                private$.beta_col_names <- colnames(self$beta_file_in_memory)
            } else {
                if (!is.null(self$beta_file)) {
                    if (endsWith(self$beta_file, ".gz")) {
                        conn <- gzfile(self$beta_file, "r")
                    } else {
                        conn <- file(self$beta_file, "r")
                    }
                } else {
                    conn <- gzfile(self$tabix_file, "r")
                }

                file_beta_col_names <- scan(conn,
                    sep = "\t",
                    what = character(),
                    nlines = 1,
                    quiet = TRUE
                )
                close(conn)

                if (!is.null(self$beta_file)) {
                    # expected first column: cpg ID
                    file_beta_col_names <- file_beta_col_names[-1]
                } else {
                    # expected first columns: #chr start end id score strand
                    file_beta_col_names <- file_beta_col_names[7:length(file_beta_col_names)]
                }

                private$.beta_col_names <- file_beta_col_names
            }

            return(private$.beta_col_names)
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

            .log_step("Validating beta file sorting by position...", level = 2)

            sorted_locs <- private$get_sorted_locs()
            beta_row_names <- self$getBetaRowNames()

            if (is.null(self$beta_file_in_memory)) {
                # Validate that file is sorted
                if (!all(beta_row_names[orderByLoc(beta_row_names,
                    genome = self$genome,
                    genomic_locs = sorted_locs
                )] == beta_row_names)) {
                    stop("Provided beta file is not sorted by position!")
                }
            } else {
                # Sort in-memory beta data
                self$beta_file_in_memory <- self$beta_file_in_memory[
                    orderByLoc(rownames(self$beta_file_in_memory),
                        genome = self$genome,
                        genomic_locs = sorted_locs
                    ), ,
                    drop = FALSE
                ]
                rownames(self$beta_file_in_memory) <- beta_row_names[
                    orderByLoc(beta_row_names,
                        genome = self$genome,
                        genomic_locs = sorted_locs
                    )
                ]
            }

            .log_success("Beta file sorting validated", level = 2)
            private$.validated <- TRUE
            invisible(self)
        },

        #' @description Get genomic locations for beta values
        #' @return Data frame of genomic locations
        getBetaLocs = function() {
            if (!is.null(private$.beta_locs)) {
                return(private$.beta_locs)
            }
            sorted_locs <- private$get_sorted_locs()
            beta_row_names <- self$getBetaRowNames()
            private$.beta_locs <- sorted_locs[beta_row_names, , drop = FALSE]
            return(private$.beta_locs)
        },

        #' @description Extract beta values for specific CpG sites and samples
        #' @param row_names Character vector of CpG IDs to extract
        #' @param col_names Character vector of sample IDs to extract (default: NULL for all)
        #' @return Matrix of beta values
        getBeta = function(row_names, col_names = NULL) {
            self$validate()

            # Use in-memory beta data if available
            if (!is.null(self$beta_file_in_memory)) {
                .log_step("Subsetting from in-memory beta data..", level = 3)
                if (is.null(col_names)) {
                    beta_subset <- self$beta_file_in_memory[row_names, , drop = FALSE]
                } else {
                    beta_subset <- self$beta_file_in_memory[row_names, col_names, drop = FALSE]
                }
            } else if (!is.null(self$beta_file)) {
                .log_step("Subsetting from beta file..", level = 3)
                beta_subset <- .subsetBetaFile(
                    self$beta_file,
                    row_names,
                    beta_row_names = private$.beta_row_names,
                    beta_col_names = private$.beta_col_names
                )
                if (!is.null(col_names)) {
                    beta_subset <- beta_subset[, col_names, drop = FALSE]
                }
            } else {
                .log_step("Subsetting from tabix file..", level = 3)
                locs <- private$get_sorted_locs()[row_names, , drop = FALSE]
                regions <- locs[, c("chr", "start", "end")]
                regions[, "chr"] <- as.character(regions[, "chr"])
                beta_subset <- bedr::tabix(
                    regions,
                    self$tabix_file,
                    check.valid = FALSE,
                    verbose = FALSE
                )
                if (is.null(beta_subset)) {
                    dir.create("debug", showWarnings = FALSE)
                    saveRDS(row_names, "debug/row_names.rds")
                    saveRDS(private$.beta_row_names, "debug/beta_row_names.rds")
                    stop("No beta values found for the requested CpG sites in the tabix file")
                }
                # Remove BED columns (chr, start, end, name, score, strand) to get only beta values
                # Columns 1-6 are BED format, 7+ are beta values
                beta_subset <- beta_subset[, 7:ncol(beta_subset), drop = FALSE]
                if (!is.null(col_names)) {
                    beta_subset <- as.data.frame(sapply(beta_subset[, col_names, drop = FALSE], as.numeric))
                } else {
                    beta_subset <- as.data.frame(sapply(beta_subset, as.numeric))
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
        get_sorted_locs = function() {
            if (is.null(self$sorted_locs)) {
                self$sorted_locs <- getSortedGenomicLocs(
                    array = self$array,
                    genome = self$genome
                )
            }
            return(self$sorted_locs)
        }
    )
)
