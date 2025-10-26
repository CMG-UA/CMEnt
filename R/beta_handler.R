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
        #' @field array Array platform, ignore for mouse genomes
        array = "450K",
        #' @field beta_row_names_file Path to row names file
        beta_row_names_file = NULL,
        #' @field sorted_locs Sorted genomic locations
        sorted_locs = NULL,
        #' @field memory_threshold_mb Memory threshold in MB
        memory_threshold_mb = 500,
        #' @field njobs Number of parallel jobs
        njobs = 1,
        #' @description Create a new BetaHandler object
        #' @param beta Path to beta values file, or a tabix, or a beta matrix
        #' @param array Array platform type
        #' @param genome Reference genome version
        #' @param beta_row_names_file Path to row names file
        #' @param memory_threshold_mb Memory threshold in MB
        #' @param njobs Number of parallel jobs
        #' @return A new BetaHandler object
        initialize = function(beta = NULL,
                              array = c("450K", "27K", "EPIC", "EPICv2"),
                              genome = c("hg19", "hg38", "mm10", "mm39"),
                              beta_row_names_file = NULL,
                              memory_threshold_mb = 500,
                              njobs = 1) {
            # Validate inputs
            if (is.null(beta)) {
                stop("Beta values must be provided")
            }
            if (!is.null(beta) && is.character(beta) && length(beta) == 1 && !file.exists(beta)) {
                stop("Provided beta file does not exist: ", beta)
            }
            array <- strex::match_arg(array, ignore_case = TRUE)
            genome <- strex::match_arg(genome, ignore_case = TRUE)

            # Set fields
            self$beta <- beta
            self$array <- array
            self$genome <- genome
            self$beta_row_names_file <- beta_row_names_file
            self$verbose <- getOption("DMRsegal.verbose", default = 0)
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
            if (!is.character(self$beta) && length(self$beta) > 0) {
                if (inherits(self$beta, "BetaHandler")) {
                    warning("Provided beta is already a BetaHandler instance. Returning it directly.")
                    return(invisible(self$beta))
                }
                private$.beta_file_in_memory <- self$beta
                self$beta <- NULL
                private$.loaded <- TRUE
                return(invisible(self))
            }
            # Check if beta file is tabix
            if (!is.null(self$beta) && is.character(self$beta) && length(self$beta) == 1) {
                if (endsWith(self$beta, ".gz") && file.exists(paste0(self$beta, ".tbi"))) {
                    .log_step("Beta file appears to be tabix-indexed. Using as tabix file...", level = 1)
                    private$.tabix_file <- self$beta
                } else {
                    private$.beta_file <- self$beta
                }
            }
            if (!is.null(private$.beta_file)) {
                # Check file size
                file_size_mb <- file.info(private$.beta_file)$size / (1024^2)

                if (file_size_mb < self$memory_threshold_mb) {
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
                if (is.null(private$.beta_file_in_memory)) {
                    sorted_locs <- private$get_sorted_locs()

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
                if (!is.null(private$.beta_file_in_memory)) {
                    private$.beta_row_names <- rownames(private$.beta_file_in_memory)
                } else if (!is.null(private$.beta_file)) {
                    private$.beta_row_names <- unlist(data.table::fread(
                        file = private$.beta_file,
                        select = 1,
                        sep = "\t",
                        header = TRUE,
                        showProgress = TRUE,
                        nThread = self$njobs
                    ))
                } else {
                    private$.beta_row_names <- unlist(data.table::fread(
                        file = private$.tabix_file,
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

            .log_step("Validating beta file sorting by position...", level = 2)

            sorted_locs <- private$get_sorted_locs()
            beta_row_names <- self$getBetaRowNames()

            if (is.null(private$.beta_file_in_memory)) {
                # Validate that file is sorted
                if (!all(beta_row_names[
                    orderByLoc(beta_row_names,
                        genome = self$genome,
                        genomic_locs = sorted_locs
                    )
                ] == beta_row_names)) {
                    stop("Provided beta file is not sorted by position!")
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
            private$.beta_locs
        },

        #' @description Extract beta values for specific CpG sites and samples
        #' @param row_names Character vector of CpG IDs to extract
        #' @param col_names Character vector of sample IDs to extract (default: NULL for all)
        #' @return Matrix of beta values
        getBeta = function(row_names = NULL, col_names = NULL) {
            self$validate()
            if (is.null(row_names)) {
                row_names <- self$getBetaRowNames()
            }
            # Use in-memory beta data if available
            if (!is.null(private$.beta_file_in_memory)) {
                .log_step("Subsetting from in-memory beta data..", level = 3)
                if (is.null(col_names)) {
                    beta_subset <- private$.beta_file_in_memory[row_names, , drop = FALSE]
                } else {
                    beta_subset <- private$.beta_file_in_memory[row_names, col_names, drop = FALSE]
                }
            } else if (!is.null(private$.beta_file)) {
                .log_step("Subsetting from beta file..", level = 3)
                beta_subset <- .subsetBetaFile(
                    private$.beta_file,
                    row_names,
                    beta_row_names = private$.beta_row_names,
                    beta_col_names = private$.beta_col_names
                )
                if (!is.null(col_names)) {
                    beta_subset <- beta_subset[, col_names, drop = FALSE]
                }
            } else {
                .log_step("Subsetting from tabix file..", level = 3)
                locs <- self$getBetaLocs()[row_names, , drop = FALSE]
                regions <- locs[, c("chr", "start", "end")]
                regions[, "chr"] <- as.character(regions[, "chr"])
                beta_subset <- bedr::tabix(
                    regions,
                    private$.tabix_file,
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
        .beta_file = NULL,
        .tabix_file = NULL,
        .beta_file_in_memory = NULL,
        get_sorted_locs = function() {
            if (is.null(self$sorted_locs)) {
                self$sorted_locs <- getSortedGenomicLocs(
                    array = self$array,
                    genome = self$genome
                )
            }
            self$sorted_locs
        }
    )
)

#' @description Create a new BetaHandler object
#' @param beta Path to beta values file, or a tabix, or a beta matrix
#' @param array Array platform type
#' @param genome Reference genome version
#' @param beta_row_names_file Path to row names file
#' @param memory_threshold_mb Memory threshold in MB
#' @param njobs Number of parallel jobs
#' @return A new BetaHandler object
#'
#' @export
getBetaHandler <- function(beta, array = c("450K", "27K", "EPIC", "EPICv2"),
                           genome = c("hg19", "hg38", "mm10", "mm39"),
                           beta_row_names_file = NULL,
                           memory_threshold_mb = 500,
                           njobs = 1) {
    BetaHandler$new(
        beta = beta,
        array = array,
        genome = genome,
        beta_row_names_file = beta_row_names_file,
        memory_threshold_mb = memory_threshold_mb,
        njobs = njobs
    )
}
