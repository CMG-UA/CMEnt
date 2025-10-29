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
        #' @param array Array platform type. Ignored if sorted_locs have been provided.
        #' @param genome Reference genome version, eg. hg19. Only human and mouse genomes are supported. Ignored if sorted_locs have been provided.
        #' @param beta_row_names_file Path to row names file
        #' @param sorted_locs Sorted genomic locations data frame. If NULL, will be retrieved automatically
        #' @param memory_threshold_mb Memory threshold in MB
        #' @param njobs Number of parallel jobs
        #' @return A new BetaHandler object
    initialize = function(beta = NULL,
                  array = c("450K", "27K", "EPIC", "EPICv2"),
                  genome = c("hg19", "hg38", "mm10", "mm39"),
                  beta_row_names_file = NULL,
                  sorted_locs = NULL,
                  memory_threshold_mb = 500,
                  njobs = 1) {
            # Validate inputs
            if (is.null(beta)) {
                stop("Beta values must be provided")
            }
            if (!is.null(beta) && is.character(beta) && length(beta) == 1 && !file.exists(beta)) {
                stop("Provided beta file does not exist: ", beta)
            }
            if(is.null(sorted_locs)){
                array <- strex::match_arg(array, ignore_case = TRUE)
                genome <- strex::match_arg(genome, ignore_case = TRUE)
                self$array <- array
                self$genome <- genome
            }


            # Set fields
            self$beta <- beta

            self$beta_row_names_file <- beta_row_names_file
            if (!is.null(sorted_locs)){
                if (is.character(sorted_locs) && length(sorted_locs) == 1 && file.exists(sorted_locs)){
                    sorted_locs <- readRDS(sorted_locs)
                    try(sorted_locs <- bigmemory::attach.big.matrix(sorted_locs), silent = TRUE)
                }
            }
            private$.sorted_locs_is_bigmatrix <- bigmemory::is.big.matrix(sorted_locs)
            self$sorted_locs <- sorted_locs
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
                if(is.null(private$.tabix_file)){
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
                regions <- as.data.frame(regions)
                regions[, "chr"] <- as.character(regions[, "chr"])
                regions[, "start"] <- as.integer(regions[, "start"])
                regions[, "end"] <- as.integer(regions[, "end"])
                beta_subset <- bedr::tabix(
                    regions,
                    private$.tabix_file,
                    check.valid = FALSE,
                    check.sort = FALSE,
                    check.chr = FALSE,
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
                    beta_subset <- as.data.frame(as.numeric(beta_subset[, col_names, drop = FALSE]))
                } else {
                    beta_subset <- as.data.frame(as.numeric(beta_subset))
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
        .sorted_locs_is_bigmatrix = FALSE
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
#' @param memory_threshold_mb Memory threshold in MB
#' @param njobs Number of parallel jobs
#' @return A new BetaHandler object
#'
#' @examples
#' \donttest{
#' if (!requireNamespace("DMRsegaldata", quietly = TRUE)) {
#'     remotes::install_github("CMG-UA/DMRsegaldata")
#' }
#' library(DMRsegaldata)
#' beta_matrix <- DMRsegaldata::beta
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
#' }
#'
#' @export
getBetaHandler <- function(beta, array = c("450K", "27K", "EPIC", "EPICv2"),
                           genome = c("hg19", "hg38", "mm10", "mm39"),
                           beta_row_names_file = NULL,
                           sorted_locs = NULL,
                           memory_threshold_mb = 500,
                           njobs = 1) {
    BetaHandler$new(
        beta = beta,
        array = array,
        genome = genome,
        beta_row_names_file = beta_row_names_file,
        memory_threshold_mb = memory_threshold_mb,
        njobs = njobs,
        sorted_locs = sorted_locs
    )
}
