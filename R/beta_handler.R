is_file <- function(file_path) {
    is.character(file_path) && length(file_path) == 1 && file.exists(file_path)
}
file_is_tabix <- function(file_path) {
    endsWith(file_path, ".gz") && file.exists(paste0(file_path, ".tbi"))
}
is_bsseq <- function(obj) {
    inherits(obj, "BSseq")
}

#' Beta Handler Class
#'
#' @description An R6 class to handle methylation beta value files efficiently,
#' with support for in-memory loading, tabix indexing, BSseq objects, and various file formats.
#'
#'
#' @importFrom R6 R6Class
#' @importFrom GenomicRanges granges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges start end
#' @importFrom bsseq sampleNames
#' @keywords internal
BetaHandler <- R6::R6Class("BetaHandler", # nolint
    public = list(
        #' @field beta Path to beta values file, or a tabix file, or in-memory beta matrix, or BSseq object
        beta = NULL,
        #' @field genome Reference genome
        genome = "hg38",
        #' @field array Array platform, ignore for mouse genomes or when sorted_locs provided
        array = "450K",
        #' @field beta_row_names_file Path to row names file
        beta_row_names_file = NULL,
        #' @field sorted_locs Sorted genomic locations
        sorted_locs = NULL,
        #' @field njobs Number of parallel jobs
        njobs = 1,
        #' @field beta_chunk_size Chunk size for subsetting beta values
        beta_chunk_size = NULL,
        #' @description Create a new BetaHandler object
        #' @param beta Path to beta values file, or a tabix, or a beta matrix, or a BSseq object
        #' @param array Array platform type. Ignored if sorted_locs, a tabix file, or a BSseq object have been provided.
        #' @param genome Reference genome version, eg. hg38 or hs1. Only human and mouse genomes are supported. Ignored if sorted_locs, a tabix file, or a BSseq object have been provided.
        #' @param beta_row_names_file Path to row names file. If NULL, row names will be read from input `beta`.
        #' @param sorted_locs Sorted genomic locations data frame. If given, the input data will be assumed already sorted. If NULL, will be retrieved automatically
        #' @param chrom_col Chromosome column name in tabix file
        #' @param start_col Start position column name in tabix file
        #' @param njobs Number of parallel jobs
        #' @return A new BetaHandler object
        initialize = function(beta = NULL,
                              array = c("450K", "27K", "EPIC", "EPICv2"),
                              genome = "hg38",
                              beta_row_names_file = NULL,
                              sorted_locs = NULL,
                              chrom_col = "#chrom",
                              start_col = "start",
                              njobs = 1) {
            # Validate inputs
            if (is.null(beta)) {
                stop("Beta values must be provided")
            }
            if (!is.null(beta) && is.character(beta) && length(beta) == 1 && !file.exists(beta)) {
                stop("Provided beta file does not exist: ", beta)
            }
            if (is.null(sorted_locs) && !(is_file(beta) && file_is_tabix(beta)) && !is_bsseq(beta)) {
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
                }
                private$.self_contained <- TRUE
            } else if (is_file(beta) && file_is_tabix(beta)) {
                .log_info("Loading genomic locations from tabix beta file...", level = 2)
                sorted_locs <- genomicLocsFromTabix(input_tabix = beta, chrom_col = chrom_col, start_col = start_col)
                private$.self_contained <- TRUE
            } else if (is_bsseq(beta)) {
                .log_step("Extracting genomic locations from BSseq object...", level = 2)
                gr <- granges(beta)
                .log_step("Constructing sorted_locs delayed data frame..", level = 3)
                sorted_locs <- DelayedDataFrame::DelayedDataFrame(
                    chr = as.character(seqnames(gr)),
                    start = start(gr),
                    end = end(gr)
                )
                .log_success("Constructed sorted_locs delayed data frame..", level = 3)
                .log_step("Creating 'name' column for sorted_locs..", level = 3)
                rownames(sorted_locs) <- paste0(sorted_locs$chr, ":", sorted_locs$start)
                .log_success("Created 'name' column for sorted_locs..", level = 3)

                .log_success("Genomic locations extracted from BSseq object: ", nrow(sorted_locs), " CpGs", level = 2)
                private$.self_contained <- TRUE
            }
            self$sorted_locs <- sorted_locs

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
            .log_info("Loading beta data...", level = 2)
            if (!is.character(self$beta) && length(self$beta) > 0) {
                if (is_bsseq(self$beta)) {
                    private$.bsseq_object <- sort(self$beta)
                    self$beta <- NULL
                    private$.loaded <- TRUE
                    return(invisible(self))
                }
                .log_info("Beta provided as in-memory object. Using as is...", level = 2)
                private$.beta_file_in_memory <- self$beta
                self$beta <- NULL
                private$.loaded <- TRUE
                .log_info("Beta data loading complete.", level = 2)
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
                if (file_size_mb < mem_thres ) {
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
                    .log_info(
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
                        self$sorted_locs <- genomicLocsFromTabix(input_tabix = private$.tabix_file, use_id_as_rownames = TRUE)
                        private$.beta_row_names <- rownames(self$sorted_locs)
                        private$.self_contained <- TRUE
                        private$.beta_file <- NULL
                    } else {
                        .log_info("Continuing with standard beta file (tabix conversion not available)")
                    }
                } else {
                    private$.beta_file <- NULL
                }
            }
            .log_info("Beta data loading complete.", level = 2)
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
                .log_info("Reading row names from beta row names file: ",
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
                .log_info("Reading row names from input...", level = 2)
                if (!is.null(private$.bsseq_object)) {
                    .log_info("Reading from BSseq object...", level = 3)
                    gr <- granges(private$.bsseq_object)
                    private$.beta_row_names <- paste0(seqnames(gr), ":", start(gr))
                } else if (!is.null(private$.beta_file_in_memory)) {
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
                    names(private$.beta_row_names) <- NULL
                } else {
                    .log_info("Reading from tabix file...", level = 3)
                    private$.beta_row_names <- do.call(
                        paste,
                        c(data.table::fread(
                            file = private$.tabix_file,
                            select = c(1, 2),
                            sep = "\t",
                            header = TRUE,
                            showProgress = TRUE,
                            nThread = self$njobs
                        ), sep = ":"
                        )
                    )
                    names(private$.beta_row_names) <- NULL
                }
                if (!private$.self_contained) {
                    sorted_locs <- self$getGenomicLocs()
                    private$.beta_row_names <- private$.beta_row_names[private$.beta_row_names %in% rownames(sorted_locs)]
                }
                if (!is.null(self$beta_row_names_file)) {
                    writeLines(
                        paste(private$.beta_row_names, collapse = "\n"),
                        self$beta_row_names_file
                    )
                }
            }
            if (length(private$.beta_row_names) == 0) {
                stop("No row names could be read from the beta input!")
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

            if (!is.null(private$.bsseq_object)) {
                private$.beta_col_names <- sampleNames(private$.bsseq_object)
            } else if (!is.null(private$.beta_file_in_memory)) {
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
                    # expected first columns: #chrom start end id score strand
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
            .log_info("Loading genomic locations for validation...", level = 2)
            sorted_locs <- self$getGenomicLocs()
            .log_info("Getting beta row names for validation...", level = 2)
            beta_row_names <- self$getBetaRowNames()
            if (is.null(private$.beta_file_in_memory)) {
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
            } else {
                # Sort in-memory beta data
                private$.beta_file_in_memory <- private$.beta_file_in_memory[
                    orderByLoc(rownames(private$.beta_file_in_memory),
                        genome = self$genome,
                        genomic_locs = sorted_locs
                    ), ,
                    drop = FALSE
                ]
                private$.beta_row_names <- rownames(private$.beta_file_in_memory)
            }
            .log_info("Beta file validated.", level = 2)
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
            if (private$.self_contained) {
                private$.beta_locs <- sorted_locs
                return(sorted_locs)
            }
            beta_row_names <- self$getBetaRowNames()
            private$.beta_locs <- sorted_locs[beta_row_names, , drop = FALSE]
            private$.beta_locs
        },
        #' @description Get the chunk size used for beta subsetting
        #' @return Integer chunk size
        getBetaChunkSize = function() {
            self$load()
            self$beta_chunk_size
        },

        #' @description Check if the beta data is array-based (i.e. does not have row names in 'chr:pos' format)
        #' @return Logical indicating if the beta data is array-based
        isArrayBased = function() {
            if (!is.null(private$.is_array_based)) {
                return(private$.is_array_based)
            }
            self$load()
            first_loc <- self$getBetaLocs()[1, , drop = FALSE]
            private$.is_array_based <- !grepl(rownames(first_loc), pattern = "^chr[A-Za-z0-9]+:\\d+$", ignore.case = TRUE)
            private$.is_array_based
        },
        #' @description Extract beta values for specific CpG sites and samples
        #' @param row_names Character vector of CpG IDs to extract. If numeric, treated as row indices.
        #' @param col_names Character vector of sample IDs to extract (default: NULL for all)
        #' @param allow_missing Logical. If TRUE, missing CpG sites will be ignored instead of throwing an error (default: FALSE)
        #' @param chr Character vector of chromosome names to extract, cannot be used along with row_names (default: NULL for all)
        #' @return Matrix of beta values
        getBeta = function(row_names = NULL, col_names = NULL, allow_missing = FALSE, chr = NULL) {
            if (!is.null(row_names) && !is.null(chr)) {
                stop("Cannot specify both row_names and chr for subsetting.")
            }
            self$validate()
            # Fast-path numeric row indexing for non in-memory backends:
            # convert indices to row names once and avoid repeated set operations downstream.
            if (!is.null(row_names) && is.numeric(row_names) && is.null(private$.beta_file_in_memory)) {
                row_idx <- as.integer(row_names)
                all_row_names <- self$getBetaRowNames()
                n_all <- length(all_row_names)
                if (allow_missing) {
                    keep <- !is.na(row_idx) & row_idx >= 1L & row_idx <= n_all
                    row_idx <- row_idx[keep]
                } else {
                    bad <- is.na(row_idx) | row_idx < 1L | row_idx > n_all
                    if (any(bad)) {
                        bad_idx <- unique(row_idx[bad])
                        bad_idx <- bad_idx[!is.na(bad_idx)]
                        stop(
                            "Requested row indices out of bounds [1,", n_all, "]: ",
                            paste(head(bad_idx, 10), collapse = ", ")
                        )
                    }
                }
                row_names <- all_row_names[row_idx]
            }
            if (!is.null(private$.beta_file_in_memory)) {
                .log_step("Subsetting from in-memory beta data..", level = 3)
                if (!is.null(row_names)) {
                    if (is.numeric(row_names)) {
                        row_idx <- as.integer(row_names)
                        n_rows <- nrow(private$.beta_file_in_memory)
                        if (allow_missing) {
                            keep <- !is.na(row_idx) & row_idx >= 1L & row_idx <= n_rows
                            row_idx <- row_idx[keep]
                        } else {
                            bad <- is.na(row_idx) | row_idx < 1L | row_idx > n_rows
                            if (any(bad)) {
                                bad_idx <- unique(row_idx[bad])
                                bad_idx <- bad_idx[!is.na(bad_idx)]
                                stop(
                                    "Requested row indices out of bounds [1,", n_rows, "]: ",
                                    paste(head(bad_idx, 10), collapse = ", ")
                                )
                            }
                        }
                        beta_subset <- private$.beta_file_in_memory[row_idx, , drop = FALSE]
                    } else {
                        rcmp <- rownames(private$.beta_file_in_memory)
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
                        beta_subset <- private$.beta_file_in_memory[row_names, , drop = FALSE]
                    }
                } else if (!is.null(chr)) {
                    .log_info("Subsetting by chromosome from in-memory beta data..", level = 3)
                    .log_info("Getting genomic locations for chromosome subsetting...", level = 4)
                    all_locs <- self$getBetaLocs()
                    .log_info("Performing chromosome subsetting...", level = 4)
                    chr_rows <- rownames(all_locs)[all_locs$chr %in% chr]
                    .log_info("Found ", length(chr_rows), " CpGs on specified chromosome(s)", level = 4)
                    chr_rows <- intersect(chr_rows, rownames(private$.beta_file_in_memory))
                    beta_subset <- private$.beta_file_in_memory[chr_rows, , drop = FALSE]
                } else {
                    beta_subset <- private$.beta_file_in_memory
                }
                if (!is.null(col_names)) {
                    beta_subset <- beta_subset[, col_names, drop = FALSE]
                }
                .log_success("Beta values subsetted: ", nrow(beta_subset),
                    " CpGs x ", ncol(beta_subset), " samples",
                    level = 3
                )
                return(beta_subset)
            }

            if (!is.null(private$.beta_file)) {
                .log_step("Subsetting from beta file..", level = 3)
                if (!is.null(chr)) {
                    .log_step("Subsetting by chromosome from beta file..", level = 3)
                    sorted_locs <- self$getGenomicLocs()
                    row_names_chr <- rownames(sorted_locs)[sorted_locs$chr %in% chr]
                    row_names_chr <- intersect(row_names_chr, self$getBetaRowNames())
                    row_names <- row_names_chr
                }
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
            }
            if (!is.null(private$.tabix_file)) {
                .log_step("Subsetting from tabix file..", level = 3)
                if (!is.null(chr)) {
                    qregions <- chr
                    regions <- data.frame(chr = base::strsplit(chr, ",")[[1]])
                } else {
                    regions <- private$.regionsFromRowNames(row_names)
                    qregions <- regions[!duplicated(regions),]
                    qregions <- qregions[str_order(paste(qregions$chr, qregions$start, ":"), numeric = TRUE), 1:3, drop = FALSE]
                }
                beta_subset <- bedr::tabix(qregions, private$.tabix_file,
                    check.valid = FALSE,
                    check.sort = FALSE, check.chr = FALSE, verbose = FALSE
                )
                if (is.null(beta_subset) || (!allow_missing && !is.null(row_names) && nrow(beta_subset) < length(row_names))) {
                    stop("Requested CpG sites not found in beta tabix file")
                }
                # bedr forces the first three columns to be named "chr", "start", "stop" .... https://github.com/cran/bedr/blob/ddf228e25c7ff2084246060a38cfc073ab56db33/R/tabix.R#L91
                if (is.null(chr)) {
                    merge_on <- c("chr", "start")
                    beta_subset <- merge(
                        regions[, c(merge_on, "name")],
                        beta_subset, by.x = merge_on, by.y = merge_on,
                        all.x = TRUE, all.y = FALSE
                    )
                }
                rownames(beta_subset) <- beta_subset$name
                # order by input row_names if provided
                if (!is.null(row_names)) {
                    beta_subset <- beta_subset[row_names, , drop = FALSE]
                }
                beta_subset <- beta_subset[, 8:ncol(beta_subset), drop = FALSE]
                beta_subset <- as.data.frame(sapply(beta_subset, as.numeric))
                if (!is.null(col_names)) {
                    beta_subset <- beta_subset[, col_names, drop = FALSE]
                }
            }
            if (!is.null(private$.bsseq_object)) {
                .log_step("Extracting beta values from BSseq object..", level = 4)
                all_row_names <- self$getBetaRowNames()
                if (!is.null(chr)) {
                    all_locs <- self$getBetaLocs()
                    chr_rows <- rownames(all_locs)[all_locs$chr %in% chr]
                    selected_row_names <- intersect(chr_rows, all_row_names)
                    row_idx <- match(selected_row_names, all_row_names)
                } else if (!is.null(row_names)) {
                    if (is.numeric(row_names)) {
                        row_idx <- row_names
                        selected_row_names <- all_row_names[row_idx]
                    } else {
                        if (allow_missing) {
                            selected_row_names <- intersect(row_names, all_row_names)
                        } else {
                            missing_rows <- setdiff(row_names, all_row_names)
                            if (length(missing_rows) > 0) {
                                stop(
                                    "Requested CpG sites not found in BSseq object: ",
                                    paste(missing_rows, collapse = ", ")
                                )
                            }
                            selected_row_names <- row_names
                        }
                        row_idx <- match(selected_row_names, all_row_names)
                    }
                } else {
                    selected_row_names <- all_row_names
                    row_idx <- seq_along(all_row_names)
                }

                all_col_names <- self$getBetaColNames()
                if (is.null(col_names)) {
                    selected_col_names <- all_col_names
                    col_idx <- seq_along(all_col_names)
                } else {
                    missing_cols <- setdiff(col_names, all_col_names)
                    if (length(missing_cols) > 0) {
                        stop(
                            "Requested samples not found in BSseq object: ",
                            paste(missing_cols, collapse = ", ")
                        )
                    }
                    selected_col_names <- col_names
                    col_idx <- match(selected_col_names, all_col_names)
                }

                m_assay <- SummarizedExperiment::assay(private$.bsseq_object, "M", withDimnames = FALSE)
                cov_assay <- SummarizedExperiment::assay(private$.bsseq_object, "Cov", withDimnames = FALSE)
                beta_subset <- as.matrix(
                    m_assay[row_idx, col_idx, drop = FALSE] / cov_assay[row_idx, col_idx, drop = FALSE]
                )
                rownames(beta_subset) <- selected_row_names
                colnames(beta_subset) <- selected_col_names
            }

            .log_success("Beta values subsetted: ", nrow(beta_subset),
                " CpGs x ", ncol(beta_subset), " samples",
                level = 3
            )
            beta_subset
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
        .is_array_based = NULL,
        .beta_file_in_memory = NULL,
        .bsseq_object = NULL,
        .self_contained = FALSE,
        .regionsFromRowNames = function(row_names) {
            if (is.null(row_names)) {
                locs <- self$getBetaLocs()
            } else {
                locs <- self$getBetaLocs()[row_names, , drop = FALSE]
            }
            regions <- data.frame(
                chr = as.character(locs[, "chr"]),
                start = as.integer(locs[, "start"]),
                end = as.integer(locs[, "start"]) + 1,
                name = rownames(locs)
            )
        }
    )
)

#' Create a BetaHandler object for efficient beta value access
#'
#' @name getBetaHandler
#' @description Create a new BetaHandler object that manages methylation beta values
#' from various input formats (files, matrices, tabix, BSseq objects) with memory-efficient access patterns.
#'
#' @param beta Path to beta values file, or a tabix, or a beta matrix, or a BSseq object
#' @param array Array platform type, **ignored** if sorted_locs or a BSseq object have been provided
#' @param genome Reference genome version, eg. hg38 or hs1. Only human and mouse genomes are supported. **Ignored** if sorted_locs or a BSseq object have been provided.
#' @param beta_row_names_file Path to row names file
#' @param sorted_locs Data frame with genomic locations containing 'chr' and 'start' and 'end' columns, sorted by genomic position. If NULL, will be retrieved automatically using genome and array information, or extracted from BSseq object.
#' @param chrom_col Chromosome column name in tabix file
#' @param start_col Start position column name in tabix file
#' @param njobs Number of parallel jobs to use when reading beta file or tabix file. Default is number of available cores minus one, up to a maximum of 8.
#' @return A new BetaHandler object
#'
#' @examples
#' beta_matrix <- loadExampleInputData("beta")
#'
#' beta_handler <- getBetaHandler(
#'     beta = beta_matrix,
#'     array = "450K",
#'     genome = "hg38"
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
                           genome = c("hg38", "hg19", "hs1", "mm10", "mm39"),
                           beta_row_names_file = NULL,
                           sorted_locs = NULL,
                           chrom_col = "#chrom",
                           start_col = "start",
                           njobs = getOption("DMRsegal.njobs", min(8, future::availableCores() - 1))) {
    if (inherits(beta, "BetaHandler")) {
        return(invisible(beta))
    }
    invisible(BetaHandler$new(
        beta = beta,
        array = array,
        genome = genome,
        beta_row_names_file = beta_row_names_file,
        njobs = njobs,
        sorted_locs = sorted_locs,
        chrom_col = chrom_col,
        start_col = start_col
    ))
}
