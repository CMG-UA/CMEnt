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
#' @param njobs Integer. Number of jobs to use after normalization. On Unix-like
#'   systems, location-registry creation and tabix conversion run in parallel
#'   when `njobs > 1`.
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
#'     sample_group = c("case", "control"),
#'     row.names = c("Sample1", "Sample2")
#' )
#'
#' if (nzchar(Sys.which("tabix")) && nzchar(Sys.which("bgzip"))) {
#'     bed_file <- tempfile(fileext = ".bed")
#'     writeLines(c(
#'         "#chrom\tstart\tSample1\tSample2",
#'         "chr1\t100\t0.2\t0.8",
#'         "chr1\t200\t0.3\t0.7"
#'     ), bed_file)
#'     result <- readCustomMethylationBedData(bed_file, pheno)
#'     result$tabix_file
#' }
#'
#' @seealso
#' \code{\link{convertBetaToTabix}} for converting standard beta files to tabix format
#' \code{\link{getBetaHandler}} for creating a BetaHandler object from processed files
#'
#' @export
readCustomMethylationBedData <- function(bed_file, pheno, genome = "hg38", chrom_col = "#chrom",
                                         start_col = "start", output_dir = NULL, chunk_size = 50000,
                                         output_prefix = NULL, njobs = 1L) {
    njobs <- suppressWarnings(as.integer(njobs)[1L])
    if (is.na(njobs) || njobs < 1L) {
        stop("njobs must be a positive integer.")
    }
    tabix_available <- all(nzchar(Sys.which(c("tabix", "bgzip"))))

    cache_dir <- .getTabixCacheDir(output_dir)
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    hash <- .getFileHash(bed_file)
    normalized_bed_file <- file.path(tempdir(), paste0("bed_", hash, ".tsv"))
    on.exit(unlink(normalized_bed_file), add = TRUE)
    fallback_beta_file <- NULL
    if (!tabix_available) {
        .log_warn("tabix/bgzip not found in PATH. Using normalized beta TSV fallback for custom BED input.")
        fallback_beta_file <- .getDerivedOutputPath(output_prefix, ".input_beta.tsv")
        if (is.null(fallback_beta_file)) {
            fallback_beta_file <- file.path(cache_dir, paste0("bed_beta_", hash, ".tsv"))
        }
    }

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
    num_rows <- sum(vapply(readLines(tmp_con), function(x) nchar(x) > 0, logical(1))) - 1
    close(tmp_con)
    .log_info("Processing BED file with ", num_rows, " rows and ", length(existing_ids), " matching sample IDs.", level = 2)


    # Read chunks of the BED file to minimize memory usage
    con <- if (endsWith(bed_file, ".gz")) gzfile(bed_file, "r") else file(bed_file, "r")
    # Write normalized header to new BED file
    norm_bed_header <- c("#chrom", "start", "end", "id", "score", "strand", existing_ids)
    writeLines(paste(norm_bed_header, collapse = "\t"), normalized_bed_file)
    if (!is.null(fallback_beta_file)) {
        writeLines(paste(c("site_id", existing_ids), collapse = "\t"), fallback_beta_file)
    }
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
        if (!is.null(fallback_beta_file)) {
            beta_subset <- data.frame(
                site_id = paste0(bed_data$chr, ":", bed_data$start),
                bed_data[, existing_ids, drop = FALSE],
                check.names = FALSE
            )
            data.table::fwrite(
                beta_subset,
                file = fallback_beta_file,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                append = TRUE
            )
        }
        count <- count + nrow(bed_data)
    }
    close(con)

    locations_h5_file <- .getDerivedOutputPath(output_prefix, ".input_beta.locations.h5")
    hash <- .getFileHash(bed_file)
    tabix_file_path <- .getDerivedOutputPath(output_prefix, ".input_beta.tabix.bed.gz")
    if (is.null(tabix_file_path)) {
        tabix_file_path <- file.path(cache_dir, paste0("bed_beta_", hash, ".bed.gz"))
    }

    build_locations <- function() {
        getRegistry(
            normalized_bed_file,
            select = c("#chrom", "start"),
            rename = c("#chrom" = "chr", start = "start"),
            derive = list(
                index = list(
                    cols = c("chr", "start"),
                    fun = function(chr, start) paste0(chr, ":", start)
                )
            ),
            indices = "index",
            chunk_size = chunk_size,
            output_h5file = locations_h5_file
        )
    }
    build_tabix <- function(njobs) {
        tabix_file <- convertBetaToTabix(
            .bed_file = normalized_bed_file,
            output_file = tabix_file_path,
            chunk_size = chunk_size,
            njobs = njobs,
            output_prefix = output_prefix
        )
        if (is.null(tabix_file)) {
            stop("Failed to create tabix-indexed BED file.")
        }
        tabix_file
    }

    if (tabix_available && njobs > 1L && identical(.Platform$OS.type, "unix")) {
        locations_job <- parallel::mcparallel(build_locations(), silent = TRUE)
        tabix_result <- tryCatch(
            build_tabix(max(1L, njobs - 1L)),
            error = identity
        )
        locations_result <- parallel::mccollect(locations_job, wait = TRUE)[[1L]]
        if (inherits(tabix_result, "error")) {
            stop(tabix_result)
        }
        if (inherits(locations_result, "try-error")) {
            stop(attr(locations_result, "condition"))
        }
        tabix_file_path <- tabix_result
        locations <- locations_result
    } else {
        locations <- build_locations()
        tabix_file_path <- if (tabix_available) build_tabix(njobs) else NULL
    }


    list(
        tabix_file = tabix_file_path,
        beta_file = if (is.null(tabix_file_path)) fallback_beta_file else tabix_file_path,
        locations = locations
    )
}
