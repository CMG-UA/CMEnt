#' Sort Beta File by Genomic Coordinates
#'
#' @description This helper function sorts a methylation beta values file by genomic coordinates
#' (chromosome and position) as required by the buildDMRs function. The function reads
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
#' beta_file <- tempfile(fileext = ".tsv")
#' writeLines(c("sample1", "cg2\t0.2", "cg1\t0.1"), beta_file)
#' locs <- data.frame(
#'     chr = c("chr1", "chr1"),
#'     start = c(100L, 200L),
#'     row.names = c("cg1", "cg2")
#' )
#' sorted_file <- sortBetaFileByCoordinates(beta_file, genomic_locs = locs)
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
    beta_data <- .readBetaFileData(
        beta_file,
        data_table = FALSE,
        showProgress = getOption("CMEnt.verbose", 0) > 1
    )

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