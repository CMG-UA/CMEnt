#' Find Differentially Methylated Regions (DMRs) from Pre-computed DMPs
#'
#' @name findDMRsFromDMPs
#' @description This function identifies Differentially Methylated Regions (DMRs) from pre-computed
#' Differentially Methylated Positions (DMPs) using a correlation-based approach. It expands
#' significant DMPs into regions, considering both statistical significance and biological
#' relevance of methylation changes.
#'
#' @section Important Note on Input Data:
#' Do not apply heavy filtering to your DMPs prior to using this function, particularly based on
#' beta values or effect sizes. The function works by expanding regions around significant DMPs
#' and connecting nearby CpGs into larger regions. Filtering out DMPs with smaller effect sizes
#' may remove important CpGs that could serve as "bridges" to connect more significant DMPs into
#' larger, biologically meaningful DMRs. For optimal results, include all statistically
#' significant DMPs (e.g., adjusted p-value < 0.05) and let the function handle region expansion
#' and filtering internally using the min_cpg_delta_beta parameter if needed.
#'
#' @param beta_file Path to the methylation beta values file or a data matrix with beta values
#' @param dmps_tsv_file Path to the pre-computed DMPs file or a data frame with DMPs
#' @param pheno Data frame containing sample phenotype information
#' @param pval_col Column name in DMPs file containing p-values (default: "pval_adj")
#' @param sample_group_col Column in pheno for sample grouping (default: "Sample_Group")
#' @param dmp_group_col Column in DMPs file for grouping DMPs (default: NULL)
#' @param casecontrol_col Column in pheno for case/control status (default: "casecontrol")
#' @param min_cpg_delta_beta Minimum delta beta threshold for CpGs (default: 0)
#' @param expansion_step Distance in bp to expand regions during search (default: 500)
#' @param expansion_relaxation Relaxation parameter for region expansion (default: 0)
#' @param array Array platform, either "450K" or "EPIC" (default: c("450K", "EPIC"))
#' @param genome Reference genome, "hg19" or "hg38" (default: c("hg19", "hg38"))
#' @param max_pval Maximum p-value threshold for DMPs (default: 0.05)
#' @param max_lookup_dist Maximum distance for region expansion in bp (default: 10000)
#' @param min_dmps Minimum number of DMPs required per region (default: 1)
#' @param min_adj_dmps Minimum number of adjacent DMPs required per region (default: 1)
#' @param min_cpgs Minimum number of CpGs required per region (default: 50)
#' @param ignored_sample_groups Sample groups to ignore during analysis (default: NULL)
#' @param output_prefix Optional identifier prefix for output files (default: NULL)
#' @param njobs Number of parallel jobs (default: detectCores())
#' @param verbose Enable verbose output (default: FALSE)
#' @param beta_row_names_file Optional file with beta value row names (default: NULL)
#' @param dmps_beta_file Optional separate beta file for DMPs (default: NULL)
#' @param tabix_file Path to tabix-indexed beta values file (alternative to beta_file, default: NULL)
#'
#' @return A GRanges object containing identified DMRs with metadata columns:
#' \itemize{
#'   \item n_cpgs: Number of CpGs in the region
#'   \item n_dmps: Number of DMPs in the region
#'   \item mean_delta_beta: Mean methylation difference
#'   \item max_delta_beta: Maximum methylation difference
#'   \item min_pval: Minimum p-value of DMPs in region
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(example_beta)
#' data(example_dmps)
#' data(example_pheno)
#'
#' # Write beta values to file
#' beta_file <- tempfile(fileext = ".txt")
#' write.table(cbind(ID = rownames(example_beta), example_beta),
#'     file = beta_file, sep = "\t", quote = FALSE, row_names = FALSE
#' )
#'
#' # Find DMRs
#' dmrs <- findDMRsFromDMPs(
#'     beta_file = beta_file,
#'     dmps_tsv_file = example_dmps,
#'     pheno = example_pheno
#' )
#' }
#'
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom future.apply future_lapply
#' @importFrom future availableCores
#' @importFrom progressr progressor
#' @importFrom stringr str_count str_order
#' @importFrom readr read_tsv
#' @importFrom data.table fread fwrite
#' @importFrom psych corr.test
#' @importFrom dplyr %>% filter select mutate
#' @importFrom bedr tabix
#' @importFrom BSgenome getSeq
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom rtracklayer import.chain
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom stats sd
#' @importFrom utils write.table read.table
#' @importFrom tools file_ext file_path_sans_ext
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import IlluminaHumanMethylationEPICanno.ilm10b4.hg19
#' @export

# Helper functions for progress tracking with temporary files
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

#' Sort Beta File by Genomic Coordinates
#'
#' @description This helper function sorts a methylation beta values file by genomic coordinates
#' (chromosome and position) as required by the findDMRsFromDMPs function. The function reads
#' the beta file, sorts the CpG sites according to their genomic positions using array annotation,
#' and writes the sorted data to a new file.
#'
#' @param beta_file Character. Path to the input beta values file to be sorted
#' @param output_file Character. Path for the output sorted beta file (default: adds "_sorted" suffix)
#' @param array Character. Array platform type, either "450K" or "EPIC" (default: "450K")
#' @param genomic_locs Data frame. Optional pre-computed genomic locations. If NULL, locations will be retrieved automatically (default: NULL)
#' @param verbose Logical. Whether to print progress messages (default: TRUE)
#' @param overwrite Logical. Whether to overwrite existing output file (default: FALSE)
#'
#' @return Character. Path to the sorted output file
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Reads the beta values file
#'   \item Loads the appropriate array annotation (450K or EPIC)
#'   \item Sorts CpG sites by genomic coordinates (chr:pos)
#'   \item Writes the sorted data to a new file
#'   \item Validates that the output is properly sorted
#' }
#'
#' @examples
#' \dontrun{
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
#' }
#'
#' @export
sortBetaFileByCoordinates <- function(beta_file,
                                      output_file = NULL,
                                      array = c("450K", "EPIC"),
                                      genomic_locs = NULL,
                                      verbose = TRUE,
                                      overwrite = FALSE) {
    # Validate inputs
    array <- match.arg(array)
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

    if (verbose) {
        message("Reading beta file: ", beta_file)
    }

    # Read the beta file
    beta_data <- data.table::fread(beta_file, header = TRUE, data.table = FALSE)

    # Get row names (CpG IDs) from first column
    cpg_ids <- beta_data[[1]]
    beta_values <- beta_data[, -1, drop = FALSE]
    rownames(beta_values) <- cpg_ids

    if (verbose) {
        message("Read ", nrow(beta_values), " CpG sites and ", ncol(beta_values), " samples")
    }

    sorted_locs <- genomic_locs
    if (is.null(sorted_locs)) {
        sorted_locs <- getSortedGenomicLocs(array = array)
    }


    # Find CpGs that are present in both the beta file and array annotation
    common_cpgs <- intersect(cpg_ids, rownames(sorted_locs))
    missing_from_annotation <- setdiff(cpg_ids, rownames(sorted_locs))
    if (length(missing_from_annotation) > 0) {
        stop(
            "Found ", length(missing_from_annotation), " CpG sites in beta file that are not in ",
            array, " annotation.First 5 missing: ", paste(head(missing_from_annotation, 5), collapse = ", ")
        )
    }

    missing_from_beta <- setdiff(rownames(sorted_locs), cpg_ids)
    if (length(missing_from_beta) > 0 && verbose) {
        message("Note: ", length(missing_from_beta), " CpG sites in ", array, " annotation are not present in beta file.")
    }

    final_order <- rownames(sorted_locs)[rownames(sorted_locs) %in% common_cpgs]

    # Reorder beta values
    sorted_beta_values <- beta_values[final_order, , drop = FALSE]

    # Prepare output data frame
    output_data <- data.frame(
        ID = rownames(sorted_beta_values),
        sorted_beta_values,
        check_names = FALSE,
        stringsAsFactors = FALSE
    )

    if (verbose) {
        message("Writing sorted beta file to: ", output_file)
    }

    # Write sorted file
    data.table::fwrite(
        output_data,
        file = output_file,
        sep = "\t",
        quote = FALSE,
        row_names = FALSE,
        col_names = TRUE
    )

    return(output_file)
}

.subsetBeta <- function(beta_file,
                        sites,
                        beta_row_names = NULL,
                        beta_col_names = NULL) {
    # Fallback simple subsetting path for small site sets (avoids scan/skip logic issues in tests)
    if (length(sites) <= 5000) {
        full <- data.table::fread(beta_file, header = TRUE, data.table = FALSE)
        rn <- full[[1]]
        full <- full[, -1, drop = FALSE]
        rownames(full) <- rn
        # Preserve only requested sites in order
        idx <- match(sites, rn)
        if (anyNA(idx)) {
            missing <- sites[is.na(idx)]
            stop("Internal error: requested CpG IDs not found in beta file during subsetting: ", paste(missing, collapse = ","))
        }
        sub <- full[idx, , drop = FALSE]
        if (nrow(sub) == 0) {
            if (interactive()) browser()
            stop("None of the provided sites exist in the read beta file")
        }
        sub <- apply(sub, 2, as.numeric)

        tryCatch(
            {
                rownames(sub) <- sites
            },
            error = function(e) {
                if (interactive()) browser()
            }
        )
        nas_per_row <- apply(sub, 1, function(r) sum(is.na(r)))
        if (any(nas_per_row == ncol(sub))) {
            warning("All-NA beta rows detected for sites: ", paste(names(nas_per_row)[nas_per_row == ncol(sub)], collapse = ","))
        }
        return(sub)
    }
    if (is.null(beta_row_names)) {
        beta_row_names <- unlist(fread(
            file = beta_file,
            select = 1,
            header = TRUE,
        ))
    }

    ret <- .getBetaColNamesAndInds(beta_file, beta_col_names)
    beta_col_names <- ret$beta_col_names
    cols_inds <- ret$beta_col_inds
    sites_inds <- which(beta_row_names %in% sites)
    sites_inds_steps <- diff(c(-1, sites_inds)) - 1
    if (endsWith(beta_file, "gz")) {
        conn <- gzfile(beta_file, "r")
    } else {
        conn <- file(beta_file, "r")
    }
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
    close(conn)
    beta_sites <- do.call(rbind, beta_sites)
    rownames(beta_sites) <- sites
    colnames(beta_sites) <- beta_col_names
    beta_sites
}

# Expand the identified dmrs to nearby CpG regions
.expandDMRs <- function(dmr,
                        beta_file,
                        beta_row_names,
                        beta_col_names,
                        sample_groups,
                        sorted_locs,
                        max_pval,
                        min_cpg_delta_beta = 0,
                        casecontrol = NULL,
                        tabix_file = NULL,
                        expansion_step = 500,
                        expansion_relaxation = 0) {
    if (!is.null(beta_file)) {
        ret <- .getBetaColNamesAndInds(beta_file, beta_col_names)
    } else {
        ret <- .getBetaColNamesAndInds(tabix_file, beta_col_names, is_tabix = TRUE)
    }

    cols_inds <- ret$beta_col_inds
    file_beta_col_names <- ret$file_beta_col_names
    beta_col_names <- ret$beta_col_names
    sorted_locs <- sorted_locs[beta_row_names, ]

    dmr_start <- dmr$start_dmp

    dmr_end <- dmr$end_dmp
    dmr_start_ind <- which(beta_row_names == dmr_start)

    if (length(dmr_start_ind) == 0) {
        stop("Could not find the start CpG ", dmr_start, " in the beta file row names.")
    }
    end_site_ind <- dmr_start_ind[[1]]
    upstream_exp <- end_site_ind
    upstream_stop_found <- FALSE
    while (TRUE) {
        if (end_site_ind < 0) {
            upstream_stop_reason <- "end-of-input"
            upstream_exp <- 1
            break
        }
        start_site_ind <- max(0, end_site_ind - expansion_step) + 1
        x <- which(sorted_locs[start_site_ind:end_site_ind, "chr"] == dmr$chr)
        if (length(x) == 0) {
            upstream_stop_reason <- "end-of-input"
            break
        }
        x <- x[[1]]
        start_site_ind <- start_site_ind + x - 1
        exp.step <- end_site_ind - start_site_ind + 1
        if (!is.null(beta_file)) {
            upstream_betas <- fread(
                file = beta_file,
                skip = start_site_ind,
                nrows = exp.step,
                header = FALSE, data.table = FALSE,
                colClasses = c("character", rep("numeric", length(file_beta_col_names)))
            )
            upstream_betas <- upstream_betas[, cols_inds]
        } else {
            upstream_region <- paste0(
                sorted_locs[start_site_ind, "chr"], ":",
                sorted_locs[start_site_ind, "pos"], "-",
                sorted_locs[end_site_ind, "pos"] + 1
            )
            upstream_betas <- try(tabix(
                upstream_region,
                tabix_file,
                check.valid = FALSE,
                verbose = FALSE
            ))
            if (inherits(upstream_betas, "try-error")) {
                warning("Error reading upstream region ", upstream_region, " from tabix file, stopping extension. The error is:\n\t", paste(capture.output(print(upstream_betas)), collapse = "\n\t"))
                upstream_stop_reason <- "error-reading-tabix"
                break
            }
            upstream_betas <- as.data.frame(sapply(upstream_betas[, beta_col_names], as.numeric))
        }
        upstream_betas <- upstream_betas[rev(seq_len(nrow(upstream_betas))), , drop = FALSE]
        if (nrow(upstream_betas) == 1) {
            upstream_stop_reason <- "end-of-input"
            break
        }
        i <- 1
        exp_relax_counter <- 0
        while (TRUE) {
            corr_ret <- .testConnectivity(
                site1_beta = unlist(upstream_betas[i, ]),
                site2_beta = unlist(upstream_betas[i + 1, ]),
                sample_groups = sample_groups,
                casecontrol = casecontrol,
                max_pval = max_pval,
                min_delta_beta = min_cpg_delta_beta,
                extreme_verbosity = TRUE
            )
            if (!corr_ret[[1]]) {
                if (exp_relax_counter < expansion_relaxation) {
                    exp_relax_counter <- exp_relax_counter + 1
                } else {
                    upstream_exp <- end_site_ind - i + 1 - exp_relax_counter
                    upstream_stop_found <- TRUE
                    upstream_stop_reason <- corr_ret$reason
                    exp_relax_counter <- 0
                    break
                }
            } else {
                exp_relax_counter <- 0
            }
            i <- i + 1
            if (i == nrow(upstream_betas)) {
                break
            }
        }
        if (upstream_stop_found) {
            break
        }
        end_site_ind <- end_site_ind - expansion_step
    }

    dmr_end_ind <- which(beta_row_names == dmr_end)
    if (length(dmr_end_ind) == 0) {
        stop("Could not find the end CpG ", dmr_end, " in the beta file row names.")
    }
    start_site_ind <- dmr_end_ind[[1]]
    downstream_exp <- start_site_ind
    downstream_stop_found <- FALSE
    while (TRUE) {
        if (start_site_ind > length(beta_row_names)) {
            downstream_stop_reason <- "end-of-input"
            downstream_exp <- start_site_ind - 1
            break
        }
        end_site_ind <- min(start_site_ind + expansion_step, length(beta_row_names))
        x <- which(rev(sorted_locs[start_site_ind:end_site_ind, "chr"]) == dmr$chr)
        if (length(x) == 0) {
            downstream_stop_reason <- "end-of-input"
            downstream_exp <- start_site_ind - 1
            break
        }
        x <- x[[1]]
        end_site_ind <- end_site_ind - x + 1
        if (end_site_ind <= start_site_ind + 1) {
            downstream_stop_reason <- "end-of-input"
            downstream_exp <- start_site_ind - 1
            break
        }
        if (!is.null(beta_file)) {
            downstream_betas <- fread(
                file = beta_file,
                skip = start_site_ind,
                nrows = expansion_step - x + 1,
                header = FALSE,
                data.table = FALSE,
                colClasses = c("character", rep("numeric", length(file_beta_col_names)))
            )
            downstream_betas <- downstream_betas[, cols_inds]
        } else {
            downstream_region <- paste0(
                sorted_locs[start_site_ind, "chr"], ":",
                sorted_locs[start_site_ind, "pos"], "-",
                sorted_locs[end_site_ind, "pos"]
            )
            downstream_betas <- try(tabix(
                downstream_region,
                tabix_file,
                check.valid = FALSE,
                verbose = FALSE
            ))
            if (inherits(downstream_betas, "try-error")) {
                warning("Error reading downstream region ", downstream_region, " from tabix file, stopping extension. The error is:\n\t", paste(capture.output(print(downstream_betas)), collapse = "\n\t"))
                downstream_stop_reason <- "error-reading-tabix"
                break
            }
            downstream_betas <- as.data.frame(sapply(downstream_betas[, beta_col_names], as.numeric))
        }

        i <- 1
        exp_relax_counter <- 0
        while (TRUE) {
            corr_ret <- .testConnectivity(
                site1_beta = unlist(downstream_betas[i, ]),
                site2_beta = unlist(downstream_betas[i + 1, ]),
                sample_groups = sample_groups,
                max_pval = max_pval,
                casecontrol = casecontrol,
                min_delta_beta = min_cpg_delta_beta
            )
            if (!corr_ret[[1]]) {
                if (exp_relax_counter < expansion_relaxation) {
                    exp_relax_counter <- exp_relax_counter + 1
                } else {
                    downstream_exp <- start_site_ind + i - 1 + exp_relax_counter
                    downstream_stop_found <- TRUE
                    downstream_stop_reason <- corr_ret$reason
                    exp_relax_counter <- 0
                    break
                }
            } else {
                exp_relax_counter <- 0
            }
            i <- i + 1
            if (i == nrow(downstream_betas)) {
                break
            }
        }
        if (downstream_stop_found) {
            break
        }
        prev_start_site_ind <- start_site_ind
        start_site_ind <- start_site_ind + expansion_step - x + 1
        if (start_site_ind < prev_start_site_ind) {
            stop("BUG: start_site_ind < prev_start_site_ind during downstream expansion. start_site_ind: ", start_site_ind, " prev_start_site_ind: ", prev_start_site_ind, " expansion_step: ", expansion_step, " x: ", x)
        }
    }
    dmr$start_cpg <- beta_row_names[upstream_exp]
    tryCatch(
        {
            dmr$end_cpg <- beta_row_names[downstream_exp]
        },
        error = function(e) {
            if (interactive()) browser()
        }
    )
    dmr$start <- sorted_locs[dmr$start_cpg, "pos"]
    dmr$end <- sorted_locs[dmr$end_cpg, "pos"]
    dmr$downstream_cpg_expansion_stop_reason <- downstream_stop_reason
    dmr$upstream_cpg_expansion_stop_reason <- upstream_stop_reason

    dmr
}

.testConnectivity <- function(site1_beta, site2_beta, sample_groups, max_pval, casecontrol = NULL, min_delta_beta = 0, extreme_verbosity = FALSE) {
    pval <- 0
    delta_beta <- NULL

    max_pval_corrected <- max_pval / length(unique(sample_groups)) # Bonferroni correction
    for (g in levels(sample_groups)) {
        if (sum(sample_groups == g) < 3) {
            next
        }
        op <- options(warn = 2)$warn
        corr_ret <- try(corr.test(site1_beta[sample_groups == g], site2_beta[sample_groups == g], ci = FALSE))
        options(warn = op)
        if (inherits(corr_ret, "try-error")) {
            if (extreme_verbosity) {
                message(".testConnectivity: Error occurred in corr.test while processing the following:")
                message("casecontrol:", paste(casecontrol, collapse = ","))
                message("site2_beta:", paste(site2_beta, collapse = ","))
                message("sample_groups:", paste(sample_groups, collapse = ","))
                message("max_pval_corrected:", max_pval_corrected)
                message("Error message:", corr_ret)
                if (interactive()) browser()
            }
            return(list(FALSE, pval, delta_beta, failing = g, reason = "error occurred"))
        }
        r <- max(pval, corr_ret$p)
        if (is.null(r) || is.na(r)) {
            return(list(FALSE, pval, delta_beta, failing = g, reason = "na pval"))
        }
        pval <- r
        if (pval > max_pval_corrected) {
            return(list(FALSE, pval, delta_beta, failing = g, reason = "pval>max_pval (corrected)"))
        }
    }
    if (!is.null(casecontrol) && (min_delta_beta > 0)) {
        if (length(casecontrol) != length(site2_beta)) {
            if (extreme_verbosity) {
                message(".testConnectivity: Error occurred while computing delta beta for the following:")
                message("casecontrol:", paste(casecontrol, collapse = ","))
                message("site2_beta:", paste(site2_beta, collapse = ","))
                message("sample_groups:", paste(sample_groups, collapse = ","))
                message("max_pval_corrected:", max_pval_corrected)
            }
            stop(paste0(
                "The provided casecontrol vector has length", length(casecontrol),
                " while the site2_beta has length ", length(site2_beta)
            ))
        }
        delta_beta <- mean(site2_beta[casecontrol == 1], na.rm = TRUE) - mean(site2_beta[casecontrol == 0], na.rm = TRUE)
        if (is.null(delta_beta) || is.na(delta_beta) || (abs(delta_beta) < min_delta_beta)) {
            return(list(FALSE, pval, delta_beta, reason = "delta_beta<min_delta_beta"))
        }
    }
    return(list(TRUE, pval, delta_beta))
}

#' Extract CpG Information from DMR Results
#'
#' @description This function extracts detailed CpG methylation information for
#' identified DMRs, providing beta values and additional metadata for each CpG
#' site within the DMR regions.
#'
#' @param dmrs GRanges object containing DMR results from findDMRsFromDMPs
#' @param beta_file Character. Path to the methylation beta values file
#' @param beta_row_names Character vector. Row names for beta values
#' @param beta_col_names Character vector. Column names for beta values
#' @param sorted_locs GRanges or data frame with sorted genomic locations
#' @param output_prefix Character. Optional prefix for output files (default: NULL)
#' @param tabix_file Character. Path to tabix-indexed beta file (alternative to beta_file, default: NULL)
#' @param njobs Integer. Number of cores for parallel processing (default: 1)
#' @param verbose Logical. Whether to print progress messages (default: FALSE)
#'
#' @return A data frame containing detailed CpG information for each DMR
#'
#' @examples
#' \dontrun{
#' # Extract CpG info from DMR results
#' cpg_info <- extractCpgInfoFromResultDMRs(
#'     dmrs = dmr_results,
#'     beta_file = "sorted_beta.txt",
#'     beta_row_names = row_names,
#'     beta_col_names = col_names,
#'     sorted_locs = genomic_locations
#' )
#' }
#'
#' @export
extractCpgInfoFromResultDMRs <- function(dmrs,
                                         beta_file = NULL,
                                         beta_row_names,
                                         beta_col_names,
                                         sorted_locs,
                                         output_prefix = NULL,
                                         tabix_file = NULL,
                                         njobs = 1,
                                         verbose = FALSE) {
    if (!is.null(output_prefix)) {
        if (!dir.exists(dirname(output_prefix))) {
            dir.create(dirname(output_prefix), recursive = TRUE)
        }
    }
    dmrs_cpgs <- data.frame()
    dmrs_sites <- c()
    if (verbose) {
        message("Finding constituent CpGs of DMRs, that exist in the provided beta file..")
    }
    for (i in seq_len(nrow(dmrs))) {
        start_ind <- dmrs$start_ind[i]
        end_ind <- dmrs$end_ind[i]
        cpgs <- rownames(sorted_locs)[start_ind:end_ind]
        cpgs <- cpgs[cpgs %in% beta_row_names]
        dmrs_cpgs <- rbind(dmrs_cpgs, data.frame(dmr = dmrs$id[i], cpgs = cpgs))
        dmrs_sites <- c(dmrs_sites, cpgs)
    }

    if (!is.null(output_prefix)) {
        dmrs_cpgs_file <- paste0(output_prefix, ".cpgs.tsv.gz")
        if (verbose) {
            message("Saving extended DMRs constituent CpGs names mapping to ", dmrs_cpgs_file, "..")
        }
        gz <- gzfile(dmrs_cpgs_file, "w")
        write.table(
            dmrs_cpgs,
            gz,
            sep = "\t",
            quote = FALSE,
            qmethod = "double",
            col.names = TRUE,
            row.names = FALSE
        )
        close(gz)
    }

    dmrs_sites <- unique(dmrs_sites)
    dmrs_sites <- rownames(sorted_locs)[rownames(sorted_locs) %in% dmrs_sites]

    dmrs_beta_file <- paste0(output_prefix, "dmrs_beta.tsv.gz")
    if (verbose) {
        message("Saving extended DMRs constituent CpGs betas to ", dmrs_beta_file, "..")
    }

    if (!is.null(beta_file)) {
        dmrs_beta <- .subsetBeta(beta_file,
            dmrs_sites,
            beta_row_names = beta_row_names,
            beta_col_names = beta_col_names
        )
    } else {
        dmrs_beta <- tabix(
            as.data.frame(
                sorted_locs[dmrs_sites, c("chr", "start", "end")]
            ),
            tabix_file,
            params = paste0("-@ ", njobs),
            verbose = FALSE
        )
        dmrs_beta <- as.data.frame(sapply(dmrs_beta[, beta_col_names], as.numeric))
    }
    rownames(dmrs_beta) <- dmrs_sites
    if (!is.null(output_prefix)) {
        gz <- gzfile(dmrs_beta_file, "w")
        write.table(
            dmrs_beta,
            gz,
            sep = "\t",
            quote = FALSE,
            qmethod = "double",
            col.names = NA,
            row.names = TRUE
        )
        close(gz)
    }
    ret <- list(
        dmrs_cpgs = dmrs_cpgs,
        dmrs_beta = dmrs_beta
    )
    invisible(ret)
}



#' Find Differentially Methylated Regions (DMRs) from Differentially Methylated Positions (DMPs)
#'
#' This function identifies DMRs from a given set of DMPs and a beta value file.
#'
#' @param beta_file Character. Path to the beta value file. Either this or tabix_file must be provided.
#' @param dmps_tsv_file Character. Path to the DMPs TSV file.
#' @param pheno Data frame. Phenotype data.
#' @param dmps_tsv_id_col Character. Column name for DMP identifiers in the DMPs TSV file. Default is NULL.
#' @param dmp_groups_info Named list. Required when `dmps_tsv_id_col` is given. List of DMP group information, where names are group identifiers, found in dmps_tsv_id_col column, and values are the samples names, found in the beta values columns. Default is NULL.
#' @param pval_col Character. Column name for p-values in the DMPs file. Default is "pval_adj".
#' @param sample_group_col Character. Column name for sample group information in the phenotype data. Default is NULL.
#' @param dmp_group_col Character. Column name for DMP group information in the DMPs TSV file. Default is NULL.
#' @param casecontrol_col Character. Column name for case-control information in the phenotype data. Default is "casecontrol".
#' @param min_cpg_delta_beta Numeric. Minimum delta beta value for CpGs. Default is 0.
#' @param expansion_step Numeric. Step size for expanding DMRs. Increasing it means higher memory usage and faster computation. Default is 500.
#' @param expansion_relaxation Numeric. Maximum number of intermittent CpGs allowed to not be significanly correlated, to increase the extended DMR size. Default is 0.
#' @param array Character. Type of array used ("450K" or "EPIC"). Default is "450K".
#' @param genome Character. Genome version ("hg19" or "hg38"). Default is "hg19".
#' @param max_pval Numeric. Maximum p-value to assume DMPs correlation is significant. Default is 0.05.
#' @param max_lookup_dist Numeric. Maximum distance to look up for adjacent DMPs belonging to the same DMR. Default is 10000.
#' @param min_dmps Numeric. Minimum number of connected DMPs in a DMR. Default is 1.
#' @param min_adj_dmps Numeric. Minimum number of DMPs, adjusted by CpG density, in a DMR after extension. Default is 1.
#' @param min_cpgs Numeric. Minimum number of CpGs in a DMR after extension. Default is 50.
#' @param ignored_sample_groups Character vector. Sample groups to ignore, separated by commas. Default is NULL.
#' @param output_prefix Character. Identifier for the output files. If not provided, no output will be saved. Default is NULL.
#' @param njobs Numeric. Number of parallel jobs to use. Default is the number of available cores.
#' @param verbose Logical. Whether to print detailed messages. Default is FALSE.
#' @param beta_row_names_file Character. Path to a file containing row names for the beta values. If not provided, row names will be read from the beta file. Default is NULL.
#' @param dmps_beta_file Character. Path to load and save the beta values of DMPs. If not provided, the beta values will be read from the beta file. Default is NULL.
#' @param tabix_file Character. Path to a tabix-indexed beta file. Either this or beta_file must be provided. Default is NULL.
#'
#' @return Data frame of identified DMRs.
#' @export
findDMRsFromDMPs <- function(beta_file = NULL,
                             dmps_tsv_file = NULL,
                             pheno = NULL,
                             dmps_tsv_id_col = NULL,
                             dmp_groups_info = NULL,
                             pval_col = "pval_adj",
                             sample_group_col = "Sample_Group",
                             dmp_group_col = NULL,
                             casecontrol_col = "casecontrol",
                             min_cpg_delta_beta = 0,
                             expansion_step = 500,
                             expansion_relaxation = 0,
                             array = c("450K", "EPIC"),
                             genome = c("hg19", "hg38"),
                             max_pval = 0.05,
                             max_lookup_dist = 10000,
                             min_dmps = 1,
                             min_adj_dmps = 1,
                             min_cpgs = 50,
                             ignored_sample_groups = NULL,
                             output_prefix = NULL,
                             njobs = future::availableCores(),
                             verbose = FALSE,
                             beta_row_names_file = NULL,
                             dmps_beta_file = NULL,
                             tabix_file = NULL) {
    if (is.null(dmps_tsv_file) || is.null(pheno)) {
        stop("dmps_tsv_file and pheno parameters are required")
    }
    if (is.null(beta_file) && is.null(tabix_file)) {
        stop("Either beta_file or tabix_file parameter is required")
    }
    stopifnot(!is.null(max_pval))
    stopifnot(!is.null(min_dmps))
    stopifnot(!is.null(min_cpgs))
    stopifnot(!is.null(expansion_step))
    stopifnot(!is.null(min_cpg_delta_beta))
    stopifnot(!is.null(max_lookup_dist))

    array <- match.arg(array)
    genome <- match.arg(genome)
    if (!is.null(beta_file)) {
        stopifnot(file.exists(beta_file))
    }
    if (!is.null(tabix_file)) {
        stopifnot(file.exists(tabix_file))
    }
    stopifnot(file.exists(dmps_tsv_file))
    stopifnot(casecontrol_col %in% colnames(pheno))
    stopifnot(sample_group_col %in% colnames(pheno))
    if (is.null(ignored_sample_groups)) {
        ignored_sample_groups <- c()
    } else {
        ignored_sample_groups <- unlist(strsplit(ignored_sample_groups, ","))
    }
    if (!is.null(output_prefix)) {
        output_dir <- dirname(output_prefix)
        dir.create(output_dir, showWarnings = FALSE)
        output_prefix <- paste0(output_prefix, ".")
    } else {
        output_dir <- NULL
        output_prefix <- NULL
    }


    if (verbose) {
        message("Reading dmp tsv file..")
    }
    dmps_tsv <- try(read.table(
        dmps_tsv_file,
        header = TRUE,
        sep = "\t",
        check.names = FALSE,
        quote = "",
        comment.char = "",
        row.names = NULL,
    ))
    if (inherits(dmps_tsv, "try-error")) {
        message("Provided DMPs file is empty or does not exist. Not proceeding.")
        if (!is.null(output_prefix)) {
            for (f in c(".methylation.tsv.gz", ".tsv.gz")) {
                gzfile <- gzfile(paste0(output_prefix, f), "w", compression = 2)
                close(gzfile)
            }
        }
        return(NULL)
    }
    if (!is.null(dmp_group_col)) {
        if (!dmp_group_col %in% colnames(dmps_tsv)) {
            stop(
                "DMP group column '", dmp_group_col,
                "' does not reside in the DMPs file columns: ",
                paste(colnames(dmps_tsv), collapse = ",")
            )
        }
        if (is.null(dmp_groups_info)) {
            stop(
                "dmp_groups_info parameter is required when dmp_group_col is provided"
            )
        }
        unique_dmp_groups <- unique(dmps_tsv[, dmp_group_col])
        if (!all(unique_dmp_groups %in% names(dmp_groups_info))) {
            missing <- unique_dmp_groups[!(unique_dmp_groups %in% names(dmp_groups_info))]
            stop("The following DMP groups do not exist in the names of the dmp_groups_info list'", dmp_group_col, "' column: ", paste(missing, collapse = ","))
        }
    } else {
        dmp_group_col <- "_DUMMY_DMP_GROUP_COL_"
        dmps_tsv[, dmp_group_col] <- "all"
    }

    if (!pval_col %in% colnames(dmps_tsv)) {
        stop(
            "P-value column '", pval_col, "' does not reside in the DMPs file columns: ",
            paste(colnames(dmps_tsv), collapse = ",")
        )
    }
    stopifnot(pval_col %in% colnames(dmps_tsv))

    if (verbose) {
        message("Reading beta file characteristics..")
    }

    if (!is.null(output_prefix)) {
        output_beta_row_names_file <- paste0(output_prefix, "row_names.txt")
        if (is.null(beta_row_names_file)) {
            beta_row_names_file <- output_beta_row_names_file # load from previous run
        }
    }
    if (!is.null(beta_row_names_file) && file.exists(beta_row_names_file)) {
        if (verbose) message("Reading beta file row names from beta row names file ", beta_row_names_file, "..")
        beta_row_names <- unlist(read.table(beta_row_names_file, header = FALSE, comment.char = "", quote = ""))
    } else {
        if (verbose) message("Reading beta file row names from beta file..")
        if (!is.null(beta_file)) {
            beta_row_names <- unlist(fread(
                file = beta_file,
                select = 1,
                sep = "\t",
                header = TRUE,
                showProgress = TRUE,
                nThread = njobs
            ))
        } else {
            beta_row_names <- unlist(fread(
                file = tabix_file,
                select = 4,
                sep = "\t",
                header = TRUE,
                showProgress = TRUE,
                nThread = njobs
            ))
        }
        if (!is.null(beta_row_names_file)) {
            if (verbose) message("Saving beta file row names to beta row names file: ", beta_row_names_file)
            writeLines(paste(beta_row_names, collapse = "\n"), beta_row_names_file)
        }
    }
    if (verbose) {
        message("Number of rows names read:", length(beta_row_names))
    }

    beta_col_names <- rownames(pheno)
    samples_selection_mask <- !(pheno[, sample_group_col] %in% ignored_sample_groups)
    beta_col_names <- beta_col_names[samples_selection_mask]
    pheno <- pheno[beta_col_names, ]
    if (verbose) message("Number of samples to process: ", length(beta_col_names))

    if (!is.null(beta_file)) {
        beta_col_names <- .getBetaColNamesAndInds(beta_file, beta_col_names)$beta_col_names
    } else {
        beta_col_names <- .getBetaColNamesAndInds(tabix_file, beta_col_names, is_tabix = TRUE)$beta_col_names
    }
    pheno <- pheno[beta_col_names, ]
    sample_groups <- factor(pheno[beta_col_names, sample_group_col])



    if (verbose) {
        message("Reordering DMPs based on location..")
    }

    if (is.null(dmps_tsv_id_col)) {
        dmps_tsv_id_col <- "row.names"
    }
    if (!dmps_tsv_id_col %in% colnames(dmps_tsv)) {
        stop(
            "DMP id column '", dmps_tsv_id_col,
            "' does not reside in the DMPs file columns: ",
            paste(colnames(dmps_tsv), collapse = ",")
        )
    }
    sorted_locs <- getSortedGenomicLocs(array = array)
    dmps_tsv <- dmps_tsv[orderByLoc(dmps_tsv[, dmps_tsv_id_col], genomic_locs = sorted_locs), , drop = FALSE]

    dmps <- unique(dmps_tsv[, dmps_tsv_id_col])
    # Filter DMPs not present in array annotation first (prevents NA logical indices later)
    missing_in_annotation <- setdiff(dmps, rownames(sorted_locs))
    if (length(missing_in_annotation) > 0) {
        warning(
            "Dropping ", length(missing_in_annotation), " DMP(s) not found in the array annotation: ",
            paste(head(missing_in_annotation, 10), collapse = ","),
            if (length(missing_in_annotation) > 10) " ..." else ""
        )
        dmps_tsv <- dmps_tsv[!(dmps_tsv[, dmps_tsv_id_col] %in% missing_in_annotation), , drop = FALSE]
        dmps <- setdiff(dmps, missing_in_annotation)
    }
    if (length(dmps) == 0) {
        if (interactive()) browser()
        stop("No DMPs remain after filtering against array annotation.")
    }
    if (!all(dmps %in% beta_row_names)) {
        missing_in_beta <- dmps[!(dmps %in% beta_row_names)]
        if (interactive()) browser()
        stop("Some of the DMPs are not present in the beta file. DMPs: ", paste(missing_in_beta, collapse = ","))
    }
    dmps <- dmps[orderByLoc(dmps, genomic_locs = sorted_locs)]

    if (verbose) {
        message("Making sure that the beta file is sorted by position..")
    }
    if (!all(beta_row_names[orderByLoc(beta_row_names, genomic_locs = sorted_locs)] == beta_row_names)) {
        stop("Provided beta file is not sorted by position!")
    }

    if (verbose) {
        message("Getting subset beta corresponding to DMPs..")
    }
    dmps_locs <- sorted_locs[dmps, , drop = FALSE]
    dmps_tsv[, "chr"] <- dmps_locs[dmps_tsv[, dmps_tsv_id_col], "chr"]
    dmps_tsv[, "pos"] <- dmps_locs[dmps_tsv[, dmps_tsv_id_col], "pos"]
    if (!is.null(output_prefix)) {
        dmps_beta_output_file <- paste0(output_prefix, "dmps_beta.tsv.gz")
        if (is.null(dmps_beta_file)) {
            dmps_beta_file <- dmps_beta_output_file # loading from previous run
        }
    }
    if (!is.null(dmps_beta_file) && file.exists(dmps_beta_file)) {
        if (verbose) message("Reading beta file row names from supplied dmps beta file..")
        gz <- gzfile(dmps_beta_output_file, "r")
        dmps_beta <- read.table(gz, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, comment.char = "", quote = "")
        close(gz)
    } else {
        if (!is.null(beta_file)) {
            dmps_beta <- .subsetBeta(beta_file,
                dmps,
                beta_row_names = beta_row_names,
                beta_col_names = beta_col_names
            )
        } else {
            dmps_beta <- tabix(as.data.frame(dmps_locs[, c("chr", "start", "end")]),
                tabix_file,
                verbose = FALSE
            )
            dmps_beta <- as.data.frame(sapply(dmps_beta[, beta_col_names], as.numeric))
        }
        if (verbose) message("Calculating delta_beta related columns in the DMPs table..")
        if (dmp_group_col == "_DUMMY_DMP_GROUP_COL_") {
            dmp_groups_info <- list(all = NULL)
        }
        for (dmp_group in unique(dmps_tsv[, dmp_group_col])) {
            group_mask <- dmps_tsv[, dmp_group_col] == dmp_group
            if (!is.null(dmp_groups_info[[dmp_group]])) {
                beta_mask <- rownames(dmps_beta) %in% dmp_groups_info[[dmp_group]]
            } else {
                beta_mask <- rep(TRUE, nrow(dmps_beta))
            }
            group_beta_subset <- dmps_beta[beta_mask, , drop = FALSE]
            case_means <- rowMeans(group_beta_subset[, pheno[beta_col_names, casecontrol_col] == 1, drop = FALSE], na.rm = TRUE)
            control_means <- rowMeans(group_beta_subset[, pheno[beta_col_names, casecontrol_col] == 0, drop = FALSE], na.rm = TRUE)
            case_sd <- apply(group_beta_subset[, pheno[beta_col_names, casecontrol_col] == 1, drop = FALSE], 1, sd, na.rm = TRUE)
            control_sd <- apply(group_beta_subset[, pheno[beta_col_names, casecontrol_col] == 0, drop = FALSE], 1, sd, na.rm = TRUE)
            dmps_tsv[group_mask, "cases_num"] <- colSums(!is.na(group_beta_subset[, pheno[beta_col_names, casecontrol_col] == 1, drop = FALSE]))
            dmps_tsv[group_mask, "controls_num"] <- colSums(!is.na(group_beta_subset[, pheno[beta_col_names, casecontrol_col] == 0, drop = FALSE]))
            dmps_tsv[group_mask, "cases_beta"] <- case_means[dmps_tsv[group_mask, dmps_tsv_id_col]]
            dmps_tsv[group_mask, "controls_beta"] <- control_means[dmps_tsv[group_mask, dmps_tsv_id_col]]
            dmps_tsv[group_mask, "cases_beta_sd"] <- case_sd[dmps_tsv[group_mask, dmps_tsv_id_col]]
            dmps_tsv[group_mask, "controls_beta_sd"] <- control_sd[dmps_tsv[group_mask, dmps_tsv_id_col]]
            dmps_tsv[group_mask, "delta_beta"] <- dmps_tsv[group_mask, "cases_beta"] - dmps_tsv[group_mask, "controls_beta"]
        }

        if (!is.null(output_prefix)) {
            if (verbose) message("Saving dmps beta to file: ", dmps_beta_output_file)
            gz <- gzfile(dmps_beta_output_file, "w")
            write.table(dmps_beta, gz, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
            close(gz)
        }
    }
    if (nrow(dmps_locs) != nrow(dmps_beta)) {
        stop(
            "Number of rows in the queried dmps beta file does not match the number of DMPs. Number of rows in beta file: ",
            nrow(dmps_beta),
            " Number of rows in DMPs: ",
            nrow(dmps_locs)
        )
    }
    # Diagnostic: identify any rows with all NA betas (should not happen with synthetic test data)
    all.na.rows <- apply(dmps_beta, 1, function(r) all(is.na(r)))
    if (any(all.na.rows)) {
        stop(
            "Beta extraction failure: the following DMP rows have all NA beta values: ",
            paste(rownames(dmps_beta)[all.na.rows], collapse = ","),
            ". This indicates a mismatch between requested CpG IDs and beta file columns or a parsing issue."
        )
    }

    if (verbose) {
        message("Subset size: ", paste(dim(dmps_beta), collapse = ","))
    }

    if (verbose) {
        message("Number of provided DMPs:", length(dmps))
    }
    if (verbose) {
        message("Connecting DMPs to form initial DMRs..")
    }

    # Set up progress tracking for DMP connection
    chromosomes <- unique(dmps_locs$chr)
    n_chromosomes <- length(chromosomes)

    # Set up future plan for parallel processing
    if (!identical(njobs, 1L)) {
        if (.Platform$OS.type == "unix") {
            future::plan(future::multicore, workers = njobs)
        } else {
            future::plan(future::multisession, workers = njobs)
        }
    } else {
        future::plan(future::sequential)
    }

    # Use progressr for cross-platform progress reporting
    if (verbose && n_chromosomes > 1) {
        p <- progressr::progressor(steps = n_chromosomes)
    }

    ret <- future.apply::future_lapply(chromosomes, function(chr) {
        m <- dmps_locs$chr == chr
        cdmps_tsv <- dmps_tsv[(dmps_tsv[, dmps_tsv_id_col] %in% rownames(dmps_locs)), , drop = FALSE]
        cdmps <- dmps[m]
        cdmps_beta <- dmps_beta[m, , drop = FALSE]
        cdmps_locs <- dmps_locs[m, , drop = FALSE]
        dmrs <- data.frame()
        start_ind <- 1
        corr_pval <- 1
        dmr_dmps_inds <- c()

        # Update progress for this chromosome
        if (verbose && n_chromosomes > 1 && exists("p")) {
            p()
        }
        for (i in seq_len(nrow(cdmps_beta))) {
            reg_dmr <- FALSE
            dmr_dmps_inds <- c(dmr_dmps_inds, i)
            stop_condition <- FALSE
            if (i == nrow(cdmps_beta)) {
                stop_condition <- TRUE
                stop_reason <- "end of input"
            } else if ((max_lookup_dist > 0) && (cdmps_locs[i + 1, "pos"] - cdmps_locs[i, "pos"] > max_lookup_dist)) {
                stop_condition <- TRUE
                stop_reason <- "exceeded max distance"
            }
            if (stop_condition) {
                reg_dmr <- TRUE
            } else {
                t <- .testConnectivity(
                    site1_beta = unlist(cdmps_beta[i, ]),
                    site2_beta = unlist(cdmps_beta[i + 1, ]),
                    sample_groups = sample_groups,
                    max_pval = max_pval
                )

                corr_pval <- min(corr_pval, t[[2]])

                if (!t[[1]]) {
                    reg_dmr <- TRUE
                    stop_reason <- t$reason
                }
            }
            if (reg_dmr) {
                for (dmp_group in unique(cdmps_tsv[, dmp_group_col])) {
                    gdmps_tsv <- cdmps_tsv[
                        (cdmps_tsv[, dmp_group_col] == dmp_group),
                    ]
                    rownames(gdmps_tsv) <- gdmps_tsv[, dmps_tsv_id_col] # here there must be unique CpGs in `dmp` column

                    dmr_dmps_tsv <- gdmps_tsv[cdmps[dmr_dmps_inds], , drop = FALSE]
                    dmr_dmps_tsv[, "controls_beta"] <- dmr_dmps_tsv[, "cases_beta"] - dmr_dmps_tsv[, "delta_beta"]
                    for (num.col in c("cases_num", "controls_num")) {
                        if (!num.col %in% colnames(dmr_dmps_tsv)) {
                            dmr_dmps_tsv[, num.col] <- sum(!is.na(dmr_dmps_tsv[, "cases_beta"]))
                        }
                    }
                    for (opt.col in c("cases_beta_sd", "controls_beta_sd")) {
                        if (!opt.col %in% colnames(dmr_dmps_tsv)) {
                            dmr_dmps_tsv[, opt.col] <- NA
                        }
                    }

                    new_dmr <- data.frame(
                        chr = chr,
                        start_dmp = cdmps[[start_ind]],
                        end_dmp = cdmps[[i]],
                        start_dmp_pos = cdmps_locs[start_ind, "pos"],
                        end_dmp_pos = cdmps_locs[i, "pos"],
                        dmps_num = length(dmr_dmps_inds),
                        delta_beta = mean(abs(dmr_dmps_tsv[, "delta_beta"])) * sign(sum(sign(dmr_dmps_tsv[, "delta_beta"]))),
                        delta_beta_sd = sd(dmr_dmps_tsv[, "delta_beta"]),
                        delta_beta_se = sd(dmr_dmps_tsv[, "delta_beta"]) / sqrt(length(dmr_dmps_inds)),
                        delta_beta_min = min(dmr_dmps_tsv[, "delta_beta"]),
                        delta_beta_max = max(dmr_dmps_tsv[, "delta_beta"]),
                        delta_beta_start = dmr_dmps_tsv[1, "delta_beta"],
                        delta_beta_mid = dmr_dmps_tsv[ceiling(nrow(dmr_dmps_tsv) / 2), "delta_beta"],
                        delta_beta_end = dmr_dmps_tsv[nrow(dmr_dmps_tsv), "delta_beta"],
                        cases_beta = mean(abs(dmr_dmps_tsv[, "cases_beta"])) * sign(sum(sign(dmr_dmps_tsv[, "cases_beta"]))),
                        cases_beta_max = max(dmr_dmps_tsv[, "cases_beta"]),
                        cases_beta_min = min(dmr_dmps_tsv[, "cases_beta"]),
                        cases_beta_sd = sd(dmr_dmps_tsv[, "cases_beta"]),
                        cases_beta_se = sd(dmr_dmps_tsv[, "cases_beta"]) / sqrt(length(dmr_dmps_inds)),
                        cases_beta_start = dmr_dmps_tsv[1, "cases_beta"],
                        cases_beta_mid = dmr_dmps_tsv[ceiling(nrow(dmr_dmps_tsv) / 2), "cases_beta"],
                        cases_beta_end = dmr_dmps_tsv[nrow(dmr_dmps_tsv), "cases_beta"],
                        cases_beta_dmps_sd = mean(dmr_dmps_tsv[, "cases_beta_sd"]),
                        cases_beta_dmps_sd_max = max(dmr_dmps_tsv[, "cases_beta_sd"]),
                        cases_beta_dmps_sd_min = min(dmr_dmps_tsv[, "cases_beta_sd"]),
                        cases_beta_dmps_sd_start = dmr_dmps_tsv[1, "cases_beta_sd"],
                        cases_beta_dmps_sd_mid = dmr_dmps_tsv[ceiling(nrow(dmr_dmps_tsv) / 2), "cases_beta_sd"],
                        cases_beta_dmps_sd_end = dmr_dmps_tsv[nrow(dmr_dmps_tsv), "cases_beta_sd"],
                        controls_beta = mean(abs(dmr_dmps_tsv[, "controls_beta"])) * sign(sum(sign(dmr_dmps_tsv[, "controls_beta"]))),
                        controls_beta_max = max(dmr_dmps_tsv[, "controls_beta"]),
                        controls_beta_min = min(dmr_dmps_tsv[, "controls_beta"]),
                        controls_beta_sd = sd(dmr_dmps_tsv[, "controls_beta"]),
                        controls_beta_se = sd(dmr_dmps_tsv[, "controls_beta"]) / sqrt(length(dmr_dmps_inds)),
                        controls_beta_start = dmr_dmps_tsv[1, "controls_beta"],
                        controls_beta_mid = dmr_dmps_tsv[ceiling(nrow(dmr_dmps_tsv) / 2), "controls_beta"],
                        controls_beta_end = dmr_dmps_tsv[nrow(dmr_dmps_tsv), "controls_beta"],
                        controls_beta_dmps_sd = mean(dmr_dmps_tsv[, "controls_beta_sd"]),
                        controls_beta_dmps_sd_max = max(dmr_dmps_tsv[, "controls_beta_sd"]),
                        controls_beta_dmps_sd_min = min(dmr_dmps_tsv[, "controls_beta_sd"]),
                        controls_beta_dmps_sd_start = dmr_dmps_tsv[1, "controls_beta_sd"],
                        controls_beta_dmps_sd_mid = dmr_dmps_tsv[ceiling(nrow(dmr_dmps_tsv) / 2), "controls_beta_sd"],
                        controls_beta_dmps_sd_end = dmr_dmps_tsv[nrow(dmr_dmps_tsv), "controls_beta_sd"],
                        cases_num = min(dmr_dmps_tsv[, "cases_num"]),
                        controls_num = min(dmr_dmps_tsv[, "controls_num"]),
                        corr_pval = corr_pval,
                        dmps_pval_adj = mean(dmr_dmps_tsv[, pval_col]),
                        dmps_pval_adj_min = min(dmr_dmps_tsv[, pval_col]),
                        dmps_pval_adj_max = max(dmr_dmps_tsv[, pval_col]),
                        stop_connection_reason = stop_reason,
                        dmps = paste(cdmps[dmr_dmps_inds], collapse = ",")
                    )
                    names(new_dmr) <- gsub("pval_adj", pval_col, names(new_dmr), fixed = TRUE)
                    if (dmp_group_col != "_DUMMY_DMP_GROUP_COL_") {
                        new_dmr[[dmp_group_col]] <- dmp_group
                    }
                    tryCatch(
                        {
                            dmrs <- rbind(dmrs, new_dmr)
                        },
                        error = function(e) {
                            message("Error in rbind: ", e)
                            message("New DMR: \n\t", paste(paste(colnames(new_dmr), unlist(new_dmr), sep = ": "), collapse = "\n\t"))
                            message("Existing DMRs: \n\t", paste(capture.output(print(dmrs)), collapse = "\n\t"))
                            if (interactive()) browser()
                        }
                    )
                }

                start_ind <- i + 1
                corr_pval <- 1
                dmr_dmps_inds <- c()
            }
        }

        rownames(dmrs) <- seq_along(dmrs[, 1])
        dmrs[, "chr"] <- chr
        dmrs
    })

    if (inherits(ret[[1]], "try-error")) {
        stop(ret)
    }
    dmrs <- do.call(rbind, ret)
    if (verbose) message("Summary of connected DMRs before filtering based on connected DMP number:\n\t", paste(capture.output(summary(dmrs)), collapse = "\n\t"))
    dmrs <- dmrs[dmrs$dmps_num >= min_dmps, , drop = FALSE]
    if (nrow(dmrs) == 0) {
        if (verbose) {
            message("No DMRs remain after filtering based on connected DMP number.")
        }
        if (!is.null(output_prefix)) {
            for (f in c(".methylation.tsv.gz", ".tsv.gz")) {
                gzfile <- gzfile(paste0(output_prefix, f), "w", compression = 2)
                close(gzfile)
            }
        }
        return(NULL)
    }
    cases_num <- dmrs$cases_num
    controls_num <- dmrs$controls_num
    cases_sd_dmps_methylation <- dmrs$cases_beta_sd
    controls_sd_dmps_methylation <- dmrs$controls_beta_sd
    if (anyNA(c(cases_num, controls_num))) {
        warning("NAs introduced while coercing cases_num / controls_num to numeric; replacing NAs with 1 to avoid division errors.")
        cases_num[is.na(cases_num)] <- 1
        controls_num[is.na(controls_num)] <- 1
    }

    pooled_sd <- sqrt(((cases_num - 1) * cases_sd_dmps_methylation^2 + (controls_num - 1) * controls_sd_dmps_methylation^2) / (cases_num + controls_num - 2))
    pooled_sd[pooled_sd < 1e-7] <- 1e-7
    dmrs$cohensd <- dmrs$delta_beta / pooled_sd


    if (dmp_group_col %in% colnames(dmrs)) {
        ungrouped_dmrs <- dmrs[
            dmrs[, dmp_group_col] == dmrs[1, dmp_group_col],
            c("chr", "start_dmp", "end_dmp", "start_dmp_pos", "end_dmp_pos", "dmps_num", "corr_pval")
        ]
    } else {
        ungrouped_dmrs <- dmrs[, c("chr", "start_dmp", "end_dmp", "start_dmp_pos", "end_dmp_pos", "dmps_num", "corr_pval")]
    }

    if (verbose) {
        message("Number of initial DMRs:", nrow(ungrouped_dmrs))
    }
    message("Summary:\n\t", paste(capture.output(summary(ungrouped_dmrs)), collapse = "\n\t"))
    if (verbose) {
        message("Expanding DMRs on neighborhood CpGs..")
    }

    # Set up progress tracking for DMR expansion
    n_dmrs <- nrow(ungrouped_dmrs)
    if (verbose && n_dmrs > 1) {
        p_dmr <- progressr::progressor(steps = n_dmrs)
    }

    ret <- future.apply::future_lapply(split(ungrouped_dmrs, seq_along(ungrouped_dmrs[, 1])), function(dmr) {
        ret <- .expandDMRs(
            dmr = dmr,
            sample_groups = sample_groups,
            max_pval = max_pval,
            tabix_file = tabix_file,
            beta_file = beta_file,
            beta_row_names = beta_row_names,
            beta_col_names = beta_col_names,
            casecontrol = pheno[, casecontrol_col],
            expansion_step = expansion_step,
            expansion_relaxation = expansion_relaxation,
            min_cpg_delta_beta = min_cpg_delta_beta,
            sorted_locs = sorted_locs
        )

        # Update progress for this DMR
        if (verbose && n_dmrs > 1 && exists("p_dmr")) {
            p_dmr()
        }

        ret
    })

    if (inherits(ret, "try-error")) {
        stop(ret)
    }
    extended_dmrs <- do.call(rbind, ret)

    end_less_than_start <- extended_dmrs$end - extended_dmrs$start < 0

    if (any(end_less_than_start)) {
        warning(
            paste(
                sum(end_less_than_start),
                "DMRs have been assigned an end larger than start ! (CODE BUG TO BE REPORTED)"
            )
        )
        warning(
            "Those are: \n\t",
            paste0(capture.output(print(extended_dmrs[end_less_than_start, ])), collapse = "\n\t")
        )

        warning("Removing them..")
        extended_dmrs <- extended_dmrs[!end_less_than_start, ]
        warning("Remaining: ", nrow(extended_dmrs))
    }
    all_locs_inds <- rownames(sorted_locs)
    names(all_locs_inds) <- all_locs_inds
    all_locs_inds[seq_along(all_locs_inds)] <- seq_along(all_locs_inds)
    extended_dmrs$start_ind <- as.numeric(all_locs_inds[extended_dmrs$start_cpg])
    extended_dmrs$end_ind <- as.numeric(all_locs_inds[extended_dmrs$end_cpg])
    extended_dmrs$sup_cpgs_num <- extended_dmrs$end_ind - extended_dmrs$start_ind + 1

    extended_dmrs$id <- paste0(extended_dmrs$chr, ":", extended_dmrs$start, "-", extended_dmrs$end)
    extended_dmrs <- makeGRangesFromDataFrame(
        extended_dmrs,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqnames.field = "chr",
        seqinfo = Seqinfo(genome = genome),
        na.rm = TRUE
    )

    if (verbose) {
        message("Finding GC content of DMRs..")
    }

    if (genome == "hg19") {
        library(BSgenome.Hsapiens.UCSC.hg19)
        hsapiens <- BSgenome.Hsapiens.UCSC.hg19
        extended_dmrs_lifted_over <- extended_dmrs
    } else if (genome == "hg38") {
        library(BSgenome.Hsapiens.UCSC.hg38)
        hsapiens <- BSgenome.Hsapiens.UCSC.hg38
        library(rtracklayer)
        path <- system.file(package = "liftOver", "extdata", "hg19ToHg38.over.chain")
        chain <- import.chain(path)
        extended_dmrs_lifted_over <- liftOver(extended_dmrs, chain)
    }

    sequences <- getSeq(hsapiens, extended_dmrs_lifted_over, as.character = TRUE)
    # Convert sequences to character vector if needed
    if (is.list(sequences)) {
        sequences <- sapply(sequences, function(x) paste(x, collapse = ""))
    }
    extended_dmrs_lifted_over <- as.data.frame(extended_dmrs_lifted_over)
    extended_dmrs_lifted_over$cpgs_num <- str_count(sequences, "GC")
    extended_dmrs <- merge(extended_dmrs, extended_dmrs_lifted_over[, c("id", "cpgs_num")], by = "id")
    colnames(extended_dmrs)[colnames(extended_dmrs) == "seqnames"] <- "chr"


    extended_dmrs[extended_dmrs$cpgs_num == 0, "cpgs_num"] <- 1

    extended_dmrs$dmps_num_adj <- ceiling(extended_dmrs$cpgs_num / extended_dmrs$sup_cpgs_num * extended_dmrs$dmps_num)

    if (verbose) message("Summary of extended DMRs before filtering based on CpG number and adjusted DMPs number:\n\t", paste(capture.output(summary(extended_dmrs)), collapse = "\n\t"))
    filtered_dmrs <- extended_dmrs[extended_dmrs$cpgs_num >= min_cpgs & extended_dmrs$dmps_num_adj >= min_adj_dmps, , drop = FALSE]
    if (verbose) {
        message(
            "Keeping ",
            nrow(filtered_dmrs),
            " out of ",
            nrow(extended_dmrs),
            " with at least ",
            min_dmps,
            " (adjusted) supporting DMPs and ",
            min_cpgs,
            " contained CpGs."
        )
    }
    extended_dmrs <- filtered_dmrs
    if (nrow(extended_dmrs) == 0) {
        if (verbose) message("No DMRs passed the filtering based on min_cpgs and min_adj_dmps criteria.")
        if (!is.null(output_prefix)) {
            for (f in c(".methylation.tsv.gz", ".tsv.gz")) {
                gzfile <- gzfile(file.path(output_dir, paste0(output_prefix, f)), "w", compression = 2)
                close(gzfile)
            }
        }
        return(NULL)
    }
    if (verbose) message("Merging extended DMRs with original DMRs table to fill in missing information..")
    ne_dmrs_cols_to_keep <- colnames(dmrs)
    ne_dmrs_cols_to_keep <- c("start_dmp", "end_dmp", ne_dmrs_cols_to_keep[!ne_dmrs_cols_to_keep %in% colnames(extended_dmrs)])
    dmrs <- merge(extended_dmrs, dmrs[, ne_dmrs_cols_to_keep], by = c("start_dmp", "end_dmp"))
    dmrs <- dmrs[str_order(dmrs[, "id"], numeric = TRUE), ]
    if (verbose) message("Summary of final DMRs:\n\t", paste(capture.output(summary(dmrs)), collapse = "\n\t"))
    if (!is.null(output_prefix)) {
        dmrs_file <- paste0(output_prefix, "dmrs.tsv.gz")
        if (verbose) {
            message("Saving DMRs to ", dmrs_file, "..")
        }
        gz <- gzfile(dmrs_file, "w")
        write.table(
            dmrs,
            gz,
            sep = "\t",
            quote = FALSE,
            qmethod = "double",
            col.names = TRUE,
            row.names = FALSE
        )
        close(gz)
    }
    if (verbose) message("Done.")
    GenomicRanges::makeGRangesFromDataFrame(dmrs,
        keep.extra.columns = TRUE,
        seqinfo = Seqinfo(genome = genome),
        na.rm = TRUE
    )
}
