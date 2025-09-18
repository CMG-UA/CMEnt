.readSamplesheet <- function(samplesheet_file,
                             samplesheet_file_sep,
                             sample_group_col,
                             sample_group_case,
                             sample_group_control,
                             target_col,
                             subset = NULL,
                             max_samples = NULL,
                             max_case_samples = NULL,
                             max_control_samples = NULL,
                             verbose = TRUE) {
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
        if (!is.null(max_control_samples) &&
            max_control_samples < nrow(control_samplesheet)) {
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
        if (!is.null(max_case_samples) &&
            max_case_samples < nrow(case_samplesheet)) {
            inds <- sort(sample(seq_len(nrow(case_samplesheet)), max_case_samples, replace = FALSE))
            case_samplesheet <- case_samplesheet[inds, ]
        }
    }

    if (!is.null(sample_group_case) ||
        !is.null(sample_group_control)) {
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

    if (!is.null(max_samples) &&
        max_samples < nrow(subset_samplesheet)) {
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
    if (verbose) {
        message(
            "Read samplesheet head:\n\t",
            paste(capture.output(print(
                head(subset_samplesheet)
            )), collapse = "\n\t")
        )
    }
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
    max_missing.ratio <- args$max_missing_cov_ratio
    if (is.null(max_missing.ratio)) {
        max_missing.ratio <- 0.3
    }
    if (!is.null(sample_group_case)) {
        sample_group_case <- strsplit(sample_group_case, ",")[[1]]
    }
    if (!is.null(sample_group_control)) {
        sample_group_control <- strsplit(sample_group_control, ",")[[1]]
    }
    max_samples <- args$max_samples
    max_case.samples <- args$max_case_samples
    max_control.samples <- args$max_control_samples

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
        max_case.samples = max_case.samples,
        max_control.samples = max_control.samples,
    )
    ret <- list(samplesheet = subset_samplesheet[, c(sample_group_col, "casecontrol")])
    ret
}



#' @export
getSortedGenomicLocs <- function(array) {
  if (array == "450K") {
    sorted_locs <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
  } else if (array == "EPIC") {
    sorted_locs <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
  } else {
    stop("Unknown array type: ", array)
  }
  sorted_locs <- sorted_locs[str_order(paste0(sorted_locs[, "chr"], ":", sorted_locs[, "pos"]), numeric = T), ]
  sorted_locs[,'start'] <- sorted_locs[,'pos']
  sorted_locs[,'end'] <- sorted_locs[,'pos'] + 1
  sorted_locs
}

#' @export
orderByLoc <- function(x, array="450K", genomic_locs=NULL) {
    if (is.null(genomic_locs)){
      genomic_locs <- getSortedGenomicLocs(array)
    }
    str_order(paste0(genomic_locs[x, "chr"], ":", genomic_locs[x, "pos"]), numeric = T)
  }