.readSamplesheet <- function(samplesheet.file,
                             samplesheet.file.sep,
                             sample_group.col,
                             sample_group.case,
                             sample_group.control,
                             target.col,
                             subset = NULL,
                             max.samples = NULL,
                             max.case.samples = NULL,
                             max.control.samples = NULL,
                             verbose = TRUE) {
    require(stringr)
    require(data.table)

    ref <- read.table(
        samplesheet.file,
        header = TRUE,
        sep = samplesheet.file.sep,
        quote = "",
        comment.char = "",
        check.names = FALSE
    )
    if (is.null(sample_group.col)) {
        return(ref)
    }
    if (is.null(sample_group.case) && (is.null(sample_group.control))) {
        stop("Both sample_group.case and sample_group.control provided are NULL, cannot continue.")
    }
    if (!(sample_group.col %in% colnames(ref))) {
        stop(
            "Sample group column ",
            sample_group.col,
            " not found in samplesheet columns: ",
            paste0(colnames(ref), collapse = ",")
        )
    }
    if (!(target.col %in% colnames(ref))) {
        stop(
            "Target column ",
            target.col,
            " not found in samplesheet columns: ",
            paste0(colnames(ref), collapse = ",")
        )
    }
    if (!is.null(subset)) {
        stopifnot(any(ref[, target.col] %in% subset))
        ref <- ref[ref[, target.col] %in% subset, ]
    }
    if (!is.null(sample_group.control)) {
        control.samplesheet <- ref[ref[, sample_group.col] %in% sample_group.control, c(target.col, sample_group.col, cov.cols),
            drop =
                FALSE
        ]
        if (nrow(control.samplesheet) == 0) {
            stop(
                paste0(
                    "The control label ",
                    sample_group.control,
                    " was not found in the provided samplesheet sample_group.col ",
                    sample_group.col,
                    " which has the following values: ",
                    paste(unique(ref[, sample_group.col]), collapse = ",")
                )
            )
        }
        if (!is.null(max.control.samples) &&
            max.control.samples < nrow(control.samplesheet)) {
            inds <- sort(sample(seq_len(nrow(control.samplesheet)), max.control.samples, replace = FALSE))
            control.samplesheet <- control.samplesheet[inds, ]
        }
    }

    if (!is.null(sample_group.case)) {
        case.samplesheet <- ref[ref[, sample_group.col] %in% sample_group.case, c(target.col, sample_group.col, cov.cols),
            drop =
                FALSE
        ]

        if (nrow(case.samplesheet) == 0) {
            stop(
                paste0(
                    "The case label ",
                    sample_group.case,
                    " was not found in the provided samplesheet sample_group.col ",
                    sample_group.col,
                    " which has the following values: ",
                    paste(unique(ref[, sample_group.col]), collapse = ",")
                )
            )
        }
        if (!is.null(max.case.samples) &&
            max.case.samples < nrow(case.samplesheet)) {
            inds <- sort(sample(seq_len(nrow(case.samplesheet)), max.case.samples, replace = FALSE))
            case.samplesheet <- case.samplesheet[inds, ]
        }
    }

    if (!is.null(sample_group.case) ||
        !is.null(sample_group.control)) {
        if (!is.null(sample_group.case) && !is.null(sample_group.control)) {
            subset.samplesheet <- rbind(control.samplesheet, case.samplesheet)
        } else {
            if (is.null(sample_group.case)) {
                subset.samplesheet <- control.samplesheet
            } else {
                subset.samplesheet <- case.samplesheet
            }
        }
    } else {
        subset.samplesheet <- ref
    }

    if (!is.null(max.samples) &&
        max.samples < nrow(subset.samplesheet)) {
        inds <- sort(sample(seq_len(nrow(subset.samplesheet)), max.samples, replace = FALSE))
        subset.samplesheet <- subset.samplesheet[inds, ]
    }

    rownames(subset.samplesheet) <- subset.samplesheet[, target.col]
    subset.samplesheet <- subset.samplesheet[, sample_group.col,
        drop =
            FALSE
    ]
    subset.samplesheet <- subset.samplesheet[
        str_order(rownames(subset.samplesheet),
            numeric = TRUE
        ), ,
        drop = FALSE
    ]
    stopifnot(ncol(subset.samplesheet) != 0)
    stopifnot(nrow(subset.samplesheet) != 0)
    if (verbose) {
        message(
            "Read samplesheet head:\n\t",
            paste(capture.output(print(
                head(subset.samplesheet)
            )), collapse = "\n\t")
        )
    }
    if (!is.null(sample_group.col)) {
        if (!is.null(sample_group.control)) {
            subset.samplesheet[, "casecontrol"] <- 1
            if (length(sample_group.control) == 1) {
                sample_group.control <- strsplit(sample_group.control, ",")[[1]]
            }
            subset.samplesheet[subset.samplesheet[, sample_group.col] %in% sample_group.control, "casecontrol"] <- 0
        } else {
            subset.samplesheet[, "casecontrol"] <- 0
            if (length(sample_group.case) == 1) {
                sample_group.case <- strsplit(sample_group.case, ",")[[1]]
            }
            subset.samplesheet[subset.samplesheet[, sample_group.col] %in% sample_group.case, "casecontrol"] <- 1
        }
    }
    subset.samplesheet
}

.processSamplesheet <- function(args,
                                subset = NULL) {
    suppressWarnings(suppressMessages({
        require(stringr)
    }))
    samplesheet.file <- args$samplesheet
    target.col <- args$target_col
    sample_group.col <- args$sample_group_col
    sample_group.control <- args$sample_group_control
    sample_group.case <- args$sample_group_case
    max.missing.ratio <- args$max_missing_cov_ratio
    if (is.null(max.missing.ratio)) {
        max.missing.ratio <- 0.3
    }
    if (!is.null(sample_group.case)) {
        sample_group.case <- strsplit(sample_group.case, ",")[[1]]
    }
    if (!is.null(sample_group.control)) {
        sample_group.control <- strsplit(sample_group.control, ",")[[1]]
    }
    max.samples <- args$max_samples
    max.case.samples <- args$max_case_samples
    max.control.samples <- args$max_control_samples

    if (endsWith(samplesheet.file, ".csv")) {
        samplesheet.file.sep <- ","
    } else if (endsWith(samplesheet.file, ".tsv")) {
        samplesheet.file.sep <- "\t"
    }
    subset.samplesheet <- .readSamplesheet(
        samplesheet.file = samplesheet.file,
        samplesheet.file.sep = samplesheet.file.sep,
        sample_group.col = sample_group.col,
        sample_group.control = sample_group.control,
        sample_group.case = sample_group.case,
        target.col = target.col,
        subset = subset,
        max.samples = max.samples,
        max.case.samples = max.case.samples,
        max.control.samples = max.control.samples,
    )
    ret <- list(samplesheet = subset.samplesheet[, c(sample_group.col, "casecontrol")])
    ret
}
