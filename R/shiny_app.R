#' Launch CMEnt Interactive Viewer
#'
#' Launches a Shiny application for interactive exploration of DMR analysis
#' results from \code{\link{findDMRsFromSeeds}}.
#'
#' @param output_prefix Character. Prefix used when saving DMR analysis results.
#'   The function will look for files: \code{{output_prefix}.dmrs.tsv.gz},
#'   \code{{output_prefix}.seeds_beta.tsv.gz}, and \code{{output_prefix}.meta.rds}.
#' @param launch_browser Logical. Whether to launch in browser (default: TRUE).
#' @param port Integer. Port number for Shiny server (default: auto-assigned).
#' @param host Character. Host interface for the Shiny server (default: \code{"127.0.0.1"}).
#'   Use \code{"0.0.0.0"} when exposing the app from a Docker container.
#' @param diagnostic Logical. Whether to enable diagnostic features for block formation visualization (default: FALSE).
#'
#' @return Invisibly returns the Shiny app object.
#'
#' @details
#' The viewer provides interactive access to all visualization functions in the package:
#' \itemize{
#'   \item \strong{Overview}: Summary statistics and filterable DMR table
#'   \item \strong{Single DMR}: Detailed view of individual DMRs with heatmap and motif logo
#'   \item \strong{Manhattan}: Genome-wide scatter plot of DMR scores
#'   \item \strong{Circos}: Circular genome plot with motif interactions
#'   \item \strong{Block Formation}: Diagnostic view of DMR block detection, if `diagnostic = TRUE`
#' }
#'
#' The function expects the following output files from \code{findDMRsFromSeeds}:
#' \itemize{
#'   \item \code{{output_prefix}.dmrs.tsv.gz} - Main DMR results (required)
#'   \item \code{{output_prefix}.seeds_beta.tsv.gz} - Beta values for seeds (required)
#'   \item \code{{output_prefix}.meta.rds} - Viewer metadata with phenotype, array, and genome information (required)
#'   \item \code{{output_prefix}dmr_interactions.tsv} - DMR interactions (optional)
#'   \item \code{{output_prefix}dmr_components.tsv} - Motif components (optional)
#' }
#'
#' @examples
#' \donttest{
#' # After running findDMRsFromSeeds with output_prefix = "my_analysis"
#' launchCMEntViewer(
#'   output_prefix = "results/my_analysis"
#' )
#' }
#'
#' @seealso \code{\link{findDMRsFromSeeds}}, \code{\link{plotDMR}}, \code{\link{plotDMRsCircos}}
#'
#' @export
launchCMEntViewer <- function(
    output_prefix,
    launch_browser = TRUE,
    port = NULL,
    host = "127.0.0.1",
    diagnostic = FALSE
) {
    validation <- .validateOutputPrefix(output_prefix)
    if (!validation$valid) {
        stop(paste(validation$errors, collapse = "\n"))
    }
    if (length(validation$warnings) > 0) {
        for (w in validation$warnings) {
            warning(w)
        }
    }

    data <- .loadCMEntData(output_prefix)

    ui <- .createViewerUI(diagnostic = diagnostic)
    server <- .createViewerServer(data)

    app <- shiny::shinyApp(ui = ui, server = server)
    run_args <- c(
        list(appDir = app, launch.browser = launch_browser, host = host),
        if (is.null(port)) list() else list(port = port)
    )

    invisible(do.call(shiny::runApp, run_args))
}

.viewer_worker_env <- new.env(parent = emptyenv())

.viewerOutputFiles <- function(output_prefix) {
    list(
        dmrs_file = paste0(output_prefix, ".dmrs.tsv.gz"),
        beta_file = paste0(output_prefix, ".seeds_beta.tsv.gz"),
        meta_file = paste0(output_prefix, ".meta.rds"),
        interactions_file = paste0(output_prefix, ".dmr_interactions.tsv"),
        components_file = paste0(output_prefix, ".dmr_components.tsv")
    )
}

.readViewerTabularCache <- function(path, label) {
    if (!file.exists(path)) {
        return(NULL)
    }

    file_size <- file.info(path)$size
    if (isTRUE(is.finite(file_size) && file_size <= 1)) {
        .log_info("Loaded empty ", label, " cache from ", path, ".", level = 1)
        return(data.frame())
    }

    tryCatch(
        data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE),
        error = function(e) {
            .log_warn("Failed to load ", label, ": ", e$message, level = 1)
            NULL
        }
    )
}

.loadViewerCircosCache <- function(output_prefix) {
    files <- .viewerOutputFiles(output_prefix)

    interactions <- NULL
    if (file.exists(files$interactions_file)) {
        .log_step("Loading interactions from ", files$interactions_file, "...", level = 1)
        interactions <- .readViewerTabularCache(files$interactions_file, "interactions")
    }

    components <- NULL
    if (file.exists(files$components_file)) {
        .log_step("Loading components from ", files$components_file, "...", level = 1)
        components <- .readViewerTabularCache(files$components_file, "components")
    }

    list(
        interactions = interactions,
        components = components
    )
}

.reloadViewerCircosCache <- function(data) {
    stopifnot(is.list(data))
    if (is.null(data$output_prefix) || !nzchar(data$output_prefix)) {
        return(data)
    }

    circos_cache <- .loadViewerCircosCache(data$output_prefix)
    data$interactions <- circos_cache$interactions
    data$components <- circos_cache$components
    data
}


.validateOutputPrefix <- function(output_prefix) {
    errors <- character()
    warnings <- character()
    files <- .viewerOutputFiles(output_prefix)

    if (!file.exists(files$dmrs_file)) {
        errors <- c(errors, paste("DMRs file not found:", files$dmrs_file))
    }
    if (!file.exists(files$beta_file)) {
        errors <- c(errors, paste("Beta values file not found:", files$beta_file))
    }
    if (!file.exists(files$meta_file)) {
        errors <- c(
            errors,
            paste(
                "Viewer metadata file not found:",
                files$meta_file,
                "- re-run findDMRsFromSeeds() with output_prefix to create it"
            )
        )
    }
    if (!file.exists(files$interactions_file)) {
        warnings <- c(
            warnings,
            paste(
                "Interactions file not found:",
                files$interactions_file,
                "- the Circos panel can precompute this cache on demand"
            )
        )
    }
    if (!file.exists(files$components_file)) {
        warnings <- c(
            warnings,
            paste(
                "Components file not found:",
                files$components_file,
                "- the Circos panel can precompute this cache on demand"
            )
        )
    }

    list(
        valid = length(errors) == 0,
        errors = errors,
        warnings = warnings
    )
}

.resolveViewerSampleGroupCol <- function(pheno, sample_group_col = NULL) {
    stopifnot(is.data.frame(pheno))
    candidate_cols <- setdiff(colnames(pheno), "__casecontrol__")
    if (!is.null(sample_group_col) && sample_group_col %in% candidate_cols) {
        return(sample_group_col)
    }
    if ("Sample_Group" %in% candidate_cols) {
        return("Sample_Group")
    }
    if (length(candidate_cols) > 0) {
        return(candidate_cols[[1]])
    }

    stop("Viewer metadata phenotype table does not contain a usable sample group column.")
}

.loadCMEntViewerMeta <- function(output_prefix) {
    meta_file <- .viewerOutputFiles(output_prefix)$meta_file

    .log_step("Loading viewer metadata from ", meta_file, "...", level = 1)
    meta <- tryCatch(
        readRDS(meta_file),
        error = function(e) {
            stop("Failed to load viewer metadata file ", meta_file, ": ", e$message, call. = FALSE)
        }
    )

    if (!is.list(meta)) {
        stop("Viewer metadata file must contain a list.")
    }
    if (!("pheno" %in% names(meta)) || !is.data.frame(meta$pheno)) {
        stop("Viewer metadata file must contain a data.frame `pheno` entry.")
    }
    if (is.null(rownames(meta$pheno)) || any(!nzchar(rownames(meta$pheno)))) {
        stop("Viewer metadata phenotype table must have non-empty row names.")
    }

    genome <- .normalizeFindDMRsGenome(meta$genome)
    if (is.null(genome)) {
        stop("Viewer metadata file must contain a valid `genome` entry.")
    }

    array <- .normalizeFindDMRsArray(meta$array)
    sample_group_col <- .resolveViewerSampleGroupCol(meta$pheno, meta$sample_group_col)
    if (!identical(meta$sample_group_col, sample_group_col)) {
        .log_info(
            "Viewer metadata sample_group_col resolved to '",
            sample_group_col,
            "'.",
            level = 2
        )
    }

    .log_success("Loaded viewer metadata for ", nrow(meta$pheno), " samples", level = 1)

    list(
        pheno = meta$pheno,
        genome = genome,
        array = array,
        sample_group_col = sample_group_col
    )
}

.loadCMEntData <- function(output_prefix) {
    files <- .viewerOutputFiles(output_prefix)
    meta <- .loadCMEntViewerMeta(output_prefix)
    pheno <- meta$pheno
    genome <- meta$genome
    array <- meta$array
    sample_group_col <- meta$sample_group_col

    .log_step("Loading DMRs from ", files$dmrs_file, "...", level = 1)
    dmrs_df <- data.table::fread(files$dmrs_file, sep = "\t", header = TRUE, data.table = FALSE)
    dmrs <- convertToGRanges(dmrs_df, genome = genome)
    .log_success("Loaded ", length(dmrs), " DMRs", level = 1)

    .log_step("Loading beta values from ", files$beta_file, "...", level = 1)
    beta <- data.table::fread(files$beta_file, sep = "\t", header = TRUE, data.table = FALSE)
    rownames(beta) <- beta[[1]]
    beta <- beta[, -1, drop = FALSE]
    beta <- as.matrix(beta)
    .log_success("Loaded beta values for ", nrow(beta), " sites and ", ncol(beta), " samples", level = 1)

    beta_handler <- getBetaHandler(
        beta = beta,
        array = array,
        genome = genome
    )

    circos_cache <- .loadViewerCircosCache(output_prefix)

    list(
        dmrs = dmrs,
        beta = beta,
        beta_handler = beta_handler,
        pheno = pheno,
        interactions = circos_cache$interactions,
        components = circos_cache$components,
        genome = genome,
        array = array,
        sample_group_col = sample_group_col,
        output_prefix = output_prefix
    )
}

.viewerAssetPrefix <- "cment-viewer-assets"

.registerViewerAssets <- function() {
    asset_dir <- system.file("shiny/CMEntViewer/www", package = "CMEnt")
    if (!nzchar(asset_dir) || !dir.exists(asset_dir)) {
        return(NULL)
    }

    current_paths <- shiny::resourcePaths()
    if (!(.viewerAssetPrefix %in% names(current_paths))) {
        shiny::addResourcePath(.viewerAssetPrefix, asset_dir)
    }

    asset_dir
}

.viewerAssetHref <- function(asset_name) {
    paste0(.viewerAssetPrefix, "/", asset_name)
}

.viewerWithSpinner <- function(ui, proxy_height = NULL) {
    shinycssloaders::withSpinner(
        ui,
        type = 4,
        color = "#2C3E50",
        size = 0.9,
        hide.ui = FALSE,
        proxy.height = proxy_height
    )
}

.viewerPageSpinnerCaption <- function(
    message = "Processing...",
    detail = NULL,
    cancelable = FALSE,
    cancel_input_id = "viewer_cancel_task"
) {
    detail <- trimws(if (is.null(detail)) "" else as.character(detail))
    shiny::tags$div(
        class = "text-center",
        shiny::tags$div(class = "fw-semibold", message),
        if (nzchar(detail)) shiny::tags$div(class = "small text-muted mt-1", detail),
        if (isTRUE(cancelable)) {
            shiny::tags$div(
                class = "mt-3",
                shiny::actionButton(
                    cancel_input_id,
                    "Cancel",
                    icon = shiny::icon("ban"),
                    class = "btn-outline-secondary btn-sm"
                )
            )
        }
    )
}

.viewerShowPageSpinner <- function(
    message = "Processing...",
    detail = NULL,
    cancelable = FALSE,
    cancel_input_id = "viewer_cancel_task"
) {
    shinycssloaders::showPageSpinner(
        type = 4,
        color = "#2C3E50",
        size = 0.9,
        caption = .viewerPageSpinnerCaption(
            message = message,
            detail = detail,
            cancelable = cancelable,
            cancel_input_id = cancel_input_id
        )
    )
}

.viewerHidePageSpinner <- function() {
    shinycssloaders::hidePageSpinner()
}

.viewerBrandUI <- function(asset_dir) {
    github_link <- shiny::tags$a(
        shiny::icon("github"),
        shiny::tags$span("GitHub"),
        href = "https://github.com/CMG-UA/CMEnt",
        target = "_blank",
        class = "cment-viewer-github"
    )

    logo <- if (!is.null(asset_dir) && file.exists(file.path(asset_dir, "logo.png"))) {
        shiny::tags$img(
            src = .viewerAssetHref("logo.png"),
            alt = "CMEnt logo",
            class = "cment-viewer-logo"
        )
    } else {
        shiny::tags$span("CMEnt", class = "cment-viewer-wordmark")
    }

    shiny::tags$div(
        class = "cment-viewer-brand",
        github_link,
        logo
    )
}

.viewerTaskMessage <- function(task_type) {
    switch(task_type,
        single_dmr_plot = list(
            message = "Generating DMR plot...",
            detail = ""
        ),
        manhattan_plot = list(
            message = "Generating Manhattan plot...",
            detail = ""
        ),
        block_plot = list(
            message = "Generating block formation plot...",
            detail = ""
        ),
        circos_plot = list(
            message = "Generating Circos plot...",
            detail = ""
        ),
        circos_cache_compute = list(
            message = "Computing interactions...",
            detail = ""
        ),
        list(
            message = "Processing...",
            detail = ""
        )
    )
}

.viewerDevPackagePath <- function() {
    ns_path <- tryCatch(getNamespaceInfo(asNamespace("CMEnt"), "path"), error = function(e) NA_character_)
    cwd_path <- tryCatch(normalizePath(getwd(), winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
    candidate_paths <- unique(c(ns_path, cwd_path))
    candidate_paths <- candidate_paths[!is.na(candidate_paths) & nzchar(candidate_paths)]

    for (path in candidate_paths) {
        desc_file <- file.path(path, "DESCRIPTION")
        if (!file.exists(desc_file)) {
            next
        }

        pkg_name <- tryCatch(unname(as.character(read.dcf(desc_file, fields = "Package")[1, 1])), error = function(e) NA_character_)
        if (!identical(pkg_name, "CMEnt")) {
            next
        }

        if (
            file.exists(file.path(path, "R", "shiny_app.R")) &&
            file.exists(file.path(path, "R", "shiny_modules.R"))
        ) {
            return(path)
        }
    }

    NULL
}

.captureViewerRecordedPlot <- function(plot_fun) {
    stopifnot(is.function(plot_fun))

    grDevices::pdf(file = NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot_fun()
    grDevices::recordPlot()
}

.viewerRunBackgroundTaskFromData <- function(task_type, data, params) {
    switch(task_type,
        single_dmr_plot = list(
            task_type = task_type,
            plot = .captureViewerRecordedPlot(function() {
                plotDMR(
                    dmrs = data$dmrs,
                    dmr_index = params$dmr_index,
                    beta = data$beta_handler,
                    pheno = data$pheno,
                    genome = data$genome,
                    array = data$array,
                    sample_group_col = data$sample_group_col,
                    max_sites = params$max_sites,
                    max_samples_per_group = params$max_samples_per_group,
                    plot_title = TRUE
                )
            })
        ),
        manhattan_plot = list(
            task_type = task_type,
            plot = .captureViewerRecordedPlot(function() {
                plotDMRsManhattan(
                    dmrs = data$dmrs,
                    region = params$region,
                    genome = data$genome,
                    point_size = params$point_size,
                    point_alpha = params$point_alpha
                )
            })
        ),
        block_plot = list(
            task_type = task_type,
            plot = .captureViewerRecordedPlot(function() {
                plotDMRBlockFormation(
                    dmrs = data$dmrs,
                    chromosome = params$chromosome,
                    genome = data$genome,
                    k_neighbors = params$k_neighbors,
                    block_gap_mode = params$block_gap_mode,
                    block_gap_fixed_bp = if (identical(params$block_gap_mode, "fixed")) params$block_gap_fixed_bp else NULL
                )
            })
        ),
        circos_plot = {
            prepared_plot_state <- NULL
            plot <- .captureViewerRecordedPlot(function() {
                prepared_plot_state <<- .runViewerCircosPlot(
                    params,
                    data,
                    prepared_plot_state = params$prepared_plot_state
                )
            })
            list(
                task_type = task_type,
                plot = plot,
                prepared_plot_state = prepared_plot_state
            )
        },
        circos_cache_compute = {
            result <- computeDMRsInteraction(
                dmrs = data$dmrs,
                genome = data$genome,
                array = data$array,
                beta_locs = data$beta_handler$getBetaLocs(),
                min_similarity = params$min_similarity,
                output_prefix = data$output_prefix
            )
            cache <- .loadViewerCircosCache(data$output_prefix)
            if (is.null(cache$interactions)) {
                cache$interactions <- result$interactions
            }
            if (is.null(cache$components)) {
                cache$components <- .serializeDMRInteractionComponentsForStorage(result$components)
            }

            list(
                task_type = task_type,
                cache = cache
            )
        },
        stop("Unknown viewer background task: ", task_type, call. = FALSE)
    )
}

.viewerRunBackgroundTask <- function(task_type, output_prefix, params) {
    data <- .loadCMEntData(output_prefix)
    .viewerRunBackgroundTaskFromData(
        task_type = task_type,
        data = data,
        params = params
    )
}

.launchViewerBackgroundTask <- function(task_type, output_prefix, params, dev_pkg_path = NULL) {
    callr::r_bg(
        func = function(task_type, output_prefix, params, dev_pkg_path) {
            if (
                !is.null(dev_pkg_path) &&
                nzchar(dev_pkg_path) &&
                file.exists(file.path(dev_pkg_path, "DESCRIPTION")) &&
                requireNamespace("pkgload", quietly = TRUE)
            ) {
                pkgload::load_all(
                    dev_pkg_path,
                    quiet = TRUE,
                    export_all = FALSE,
                    helpers = FALSE,
                    attach_testthat = FALSE
                )
            }

            .viewerRunBackgroundTask(
                task_type = task_type,
                output_prefix = output_prefix,
                params = params
            )
        },
        args = list(
            task_type = task_type,
            output_prefix = output_prefix,
            params = params,
            dev_pkg_path = dev_pkg_path
        ),
        supervise = TRUE,
        # Avoid pipe-buffer deadlocks from verbose worker logging.
        stdout = NULL,
        stderr = NULL
    )
}

.createViewerWorkerSession <- function(output_prefix, dev_pkg_path = NULL) {
    worker <- callr::r_session$new(wait = TRUE)
    bootstrap_error <- tryCatch(
        {
            worker$run(
                func = function(output_prefix, dev_pkg_path) {
                    if (
                        !is.null(dev_pkg_path) &&
                        nzchar(dev_pkg_path) &&
                        file.exists(file.path(dev_pkg_path, "DESCRIPTION")) &&
                        requireNamespace("pkgload", quietly = TRUE)
                    ) {
                        pkgload::load_all(
                            dev_pkg_path,
                            quiet = TRUE,
                            export_all = FALSE,
                            helpers = FALSE,
                            attach_testthat = FALSE
                        )
                    }

                    assign(
                        ".cment_viewer_worker_data",
                        suppressMessages(
                            .loadCMEntData(output_prefix)
                        ),
                        envir = .viewer_worker_env
                    )
                    TRUE
                },
                args = list(
                    output_prefix = output_prefix,
                    dev_pkg_path = dev_pkg_path
                )
            )
            NULL
        },
        error = function(e) e
    )

    if (inherits(bootstrap_error, "error")) {
        try(worker$close(), silent = TRUE)
        stop(bootstrap_error)
    }

    worker
}

.readViewerWorkerTask <- function(worker) {
    task_error <- NULL
    task_result <- NULL
    close_worker <- FALSE
    task_done <- FALSE

    repeat {
        next_event <- tryCatch(worker$read(), error = function(e) e)
        if (inherits(next_event, "error")) {
            task_error <- next_event
            close_worker <- TRUE
            task_done <- TRUE
            break
        }
        if (is.null(next_event)) {
            break
        }

        if (identical(next_event$code, worker$status$MSG)) {
            next
        }

        if (identical(next_event$code, worker$status$DONE)) {
            task_result <- next_event$result
            task_error <- next_event$error
            task_done <- TRUE
            break
        }

        if (next_event$code %in% c(worker$status$EXITED, worker$status$CRASHED, worker$status$CLOSED)) {
            task_error <- simpleError("Viewer background worker stopped unexpectedly.")
            close_worker <- TRUE
            task_done <- TRUE
            break
        }
    }

    if (!isTRUE(task_done)) {
        worker_state <- tryCatch(worker$get_state(), error = function(e) "finished")
        worker_alive <- isTRUE(tryCatch(worker$is_alive(), error = function(e) FALSE))

        if (worker_state %in% c("starting", "busy") && worker_alive) {
            return(list(
                done = FALSE,
                result = NULL,
                error = NULL,
                close_worker = FALSE
            ))
        }

        if (!worker_alive || worker_state %in% c("finished", "closed")) {
            return(list(
                done = TRUE,
                result = NULL,
                error = simpleError("Viewer background worker stopped unexpectedly."),
                close_worker = TRUE
            ))
        }

        if (identical(worker_state, "idle")) {
            return(list(
                done = TRUE,
                result = NULL,
                error = simpleError("Viewer background task finished without returning a result."),
                close_worker = FALSE
            ))
        }
    }

    list(
        done = TRUE,
        result = task_result,
        error = task_error,
        close_worker = close_worker
    )
}

.createViewerTaskController <- function(session, data) {
    active_task <- shiny::reactiveVal(NULL)
    task_state <- shiny::reactiveVal(list(
        active = FALSE,
        task_type = NULL,
        message = NULL,
        detail = NULL,
        cancelable = FALSE
    ))
    worker_session <- NULL
    worker_dev_pkg_path <- .viewerDevPackagePath()
    task_counter <- 0L

    set_task_state <- function(
        active = FALSE,
        task_type = NULL,
        message = NULL,
        detail = NULL,
        cancelable = FALSE
    ) {
        task_state(list(
            active = active,
            task_type = task_type,
            message = message,
            detail = detail,
            cancelable = cancelable
        ))
    }

    get_active_task <- function() {
        shiny::isolate(active_task())
    }

    close_worker <- function() {
        if (!is.null(worker_session)) {
            try(worker_session$close(), silent = TRUE)
            worker_session <<- NULL
        }
        invisible(TRUE)
    }

    ensure_worker <- function() {
        if (!is.null(worker_session)) {
            worker_state <- tryCatch(worker_session$get_state(), error = function(e) "finished")
            if (
                identical(worker_state, "idle") &&
                isTRUE(tryCatch(worker_session$is_alive(), error = function(e) FALSE))
            ) {
                return(worker_session)
            }
            close_worker()
        }

        worker_session <<- .createViewerWorkerSession(
            output_prefix = data$output_prefix,
            dev_pkg_path = worker_dev_pkg_path
        )
        worker_session
    }

    cancel <- function() {
        current <- get_active_task()
        if (is.null(current)) {
            return(FALSE)
        }

        callback <- current$on_cancel
        set_task_state(
            active = TRUE,
            task_type = current$task_type,
            message = current$message,
            detail = "Stopping background process...",
            cancelable = FALSE
        )
        close_worker()
        active_task(NULL)
        set_task_state()
        if (is.function(callback)) {
            callback()
        }
        TRUE
    }

    shutdown <- function() {
        active_task(NULL)
        set_task_state()
        close_worker()
        TRUE
    }

    start <- function(
        task_type,
        params,
        on_success = NULL,
        on_error = NULL,
        on_cancel = NULL
    ) {
        current <- get_active_task()
        if (!is.null(current)) {
            shiny::showNotification(
                "A viewer task is already running. Cancel it before starting a new one.",
                type = "warning",
                duration = 5
            )
            return(FALSE)
        }

        descriptor <- .viewerTaskMessage(task_type)
        task_counter <<- task_counter + 1L
        task_id <- task_counter
        pending_detail <- "Starting background worker..."

        active_task(list(
            backend = "pending",
            worker = NULL,
            task_id = task_id,
            task_type = task_type,
            message = descriptor$message,
            detail = pending_detail,
            cancelable = TRUE,
            on_success = on_success,
            on_error = on_error,
            on_cancel = on_cancel
        ))
        set_task_state(
            active = TRUE,
            task_type = task_type,
            message = descriptor$message,
            detail = pending_detail,
            cancelable = TRUE
        )

        session$onFlushed(function() {
            current <- get_active_task()
            if (is.null(current) || !identical(current$task_id, task_id)) {
                return()
            }

            worker <- tryCatch(
                ensure_worker(),
                error = function(e) e
            )

            if (inherits(worker, "error")) {
                active_task(NULL)
                set_task_state()
                if (is.function(on_error)) {
                    on_error(worker)
                } else {
                    shiny::showNotification(
                        paste0("Unable to start viewer task: ", conditionMessage(worker)),
                        type = "error",
                        duration = NULL
                    )
                }
                return()
            }

            launch_error <- tryCatch(
                {
                    worker$call(
                        func = function(task_type, params) {
                            suppressMessages({
                                if (!exists(".cment_viewer_worker_data", envir = .viewer_worker_env, inherits = FALSE)) {
                                    stop("Viewer worker data is not initialized.", call. = FALSE)
                                }

                                data <- get(".cment_viewer_worker_data", envir = .viewer_worker_env, inherits = FALSE)
                                result <- .viewerRunBackgroundTaskFromData(
                                    task_type = task_type,
                                    data = data,
                                    params = params
                                )

                                if (identical(task_type, "circos_cache_compute") && is.list(result$cache)) {
                                    data$interactions <- result$cache$interactions
                                    data$components <- result$cache$components
                                    assign(".cment_viewer_worker_data", data, envir = .viewer_worker_env)
                                }

                                result
                            })
                        },
                        args = list(
                            task_type = task_type,
                            params = params
                        )
                    )
                    NULL
                },
                error = function(e) e
            )

            if (inherits(launch_error, "error")) {
                close_worker()
                active_task(NULL)
                set_task_state()
                if (is.function(on_error)) {
                    on_error(launch_error)
                } else {
                    shiny::showNotification(
                        paste0("Unable to start viewer task: ", conditionMessage(launch_error)),
                        type = "error",
                        duration = NULL
                    )
                }
                return()
            }

            current <- get_active_task()
            if (is.null(current) || !identical(current$task_id, task_id)) {
                close_worker()
                return()
            }

            current$backend <- "session"
            current$worker <- worker
            current$detail <- descriptor$detail
            active_task(current)
            set_task_state(
                active = TRUE,
                task_type = task_type,
                message = descriptor$message,
                detail = descriptor$detail,
                cancelable = TRUE
            )
        }, once = TRUE)

        TRUE
    }

    shiny::observe({
        current <- active_task()
        if (is.null(current)) {
            return()
        }

        shiny::invalidateLater(250, session)
        worker <- current$worker
        if (is.null(worker)) {
            return()
        }
        task_status <- .readViewerWorkerTask(worker)
        if (!isTRUE(task_status$done)) {
            return()
        }

        if (isTRUE(task_status$close_worker)) {
            close_worker()
        }

        active_task(NULL)
        set_task_state()

        if (!is.null(task_status$error)) {
            if (is.function(current$on_error)) {
                current$on_error(task_status$error)
            } else {
                shiny::showNotification(
                    paste0("Viewer task failed: ", task_status$error$message),
                    type = "error",
                    duration = NULL
                )
            }
            return()
        }

        if (is.function(current$on_success)) {
            current$on_success(task_status$result)
        }
    })

    list(
        initialize = function() {
            invisible(try(ensure_worker(), silent = TRUE))
        },
        start = start,
        cancel = cancel,
        shutdown = shutdown,
        state = shiny::reactive({
            task_state()
        })
    )
}

.createViewerUI <- function(diagnostic = FALSE) {
    asset_dir <- .registerViewerAssets()
    asset_head <- if (!is.null(asset_dir)) {
        shiny::tags$head(
            shiny::tags$link(
                rel = "stylesheet",
                type = "text/css",
                href = .viewerAssetHref("custom.css")
            )
        )
    }
    panels_list <- list(
        bslib::nav_panel(
            title = "Overview",
            icon = shiny::icon("table"),
            .overviewUI("overview")
        ),
        bslib::nav_panel(
            title = "Single DMR",
            icon = shiny::icon("dna"),
            .plotDMRUI("single_dmr")
        ),
        bslib::nav_panel(
            title = "Manhattan",
            icon = shiny::icon("chart-line"),
            .manhattanUI("manhattan")
        ),
        bslib::nav_panel(
            title = "Circos",
            icon = shiny::icon("circle-notch"),
            .circosUI("circos")
        )
    )
    if (diagnostic) {
        panels_list <- c(panels_list, list(
            bslib::nav_panel(
                title = "Block Formation",
                icon = shiny::icon("layer-group"),
                .blockFormationUI("block_formation")
            )
        ))
    }
    navbar <- do.call(bslib::page_navbar, c(
        list(
            title = "CMEnt Viewer",
            id = "navbar",
            theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
            fillable = TRUE
        ),
        panels_list,
        list(
            bslib::nav_spacer(),
            bslib::nav_item(.viewerBrandUI(asset_dir))
        )
    ))

    shiny::tagList(
        asset_head,
        navbar
    )
}

.createViewerServer <- function(data) {
    function(input, output, session) {
        task_controller <- .createViewerTaskController(session, data)
        page_spinner_visible <- shiny::reactiveVal(FALSE)
        page_spinner_state_key <- shiny::reactiveVal(NULL)

        shiny::observeEvent(input$viewer_cancel_task, {
            task_controller$cancel()
        }, ignoreInit = TRUE)

        shiny::observe({
            state <- task_controller$state()
            spinner_visible <- page_spinner_visible()
            spinner_state_key <- paste(
                isTRUE(state$active),
                if (is.null(state$task_type)) "" else state$task_type,
                if (is.null(state$message)) "" else state$message,
                if (is.null(state$detail)) "" else state$detail,
                isTRUE(state$cancelable),
                sep = "\r"
            )

            if (
                isTRUE(state$active) &&
                (!spinner_visible || !identical(spinner_state_key, page_spinner_state_key()))
            ) {
                .viewerShowPageSpinner(
                    message = if (is.null(state$message)) "Processing..." else state$message,
                    detail = state$detail,
                    cancelable = isTRUE(state$cancelable),
                    cancel_input_id = "viewer_cancel_task"
                )
                page_spinner_visible(TRUE)
                page_spinner_state_key(spinner_state_key)
                return()
            }

            if (!isTRUE(state$active) && spinner_visible) {
                .viewerHidePageSpinner()
                page_spinner_visible(FALSE)
                page_spinner_state_key(NULL)
            }
        })

        selected_dmr_id <- .overviewServer("overview", data)
        .plotDMRServer("single_dmr", data, task_controller)
        .manhattanServer("manhattan", data, task_controller)
        .blockFormationServer("block_formation", data, task_controller)
        .circosServer("circos", data, task_controller)

        shiny::observe({
            dmr_id <- selected_dmr_id()
            if (!is.null(dmr_id) && dmr_id %in% data$dmrs$id) {
                dmr_index <- which(data$dmrs$id == dmr_id)
                if (length(dmr_index) == 1) {
                    shiny::updateNavbarPage(session, "navbar", selected = "Single DMR")
                    shiny::updateNumericInput(session, "single_dmr-dmr_index", value = dmr_index)
                }
            }
        })

        session$onSessionEnded(function() {
            if (isTRUE(shiny::isolate(page_spinner_visible()))) {
                try(.viewerHidePageSpinner(), silent = TRUE)
            }
            try(task_controller$shutdown(), silent = TRUE)
            .log_info("CMEnt Viewer session ended.", level = 1)
        })
    }
}
