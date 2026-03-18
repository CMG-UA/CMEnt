#' Launch DMRsegal Interactive Viewer
#'
#' Launches a Shiny application for interactive exploration of DMR analysis
#' results from \code{\link{findDMRsFromSeeds}}.
#'
#' @param output_prefix Character. Prefix used when saving DMR analysis results.
#'   The function will look for files: \code{{output_prefix}.dmrs.tsv.gz},
#'   \code{{output_prefix}.seeds_beta.tsv.gz}, and \code{{output_prefix}_meta.rds}.
#' @param launch_browser Logical. Whether to launch in browser (default: TRUE).
#' @param port Integer. Port number for Shiny server (default: auto-assigned).
#'
#' @return Invisibly returns the Shiny app object.
#'
#' @details
#' The viewer provides interactive access to all visualization functions in the package:
#' \itemize{
#'   \item \strong{Overview}: Summary statistics and filterable DMR table
#'   \item \strong{Single DMR}: Detailed view of individual DMRs with heatmap and motif logo
#'   \item \strong{Multi DMR}: Grid view of multiple DMRs
#'   \item \strong{Manhattan}: Genome-wide scatter plot of DMR scores
#'   \item \strong{Block Formation}: Diagnostic view of DMR block detection
#'   \item \strong{Circos}: Circular genome plot with motif interactions
#' }
#'
#' The function expects the following output files from \code{findDMRsFromSeeds}:
#' \itemize{
#'   \item \code{{output_prefix}.dmrs.tsv.gz} - Main DMR results (required)
#'   \item \code{{output_prefix}.seeds_beta.tsv.gz} - Beta values for seeds (required)
#'   \item \code{{output_prefix}_meta.rds} - Viewer metadata with phenotype, array, and genome information (required)
#'   \item \code{{output_prefix}dmr_interactions.csv} - DMR interactions (optional)
#'   \item \code{{output_prefix}dmr_components.csv} - Motif components (optional)
#' }
#'
#' @examples
#' \dontrun{
#' # After running findDMRsFromSeeds with output_prefix = "my_analysis"
#' launchDMRsegalViewer(
#'   output_prefix = "results/my_analysis"
#' )
#' }
#'
#' @seealso \code{\link{findDMRsFromSeeds}}, \code{\link{plotDMR}}, \code{\link{plotDMRsCircos}}
#'
#' @export
launchDMRsegalViewer <- function(
    output_prefix,
    launch_browser = TRUE,
    port = NULL
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

    data <- .loadDMRsegalData(output_prefix)

    ui <- .createViewerUI()
    server <- .createViewerServer(data)

    app <- shiny::shinyApp(ui = ui, server = server)
    run_args <- c(
        list(appDir = app, launch.browser = launch_browser),
        if (is.null(port)) list() else list(port = port)
    )

    invisible(do.call(shiny::runApp, run_args))
}

.viewerOutputFiles <- function(output_prefix) {
    list(
        dmrs_file = paste0(output_prefix, ".dmrs.tsv.gz"),
        beta_file = paste0(output_prefix, ".seeds_beta.tsv.gz"),
        meta_file = paste0(output_prefix, "_meta.rds"),
        interactions_file = paste0(output_prefix, "dmr_interactions.csv"),
        components_file = paste0(output_prefix, "dmr_components.csv")
    )
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
                "- Circos interactions will be computed on-the-fly"
            )
        )
    }
    if (!file.exists(files$components_file)) {
        warnings <- c(
            warnings,
            paste(
                "Components file not found:",
                files$components_file,
                "- Circos components will be computed on-the-fly"
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

.loadDMRsegalViewerMeta <- function(output_prefix) {
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

.loadDMRsegalData <- function(output_prefix) {
    files <- .viewerOutputFiles(output_prefix)
    meta <- .loadDMRsegalViewerMeta(output_prefix)
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
    .log_success("Loaded beta values for ", nrow(beta), " CpGs and ", ncol(beta), " samples", level = 1)

    beta_handler <- getBetaHandler(
        beta = beta,
        array = array,
        genome = genome
    )

    interactions <- NULL
    if (file.exists(files$interactions_file)) {
        .log_step("Loading interactions from ", files$interactions_file, "...", level = 1)
        interactions <- tryCatch(
            data.table::fread(files$interactions_file, sep = "\t", header = TRUE, data.table = FALSE),
            error = function(e) {
                .log_warn("Failed to load interactions: ", e$message, level = 1)
                NULL
            }
        )
    }

    components <- NULL
    if (file.exists(files$components_file)) {
        .log_step("Loading components from ", files$components_file, "...", level = 1)
        components <- tryCatch(
            data.table::fread(files$components_file, sep = "\t", header = TRUE, data.table = FALSE),
            error = function(e) {
                .log_warn("Failed to load components: ", e$message, level = 1)
                NULL
            }
        )
    }

    list(
        dmrs = dmrs,
        beta = beta,
        beta_handler = beta_handler,
        pheno = pheno,
        interactions = interactions,
        components = components,
        genome = genome,
        array = array,
        sample_group_col = sample_group_col,
        output_prefix = output_prefix
    )
}

.viewerAssetPrefix <- "dmrsegal-viewer-assets"

.registerViewerAssets <- function() {
    asset_dir <- system.file("shiny/DMRsegalViewer/www", package = "DMRsegal")
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

.viewerBrandUI <- function(asset_dir) {
    github_link <- shiny::tags$a(
        shiny::icon("github"),
        shiny::tags$span("GitHub"),
        href = "https://github.com/CMG-UA/DMRsegal",
        target = "_blank",
        class = "dmrsegal-viewer-github"
    )

    logo <- if (!is.null(asset_dir) && file.exists(file.path(asset_dir, "logo.png"))) {
        shiny::tags$img(
            src = .viewerAssetHref("logo.png"),
            alt = "DMRsegal logo",
            class = "dmrsegal-viewer-logo"
        )
    } else {
        shiny::tags$span("DMRsegal", class = "dmrsegal-viewer-wordmark")
    }

    shiny::tags$div(
        class = "dmrsegal-viewer-brand",
        github_link,
        logo
    )
}

.createViewerUI <- function() {
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

    shiny::tagList(
        asset_head,
        bslib::page_navbar(
            title = "DMRsegal Viewer",
            id = "navbar",
            theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
            fillable = TRUE,
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
                title = "Multi DMR",
                icon = shiny::icon("grip"),
                .plotDMRsUI("multi_dmr")
            ),
            bslib::nav_panel(
                title = "Manhattan",
                icon = shiny::icon("chart-line"),
                .manhattanUI("manhattan")
            ),
            bslib::nav_panel(
                title = "Block Formation",
                icon = shiny::icon("layer-group"),
                .blockFormationUI("block_formation")
            ),
            bslib::nav_panel(
                title = "Circos",
                icon = shiny::icon("circle-notch"),
                .circosUI("circos")
            ),
            bslib::nav_spacer(),
            bslib::nav_item(.viewerBrandUI(asset_dir))
        )
    )
}

.createViewerServer <- function(data) {
    function(input, output, session) {
        selected_dmr_id <- .overviewServer("overview", data)
        .plotDMRServer("single_dmr", data)
        .plotDMRsServer("multi_dmr", data)
        .manhattanServer("manhattan", data)
        .blockFormationServer("block_formation", data)
        .circosServer("circos", data)

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
            .log_info("DMRsegal Viewer session ended.", level = 1)
        })
    }
}
