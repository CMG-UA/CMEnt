#' Launch DMRsegal Interactive Viewer
#'
#' Launches a Shiny application for interactive exploration of DMR analysis
#' results from \code{\link{findDMRsFromSeeds}}.
#'
#' @param output_prefix Character. Prefix used when saving DMR analysis results.
#'   The function will look for files: \code{{output_prefix}.dmrs.tsv.gz},
#'   \code{{output_prefix}.seeds_beta.tsv.gz}, and optionally interaction/component files.
#' @param pheno Data frame or path to phenotype file. Required for heatmap visualizations.
#' @param genome Character. Genome assembly (default: "hg38").
#' @param array Character. Array platform: "450K", "27K", "EPIC", or "EPICv2" (default: "450K").
#' @param sample_group_col Character. Column name in pheno for sample grouping (default: "Sample_Group").
#' @param launch_browser Logical. Whether to launch in browser (default: TRUE).
#' @param port Integer. Port number for Shiny server (default: auto-assigned).
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
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
#'   \item \code{{output_prefix}dmr_interactions.csv} - DMR interactions (optional)
#'   \item \code{{output_prefix}dmr_components.csv} - Motif components (optional)
#' }
#'
#' @examples
#' \dontrun{
#' # After running findDMRsFromSeeds with output_prefix = "my_analysis"
#' launchDMRsegalViewer(
#'   output_prefix = "results/my_analysis",
#'   pheno = pheno_df,
#'   genome = "hg38",
#'   array = "EPIC"
#' )
#' }
#'
#' @seealso \code{\link{findDMRsFromSeeds}}, \code{\link{plotDMR}}, \code{\link{plotDMRsCircos}}
#'
#' @export
launchDMRsegalViewer <- function(
    output_prefix,
    pheno = NULL,
    genome = "hg38",
    array = c("450K", "27K", "EPIC", "EPICv2"),
    sample_group_col = "Sample_Group",
    launch_browser = TRUE,
    port = NULL,
    ...
) {
    array <- strex::match_arg(array, ignore_case = TRUE)
    validation <- .validateOutputPrefix(output_prefix)
    if (!validation$valid) {
        stop(paste(validation$errors, collapse = "\n"))
    }
    if (length(validation$warnings) > 0) {
        for (w in validation$warnings) {
            warning(w)
        }
    }

    if (is.null(pheno)) {
        stop("pheno is required for the viewer. Please provide a data.frame or path to the phenotype file.")
    }
    if (is.character(pheno) && length(pheno) == 1) {
        if (!file.exists(pheno)) {
            stop("Phenotype file not found: ", pheno)
        }
        pheno <- .readSamplesheet(pheno)
    }

    data <- .loadDMRsegalData(output_prefix, pheno, genome, array, sample_group_col)

    ui <- .createViewerUI()
    server <- .createViewerServer(data)

    app <- shiny::shinyApp(ui = ui, server = server)

    if (!is.null(port)) {
        shiny::runApp(app, launch.browser = launch_browser, port = port, ...)
    } else {
        shiny::runApp(app, launch.browser = launch_browser, ...)
    }
}

.validateOutputPrefix <- function(output_prefix) {
    errors <- character()
    warnings <- character()

    dmrs_file <- paste0(output_prefix, ".dmrs.tsv.gz")
    beta_file <- paste0(output_prefix, ".seeds_beta.tsv.gz")
    interactions_file <- paste0(output_prefix, "dmr_interactions.csv")
    components_file <- paste0(output_prefix, "dmr_components.csv")

    if (!file.exists(dmrs_file)) {
        errors <- c(errors, paste("DMRs file not found:", dmrs_file))
    }
    if (!file.exists(beta_file)) {
        errors <- c(errors, paste("Beta values file not found:", beta_file))
    }
    if (!file.exists(interactions_file)) {
        warnings <- c(warnings, paste("Interactions file not found:", interactions_file, "- Circos interactions will be computed on-the-fly"))
    }
    if (!file.exists(components_file)) {
        warnings <- c(warnings, paste("Components file not found:", components_file, "- Circos components will be computed on-the-fly"))
    }

    list(
        valid = length(errors) == 0,
        errors = errors,
        warnings = warnings
    )
}

.loadDMRsegalData <- function(output_prefix, pheno, genome, array, sample_group_col) {
    dmrs_file <- paste0(output_prefix, ".dmrs.tsv.gz")
    beta_file <- paste0(output_prefix, ".seeds_beta.tsv.gz")
    interactions_file <- paste0(output_prefix, "dmr_interactions.csv")
    components_file <- paste0(output_prefix, "dmr_components.csv")

    .log_step("Loading DMRs from ", dmrs_file, "...", level = 1)
    dmrs_df <- data.table::fread(dmrs_file, sep = "\t", header = TRUE, data.table = FALSE)
    dmrs <- convertToGRanges(dmrs_df, genome = genome)
    .log_success("Loaded ", length(dmrs), " DMRs", level = 1)

    .log_step("Loading beta values from ", beta_file, "...", level = 1)
    beta <- data.table::fread(beta_file, sep = "\t", header = TRUE, data.table = FALSE)
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
    if (file.exists(interactions_file)) {
        .log_step("Loading interactions from ", interactions_file, "...", level = 1)
        interactions <- tryCatch({
            data.table::fread(interactions_file, sep = "\t", header = TRUE, data.table = FALSE)
        }, error = function(e) {
            .log_warn("Failed to load interactions: ", e$message, level = 1)
            NULL
        })
    }

    components <- NULL
    if (file.exists(components_file)) {
        .log_step("Loading components from ", components_file, "...", level = 1)
        components <- tryCatch({
            data.table::fread(components_file, sep = "\t", header = TRUE, data.table = FALSE)
        }, error = function(e) {
            .log_warn("Failed to load components: ", e$message, level = 1)
            NULL
        })
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

.createViewerUI <- function() {
    bslib::page_navbar(
        title = "DMRsegal Viewer",
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
        bslib::nav_item(
            shiny::tags$a(
                shiny::icon("github"),
                href = "https://github.com/CMG-UA/DMRsegal",
                target = "_blank",
                class = "nav-link"
            )
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
