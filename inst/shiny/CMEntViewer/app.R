# CMEnt Viewer - Standalone Shiny Application
# This file allows running the viewer as a standalone app
# Usage: shiny::runApp(system.file("shiny/CMEntViewer", package = "CMEnt"))

suppressPackageStartupMessages({
    library(CMEnt)
    library(shiny)
    library(bslib)
    library(DT)
})

output_prefix <- Sys.getenv("CMENT_OUTPUT_PREFIX", unset = "")

if (output_prefix == "") {
    ui <- fluidPage(
        theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
        br(),
        div(
            class = "container",
            div(
                class = "alert alert-info",
                h4("CMEnt Viewer Configuration"),
                p("Please set the following environment variable before running this app:"),
                tags$ul(
                    tags$li(tags$code("CMENT_OUTPUT_PREFIX"), " - Path prefix for output files created by findDMRsFromSeeds()")
                ),
                hr(),
                p("The viewer loads phenotype, array, and genome information from ", tags$code(".meta.rds"), " automatically."),
                p("Alternatively, use the ", tags$code("launchCMEntViewer()"), " function from R:"),
                pre(
                    "launchCMEntViewer(\n",
                    "  output_prefix = 'path/to/output'\n",
                    ")"
                )
            )
        )
    )
    server <- function(input, output, session) {}
} else {
    validation <- CMEnt:::.validateOutputPrefix(output_prefix)
    if (!validation$valid) {
        ui <- fluidPage(
            theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
            br(),
            div(
                class = "container",
                div(
                    class = "alert alert-danger",
                    h4("CMEnt Viewer Configuration Error"),
                    p("The configured output prefix is missing required viewer files:"),
                    tags$ul(lapply(validation$errors, tags$li)),
                    hr(),
                    p("Re-run ", tags$code("findDMRsFromSeeds()"), " with ", tags$code("output_prefix"), " to regenerate the viewer outputs.")
                )
            )
        )
        server <- function(input, output, session) {}
    } else {
        if (length(validation$warnings) > 0) {
            for (w in validation$warnings) {
                warning(w)
            }
        }
        data <- CMEnt:::.loadCMEntData(output_prefix)
        ui <- CMEnt:::.createViewerUI()
        server <- CMEnt:::.createViewerServer(data)
    }
}

shinyApp(ui = ui, server = server)
