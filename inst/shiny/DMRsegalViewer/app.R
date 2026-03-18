# DMRsegal Viewer - Standalone Shiny Application
# This file allows running the viewer as a standalone app
# Usage: shiny::runApp(system.file("shiny/DMRsegalViewer", package = "DMRsegal"))

library(DMRsegal)
library(shiny)
library(bslib)
library(DT)

output_prefix <- Sys.getenv("DMRSEGAL_OUTPUT_PREFIX", unset = "")
pheno_file <- Sys.getenv("DMRSEGAL_PHENO_FILE", unset = "")
genome <- Sys.getenv("DMRSEGAL_GENOME", unset = "hg38")
array <- Sys.getenv("DMRSEGAL_ARRAY", unset = "450K")
sample_group_col <- Sys.getenv("DMRSEGAL_SAMPLE_GROUP_COL", unset = "Sample_Group")

if (output_prefix == "" || pheno_file == "") {
    ui <- fluidPage(
        theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
        br(),
        div(
            class = "container",
            div(
                class = "alert alert-info",
                h4("DMRsegal Viewer Configuration"),
                p("Please set the following environment variables before running this app:"),
                tags$ul(
                    tags$li(tags$code("DMRSEGAL_OUTPUT_PREFIX"), " - Path prefix for output files (required)"),
                    tags$li(tags$code("DMRSEGAL_PHENO_FILE"), " - Path to phenotype file (required)"),
                    tags$li(tags$code("DMRSEGAL_GENOME"), " - Genome assembly (default: hg38)"),
                    tags$li(tags$code("DMRSEGAL_ARRAY"), " - Array type: 450K, EPIC, EPICv2 (default: 450K)"),
                    tags$li(tags$code("DMRSEGAL_SAMPLE_GROUP_COL"), " - Sample group column name (default: Sample_Group)")
                ),
                hr(),
                p("Alternatively, use the ", tags$code("launchDMRsegalViewer()"), " function from R:"),
                pre(
                    "launchDMRsegalViewer(\n",
                    "  output_prefix = 'path/to/output',\n",
                    "  pheno = pheno_df,\n",
                    "  genome = 'hg38',\n",
                    "  array = 'EPIC'\n",
                    ")"
                )
            )
        )
    )
    server <- function(input, output, session) {}
} else {
    pheno <- DMRsegal:::.readSamplesheet(pheno_file)
    data <- DMRsegal:::.loadDMRsegalData(output_prefix, pheno, genome, array, sample_group_col)
    ui <- DMRsegal:::.createViewerUI()
    server <- DMRsegal:::.createViewerServer(data)
}

shinyApp(ui = ui, server = server)
