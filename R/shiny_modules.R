#' @importFrom bslib card card_header card_body value_box page_navbar nav_panel
#' @importFrom bslib nav_spacer nav_item bs_theme
#' @importFrom DT DTOutput renderDT datatable
NULL

.overviewUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(3, bslib::value_box(
                title = "Total DMRs",
                value = shiny::textOutput(ns("n_dmrs")),
                showcase = shiny::icon("dna"),
                theme = "primary"
            )),
            shiny::column(3, bslib::value_box(
                title = "Total Samples",
                value = shiny::textOutput(ns("n_samples")),
                showcase = shiny::icon("users"),
                theme = "secondary"
            )),
            shiny::column(3, bslib::value_box(
                title = "Genome",
                value = shiny::textOutput(ns("genome_info")),
                showcase = shiny::icon("globe"),
                theme = "info"
            )),
            shiny::column(3, bslib::value_box(
                title = "Array",
                value = shiny::textOutput(ns("array_info")),
                showcase = shiny::icon("microchip"),
                theme = "success"
            ))
        ),
        shiny::br(),
        bslib::card(
            bslib::card_header("DMR Summary Table"),
            bslib::card_body(
                shiny::fluidRow(
                    shiny::column(4, shiny::selectInput(
                        ns("filter_chr"), "Filter by Chromosome",
                        choices = c("All" = ""), multiple = FALSE
                    )),
                    shiny::column(4, shiny::sliderInput(
                        ns("filter_delta_beta"), "Min |Delta Beta|",
                        min = 0, max = 1, value = 0, step = 0.01
                    )),
                    shiny::column(4, shiny::numericInput(
                        ns("filter_min_cpgs"), "Min CpGs",
                        value = 1, min = 1, step = 1
                    ))
                ),
                DT::DTOutput(ns("dmr_table"))
            )
        )
    )
}

.overviewServer <- function(id, data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$n_dmrs <- shiny::renderText({
            shiny::req(data$dmrs)
            as.character(length(data$dmrs))
        })

        output$n_samples <- shiny::renderText({
            shiny::req(data$beta)
            as.character(ncol(data$beta) - 1)
        })

        output$genome_info <- shiny::renderText({
            data$genome
        })

        output$array_info <- shiny::renderText({
            data$array
        })

        shiny::observe({
            shiny::req(data$dmrs)
            chrs <- unique(as.character(GenomicRanges::seqnames(data$dmrs)))
            chrs <- chrs[order(as.numeric(gsub("chr", "", chrs, ignore.case = TRUE)))]
            shiny::updateSelectInput(session, "filter_chr",
                choices = c("All" = "", chrs)
            )
        })

        filtered_dmrs <- shiny::reactive({
            shiny::req(data$dmrs)
            dmrs_df <- convertToDataFrame(data$dmrs)
            if (!is.null(input$filter_chr) && input$filter_chr != "") {
                dmrs_df <- dmrs_df[dmrs_df$chr == input$filter_chr, , drop = FALSE]
            }
            if (!is.null(input$filter_delta_beta) && input$filter_delta_beta > 0) {
                dmrs_df <- dmrs_df[abs(dmrs_df$delta_beta) >= input$filter_delta_beta, , drop = FALSE]
            }
            if (!is.null(input$filter_min_cpgs) && input$filter_min_cpgs > 1) {
                dmrs_df <- dmrs_df[dmrs_df$cpgs_num >= input$filter_min_cpgs, , drop = FALSE]
            }
            dmrs_df
        })

        output$dmr_table <- DT::renderDT({
            shiny::req(filtered_dmrs())
            dmrs_df <- filtered_dmrs()
            display_cols <- intersect(
                c("id", "chr", "start", "end", "width", "cpgs_num", "seeds_num",
                  "delta_beta", "cases_beta", "controls_beta", "rank", "block_id",
                  "in_promoter_of", "in_gene_body_of"),
                colnames(dmrs_df)
            )
            DT::datatable(
                dmrs_df[, display_cols, drop = FALSE],
                options = list(
                    pageLength = 15,
                    scrollX = TRUE,
                    dom = "frtip"
                ),
                selection = "single",
                rownames = FALSE
            )
        })

        return(shiny::reactive({
            sel <- input$dmr_table_rows_selected
            if (is.null(sel) || length(sel) == 0) return(NULL)
            filtered_dmrs()$id[sel]
        }))
    })
}

.plotDMRUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(3,
                bslib::card(
                    bslib::card_header("DMR Selection"),
                    bslib::card_body(
                        shiny::numericInput(ns("dmr_index"), "DMR Index", value = 1, min = 1, step = 1),
                        shiny::hr(),
                        shiny::h6("Plot Options"),
                        shiny::numericInput(ns("max_cpgs"), "Max CpGs", value = 100, min = 10, max = 500, step = 10),
                        shiny::numericInput(ns("max_samples_per_group"), "Max Samples/Group", value = 10, min = 1, max = 50, step = 1),
                        shiny::checkboxInput(ns("plot_motif"), "Show Motif Logo", value = TRUE),
                        shiny::numericInput(ns("motif_cpg_flank_size"), "Motif Flank Size", value = 5, min = 1, max = 20, step = 1),
                        shiny::hr(),
                        shiny::actionButton(ns("generate_plot"), "Generate Plot", class = "btn-primary w-100"),
                        shiny::br(), shiny::br(),
                        shiny::downloadButton(ns("download_plot"), "Download PDF", class = "w-100")
                    )
                )
            ),
            shiny::column(9,
                bslib::card(
                    bslib::card_header("Single DMR Visualization"),
                    bslib::card_body(
                        shiny::conditionalPanel(
                            condition = "output.has_plot == false",
                            ns = ns,
                            shiny::div(
                                class = "text-center text-muted py-5",
                                shiny::icon("chart-line", class = "fa-3x mb-3"),
                                shiny::p("Click 'Generate Plot' to visualize a DMR")
                            )
                        ),
                        shiny::plotOutput(ns("dmr_plot"), height = "700px")
                    )
                )
            )
        )
    )
}

.plotDMRServer <- function(id, data) {
    shiny::moduleServer(id, function(input, output, session) {
        shiny::observe({
            shiny::req(data$dmrs)
            shiny::updateNumericInput(session, "dmr_index",
                max = length(data$dmrs),
                value = min(1, length(data$dmrs))
            )
        })

        plot_data <- shiny::eventReactive(input$generate_plot, {
            shiny::req(data$dmrs, input$dmr_index)
            shiny::validate(
                shiny::need(input$dmr_index >= 1 && input$dmr_index <= length(data$dmrs),
                    "Invalid DMR index")
            )
            list(
                dmr_index = input$dmr_index,
                max_cpgs = input$max_cpgs,
                max_samples_per_group = input$max_samples_per_group,
                plot_motif = input$plot_motif,
                motif_cpg_flank_size = input$motif_cpg_flank_size
            )
        })

        output$has_plot <- shiny::reactive({
            !is.null(plot_data())
        })
        shiny::outputOptions(output, "has_plot", suspendWhenHidden = FALSE)

        output$dmr_plot <- shiny::renderPlot({
            shiny::req(plot_data())
            params <- plot_data()
            shiny::withProgress(message = "Generating DMR plot...", value = 0.5, {
                plotDMR(
                    dmrs = data$dmrs,
                    dmr_index = params$dmr_index,
                    beta = data$beta_handler,
                    pheno = data$pheno,
                    genome = data$genome,
                    array = data$array,
                    sample_group_col = data$sample_group_col,
                    max_cpgs = params$max_cpgs,
                    max_samples_per_group = params$max_samples_per_group,
                    plot_motif = params$plot_motif,
                    motif_cpg_flank_size = params$motif_cpg_flank_size,
                    plot_title = TRUE
                )
            })
        }, res = 100)

        output$download_plot <- shiny::downloadHandler(
            filename = function() {
                paste0("DMR_", input$dmr_index, "_", Sys.Date(), ".pdf")
            },
            content = function(file) {
                params <- plot_data()
                shiny::req(params)
                plotDMR(
                    dmrs = data$dmrs,
                    dmr_index = params$dmr_index,
                    beta = data$beta_handler,
                    pheno = data$pheno,
                    genome = data$genome,
                    array = data$array,
                    sample_group_col = data$sample_group_col,
                    max_cpgs = params$max_cpgs,
                    max_samples_per_group = params$max_samples_per_group,
                    plot_motif = params$plot_motif,
                    motif_cpg_flank_size = params$motif_cpg_flank_size,
                    output_file = file
                )
            }
        )
    })
}

.plotDMRsUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(3,
                bslib::card(
                    bslib::card_header("Grid Settings"),
                    bslib::card_body(
                        shiny::textInput(ns("dmr_indices"), "DMR Indices (comma-separated)", placeholder = "e.g., 1,2,3,4"),
                        shiny::numericInput(ns("top_n"), "Top N DMRs (if indices empty)", value = 4, min = 1, max = 20, step = 1),
                        shiny::numericInput(ns("ncol"), "Number of Columns", value = 2, min = 1, max = 4, step = 1),
                        shiny::hr(),
                        shiny::actionButton(ns("generate_plot"), "Generate Grid", class = "btn-primary w-100"),
                        shiny::br(), shiny::br(),
                        shiny::downloadButton(ns("download_plot"), "Download PDF", class = "w-100")
                    )
                )
            ),
            shiny::column(9,
                bslib::card(
                    bslib::card_header("Multi-DMR Grid"),
                    bslib::card_body(
                        shiny::plotOutput(ns("dmrs_plot"), height = "800px")
                    )
                )
            )
        )
    )
}

.plotDMRsServer <- function(id, data) {
    shiny::moduleServer(id, function(input, output, session) {
        parsed_indices <- shiny::reactive({
            idx_str <- trimws(input$dmr_indices)
            if (is.null(idx_str) || idx_str == "") return(NULL)
            idx <- tryCatch({
                as.integer(strsplit(idx_str, ",")[[1]])
            }, error = function(e) NULL)
            if (is.null(idx) || any(is.na(idx))) return(NULL)
            idx[idx >= 1 & idx <= length(data$dmrs)]
        })

        plot_data <- shiny::eventReactive(input$generate_plot, {
            shiny::req(data$dmrs)
            idx <- parsed_indices()
            list(
                dmr_indices = idx,
                top_n = input$top_n,
                ncol = input$ncol
            )
        })

        output$dmrs_plot <- shiny::renderPlot({
            shiny::req(plot_data())
            params <- plot_data()
            shiny::withProgress(message = "Generating DMR grid...", value = 0.5, {
                plotDMRs(
                    dmrs = data$dmrs,
                    dmr_indices = params$dmr_indices,
                    top_n = params$top_n,
                    beta = data$beta_handler,
                    pheno = data$pheno,
                    sample_group_col = data$sample_group_col,
                    genome = data$genome,
                    array = data$array,
                    ncol = params$ncol
                )
            })
        }, res = 100)

        output$download_plot <- shiny::downloadHandler(
            filename = function() {
                paste0("DMRs_grid_", Sys.Date(), ".pdf")
            },
            content = function(file) {
                params <- plot_data()
                shiny::req(params)
                plotDMRs(
                    dmrs = data$dmrs,
                    dmr_indices = params$dmr_indices,
                    top_n = params$top_n,
                    beta = data$beta_handler,
                    pheno = data$pheno,
                    sample_group_col = data$sample_group_col,
                    genome = data$genome,
                    array = data$array,
                    ncol = params$ncol,
                    output_file = file,
                    width = 10,
                    height = 14
                )
            }
        )
    })
}

.manhattanUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(3,
                bslib::card(
                    bslib::card_header("Manhattan Plot Settings"),
                    bslib::card_body(
                        shiny::selectInput(ns("score_col"), "Score Column", choices = c("score", "delta_beta")),
                        shiny::checkboxInput(ns("show_blocks"), "Show Blocks", value = TRUE),
                        shiny::sliderInput(ns("point_size"), "Point Size", min = 0.5, max = 3, value = 1.1, step = 0.1),
                        shiny::sliderInput(ns("point_alpha"), "Point Transparency", min = 0.1, max = 1, value = 0.75, step = 0.05),
                        shiny::sliderInput(ns("block_alpha"), "Block Transparency", min = 0.05, max = 0.5, value = 0.12, step = 0.01),
                        shiny::hr(),
                        shiny::actionButton(ns("generate_plot"), "Generate Plot", class = "btn-primary w-100"),
                        shiny::br(), shiny::br(),
                        shiny::downloadButton(ns("download_plot"), "Download PDF", class = "w-100")
                    )
                )
            ),
            shiny::column(9,
                bslib::card(
                    bslib::card_header("Manhattan Plot"),
                    bslib::card_body(
                        shiny::plotOutput(ns("manhattan_plot"), height = "500px")
                    )
                )
            )
        )
    )
}

.manhattanServer <- function(id, data) {
    shiny::moduleServer(id, function(input, output, session) {
        shiny::observe({
            shiny::req(data$dmrs)
            cols <- colnames(S4Vectors::mcols(data$dmrs))
            score_cols <- intersect(c("score", "rank", "delta_beta"), cols)
            shiny::updateSelectInput(session, "score_col", choices = score_cols)
        })

        plot_data <- shiny::eventReactive(input$generate_plot, {
            shiny::req(data$dmrs)
            list(
                score_col = input$score_col,
                show_blocks = input$show_blocks,
                point_size = input$point_size,
                point_alpha = input$point_alpha,
                block_alpha = input$block_alpha
            )
        })

        output$manhattan_plot <- shiny::renderPlot({
            shiny::req(plot_data())
            params <- plot_data()
            shiny::withProgress(message = "Generating Manhattan plot...", value = 0.5, {
                plotDMRsManhattan(
                    dmrs = data$dmrs,
                    score_col = params$score_col,
                    genome = data$genome,
                    show_blocks = params$show_blocks,
                    point_size = params$point_size,
                    point_alpha = params$point_alpha,
                    block_alpha = params$block_alpha
                )
            })
        }, res = 100)

        output$download_plot <- shiny::downloadHandler(
            filename = function() {
                paste0("Manhattan_", Sys.Date(), ".pdf")
            },
            content = function(file) {
                params <- plot_data()
                shiny::req(params)
                plotDMRsManhattan(
                    dmrs = data$dmrs,
                    score_col = params$score_col,
                    genome = data$genome,
                    show_blocks = params$show_blocks,
                    point_size = params$point_size,
                    point_alpha = params$point_alpha,
                    block_alpha = params$block_alpha,
                    output_file = file,
                    width = 14,
                    height = 6
                )
            }
        )
    })
}

.blockFormationUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(3,
                bslib::card(
                    bslib::card_header("Block Formation Settings"),
                    bslib::card_body(
                        shiny::selectInput(ns("chromosome"), "Chromosome", choices = c()),
                        shiny::selectInput(ns("score_col"), "Score Column", choices = c("score")),
                        shiny::hr(),
                        shiny::h6("Block Gap Mode"),
                        shiny::radioButtons(ns("block_gap_mode"), NULL,
                            choices = c("Adaptive" = "adaptive", "Fixed" = "fixed", "None" = "none"),
                            selected = "adaptive", inline = TRUE
                        ),
                        shiny::conditionalPanel(
                            condition = "input.block_gap_mode == 'fixed'",
                            ns = ns,
                            shiny::numericInput(ns("block_gap_fixed_bp"), "Fixed Gap (bp)", value = 500000, min = 1000, step = 10000)
                        ),
                        shiny::numericInput(ns("k_neighbors"), "K Neighbors", value = 5, min = 1, max = 20, step = 1),
                        shiny::hr(),
                        shiny::actionButton(ns("generate_plot"), "Generate Plot", class = "btn-primary w-100")
                    )
                )
            ),
            shiny::column(9,
                bslib::card(
                    bslib::card_header("Block Formation Diagnostic"),
                    bslib::card_body(
                        shiny::plotOutput(ns("block_plot"), height = "600px")
                    )
                )
            )
        )
    )
}

.blockFormationServer <- function(id, data) {
    shiny::moduleServer(id, function(input, output, session) {
        shiny::observe({
            shiny::req(data$dmrs)
            chrs <- unique(as.character(GenomicRanges::seqnames(data$dmrs)))
            chrs <- chrs[order(suppressWarnings(as.numeric(gsub("chr", "", chrs, ignore.case = TRUE))))]
            shiny::updateSelectInput(session, "chromosome", choices = chrs)

            cols <- colnames(S4Vectors::mcols(data$dmrs))
            score_cols <- intersect(c("score", "rank", "delta_beta"), cols)
            shiny::updateSelectInput(session, "score_col", choices = score_cols)
        })

        plot_data <- shiny::eventReactive(input$generate_plot, {
            shiny::req(data$dmrs, input$chromosome)
            list(
                chromosome = input$chromosome,
                score_col = input$score_col,
                block_gap_mode = input$block_gap_mode,
                block_gap_fixed_bp = input$block_gap_fixed_bp,
                k_neighbors = input$k_neighbors
            )
        })

        output$block_plot <- shiny::renderPlot({
            shiny::req(plot_data())
            params <- plot_data()
            shiny::withProgress(message = "Generating block formation plot...", value = 0.5, {
                plotDMRBlockFormation(
                    dmrs = data$dmrs,
                    chromosome = params$chromosome,
                    score_col = params$score_col,
                    genome = data$genome,
                    k_neighbors = params$k_neighbors,
                    block_gap_mode = params$block_gap_mode,
                    block_gap_fixed_bp = if (params$block_gap_mode == "fixed") params$block_gap_fixed_bp else NULL
                )
            })
        }, res = 100)
    })
}

.circosUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(3,
                bslib::card(
                    bslib::card_header("Circos Plot Settings"),
                    bslib::card_body(
                        shiny::radioButtons(ns("mode"), "Mode",
                            choices = c("Auto Selection" = "auto", "Manual Region" = "manual"),
                            selected = "auto", inline = TRUE
                        ),
                        shiny::hr(),
                        shiny::conditionalPanel(
                            condition = "input.mode == 'auto'",
                            ns = ns,
                            shiny::selectInput(ns("method"), "Selection Method",
                                choices = c("blocks", "components", "hybrid", "quick"),
                                selected = "blocks"
                            ),
                            shiny::numericInput(ns("n_regions"), "Number of Regions", value = 6, min = 1, max = 20, step = 1)
                        ),
                        shiny::conditionalPanel(
                            condition = "input.mode == 'manual'",
                            ns = ns,
                            shiny::textInput(ns("region"), "Region (chr:start-end)", placeholder = "e.g., chr7:100000-2000000"),
                            shiny::selectizeInput(ns("chromosomes"), "Chromosomes",
                                choices = c(), multiple = TRUE, options = list(placeholder = "Select chromosomes")
                            )
                        ),
                        shiny::hr(),
                        shiny::h6("Common Settings"),
                        shiny::numericInput(ns("max_dmrs_per_chr"), "Max DMRs/Chromosome", value = 10, min = 1, max = 50, step = 1),
                        shiny::numericInput(ns("max_cpgs_per_dmr"), "Max CpGs/DMR", value = 5, min = 1, max = 20, step = 1),
                        shiny::sliderInput(ns("min_similarity"), "Min Motif Similarity", min = 0.5, max = 1, value = 0.8, step = 0.05),
                        shiny::hr(),
                        shiny::actionButton(ns("generate_plot"), "Generate Plot", class = "btn-primary w-100"),
                        shiny::br(), shiny::br(),
                        shiny::downloadButton(ns("download_plot"), "Download PDF", class = "w-100")
                    )
                )
            ),
            shiny::column(9,
                bslib::card(
                    bslib::card_header("Circos Plot"),
                    bslib::card_body(
                        shiny::div(
                            style = "display: flex; justify-content: center;",
                            shiny::plotOutput(ns("circos_plot"), height = "700px", width = "700px")
                        )
                    )
                )
            )
        )
    )
}

.circosServer <- function(id, data) {
    shiny::moduleServer(id, function(input, output, session) {
        shiny::observe({
            shiny::req(data$dmrs)
            chrs <- unique(as.character(GenomicRanges::seqnames(data$dmrs)))
            chrs <- chrs[order(suppressWarnings(as.numeric(gsub("chr", "", chrs, ignore.case = TRUE))))]
            shiny::updateSelectizeInput(session, "chromosomes", choices = chrs, server = TRUE)
        })

        plot_data <- shiny::eventReactive(input$generate_plot, {
            shiny::req(data$dmrs, data$beta_handler, data$pheno)
            list(
                mode = input$mode,
                method = input$method,
                n_regions = input$n_regions,
                region = if (input$mode == "manual" && nzchar(trimws(input$region))) input$region else NULL,
                chromosomes = if (input$mode == "manual" && length(input$chromosomes) > 0) input$chromosomes else NULL,
                max_dmrs_per_chr = input$max_dmrs_per_chr,
                max_cpgs_per_dmr = input$max_cpgs_per_dmr,
                min_similarity = input$min_similarity
            )
        })

        output$circos_plot <- shiny::renderPlot({
            shiny::req(plot_data())
            params <- plot_data()
            shiny::withProgress(message = "Generating Circos plot...", value = 0.5, {
                if (params$mode == "auto") {
                    plotAutoDMRsCircos(
                        dmrs = data$dmrs,
                        beta = data$beta_handler,
                        pheno = data$pheno,
                        method = params$method,
                        n_regions = params$n_regions,
                        genome = data$genome,
                        array = data$array,
                        components = data$components,
                        interactions = data$interactions,
                        min_similarity = params$min_similarity,
                        max_dmrs_per_chr = params$max_dmrs_per_chr,
                        max_cpgs_per_dmr = params$max_cpgs_per_dmr,
                        sample_group_col = data$sample_group_col
                    )
                } else {
                    plotDMRsCircos(
                        dmrs = data$dmrs,
                        beta = data$beta_handler,
                        pheno = data$pheno,
                        genome = data$genome,
                        array = data$array,
                        region = params$region,
                        chromosomes = params$chromosomes,
                        components = data$components,
                        interactions = data$interactions,
                        min_similarity = params$min_similarity,
                        max_dmrs_per_chr = params$max_dmrs_per_chr,
                        max_cpgs_per_dmr = params$max_cpgs_per_dmr,
                        sample_group_col = data$sample_group_col
                    )
                }
            })
        }, res = 100)

        output$download_plot <- shiny::downloadHandler(
            filename = function() {
                paste0("Circos_", Sys.Date(), ".pdf")
            },
            content = function(file) {
                params <- plot_data()
                shiny::req(params)
                if (params$mode == "auto") {
                    plotAutoDMRsCircos(
                        dmrs = data$dmrs,
                        beta = data$beta_handler,
                        pheno = data$pheno,
                        method = params$method,
                        n_regions = params$n_regions,
                        genome = data$genome,
                        array = data$array,
                        components = data$components,
                        interactions = data$interactions,
                        min_similarity = params$min_similarity,
                        max_dmrs_per_chr = params$max_dmrs_per_chr,
                        max_cpgs_per_dmr = params$max_cpgs_per_dmr,
                        sample_group_col = data$sample_group_col,
                        output_file = file
                    )
                } else {
                    plotDMRsCircos(
                        dmrs = data$dmrs,
                        beta = data$beta_handler,
                        pheno = data$pheno,
                        genome = data$genome,
                        array = data$array,
                        region = params$region,
                        chromosomes = params$chromosomes,
                        components = data$components,
                        interactions = data$interactions,
                        min_similarity = params$min_similarity,
                        max_dmrs_per_chr = params$max_dmrs_per_chr,
                        max_cpgs_per_dmr = params$max_cpgs_per_dmr,
                        sample_group_col = data$sample_group_col,
                        output_file = file
                    )
                }
            }
        )
    })
}
