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
                .viewerWithSpinner(
                    DT::DTOutput(ns("dmr_table")),
                    proxy_height = "420px"
                )
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
            as.character(ncol(data$beta))
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
                  "delta_beta", "cases_beta", "controls_beta", "score", "block_id",
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
                        shiny::hr(),
                        shiny::actionButton(ns("generate_plot"), "Generate Plot", class = "btn-primary w-100"),
                        shiny::br(), shiny::br(),
                        shiny::conditionalPanel(
                            condition = "output.has_plot",
                            ns = ns,
                            shiny::downloadButton(ns("download_plot"), "Download PDF", class = "w-100")
                        )
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
                        .viewerWithSpinner(
                            shiny::plotOutput(ns("dmr_plot"), height = "700px")
                        )
                    )
                )
            )
        )
    )
}

.plotDMRServer <- function(id, data, task_controller) {
    shiny::moduleServer(id, function(input, output, session) {
        plot_request <- shiny::reactiveVal(0L)
        plot_params <- shiny::reactiveVal(NULL)
        requested_params <- shiny::reactiveVal(NULL)

        shiny::observe({
            shiny::req(data$dmrs)
            shiny::updateNumericInput(session, "dmr_index",
                max = length(data$dmrs),
                value = min(1, length(data$dmrs))
            )
        })

        shiny::observeEvent(input$generate_plot, {
            shiny::req(data$dmrs, input$dmr_index)
            shiny::validate(
                shiny::need(input$dmr_index >= 1 && input$dmr_index <= length(data$dmrs),
                    "Invalid DMR index")
            )
            params <- list(
                dmr_index = input$dmr_index,
                max_cpgs = input$max_cpgs,
                max_samples_per_group = input$max_samples_per_group
            )
            requested_params(params)
            plot_request(plot_request() + 1L)
        }, ignoreInit = TRUE)

        output$has_plot <- shiny::reactive({
            !is.null(plot_params())
        })
        shiny::outputOptions(output, "has_plot", suspendWhenHidden = FALSE)

        output$dmr_plot <- shiny::renderPlot({
            shiny::req(plot_request() > 0L, requested_params())
            params <- requested_params()
            tryCatch(
                {
                    plot_obj <- plotDMR(
                        dmrs = data$dmrs,
                        dmr_index = params$dmr_index,
                        beta = data$beta_handler,
                        pheno = data$pheno,
                        genome = data$genome,
                        array = data$array,
                        sample_group_col = data$sample_group_col,
                        max_cpgs = params$max_cpgs,
                        max_samples_per_group = params$max_samples_per_group,
                        plot_title = TRUE
                    )
                    shiny::isolate(plot_params(params))
                    plot_obj
                },
                error = function(e) {
                    shiny::showNotification(
                        paste0("DMR plot failed: ", conditionMessage(e)),
                        type = "error",
                        duration = NULL
                    )
                    stop(e)
                }
            )
        }, res = 100, execOnResize = FALSE)

        output$download_plot <- shiny::downloadHandler(
            filename = function() {
                paste0("DMR_", input$dmr_index, "_", Sys.Date(), ".pdf")
            },
            content = function(file) {
                params <- plot_params()
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
                        .viewerWithSpinner(
                            shiny::plotOutput(ns("dmrs_plot"), height = "800px")
                        )
                    )
                )
            )
        )
    )
}

.viewerPlotlyPanelRange <- function(plot_build, axis = c("x", "y")) {
    axis <- match.arg(axis)
    panel_params <- plot_build$layout$panel_params[[1]]
    range_name <- paste0(axis, ".range")

    if (!is.null(panel_params[[range_name]])) {
        return(as.numeric(panel_params[[range_name]]))
    }

    axis_params <- panel_params[[axis]]
    if (is.null(axis_params)) {
        return(NULL)
    }

    if (!is.null(axis_params$continuous_range)) {
        return(as.numeric(axis_params$continuous_range))
    }
    if (!is.null(axis_params$range$range)) {
        return(as.numeric(axis_params$range$range))
    }
    if (!is.null(axis_params$range)) {
        return(as.numeric(axis_params$range))
    }

    NULL
}

.prepareViewerPlotForPlotly <- function(plot_obj) {
    plot_build <- ggplot2::ggplot_build(plot_obj)
    x_range <- .viewerPlotlyPanelRange(plot_build, axis = "x")
    y_range <- .viewerPlotlyPanelRange(plot_build, axis = "y")

    if (is.null(x_range) && is.null(y_range)) {
        return(plot_obj)
    }

    prepared_plot <- plot_obj
    for (i in seq_along(prepared_plot$layers)) {
        if (!inherits(prepared_plot$layers[[i]]$geom, "GeomRect")) {
            next
        }

        layer_data <- prepared_plot$layers[[i]]$data
        if (!is.data.frame(layer_data) || nrow(layer_data) == 0) {
            next
        }

        if (!is.null(x_range)) {
            if ("xmin" %in% colnames(layer_data)) {
                layer_data$xmin[!is.finite(layer_data$xmin)] <- x_range[1]
            }
            if ("xmax" %in% colnames(layer_data)) {
                layer_data$xmax[!is.finite(layer_data$xmax)] <- x_range[2]
            }
        }

        if (!is.null(y_range)) {
            if ("ymin" %in% colnames(layer_data)) {
                layer_data$ymin[!is.finite(layer_data$ymin)] <- y_range[1]
            }
            if ("ymax" %in% colnames(layer_data)) {
                layer_data$ymax[!is.finite(layer_data$ymax)] <- y_range[2]
            }
        }

        prepared_plot$layers[[i]]$data <- layer_data
    }

    prepared_plot
}

.repairViewerPlotlyFilledTraces <- function(widget) {
    x_range <- suppressWarnings(as.numeric(widget$x$layout$xaxis$range))
    y_range <- suppressWarnings(as.numeric(widget$x$layout$yaxis$range))

    for (i in seq_along(widget$x$data)) {
        trace <- widget$x$data[[i]]
        if (!identical(trace$fill, "toself")) {
            next
        }

        if (!is.null(trace$x) && length(x_range) == 2L) {
            x_vals <- suppressWarnings(as.numeric(trace$x))
            if (any(!is.finite(x_vals))) {
                x_vals[is.infinite(x_vals) & x_vals < 0] <- x_range[1]
                x_vals[is.infinite(x_vals) & x_vals > 0] <- x_range[2]
                trace$x <- x_vals
            }
        }

        if (!is.null(trace$y) && length(y_range) == 2L) {
            y_vals <- suppressWarnings(as.numeric(trace$y))
            if (any(!is.finite(y_vals))) {
                y_vals[is.infinite(y_vals) & y_vals < 0] <- y_range[1]
                y_vals[is.infinite(y_vals) & y_vals > 0] <- y_range[2]
                trace$y <- y_vals
            }
        }

        widget$x$data[[i]] <- trace
    }

    widget
}

.configureViewerPlotly <- function(plot_obj) {
    plot_obj <- .prepareViewerPlotForPlotly(plot_obj)
    widget <- plotly::ggplotly(
        plot_obj,
        tooltip = "text",
        dynamicTicks = FALSE
    )
    widget <- .repairViewerPlotlyFilledTraces(widget)
    widget <- plotly::layout(
        widget,
        dragmode = "pan",
        hovermode = "closest",
        legend = list(orientation = "h", x = 0, y = -0.18),
        margin = list(l = 60, r = 20, t = 70, b = 80)
    )
    plotly::config(
        widget,
        displaylogo = FALSE,
        responsive = TRUE,
        scrollZoom = TRUE
    )
}

.manhattanUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(3,
                bslib::card(
                    bslib::card_header("Manhattan Plot Settings"),
                    bslib::card_body(
                        shiny::textInput(
                            ns("region"),
                            "Plotting Region (optional)",
                            placeholder = "e.g., chr7:100000-2000000 or chr5:1-5e7,chr11:1-6e7"
                        ),
                        shiny::helpText("Leave empty to plot the full chromosome span for every chromosome with DMRs."),
                        shiny::sliderInput(ns("point_size"), "Point Size", min = 0.5, max = 3, value = 1.1, step = 0.1),
                        shiny::sliderInput(ns("point_alpha"), "Point Transparency", min = 0.1, max = 1, value = 0.75, step = 0.05),
                        shiny::hr(),
                        shiny::actionButton(ns("generate_plot"), "Generate Plot", class = "btn-primary w-100"),
                        shiny::br(), shiny::br(),
                        shiny::conditionalPanel(
                            condition = "output.has_plot",
                            ns = ns,
                            shiny::downloadButton(ns("download_plot"), "Download PDF", class = "w-100")
                        )
                    )
                )
            ),
            shiny::column(9,
                bslib::card(
                    bslib::card_header("Manhattan Plot"),
                    bslib::card_body(
                        shiny::div(
                            class = "small text-muted mb-2",
                            "Scroll to zoom, drag to pan, and hover over points or shaded blocks for DMR details."
                        ),
                        .viewerWithSpinner(
                            plotly::plotlyOutput(ns("manhattan_plot"), height = "500px")
                        )
                    )
                )
            )
        )
    )
}

.manhattanServer <- function(id, data, task_controller) {
    shiny::moduleServer(id, function(input, output, session) {
        plot_request <- shiny::reactiveVal(0L)
        plot_params <- shiny::reactiveVal(NULL)
        requested_params <- shiny::reactiveVal(NULL)

        shiny::observe({
            shiny::req(data$dmrs)
        })

        shiny::observeEvent(input$generate_plot, {
            shiny::req(data$dmrs)
            params <- list(
                region = if (nzchar(trimws(input$region))) trimws(input$region) else NULL,
                point_size = input$point_size,
                point_alpha = input$point_alpha
            )
            requested_params(params)
            plot_request(plot_request() + 1L)
        }, ignoreInit = TRUE)

        output$has_plot <- shiny::reactive({
            !is.null(plot_params())
        })
        shiny::outputOptions(output, "has_plot", suspendWhenHidden = FALSE)

        output$manhattan_plot <- plotly::renderPlotly({
            shiny::req(plot_request() > 0L, requested_params())
            params <- requested_params()
            tryCatch(
                {
                    .log_info("Generating Manhattan plot..", level = 3)
                    plot_obj <- plotDMRsManhattan(
                        dmrs = data$dmrs,
                        region = params$region,
                        genome = data$genome,
                        point_size = params$point_size,
                        point_alpha = params$point_alpha
                    )
                    .log_info("Manhattan plot generated successfully.", level = 3)
                    .log_info("Configuring plotly interactivity for Manhattan plot..", level = 3)
                    shiny::isolate(plot_params(params))
                    ret <- .configureViewerPlotly(plot_obj)
                    .log_info("Plotly interactivity configured successfully.", level = 3)
                    ret
                },
                error = function(e) {
                    shiny::showNotification(
                        paste0("Manhattan plot failed: ", conditionMessage(e)),
                        type = "error",
                        duration = NULL
                    )
                    stop(e)
                }
            )
        })

        output$download_plot <- shiny::downloadHandler(
            filename = function() {
                paste0("Manhattan_", Sys.Date(), ".pdf")
            },
            content = function(file) {
                params <- plot_params()
                shiny::req(params)
                plotDMRsManhattan(
                    dmrs = data$dmrs,
                    region = params$region,
                    genome = data$genome,
                    point_size = params$point_size,
                    point_alpha = params$point_alpha,
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
                        shiny::div(
                            class = "small text-muted mb-2",
                            "Use the plotly controls to zoom into local score transitions, inspect candidate spans, and follow accepted blocks."
                        ),
                        .viewerWithSpinner(
                            plotly::plotlyOutput(ns("block_plot"), height = "600px")
                        )
                    )
                )
            )
        )
    )
}

.blockFormationServer <- function(id, data, task_controller) {
    shiny::moduleServer(id, function(input, output, session) {
        plot_request <- shiny::reactiveVal(0L)
        plot_params <- shiny::reactiveVal(NULL)

        shiny::observe({
            shiny::req(data$dmrs)
            chrs <- unique(as.character(GenomicRanges::seqnames(data$dmrs)))
            chrs <- chrs[order(suppressWarnings(as.numeric(gsub("chr", "", chrs, ignore.case = TRUE))))]
            shiny::updateSelectInput(session, "chromosome", choices = chrs)
        })

        shiny::observeEvent(input$generate_plot, {
            shiny::req(data$dmrs, input$chromosome)
            params <- list(
                chromosome = input$chromosome,
                block_gap_mode = input$block_gap_mode,
                block_gap_fixed_bp = input$block_gap_fixed_bp,
                k_neighbors = input$k_neighbors
            )
            plot_params(params)
            plot_request(plot_request() + 1L)
        }, ignoreInit = TRUE)

        output$block_plot <- plotly::renderPlotly({
            shiny::req(plot_request() > 0L, plot_params())
            params <- plot_params()
            tryCatch(
                {
                    .log_info("Generating block formation plot..", params$chromosome, params$block_gap_mode, level = 3)
                    plot_obj <- plotDMRBlockFormation(
                        dmrs = data$dmrs,
                        chromosome = params$chromosome,
                        genome = data$genome,
                        k_neighbors = params$k_neighbors,
                        block_gap_mode = params$block_gap_mode,
                        block_gap_fixed_bp = if (identical(params$block_gap_mode, "fixed")) params$block_gap_fixed_bp else NULL
                    )
                    .log_info("Block formation plot generated successfully.", params$chromosome, level = 3)
                    .log_info("Configuring plotly interactivity for block formation plot..", params$chromosome, level = 3)
                    ret <- .configureViewerPlotly(plot_obj)
                    .log_info("Plotly interactivity configured successfully.", params$chromosome, level = 3)
                    ret
                },
                error = function(e) {
                    shiny::showNotification(
                        paste0("Block formation plot failed: ", conditionMessage(e)),
                        type = "error",
                        duration = NULL
                    )
                    stop(e)
                }
            )
        })
    })
}

.viewerHasUsableCircosComponents <- function(components) {
    is.data.frame(components) &&
        nrow(components) > 0 &&
        "component_id" %in% colnames(components) &&
        "indices" %in% colnames(components) &&
        any(vapply(seq_len(nrow(components)), function(i) {
            length(.parseDMRComponentIndices(components$indices[[i]])) > 0
        }, logical(1)))
}

.resolveViewerModuleData <- function(data) {
    if (is.function(data)) {
        data()
    } else {
        data
    }
}

.viewerHasPrecomputedCircosCache <- function(data) {
    data <- .resolveViewerModuleData(data)
    is.list(data) &&
        is.data.frame(data$components) &&
        is.data.frame(data$interactions)
}

.viewerHasPrecomputedCircosInteractions <- function(data) {
    data <- .resolveViewerModuleData(data)
    is.list(data) &&
        .viewerHasUsableCircosComponents(data$components) &&
        is.data.frame(data$interactions) &&
        nrow(data$interactions) > 0
}

.mergeViewerCircosCache <- function(data, circos_cache = NULL) {
    if (!is.list(data) || is.null(circos_cache)) {
        return(data)
    }
    if (!is.null(circos_cache$interactions)) {
        data$interactions <- circos_cache$interactions
    }
    if (!is.null(circos_cache$components)) {
        data$components <- circos_cache$components
    }
    data
}

.circosInteractionStatusUI <- function(
    ns,
    has_precomputed_cache = FALSE,
    has_precomputed_interactions = FALSE,
    is_computing_interactions = FALSE
) {
    if (isTRUE(has_precomputed_interactions)) {
        return(
            shiny::div(
                class = "alert alert-success py-2 small",
                shiny::tags$strong("Precomputed Circos cache available. "),
                "Component-aware and hybrid auto-selection are ready to use."
            )
        )
    }

    button_label <- if (isTRUE(has_precomputed_cache)) {
        "Recompute Interactions"
    } else {
        "Compute Interactions"
    }
    waiting_message <- if (isTRUE(is_computing_interactions)) {
        shiny::tags$div(
            class = "mt-2 text-muted",
            "Computing interactions now. The cache will reload automatically when ready."
        )
    } else {
        shiny::tags$div(
            class = "mt-2",
            shiny::actionButton(
                ns("compute_interactions"),
                button_label,
                class = "btn-outline-primary btn-sm w-100"
            )
        )
    }

    if (isTRUE(has_precomputed_cache)) {
        return(
            shiny::div(
                class = "alert alert-info py-2 small",
                shiny::tags$strong("Precomputed Circos cache loaded, but no interaction links are available. "),
                "The viewer is staying in lightweight mode until a usable interaction cache is available.",
                waiting_message
            )
        )
    }

    shiny::div(
        class = "alert alert-secondary py-2 small",
        shiny::tags$strong("Precomputed Circos cache not found or is incomplete. "),
        "The viewer is staying in lightweight mode, so only block-based and quick auto-selection are available.",
        waiting_message
    )
}

.viewerCircosInteractionArgs <- function(data, has_precomputed_interactions = TRUE) {
    if (isTRUE(has_precomputed_interactions)) {
        list(components = data$components, interactions = data$interactions)
    } else {
        list(components = data.frame(), interactions = data.frame())
    }
}

.circosViewerAutoMethodChoices <- function(has_precomputed_interactions = TRUE) {
    if (isTRUE(has_precomputed_interactions)) {
        c("blocks" = "blocks", "components" = "components", "hybrid" = "hybrid", "quick" = "quick")
    } else {
        c("blocks" = "blocks", "quick" = "quick")
    }
}

.sanitizeViewerCircosParams <- function(params, has_precomputed_interactions = TRUE) {
    ret <- params
    ret$use_precomputed_interactions <- isTRUE(has_precomputed_interactions)
    ret$query_components_with_jaspar <- isTRUE(has_precomputed_interactions)
    ret$method_fallback <- FALSE
    if (
        !ret$use_precomputed_interactions &&
        identical(ret$mode, "auto") &&
        ret$method %in% c("components", "hybrid")
    ) {
        ret$method <- "blocks"
        ret$method_fallback <- TRUE
    }
    ret
}

.viewerCircosPreparedStateKey <- function(params) {
    paste(
        c(
            paste0("precomputed=", isTRUE(params$use_precomputed_interactions)),
            paste0("jaspar=", isTRUE(params$query_components_with_jaspar)),
            paste0("min_similarity=", format(params$min_similarity, trim = TRUE, scientific = FALSE)),
            paste0("max_dmrs_per_chr=", params$max_dmrs_per_chr),
            paste0("max_cpgs_per_dmr=", params$max_cpgs_per_dmr)
        ),
        collapse = "::"
    )
}

.buildViewerCircosPreparedState <- function(params, data) {
    interaction_args <- .viewerCircosInteractionArgs(
        data,
        has_precomputed_interactions = params$use_precomputed_interactions
    )

    do.call(
        .prepareCircosPlotState,
        c(
            list(
                dmrs = data$dmrs,
                beta = data$beta_handler,
                pheno = data$pheno,
                genome = data$genome,
                array = data$array,
                sample_group_col = data$sample_group_col,
                min_similarity = params$min_similarity,
                max_dmrs_per_chr = params$max_dmrs_per_chr,
                max_cpgs_per_dmr = params$max_cpgs_per_dmr,
                query_components_with_jaspar = params$query_components_with_jaspar
            ),
            interaction_args
        )
    )
}

.runViewerCircosPlot <- function(params, data, output_file = NULL, prepared_plot_state = NULL) {
    if (is.null(prepared_plot_state)) {
        prepared_plot_state <- .buildViewerCircosPreparedState(params, data)
    }

    if (params$mode == "auto") {
        region_df <- .selectCircosRegions(
            dmrs = prepared_plot_state$dmrs,
            method = params$method,
            n_regions = params$n_regions,
            components = prepared_plot_state$components
        )
        .log_info(
            "Automatically selected Circos regions (", params$method, "): ",
            .formatCircosRegionSelection(region_df),
            level = 2
        )
        .plotDMRsCircosFromPreparedState(
            prepared_state = prepared_plot_state,
            region = region_df,
            output_file = output_file
        )
    } else {
        .plotDMRsCircosFromPreparedState(
            prepared_state = prepared_plot_state,
            region = params$region,
            chromosomes = params$chromosomes,
            output_file = output_file
        )
    }

    invisible(prepared_plot_state)
}

.circosUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(3,
                bslib::card(
                    bslib::card_header("Circos Plot Settings"),
                    bslib::card_body(
                        .viewerWithSpinner(
                            shiny::uiOutput(ns("interaction_status")),
                            proxy_height = "110px"
                        ),
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
                        shiny::conditionalPanel(
                            condition = "output.has_plot",
                            ns = ns,
                            shiny::downloadButton(ns("download_plot"), "Download PDF", class = "w-100")
                        )
                    )
                )
            ),
            shiny::column(9,
                bslib::card(
                    bslib::card_header("Circos Plot"),
                    bslib::card_body(
                        shiny::div(
                            class = "small text-muted mb-2",
                            "This viewer uses the static circlize plot in the app. PDF export remains available below."
                        ),
                        .viewerWithSpinner(
                            shiny::div(
                                style = "display: flex; justify-content: center;",
                                shiny::plotOutput(ns("circos_plot"), height = "700px", width = "700px")
                            )
                        )
                    )
                )
            )
        )
    )
}

.circosServer <- function(id, data, task_controller) {
    shiny::moduleServer(id, function(input, output, session) {
        circos_cache_override <- shiny::reactiveVal(NULL)
        circos_prepared_state <- shiny::reactiveVal(NULL)
        plot_request <- shiny::reactiveVal(0L)
        plot_params <- shiny::reactiveVal(NULL)
        requested_params <- shiny::reactiveVal(NULL)
        viewer_data <- shiny::reactive({
            .mergeViewerCircosCache(
                .resolveViewerModuleData(data),
                circos_cache_override()
            )
        })
        is_computing_interactions <- shiny::reactive({
            state <- task_controller$state()
            isTRUE(state$active) && identical(state$task_type, "circos_cache_compute")
        })
        has_precomputed_cache <- shiny::reactive({
            .viewerHasPrecomputedCircosCache(viewer_data())
        })
        has_precomputed_interactions <- shiny::reactive({
            .viewerHasPrecomputedCircosInteractions(viewer_data())
        })

        shiny::observe({
            shiny::req(viewer_data()$dmrs)
            chrs <- unique(as.character(GenomicRanges::seqnames(viewer_data()$dmrs)))
            chrs <- chrs[order(suppressWarnings(as.numeric(gsub("chr", "", chrs, ignore.case = TRUE))))]
            shiny::updateSelectizeInput(session, "chromosomes", choices = chrs, server = TRUE)
        })

        output$interaction_status <- shiny::renderUI({
            .circosInteractionStatusUI(
                ns = session$ns,
                has_precomputed_cache = has_precomputed_cache(),
                has_precomputed_interactions = has_precomputed_interactions(),
                is_computing_interactions = is_computing_interactions()
            )
        })

        shiny::observe({
            choices <- .circosViewerAutoMethodChoices(has_precomputed_interactions())
            selected <- if (!is.null(input$method) && input$method %in% choices) input$method else "blocks"
            shiny::updateSelectInput(session, "method", choices = choices, selected = selected)
        })

        shiny::observeEvent(input$compute_interactions, {
            current_data <- viewer_data()
            shiny::req(
                current_data$dmrs,
                current_data$beta_handler,
                current_data$output_prefix
            )

            shiny::removeNotification(id = session$ns("circos_cache_status"))
            task_controller$start(
                task_type = "circos_cache_compute",
                params = list(min_similarity = input$min_similarity),
                on_success = function(result) {
                    computed_cache <- result$cache
                    circos_cache_override(computed_cache)
                    circos_prepared_state(NULL)
                    updated_data <- .mergeViewerCircosCache(current_data, computed_cache)
                    notification_text <- if (.viewerHasPrecomputedCircosInteractions(updated_data)) {
                        "Precomputed Circos cache reloaded."
                    } else {
                        "Precomputed Circos cache reloaded, but no usable interaction links were found for the current settings."
                    }
                    notification_type <- if (.viewerHasPrecomputedCircosInteractions(updated_data)) {
                        "message"
                    } else {
                        "warning"
                    }
                    shiny::showNotification(
                        notification_text,
                        type = notification_type,
                        duration = 6,
                        id = session$ns("circos_cache_status")
                    )
                },
                on_error = function(e) {
                    .log_warn("Failed to compute precomputed Circos cache in viewer: ", conditionMessage(e), level = 1)
                    shiny::showNotification(
                        paste0("Computing interactions failed: ", conditionMessage(e)),
                        type = "error",
                        duration = NULL,
                        id = session$ns("circos_cache_status")
                    )
                },
                on_cancel = function() {
                    shiny::showNotification(
                        "Computing interactions canceled.",
                        type = "message",
                        duration = 4,
                        id = session$ns("circos_cache_status")
                    )
                }
            )
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$generate_plot, {
            shiny::req(viewer_data()$dmrs, viewer_data()$beta_handler, viewer_data()$pheno)
            params <- .sanitizeViewerCircosParams(
                list(
                    mode = input$mode,
                    method = input$method,
                    n_regions = input$n_regions,
                    region = if (input$mode == "manual" && nzchar(trimws(input$region))) input$region else NULL,
                    chromosomes = if (input$mode == "manual" && length(input$chromosomes) > 0) input$chromosomes else NULL,
                    max_dmrs_per_chr = input$max_dmrs_per_chr,
                    max_cpgs_per_dmr = input$max_cpgs_per_dmr,
                    min_similarity = input$min_similarity
                ),
                has_precomputed_interactions = has_precomputed_interactions()
            )
            prepared_state_key <- .viewerCircosPreparedStateKey(params)
            cached_prepared_state <- circos_prepared_state()
            if (
                is.list(cached_prepared_state) &&
                identical(cached_prepared_state$key, prepared_state_key) &&
                !is.null(cached_prepared_state$state)
            ) {
                params$prepared_plot_state <- cached_prepared_state$state
            }

            shiny::removeNotification(id = session$ns("circos_error"))
            requested_params(params)
            plot_request(plot_request() + 1L)
        }, ignoreInit = TRUE)

        output$has_plot <- shiny::reactive({
            !is.null(plot_params())
        })
        shiny::outputOptions(output, "has_plot", suspendWhenHidden = FALSE)

        output$circos_plot <- shiny::renderPlot({
            shiny::req(plot_request() > 0L, requested_params())
            params <- requested_params()
            prepared_state_key <- .viewerCircosPreparedStateKey(params)
            tryCatch(
                {
                    prepared_plot_state <- params$prepared_plot_state
                    if (is.null(prepared_plot_state)) {
                        prepared_plot_state <- .runViewerCircosPlot(
                            params,
                            viewer_data(),
                            prepared_plot_state = params$prepared_plot_state
                        )
                    }
                    shiny::isolate(plot_params(params))
                    if (!is.null(prepared_plot_state)) {
                        shiny::isolate(circos_prepared_state(list(
                            key = prepared_state_key,
                            state = prepared_plot_state
                        )))
                    }
                    shiny::removeNotification(id = session$ns("circos_error"))
                },
                error = function(e) {
                    .log_warn("Failed to generate Circos plot in viewer: ", conditionMessage(e), level = 1)
                    shiny::showNotification(
                        paste0("Circos plot failed: ", conditionMessage(e)),
                        type = "error",
                        duration = NULL,
                        id = session$ns("circos_error")
                    )
                    stop(e)
                }
            )
        }, res = 100, execOnResize = FALSE)

        output$download_plot <- shiny::downloadHandler(
            filename = function() {
                paste0("Circos_", Sys.Date(), ".pdf")
            },
            content = function(file) {
                params <- plot_params()
                shiny::req(params)
                cached_prepared_state <- circos_prepared_state()
                prepared_state <- NULL
                if (
                    is.list(cached_prepared_state) &&
                    identical(cached_prepared_state$key, .viewerCircosPreparedStateKey(params))
                ) {
                    prepared_state <- cached_prepared_state$state
                }
                .runViewerCircosPlot(
                    params,
                    viewer_data(),
                    output_file = file,
                    prepared_plot_state = prepared_state
                )
            }
        )
    })
}
