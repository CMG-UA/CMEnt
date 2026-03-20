# Avoid NSE warnings from R CMD check for ggplot2 aes()
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "Sample", "Beta", "Position", "x", "xend", "y", "yend", "start", "position",
        "score", "region_class", "chr", "target_x", "target_y", "label_x", "label_y",
        "label", "Group", "block_id", "xmin", "xmax", "ymin", "ymax", "midpoint",
        "score_raw", "score_smoothed", "end_bp", "start_bp", "right_bp", "slope"
    ))
}

.resolveShowtextDpi <- function(default = 300) {
    tryCatch({
        px <- grDevices::dev.size("px")
        inches <- grDevices::dev.size("in")

        if (
            length(px) != 2 ||
            length(inches) != 2 ||
            any(!is.finite(px)) ||
            any(!is.finite(inches)) ||
            any(inches <= 0)
        ) {
            return(default)
        }

        dpi <- stats::median(px / inches)
        if (!is.finite(dpi) || dpi <= 0) {
            return(default)
        }

        dpi
    }, error = function(e) default)
}

.isShowtextAutoEnabled <- function() {
    hooks <- getHook("plot.new")
    any(vapply(hooks, inherits, logical(1), "showtext_hook"))
}

#' Plot DMR Structure with seeds and Extended CpGs
#'
#' @description Visualizes the structure of Differentially Methylated Regions (DMRs)
#' identified by findDMRsFromSeeds, showing the underlying seeds as stem plots connected
#' by horizontal lines to form DMRs, with extended CpG regions shown as vertical lines.
#' The plot distinguishes between seeds (differentially methylated positions), supporting
#' CpGs that extend the DMR, and non-supporting CpGs in the surrounding region.
#'
#' @param dmrs GRanges object. Output from findDMRsFromSeeds containing DMR information.
#' @param dmr_index Integer. Which DMR to plot (default: 1).
#' @param array Character. Array platform type: "450K", "27K", "EPIC", or "EPICv2" (default: "450K"). Ignored if beta_locs is provided.
#' @param genome Character. Genome version (default: "hg38"). Ignored if beta_locs is provided.
#' @param beta_locs Data frame. Genomic locations sorted by position (optional). If NULL, will be fetched based on array and genome.
#' @param extend_by_dmr_size_ratio Numeric. Ratio of the DMR width to extend the plot region outside of the DMR on both sides (default: 0.2).
#' @param min_extension_bp Integer. Minimum extension in base pairs for the plot region (default: 50).
#' @param plot_title Logical. Whether to display the title on the plot. If FALSE, the title is logged instead (default: TRUE).
#' @param .ret_details Logical. Internal parameter to return additional details (breaks, labels, chromosome, locations) for use by plotDMR (default: FALSE).
#'
#' @return A ggplot2 object showing the DMR structure. If .ret_details is TRUE, returns a list containing the plot and additional information.
.plotDMRStructure <- function(dmrs,
                              dmr_index = 1,
                              array = c("450K", "27K", "EPIC", "EPICv2"),
                              genome = "hg38",
                              beta_locs = NULL,
                              extend_by_dmr_size_ratio = 0.2,
                              min_extension_bp = 50,
                              plot_title = TRUE,
                              .ret_details = FALSE) {
    # Validate input
    if (dmr_index < 1 || dmr_index > length(dmrs)) {
        stop("dmr_index must be between 1 and ", length(dmrs))
    }

    # Extract the specific DMR
    dmr <- dmrs[dmr_index]
    dmr_data <- S4Vectors::mcols(dmr)

    # Get genomic locations if not provided
    if (is.null(beta_locs)) {
        array <- strex::match_arg(array, ignore_case = TRUE)
        beta_locs <- getSortedGenomicLocs(array = array, genome = genome)
    }

    # Extract DMR information
    seeds <- strsplit(dmr_data$seeds, split = ",")[[1]]
    cpgs <- strsplit(dmr_data$cpgs, split = ",")[[1]]
    downstream_sup_cpgs <- strsplit(dmr_data$downstream_cpgs, split = ",")[[1]]
    downstream_sup_cpgs <- setdiff(downstream_sup_cpgs, seeds)
    downstream_sup_cpgs_locs <- as.data.frame(beta_locs[downstream_sup_cpgs, , drop = FALSE])
    upstream_sup_cpgs <- strsplit(dmr_data$upstream_cpgs, split = ",")[[1]]
    upstream_sup_cpgs <- setdiff(upstream_sup_cpgs, seeds)
    upstream_sup_cpgs_locs <- as.data.frame(beta_locs[upstream_sup_cpgs, , drop = FALSE])
    if (length(upstream_sup_cpgs) == 0) {
        start_cpg <- seeds[[1]]
    } else {
        start_cpg <- upstream_sup_cpgs[[1]]
    }
    if (length(downstream_sup_cpgs) == 0) {
        end_cpg <- seeds[[length(seeds)]]
    } else {
        end_cpg <- downstream_sup_cpgs[[length(downstream_sup_cpgs)]]
    }
    beta_locs_rownames <- rownames(beta_locs)

    dmr_locs <- beta_locs[match(start_cpg, beta_locs_rownames):match(end_cpg, beta_locs_rownames), , drop = FALSE]

    nsup_cpgs <- setdiff(rownames(dmr_locs), cpgs)
    nsup_cpgs_locs <- as.data.frame(dmr_locs[nsup_cpgs, , drop = FALSE])

    chr <- as.character(GenomicRanges::seqnames(dmr))
    dmr_start <- GenomicRanges::start(dmr)
    dmr_end <- GenomicRanges::end(dmr)


    # Get positions

    start_cpg_pos <- as.integer(dmr_locs[start_cpg, "start"])
    end_cpg_pos <- as.integer(dmr_locs[end_cpg, "start"])
    seed_positions <- as.integer(dmr_locs[seeds, "start"])
    start_seed_pos <- dmr_data$start_seed_pos
    end_seed_pos <- dmr_data$end_seed_pos

    if (extend_by_dmr_size_ratio > 0) {
        dmr_size <- dmr_end - dmr_start
        ext <- round(dmr_size * extend_by_dmr_size_ratio)
    } else {
        ext <- 0
    }
    ext <- max(ext, min_extension_bp)
    plot_start <- max(1, start_cpg_pos - ext)
    plot_end <- end_cpg_pos + ext

    downstream_nsup_cpgs_locs <- as.data.frame(dmr_locs[which(dmr_locs$start > dmr_end & dmr_locs$start <= plot_end), , drop = FALSE])
    upstream_nsup_cpgs_locs <- as.data.frame(dmr_locs[which(dmr_locs$start < dmr_start & dmr_locs$start >= plot_start), , drop = FALSE])
    extended_nsup_cpgs_locs <- rbind(
        nsup_cpgs_locs,
        upstream_nsup_cpgs_locs,
        downstream_nsup_cpgs_locs
    )


    # Create plotting data frame
    # 1. seeds (stem plots at y=1)
    seeds_df <- data.frame(
        start = seed_positions,
        y = 1,
        type = "seed",
        stringsAsFactors = FALSE
    )

    # 2. DMR connection (horizontal line connecting seeds at y=1)
    dmr_line <- data.frame(
        x = start_seed_pos,
        xend = end_seed_pos,
        y = 1,
        yend = 1,
        type = "DMR",
        stringsAsFactors = FALSE
    )

    # 3. Extended supporting CpGs (vertical lines at y=0.5)
    if (start_cpg_pos != start_seed_pos) {
        dmr_upstream_line <- data.frame(
            x = start_cpg_pos,
            xend = start_seed_pos,
            y = 0.5,
            yend = 1,
            type = "DMR_Extension",
            stringsAsFactors = FALSE
        )
    } else {
        dmr_upstream_line <- data.frame(
            x = numeric(0),
            xend = numeric(0),
            y = numeric(0),
            yend = numeric(0),
            type = character(0),
            stringsAsFactors = FALSE
        )
    }
    if (end_cpg_pos != end_seed_pos) {
        dmr_downstream_line <- data.frame(
            x = end_seed_pos,
            xend = end_cpg_pos,
            y = 1,
            yend = 0.5,
            type = "DMR_Extension",
            stringsAsFactors = FALSE
        )
    } else {
        dmr_downstream_line <- data.frame(
            x = numeric(0),
            xend = numeric(0),
            y = numeric(0),
            yend = numeric(0),
            type = character(0),
            stringsAsFactors = FALSE
        )
    }

    extended_sup_cpgs_locs <- rbind(
        upstream_sup_cpgs_locs,
        downstream_sup_cpgs_locs
    )
    if (nrow(extended_sup_cpgs_locs) > 0) {
        extended_sup_cpgs_df <- data.frame(
            start = extended_sup_cpgs_locs$start,
            y = 0.5,
            type = "Extended_CpG",
            stringsAsFactors = FALSE
        )
    } else {
        extended_sup_cpgs_df <- data.frame(
            start = numeric(0),
            y = numeric(0),
            type = character(0),
            stringsAsFactors = FALSE
        )
    }
    # 4. Non-supporting CpGs in extended region
    if (nrow(extended_nsup_cpgs_locs) > 0) {
        extended_nsup_cpgs_df <- data.frame(
            start = extended_nsup_cpgs_locs$start,
            y = 0.5,
            type = "Extended_CpG",
            stringsAsFactors = FALSE
        )
    } else {
        extended_nsup_cpgs_df <- data.frame(
            start = numeric(0),
            y = numeric(0),
            type = character(0),
            stringsAsFactors = FALSE
        )
    }

    # Create the plot
    p <- ggplot2::ggplot()


    # Plot DMR connection line
    p <- p + ggplot2::geom_segment(
        data = dmr_line,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        color = "#E41A1C",
        linewidth = 1,
        alpha = 1
    )

    # Add DMR region shading
    p <- p + ggplot2::annotate(
        "rect",
        xmin = start_seed_pos,
        xmax = end_seed_pos,
        ymin = 0,
        ymax = 1,
        alpha = 0.1,
        fill = "#E41A1C"
    )
    # if upstream extended CpGs exist add shading in the form of a trapezoid
    if (nrow(upstream_sup_cpgs_locs) > 0) {
        p <- p + ggplot2::geom_segment(
            data = dmr_upstream_line,
            ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
            color = "#E41A1C",
            linewidth = 1,
            alpha = 0.8
        )
        p <- p + ggplot2::annotate(
            "polygon",
            x = c(min(upstream_sup_cpgs_locs$start), start_seed_pos, start_seed_pos, min(upstream_sup_cpgs_locs$start)),
            y = c(0, 0, 1, 0.5),
            alpha = 0.1,
            fill = "#E41A1C"
        )
    }

    # if downstream extended CpGs exist add shading in the form of a trapezoid
    if (nrow(downstream_sup_cpgs_locs) > 0) {
        p <- p + ggplot2::geom_segment(
            data = dmr_downstream_line,
            ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
            color = "#E41A1C",
            linewidth = 1,
            alpha = 0.8
        )
        p <- p + ggplot2::annotate(
            "polygon",
            x = c(end_seed_pos, max(downstream_sup_cpgs_locs$start), max(downstream_sup_cpgs_locs$start), end_seed_pos),
            y = c(0, 0, 0.5, 1),
            alpha = 0.1,
            fill = "#E41A1C"
        )
    }


    # Plot seed stems
    p <- p + ggplot2::geom_segment(
        data = seeds_df,
        ggplot2::aes(x = start, xend = start, y = 0, yend = y),
        color = "#377EB8",
        linewidth = 0.8,
        alpha = 1
    )

    # Plot seed points
    p <- p + ggplot2::geom_point(
        data = seeds_df,
        ggplot2::aes(x = start, y = y),
        color = "#377EB8",
        size = 2,
        shape = 16,
        alpha = 1
    )

    # Plot supported CpG stems
    if (nrow(extended_sup_cpgs_df) > 0) {
        p <- p + ggplot2::geom_segment(
            data = extended_sup_cpgs_df,
            ggplot2::aes(x = start, xend = start, y = 0, yend = y),
            color = "#377EB8",
            linewidth = 0.8,
            alpha = 1
        )
    }

    # Plot supporting CpG points
    if (nrow(extended_sup_cpgs_df) > 0) {
        p <- p + ggplot2::geom_point(
            data = extended_sup_cpgs_df,
            ggplot2::aes(x = start, y = y),
            color = "#377EB8",
            size = 2,
            shape = 16,
            alpha = 1
        )
    }

    # Plot non-supported CpG stems
    if (nrow(extended_nsup_cpgs_df) > 0) {
        p <- p + ggplot2::geom_segment(
            data = extended_nsup_cpgs_df,
            ggplot2::aes(x = start, xend = start, y = 0, yend = y),
            color = "gray50",
            linewidth = 0.8,
            alpha = 0.5
        )
    }

    # Plot non-supporting CpG points
    if (nrow(extended_nsup_cpgs_df) > 0) {
        p <- p + ggplot2::geom_point(
            data = extended_nsup_cpgs_df,
            ggplot2::aes(x = start, y = y),
            color = "gray50",
            size = 2,
            shape = 16,
            alpha = 0.5
        )
    }


    # Add labels for DMR extensions if they exist
    extension_df <- list()
    if (nrow(upstream_sup_cpgs_locs) > 0) {
        upstream_mid_x <- (dmr_upstream_line$x + dmr_upstream_line$xend) / 2
        upstream_mid_y <- (dmr_upstream_line$y + dmr_upstream_line$yend) / 2
        upstream_label_df <- data.frame(
            label_x = upstream_mid_x + (dmr_upstream_line$xend - dmr_upstream_line$x) * 0.14,
            label_y = pmin(1.11, upstream_mid_y + 0.28),
            target_x = upstream_mid_x,
            target_y = upstream_mid_y,
            label = "Upstream Extension",
            stringsAsFactors = FALSE
        )
        extension_df <- c(extension_df, list(upstream_label_df))
    }
    if (nrow(downstream_sup_cpgs_locs) > 0) {
        downstream_mid_x <- (dmr_downstream_line$x + dmr_downstream_line$xend) / 2
        downstream_mid_y <- (dmr_downstream_line$y + dmr_downstream_line$yend) / 2
        downstream_label_df <- data.frame(
            label_x = downstream_mid_x - (dmr_downstream_line$xend - dmr_downstream_line$x) * 0.12,
            label_y = pmin(1.11, downstream_mid_y + 0.28),
            target_x = downstream_mid_x,
            target_y = downstream_mid_y,
            label = "Downstream Extension",
            stringsAsFactors = FALSE
        )
        extension_df <- c(extension_df, list(downstream_label_df))
    }
    if (length(extension_df) > 0) {
        p <- p + ggplot2::geom_segment(
            data = do.call(rbind, extension_df),
            ggplot2::aes(x = label_x, y = label_y - 0.02, xend = target_x, yend = target_y),
            linewidth = 0.2,
            color = "#555555"
        ) +
            ggplot2::geom_text(
                data = do.call(rbind, extension_df),
                ggplot2::aes(x = label_x, y = label_y, label = label),
                vjust = 0,
                family = "sans",
                color = "#4A4A4A",
                size = 2
            )
    }

    # Create title if not provided
    title <- sprintf(
        "DMR #%d: %s:%s-%s\nScore: %.2f\nCV Accuracy: %.2f\n%d seeds (\u0394\u03b2=%.3f)",
        dmr_index,
        chr,
        format(dmr_start, big.mark = ",", scientific = FALSE),
        format(dmr_end, big.mark = ",", scientific = FALSE),
        dmr_data$score,
        dmr_data$cv_accuracy,
        dmr_data$seeds_num,
        dmr_data$delta_beta
    )

    if (!"in_promoter_of" %in% colnames(mcols(dmrs))) {
        ret <- try(annotateDMRsWithGenes(dmr, genome = genome))
        if (inherits(ret, "try-error")) {
            warning("Failed to annotate DMR with gene information.")
        } else {
            dmr <- ret
            in_promoter_of <- dmr$in_promoter_of
            if (length(in_promoter_of) > 0 && !all(is.na(in_promoter_of))) {
                gene_labels <- paste(in_promoter_of, collapse = ", ")
                title <- paste0(title, "\n Overlapping Promoters: ", gene_labels)
            }
            body_genes <- dmr$in_gene_body_of
            if (length(body_genes) > 0 && !all(is.na(body_genes))) {
                gene_labels <- paste(body_genes, collapse = ", ")
                title <- paste0(title, "\n Overlapping Gene Bodies: ", gene_labels)
            }
        }
    }

    # Styling
    if (nrow(extended_nsup_cpgs_locs) > 0 || nrow(extended_sup_cpgs_locs) > 0) {
        p <- p +
            ggplot2::scale_y_continuous(
                breaks = c(0.5, 1),
                labels = c("Array CpGs", "Seeds\n(Used in Motif Analysis)"),
                limits = c(-0.1, 1.15)
            )
    } else {
        p <- p +
            ggplot2::scale_y_continuous(
                breaks = c(1),
                labels = c("Seeds\n(Used in Motif Analysis)"),
                limits = c(-0.1, 1.15)
            )
    }
    p <- p +
        ggplot2::labs(
            title = if (plot_title) title else NULL,
            y = ""
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 11),
            axis.text.y = ggplot2::element_text(hjust = 0.5),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank()
        )

    # Add ticks for the seeds on the x-axis
    start_cpg_ind <- which(beta_locs_rownames == start_cpg)
    end_cpg_ind <- which(beta_locs_rownames == end_cpg)
    breaks <- c(plot_start, as.integer(beta_locs[start_cpg_ind:end_cpg_ind, "start"]), plot_end)
    cpg_positions <- as.integer(beta_locs[start_cpg_ind:end_cpg_ind, "start"])
    cpg_ids <- rownames(beta_locs[start_cpg_ind:end_cpg_ind, , drop = FALSE])
    cpgs_labs <- paste0(
        format(cpg_positions, big.mark = ",", scientific = FALSE),
        " (", cpg_ids, ")"
    )
    breaks_labels <- c(
        format(plot_start, big.mark = ",", scientific = FALSE), cpgs_labs,
        format(plot_end, big.mark = ",", scientific = FALSE)
    )
    p <- p + ggplot2::scale_x_continuous(
        breaks = breaks,
        labels = breaks_labels
    )
    p <- p + ggplot2::coord_cartesian(xlim = c(breaks[1], breaks[length(breaks)]))

    if (!plot_title) {
        .log_info("Title of the generated plot:\n", title)
    }

    if (.ret_details) {
        nsup_df <- extended_nsup_cpgs_locs
        if (nrow(nsup_df) > 0) {
            nsup_df$cpg_id <- rownames(nsup_df)
        } else {
            nsup_df$cpg_id <- character(0)
        }
        sup_df <- extended_sup_cpgs_locs
        if (nrow(sup_df) > 0) {
            sup_df$cpg_id <- rownames(sup_df)
        } else {
            sup_df$cpg_id <- character(0)
        }
        seed_df <- as.data.frame(dmr_locs[seeds, , drop = FALSE])
        seed_df$cpg_id <- rownames(seed_df)
        total_shown_positions <- rbind(nsup_df, sup_df, seed_df)
        total_shown_positions <- total_shown_positions[!duplicated(total_shown_positions$cpg_id), , drop = FALSE]
        rownames(total_shown_positions) <- total_shown_positions$cpg_id
        total_shown_positions <- total_shown_positions[order(total_shown_positions$start), , drop = FALSE]
        total_shown_positions$cpg_id <- NULL
        return(invisible(list(structure_plot = p, breaks = breaks, breaks_labels = breaks_labels, chr = chr, total_locs = total_shown_positions)))
    }
    invisible(p)
}


# Create beta heatmap plot
.plotBetaHeatmap <- function(dmr_data, beta_data, total_shown_positions, pheno = NULL, max_samples_per_group = 10, sample_group_col = "Sample_Group") {
    cpg_ids <- rownames(total_shown_positions)
    cpg_locs <- total_shown_positions[, c("chr", "start")]


    # Mark seeds
    seed_ids <- unlist(strsplit(as.character(dmr_data$seeds), ","))
    is_seed <- cpg_ids %in% seed_ids

    # Create heatmap
    # Prepare data
    beta_data <- as.data.frame(beta_data)
    beta_data[, "CpG"] <- rownames(beta_data)

    # if there are more than max_samples_per_group samples in any group, limit to max_samples_per_group samples per group for plotting
    selected_samples <- NULL
    if (!is.null(pheno) && !is.null(sample_group_col)) {
        group_counts <- table(pheno[[sample_group_col]])
        if (any(group_counts > max_samples_per_group)) {
            .log_info("Limiting to ", max_samples_per_group, " samples per group for plotting. Original group counts:\n", paste(names(group_counts), group_counts, sep = ": ", collapse = "\n"))
            set.seed(123) # for reproducibility
            selected_samples <- unlist(
                lapply(
                    names(group_counts),
                    function(g) {
                        group_samples <- rownames(pheno)[pheno[[sample_group_col]] == g]
                        if (length(group_samples) > max_samples_per_group) {
                            sample(group_samples, max_samples_per_group)
                        } else {
                            group_samples
                        }
                    }
                )
            )

        }
    } else {
        # if no sample grouping information is provided, limit to max_samples_per_group samples total for plotting
        if (ncol(beta_data) - 1 > max_samples_per_group) {
            .log_info("Limiting to ", max_samples_per_group, " samples for plotting. Original number of samples: ", ncol(beta_data) - 1)
            set.seed(123) # for reproducibility
            sample_cols <- setdiff(colnames(beta_data), "CpG")
            selected_samples <- sample(sample_cols, max_samples_per_group)
        }
    }
    if (!is.null(selected_samples)) {
        beta_data <- beta_data[, c("CpG", selected_samples), drop = FALSE]
        if (!is.null(pheno) && !is.null(sample_group_col)) {
            pheno <- pheno[selected_samples, , drop = FALSE]
        }
    }

    beta_melted <- suppressWarnings(suppressMessages(reshape2::melt(beta_data, id_vars = "CpG")))
    colnames(beta_melted) <- c("CpG", "Sample", "Beta")
    beta_melted$Position <- cpg_locs[as.character(beta_melted$CpG), "start"]
    beta_melted$is_seed <- is_seed[match(beta_melted$CpG, cpg_ids)]
    sample_order <- unique(as.character(beta_melted$Sample))
    sample_label_colors <- rep("#222222", length(sample_order))
    names(sample_label_colors) <- sample_order
    if (!is.null(pheno) && !is.null(sample_group_col)) {
        beta_melted$Group <- pheno[as.character(beta_melted$Sample), sample_group_col]
        # Order samples by group
        sample_order <- rownames(pheno)[order(pheno[[sample_group_col]])]
        beta_melted$Sample <- factor(beta_melted$Sample, levels = sample_order)
        group_per_sample <- as.character(pheno[sample_order, sample_group_col])
        unique_groups <- unique(group_per_sample)
        group_colors <- colorspace::qualitative_hcl(length(unique_groups), palette = "Dark 3")
        names(group_colors) <- unique_groups
        sample_label_colors <- unname(group_colors[group_per_sample])
        names(sample_label_colors) <- sample_order
    } else {
        beta_melted$Sample <- factor(beta_melted$Sample, levels = sample_order)
    }
    valid_beta <- beta_melted$Beta[is.finite(beta_melted$Beta)]
    beta_limits <- range(valid_beta, na.rm = TRUE)
    q <- stats::quantile(valid_beta, probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE, type = 8)
    if (beta_limits[1] < 0.5 && beta_limits[2] > 0.5) {
        q <- sort(c(q[1], 0.5, q[2]))
        coloring <- ggplot2::scale_fill_gradientn(
            colours = c("#2b83ba", "#f7f7f7", "#d7191c"),
            breaks = signif(q, digits = 2),
            limits = beta_limits,
            name = "\u03b2-values"
        )
    } else if (beta_limits[2] <= 0.5) {
        coloring <- ggplot2::scale_fill_gradient(
            low = "#2b83ba",
            high = "#f7f7f7",
            breaks = signif(q, digits = 2),
            limits = beta_limits,
            name = "\u03b2-values"
        )
    } else {
        coloring <- ggplot2::scale_fill_gradient(
            low = "#f7f7f7",
            high = "#d7191c",
            breaks = signif(q, digits = 2),
            limits = beta_limits,
            name = "\u03b2-values"
        )
    }
    heatmap_plot <- ggplot2::ggplot(beta_melted) +
        ggplot2::geom_tile(ggplot2::aes(x = Position, y = Sample, fill = Beta)) +
        coloring +
        ggplot2::labs(
            y = "Sample"
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7, color = sample_label_colors[levels(beta_melted$Sample)]),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank()
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
        )
    heatmap_plot
}


minmaxscale <- function(x) {
    (x - min(x)) / max(max(x) - min(x), 1e-10)
}

.summarizeMotifContext <- function(dmrs, dmr_index, genome, array, beta_locs, motif_cpg_flank_size = 5) {
    ret <- list(top_interactions = data.frame(), jaspar = data.frame())
    if (length(dmrs) < 2) {
        return(ret)
    }

    top_n <- as.integer(getOption("DMRsegal.plotDMR_top_motif_interactions", 5L))
    pool_size <- as.integer(getOption("DMRsegal.plotDMR_interaction_pool_size", 300L))
    top_n <- ifelse(is.na(top_n) || top_n < 1, 5L, top_n)
    pool_size <- ifelse(is.na(pool_size) || pool_size < 2, min(300L, length(dmrs)), pool_size)
    pool_size <- min(pool_size, length(dmrs))

    if ("score" %in% colnames(S4Vectors::mcols(dmrs))) {
        ord <- order(S4Vectors::mcols(dmrs)$score, na.last = TRUE, decreasing = TRUE)
    } else if ("delta_beta" %in% colnames(S4Vectors::mcols(dmrs))) {
        ord <- order(abs(S4Vectors::mcols(dmrs)$delta_beta), decreasing = TRUE, na.last = TRUE)
    } else {
        ord <- seq_along(dmrs)
    }
    candidate_inds <- unique(c(ord[seq_len(pool_size)], dmr_index))
    cand_dmrs <- dmrs[candidate_inds]

    if (!"pwm" %in% colnames(S4Vectors::mcols(cand_dmrs))) {
        cand_dmrs <- extractDMRMotifs(
            cand_dmrs,
            genome = genome,
            array = array,
            beta_locs = beta_locs,
            motif_cpg_flank_size = motif_cpg_flank_size
        )
    }
    local_idx <- match(dmr_index, candidate_inds)
    if (is.na(local_idx) || !"pwm" %in% colnames(S4Vectors::mcols(cand_dmrs))) {
        return(ret)
    }
    target_pwm <- S4Vectors::mcols(cand_dmrs)$pwm[[local_idx]]
    if (is.null(target_pwm) || !is.matrix(target_pwm)) {
        return(ret)
    }

    sim_matrix <- .extractMotifsSimilarity(cand_dmrs, motif_cpg_flank_size = motif_cpg_flank_size)
    sim_vec <- sim_matrix[local_idx, ]
    sim_vec[local_idx] <- NA_real_
    valid <- which(is.finite(sim_vec))
    if (length(valid) > 0) {
        ord_sim <- valid[order(sim_vec[valid], decreasing = TRUE)]
        ord_sim <- head(ord_sim, top_n)
        top_tbl <- data.frame(
            dmr_index = candidate_inds[ord_sim],
            chr = as.character(GenomicRanges::seqnames(cand_dmrs[ord_sim])),
            start = GenomicRanges::start(cand_dmrs[ord_sim]),
            end = GenomicRanges::end(cand_dmrs[ord_sim]),
            sim = sim_vec[ord_sim],
            stringsAsFactors = FALSE
        )
        if ("score" %in% colnames(S4Vectors::mcols(dmrs))) {
            top_tbl$score <- S4Vectors::mcols(dmrs)$score[top_tbl$dmr_index]
        }
        ret$top_interactions <- top_tbl
    }

    ret$jaspar <- tryCatch(
        comparePWMToJaspar(list(target_pwm)),
        error = function(e) data.frame()
    )
    ret
}

.buildMotifContextLines <- function(motif_context, top_n = 5) {
    lines <- c("Motif Context", "")
    top_tbl <- motif_context$top_interactions
    if (nrow(top_tbl) > 0) {
        lines <- c(lines, "Top motif interactions:")
        for (i in seq_len(min(top_n, nrow(top_tbl)))) {
            score_txt <- if ("score" %in% colnames(top_tbl) && !is.na(top_tbl$score[i])) {
                paste0(" score=", top_tbl$score[i])
            } else {
                ""
            }
            lines <- c(lines, sprintf(
                "%d) #%d %s:%s-%s sim=%.3f%s",
                i,
                top_tbl$dmr_index[i],
                top_tbl$chr[i],
                format(top_tbl$start[i], big.mark = ",", scientific = FALSE),
                format(top_tbl$end[i], big.mark = ",", scientific = FALSE),
                top_tbl$sim[i],
                score_txt
            ))
        }
    } else {
        lines <- c(lines, "Top motif interactions:", "none")
    }

    jas_tbl <- motif_context$jaspar
    lines <- c(lines, "", "JASPAR matches:")
    if (nrow(jas_tbl) > 0 && "jaspar_names" %in% colnames(jas_tbl) && !is.na(jas_tbl$jaspar_names[1])) {
        names <- strsplit(jas_tbl$jaspar_names[1], ",", fixed = TRUE)[[1]]
        cors <- if ("jaspar_corr" %in% colnames(jas_tbl) && !is.na(jas_tbl$jaspar_corr[1])) {
            strsplit(jas_tbl$jaspar_corr[1], ",", fixed = TRUE)[[1]]
        } else {
            rep("", length(names))
        }
        n_show <- min(5L, length(names))
        for (i in seq_len(n_show)) {
            corr_txt <- if (i <= length(cors) && nzchar(cors[i])) paste0(" (", cors[i], ")") else ""
            lines <- c(lines, paste0("- ", trimws(names[i]), corr_txt))
        }
        if (length(names) > n_show) {
            lines <- c(lines, sprintf("... +%d more", length(names) - n_show))
        }
    } else {
        lines <- c(lines, "none")
    }
    lines
}

.plotMotifContext <- function(lines) {
    df <- data.frame(
        x = 0,
        y = rev(seq_along(lines)),
        label = lines,
        stringsAsFactors = FALSE
    )
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, label = label)) +
        ggplot2::geom_text(
            hjust = 0,
            vjust = 1,
            size = 3,
            family = "mono",
            lineheight = 0.95
        ) +
        ggplot2::xlim(0, 1) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.margin = ggplot2::margin(4, 4, 4, 4),
            panel.background = ggplot2::element_rect(fill = "white", colour = NA)
        )
}

.plotPWM <- function(dmr, genome, array, beta_locs, motif_cpg_flank_size = 5) {
    # Extract DMR motifs if not already present
    if (!"pwm" %in% colnames(S4Vectors::mcols(dmr))) {
        dmr <- extractDMRMotifs(dmr,
            genome = genome, array = array,
            beta_locs = beta_locs, motif_cpg_flank_size = motif_cpg_flank_size
        )
    }
    pwm <- mcols(dmr)$pwm[[1]]
    consensus_seq <- mcols(dmr)$consensus_seq[[1]]

    if (is.null(pwm) || !is.matrix(pwm)) {
        return(NULL)
    }

    n_positions <- ncol(pwm)
    if (n_positions == 0) {
        return(NULL)
    }

    seed_count <- NA_integer_
    if ("seeds" %in% colnames(S4Vectors::mcols(dmr))) {
        seed_ids <- strsplit(as.character(S4Vectors::mcols(dmr)$seeds[[1]]), ",", fixed = TRUE)[[1]]
        seed_count <- sum(nzchar(seed_ids))
    }
    consensus_only <- all(colSums(pwm > 0) <= 1)
    pwm_subtitle <- NULL
    if (consensus_only) {
        seed_label <- if (!is.na(seed_count)) {
            paste0(seed_count, " seed window", if (seed_count == 1) "" else "s")
        } else {
            "seed windows"
        }
        seed_verb <- if (!is.na(seed_count) && seed_count == 1) "contributes" else "contribute"
        pwm_subtitle <- paste0(
            "Consensus-only logo: ",
            seed_label,
            " ",
            seed_verb,
            " the same base at every position"
        )
    } else if (!is.na(seed_count)) {
        pwm_subtitle <- paste0(
            "Built from ",
            seed_count,
            " seed window",
            if (seed_count == 1) "" else "s"
        )
    }

    position_labels <- c(seq(-motif_cpg_flank_size, 0), seq(0, motif_cpg_flank_size))

    rownames(pwm) <- Biostrings::DNA_BASES

    suppressWarnings(suppressMessages({
        pwm_plot <- ggseqlogo::ggseqlogo(pwm, method = "custom", seq_type = "dna") +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                panel.grid.major.x = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_text(size = 8),
                axis.text.y = ggplot2::element_text(size = 8),
                axis.title = ggplot2::element_text(size = 9, face = "bold"),
                plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0),
                plot.margin = ggplot2::margin(4, 4, 4, 4),
                aspect.ratio = 0.22
            ) +
            ggplot2::labs(
                x = "Position Relative to CpG",
                y = "Relative base weight",
                subtitle = pwm_subtitle,
                title = paste0("Motif PWM (consensus: ", consensus_seq, ")")
            ) +
            ggplot2::scale_x_continuous(breaks = 1:n_positions, labels = as.character(position_labels))
    }))

    pwm_plot
}

#' Plot Multiple DMRs in a Grid
#'
#' @description Creates a grid of DMR plots for multiple regions. Can optionally
#' include beta value heatmaps if beta and pheno data are provided.
#'
#' @param dmrs GRanges object. Output from findDMRsFromSeeds.
#' @param dmr_indices Integer vector. Which DMRs to plot. If NULL, selects top_n DMRs based on delta_beta and p-value score.
#' @param top_n Integer. Number of top DMRs to plot when dmr_indices is NULL (default: 4).
#' @param beta BetaHandler object, character path to beta file, or beta values matrix (optional). If provided, creates plots with heatmaps.
#' @param pheno Data frame or character path to phenotype file (optional). Required when beta is provided.
#' @param sample_group_col Character. Column in pheno for sample grouping (default: "Sample_Group").
#' @param genome Character. Genome version (default: "hg38").
#' @param array Character. Array platform type (default: "450K"). Ignored if beta_locs is provided.
#' @param beta_locs Data frame. Genomic locations sorted by position (optional). If NULL, will be fetched based on array and genome.
#' @param ncol Integer. Number of columns in the grid (default: 1).
#' @param output_file Character. Path to save the combined plot as PDF (optional). If NULL, the plot is not saved to file.
#' @param width Numeric. Width of the output PDF in inches (default: 8).
#' @param height Numeric. Height of the output PDF in inches (default: 12).
#' @param ... Additional arguments passed to plotDMR.
#'
#' @return If beta is NULL: A gtable object.
#'   If beta is provided: A list of combined plot objects with structure and heatmap.
#'
#' @examples
#' # Plot structure only
#' dmrs <- readRDS("dmrs.rds")
#' plotDMRs(dmrs, dmr_indices = 1:6, ncol = 3)
#'
#' # Plot with beta values heatmap
#' plotDMRs(dmrs, top_n = 4, beta = "beta.txt", pheno = pheno_df)
#'
#' @export
plotDMRs <- function(dmrs,
                     dmr_indices = NULL,
                     top_n = 4,
                     beta = NULL,
                     pheno = NULL,
                     sample_group_col = "Sample_Group",
                     genome = "hg38",
                     array = c("450K", "27K", "EPIC", "EPICv2"),
                     beta_locs = NULL,
                     ncol = 1,
                     output_file = NULL,
                     width = 8,
                     height = 12,
                     ...) {
    if (!is.null(output_file)) {
        if (!grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
            stop("Output file must have a .pdf extension.")
        }
    }
    if (!is.null(array)) {
        if (length(array) > 1) {
            array <- array[[1]]
        }
        array <- strex::match_arg(array, ignore_case = TRUE)
    }
    if (is.null(beta_locs)) {
        array <- strex::match_arg(array, ignore_case = TRUE)
    }
    dmrs <- convertToGRanges(dmrs, genome)
    if (is.null(dmr_indices)) {
        score <- minmaxscale(abs(mcols(dmrs)$delta_beta))
        ord <- order(score, decreasing = TRUE)
        dmr_indices <- ord[seq_len(min(top_n, length(dmrs)))]
    }

    # Create individual plots
    if (!is.null(beta)) {
        beta <- getBetaHandler(
            beta = beta,
            array = array,
            genome = genome,
            sorted_locs = beta_locs
        )
        beta_locs <- beta$getBetaLocs()
    }
    if (!is.null(output_file)) {
        grDevices::cairo_pdf(output_file, width = width, height = height)
    }
    if (.Device == "null device") {
        grDevices::cairo_pdf(width = width, height = height)
    }
    plot_list <- lapply(seq_along(dmr_indices), function(idx) {
        i <- dmr_indices[idx]
        plotDMR(
            dmrs = dmrs,
            dmr_index = i,
            beta = beta,
            pheno = pheno,
            genome = genome,
            array = array,
            beta_locs = beta_locs,
            sample_group_col = sample_group_col,
            ...
        )
    })
    if (!is.null(output_file)) {
        grDevices::dev.off()
    }
    invisible(plot_list)
}


#' Plot DMR
#'
#' @description Creates a detailed DMR plot with an integrated heatmap showing
#' beta values across samples for seeds and surrounding CpGs. The plot consists of
#' two panels: the top panel shows the DMR structure with seeds and extended CpGs,
#' and the bottom panel displays a heatmap of beta values for all samples, if beta values are provided.
#' Additionally, if motif information is available or can be extracted, a sequence logo
#' plot is added showing the nucleotide composition and information content around CpG sites in the DMR.
#'
#' @param dmrs GRanges object. Output from findDMRsFromSeeds.
#' @param dmr_index Integer. Which DMR to plot.
#' @param beta BetaHandler object, character path to beta file, or beta values matrix.
#'   If a character path or matrix is provided, a BetaHandler will be created automatically.
#' @param pheno Data frame or character path to phenotype file. Sample information with rownames matching beta column names (required).
#' @param genome Character. Genome version (default: "hg38").
#' @param array Character. Array platform type. Must be NULL if input is not array-based. Ignored if beta_locs is provided. (default: "450K")
#' @param beta_locs Data frame. Genomic locations sorted by position (optional).
#' @param sample_group_col Character. Column in pheno for sample grouping (default: "Sample_Group").
#' @param extend_by_dmr_size_ratio Numeric. Ratio of the DMR width to extend the plot region outside of the DMR in both sides (default: 0.2).
#' @param min_extension_bp Integer. Minimum extension in base pairs (default: 50).
#' @param max_cpgs Integer. Maximum number of CpGs to show in heatmap (default: 100).
#' @param max_samples_per_group Integer. Maximum number of samples to show per group in heatmap (default: 10).
#' @param plot_motif Logical. Whether to plot the sequence logo motif (default: TRUE).
#' @param motif_cpg_flank_size Integer. Number of base pairs to include as flanking regions around each CpG site for motif extraction (default: 5).
#' @param plot_title Logical. Whether to display the title on the plot. If FALSE, the title is shown in the logs (default: TRUE).
#' @param output_file Character. If provided, saves the plot to the specified file path (PDF format).
#'
#' @return A combined plot object (gridExtra) containing the DMR structure plot, beta values heatmap (if beta is provided),
#'   and sequence logo motif plot (if motif information is available and plot_motif is TRUE).
#'
#' @examples
#' # Using BetaHandler
#' beta_handler <- getBetaHandler(beta = "beta.txt", array = "450K", genome = "hg38")
#' plotDMR(dmrs, 1, beta = beta_handler, pheno = pheno_df)
#'
#' # Using a file path (handler created automatically)
#' plotDMR(dmrs, 1, beta = "beta.txt", pheno = pheno_df)
#'
#' # Using a beta matrix
#' plotDMR(dmrs, 1, beta = beta_matrix, pheno = pheno_df)
#'
#' # With custom flank size for motif extraction
#' plotDMR(dmrs, 1, beta = beta_matrix, pheno = pheno_df, motif_cpg_flank_size = 10)
#'
#' # Without motif plot
#' plotDMR(dmrs, 1, beta = beta_matrix, pheno = pheno_df, plot_motif = FALSE)
#'
#' @export
plotDMR <- function(dmrs,
                    dmr_index,
                    beta = NULL,
                    pheno = NULL,
                    genome = "hg38",
                    array = c("450K", "27K", "EPIC", "EPICv2"),
                    beta_locs = NULL,
                    sample_group_col = "Sample_Group",
                    extend_by_dmr_size_ratio = 0.2,
                    min_extension_bp = 50,
                    max_cpgs = 100,
                    max_samples_per_group = 10,
                    plot_motif = TRUE,
                    motif_cpg_flank_size = 5,
                    plot_title = TRUE,
                    output_file = NULL,
                    width = 8,
                    height = 12) {
    if (!is.null(output_file)) {
        if (!grepl("\\.pdf$", output_file, ignore.case = TRUE)) {
            stop("Output file must have a .pdf extension.")
        }
    }
    dmrs <- convertToGRanges(dmrs, genome)
    if (!is.null(array)) {
        if (length(array) > 1) {
            array <- array[[1]]
        }
        array <- strex::match_arg(array, ignore_case = TRUE)
    }

    old_showtext_opts <- showtext::showtext_opts()
    old_showtext_auto <- .isShowtextAutoEnabled()
    on.exit({
        showtext::showtext_opts(old_showtext_opts)
        if (!old_showtext_auto) {
            showtext::showtext_auto(enable = FALSE)
        }
    }, add = TRUE)
    showtext::showtext_auto(enable = TRUE)
    showtext::showtext_opts(dpi = if (!is.null(output_file) || .Device == "null device") 300 else .resolveShowtextDpi(300))
    if (!is.null(output_file)) {
        grDevices::cairo_pdf(output_file, width = width, height = height)
    }
    if (.Device == "null device") {
        grDevices::cairo_pdf(width = width, height = height)
    }

    # Create BetaHandler if a file path or matrix was provided
    if (!is.null(beta)) {
        if ((is.character(beta) && length(beta) == 1 && file.exists(beta)) || is.matrix(beta) || is.data.frame(beta)) {
            beta_handler <- getBetaHandler(
                beta = beta,
                array = array,
                sorted_locs = beta_locs,
                genome = genome
            )
        } else if (!"BetaHandler" %in% class(beta)) {
            stop("beta_handler must be either a file path (character) or a BetaHandler object")
        } else {
            beta_handler <- beta
        }
        beta_locs <- beta_handler$getBetaLocs()
    }

    if (!is.null(pheno) && !(sample_group_col %in% colnames(pheno))) {
        stop(sprintf("sample_group_col '%s' not found in pheno data frame", sample_group_col))
    }
    if (dmr_index < 1 || dmr_index > length(dmrs)) {
        stop("dmr_index (", dmr_index, ") is out of bounds. There are only ", length(dmrs), " DMRs available.")
    }
    dmr <- dmrs[dmr_index]
    dmr_data <- S4Vectors::mcols(dmr)

    if (!is.null(pheno)) {
        .log_info("Processing phenotype data...", level = 3)
        if (is.character(pheno) && length(pheno) == 1 && file.exists(pheno)) {
            pheno <- read.table(pheno, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
        }
        if (!is.data.frame(pheno)) {
            stop("pheno must be a data frame or a valid file path to a phenotype table")
        }
        pheno <- pheno[rownames(pheno) %in% beta_handler$getBetaColNames(), , drop = FALSE]
        if (nrow(pheno) == 0) {
            stop("No samples in pheno match the samples in beta values")
        }
    }

    .log_info(sprintf("Generating structure plot...", dmr_index), level = 3)

    # Create structure plot
    ret <- .plotDMRStructure(
        dmrs = dmrs,
        dmr_index = dmr_index,
        array = array,
        genome = genome,
        beta_locs = beta_locs,
        plot_title = plot_title,
        extend_by_dmr_size_ratio = extend_by_dmr_size_ratio,
        min_extension_bp = min_extension_bp,
        .ret_details = TRUE
    )
    structure_plot <- ret$structure_plot
    breaks <- ret$breaks
    breaks_labels <- ret$breaks_labels
    total_shown_positions <- ret$total_locs


    if (!is.null(beta)) {
        # Read beta values using BetaHandler
        .log_info("Fetching beta values for heatmap...", level = 3)
        beta_data <- beta_handler$getBeta(
            row_names = rownames(total_shown_positions),
            col_names = if (is.null(pheno)) NULL else rownames(pheno)
        )
        .log_info("Generating beta values heatmap...", level = 3)
        heatmap_plot <- .plotBetaHeatmap(
            dmr_data = dmr_data,
            beta_data = beta_data,
            pheno = pheno,
            max_samples_per_group = max_samples_per_group,
            sample_group_col = sample_group_col,
            total_shown_positions = total_shown_positions
        )
        structure_plot <- structure_plot +
            ggplot2::labs(x = NULL) +
            ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                plot.margin = ggplot2::margin(t = 5, r = 5, b = 1, l = 5)
            )
        heatmap_plot <- heatmap_plot +
            ggplot2::scale_x_continuous(
                breaks = breaks,
                labels = breaks_labels
            ) +
            ggplot2::coord_cartesian(xlim = c(breaks[1], breaks[length(breaks)])) +
            ggplot2::labs(
                x = sprintf("Genomic Position on %s (bp)", ret$chr)
            ) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
                axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
                plot.margin = ggplot2::margin(t = 0, r = 5, b = 5, l = 5)
            )

        g1 <- ggplot2::ggplotGrob(structure_plot)
        g2 <- ggplot2::ggplotGrob(heatmap_plot)
        grobs <- list(g1, g2)
    } else {
        structure_plot <- structure_plot + ggplot2::labs(
            x = sprintf("Genomic Position on %s (bp)", ret$chr)
        )
        .log_info("Beta values not provided. Only the DMR structure plot will be returned.", level = 2)
        grobs <- list(ggplot2::ggplotGrob(structure_plot))
    }

    if (plot_motif) {
        .log_info("Generating motif PWM plot...", level = 3)
        pwm_plot <- .plotPWM(dmr, genome = genome, array = array, beta_locs = beta_locs, motif_cpg_flank_size = motif_cpg_flank_size)
        if (!is.null(pwm_plot)) {
            motif_context <- .summarizeMotifContext(
                dmrs = dmrs,
                dmr_index = dmr_index,
                genome = genome,
                array = array,
                beta_locs = beta_locs,
                motif_cpg_flank_size = motif_cpg_flank_size
            )
            motif_lines <- .buildMotifContextLines(motif_context)
            motif_context_plot <- .plotMotifContext(motif_lines)
            motif_panel <- gridExtra::arrangeGrob(
                ggplot2::ggplotGrob(pwm_plot),
                ggplot2::ggplotGrob(motif_context_plot),
                ncol = 2,
                widths = c(0.58, 0.42)
            )
            grobs <- c(grobs, list(motif_panel))
        }
    }

    # Keep genomic schematic and beta heatmap perfectly aligned on x,
    # but avoid forcing motif/context panel into the same gtable geometry.
    if (!is.null(beta) && length(grobs) >= 2) {
        max_width <- grid::unit.pmax(grobs[[1]]$widths, grobs[[2]]$widths)
        grobs[[1]]$widths <- max_width
        grobs[[2]]$widths <- max_width
    }

    if (length(grobs) == 3) {
        combined <- gridExtra::arrangeGrob(
            grobs = grobs,
            ncol = 1,
            heights = c(1, 1, 0.55)
        )
    } else if (length(grobs) == 2 && plot_motif && is.null(beta)) {
        combined <- gridExtra::arrangeGrob(
            grobs = grobs,
            ncol = 1,
            heights = c(1, 0.6)
        )
    } else {
        combined <- do.call(gridExtra::gtable_rbind, grobs)
    }
    grid::grid.draw(combined)
    if (!is.null(output_file)) {
        grDevices::dev.off()
    }
    invisible(combined)
}

.classifyDMRRegionForManhattan <- function(dmrs, promoter_col = "in_promoter_of", gene_body_col = "in_gene_body_of") {
    mcols_df <- S4Vectors::mcols(dmrs)
    has_promoter <- promoter_col %in% colnames(mcols_df)
    has_gene_body <- gene_body_col %in% colnames(mcols_df)

    promoter_mask <- if (has_promoter) .hasNonEmptyString(mcols_df[[promoter_col]]) else rep(FALSE, length(dmrs))
    gene_body_mask <- if (has_gene_body) .hasNonEmptyString(mcols_df[[gene_body_col]]) else rep(FALSE, length(dmrs))

    region_class <- rep("Intergenic", length(dmrs))
    region_class[gene_body_mask] <- "Gene Body"
    region_class[promoter_mask] <- "Promoter"
    region_class[promoter_mask & gene_body_mask] <- NA

    factor(region_class, levels = c("Promoter", "Gene Body", "Intergenic"))
}

.buildManhattanBlockRects <- function(dmr_df, block_col = "block_id") {
    if (!(block_col %in% colnames(dmr_df))) {
        return(data.frame())
    }
    block_vals <- as.character(dmr_df[[block_col]])
    keep <- !is.na(block_vals) & nzchar(block_vals)
    if (!any(keep)) {
        return(data.frame())
    }
    block_df <- dmr_df[keep, c("position", "score", block_col), drop = FALSE]
    colnames(block_df)[3] <- "block_id"
    by_block <- split(block_df, block_df$block_id)
    if (length(by_block) == 0) {
        return(data.frame())
    }

    y_rng <- range(dmr_df$score, na.rm = TRUE)
    y_span <- diff(y_rng)
    y_pad <- max(1e-6, y_span * 0.03)
    x_unique <- sort(unique(dmr_df$position))
    if (length(x_unique) >= 2L) {
        x_pad <- max(1, min(diff(x_unique), na.rm = TRUE) * 0.35)
    } else {
        x_pad <- 1
    }

    rects <- do.call(rbind, lapply(names(by_block), function(id) {
        d <- by_block[[id]]
        xmin <- min(d$position, na.rm = TRUE)
        xmax <- max(d$position, na.rm = TRUE)
        if (xmin == xmax) {
            xmin <- xmin - x_pad
            xmax <- xmax + x_pad
        }
        data.frame(
            block_id = id,
            xmin = xmin,
            xmax = xmax,
            ymin = min(d$score, na.rm = TRUE) - y_pad,
            ymax = max(d$score, na.rm = TRUE) + y_pad,
            stringsAsFactors = FALSE
        )
    }))
    rects <- rects[order(rects$xmin, rects$xmax), , drop = FALSE]
    rects$fill <- colorspace::qualitative_hcl(max(3L, nrow(rects)), palette = "Set 2")[seq_len(nrow(rects))]
    rownames(rects) <- NULL
    rects
}

#' Plot Intermediate DMR Block Formation Diagnostics
#'
#' @description Visualizes how DMR blocks are formed on a single chromosome by
#' layering raw scores, smoothed scores, piecewise-linear segments, slope-based
#' candidate blocks, distance-based split boundaries, and final accepted blocks.
#'
#' @param dmrs GRanges object or data frame. DMR results from `findDMRsFromSeeds`
#' or `scoreDMRs`.
#' @param chromosome Character. Chromosome to inspect (e.g., `"chr7"` or `"7"`).
#' @param genome Character. Genome version passed to `convertToGRanges`
#' (default: `"hg38"`).
#' @param k_neighbors Integer. Number of nearest neighbors used in adaptive Gaussian
#' smoothing (default: `5`).
#' @param min_segment_size Integer. Minimum size of linear segments in PELT
#' segmentation (default: `2`).
#' @param block_gap_mode Character. Gap rule for block splitting:
#' `"adaptive"` (default), `"fixed"`, or `"none"`.
#' @param block_gap_fixed_bp Numeric. Maximum allowed midpoint gap (bp) when
#' `block_gap_mode = "fixed"`. Ignored otherwise.
#' @param block_gap_quantile Numeric in `(0, 1)`. Gap quantile used for adaptive
#' thresholding (default: `0.95`).
#' @param block_gap_multiplier Numeric > 0. Multiplier for adaptive gap threshold
#' (default: `1.5`).
#' @param block_gap_min_bp Numeric >= 0. Lower clamp for adaptive gap threshold (bp).
#' Default is `250000`.
#' @param block_gap_max_bp Numeric >= `block_gap_min_bp`. Upper clamp for adaptive
#' gap threshold (bp). Default is `5000000`.
#' @param point_alpha Numeric. Alpha for raw score points in `[0, 1]`
#' (default: `0.7`).
#' @param point_size Numeric. Point size for raw scores (default: `1.0`).
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' dmrs <- readRDS("dmrs.rds")
#' p <- plotDMRBlockFormation(dmrs, chromosome = "chr7")
#' print(p)
#' }
#' @export
plotDMRBlockFormation <- function(dmrs,
                                  chromosome,
                                  genome = "hg38",
                                  k_neighbors = 5L,
                                  min_segment_size = 2L,
                                  block_gap_mode = "adaptive",
                                  block_gap_fixed_bp = NULL,
                                  block_gap_quantile = 0.95,
                                  block_gap_multiplier = 1.5,
                                  block_gap_min_bp = 250000,
                                  block_gap_max_bp = 5000000,
                                  point_alpha = 0.7,
                                  point_size = 1.0) {
    dmrs <- convertToGRanges(dmrs, genome = genome)
    if (!( "score" %in% colnames(S4Vectors::mcols(dmrs)))) {
        stop("Column 'score' not found in DMR metadata.")
    }
    if (!is.numeric(point_alpha) || length(point_alpha) != 1 || is.na(point_alpha) || point_alpha < 0 || point_alpha > 1) {
        stop("point_alpha must be a single numeric value in [0, 1].")
    }
    if (!is.numeric(point_size) || length(point_size) != 1 || is.na(point_size) || point_size <= 0) {
        stop("point_size must be a single positive numeric value.")
    }

    chromosome <- as.character(chromosome)[1]
    if (!grepl("^chr", chromosome)) {
        chromosome <- paste0("chr", chromosome)
    }
    chr_values <- as.character(GenomicRanges::seqnames(dmrs))
    if (!(chromosome %in% chr_values)) {
        stop(
            "Chromosome '", chromosome, "' not found. Available chromosomes include: ",
            paste(head(unique(chr_values), 8), collapse = ","),
            if (length(unique(chr_values)) > 8) " ..." else ""
        )
    }
    score_values <- suppressWarnings(as.numeric(S4Vectors::mcols(dmrs)[["score"]]))
    keep <- chr_values == chromosome & is.finite(score_values)
    if (sum(keep) < 3L) {
        stop("At least 3 finite DMR scores are required on ", chromosome, " to build diagnostics.")
    }

    chr_dmrs <- dmrs[keep]

    details_result <- .assignDMRBlocksFromScores(
        dmrs = chr_dmrs,
        k_neighbors = k_neighbors,
        min_segment_size = min_segment_size,
        block_gap_mode = block_gap_mode,
        block_gap_fixed_bp = block_gap_fixed_bp,
        block_gap_quantile = block_gap_quantile,
        block_gap_multiplier = block_gap_multiplier,
        block_gap_min_bp = block_gap_min_bp,
        block_gap_max_bp = block_gap_max_bp,
        return_details = TRUE
    )

    details <- details_result$details[[chromosome]]
    if (is.null(details)) {
        stop("No block-formation details were produced for chromosome ", chromosome, ".")
    }

    dmr_df <- details$dmr_df
    seg_df <- details$segments_df
    candidate_df <- details$candidates_df
    split_df <- details$split_events_df
    blocks_df <- details$blocks_df

    p <- ggplot2::ggplot(dmr_df, ggplot2::aes(x = midpoint, y = score_raw)) +
        ggplot2::geom_point(color = "#4D4D4D", alpha = point_alpha, size = point_size) +
        ggplot2::geom_line(ggplot2::aes(y = score_smoothed), color = "#1F78B4", linewidth = 0.6)

    if (nrow(candidate_df) > 0) {
        p <- p + ggplot2::geom_rect(
            data = candidate_df,
            ggplot2::aes(xmin = start_bp, xmax = end_bp, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "#FDB863",
            alpha = 0.08,
            color = NA
        )
    }

    if (nrow(blocks_df) > 0) {
        p <- p + ggplot2::geom_rect(
            data = blocks_df,
            ggplot2::aes(xmin = start_bp, xmax = end_bp, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "#B2DF8A",
            color = "#33A02C",
            linewidth = 0.35,
            alpha = 0.14
        )
    }

    if (nrow(seg_df) > 1) {
        p <- p + ggplot2::geom_vline(
            data = seg_df[-nrow(seg_df), , drop = FALSE],
            ggplot2::aes(xintercept = end_bp),
            inherit.aes = FALSE,
            color = "#6A3D9A",
            linetype = "dotted",
            linewidth = 0.25,
            alpha = 0.8
        )
    }

    if (nrow(split_df) > 0) {
        p <- p + ggplot2::geom_vline(
            data = split_df,
            ggplot2::aes(xintercept = right_bp),
            inherit.aes = FALSE,
            color = "#E31A1C",
            linetype = "dashed",
            linewidth = 0.35,
            alpha = 0.9
        )
    }

    subtitle <- paste0(
        "gap_mode=", block_gap_mode,
        ", gap_threshold=", format(round(details$gap_threshold_bp), big.mark = ",", scientific = FALSE), " bp",
        ", slope_d=", signif(details$slope_threshold, digits = 3),
        ", candidates=", nrow(candidate_df),
        ", blocks=", nrow(blocks_df)
    )

    p <- p +
        ggplot2::labs(
            title = paste0("DMR Block Formation Diagnostics (", chromosome, ")"),
            subtitle = subtitle,
            x = "Genomic Midpoint (bp)",
            y = "score"
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 9)
        ) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.05)))

    p
}

#' Plot Manhattan-Style View of DMR Scores
#'
#' @description Creates a Manhattan-style genome-wide scatter plot of DMR scores.
#' DMRs are colored by their dominant genomic region class inferred from DMR annotations.
#'
#' @param dmrs GRanges object or data frame. DMR results from findDMRsFromSeeds.
#' @param genome Character. Genome version passed to `convertToGRanges` (default: `"hg38"`).
#' @param promoter_col Character. Metadata column indicating promoter overlap (default: `"in_promoter_of"`).
#' @param gene_body_col Character. Metadata column indicating gene-body overlap (default: `"in_gene_body_of"`).
#' @param point_size Numeric. Point size (default: `1.1`).
#' @param point_alpha Numeric. Point alpha in `[0, 1]` (default: `0.75`).
#' @param block_col Character. Metadata column containing block IDs (default: `"block_id"`).
#' @param show_blocks Logical. If TRUE, draw translucent rectangles for identified blocks (default: `TRUE`).
#' @param block_alpha Numeric. Alpha for block rectangles in `[0, 1]` (default: `0.12`).
#' @param block_linewidth Numeric. Line width for block rectangle borders (default: `0.25`).
#' @param output_file Character or NULL. If non-NULL, path to save the plot as a PDF (default: `NULL`).
#' @param width Numeric. Width of the output PDF in inches (default: `12`).
#' @param height Numeric. Height of the output PDF in inches (default: `6`).
#' 
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' dmrs <- readRDS("dmrs.rds")
#' p <- plotDMRsManhattan(dmrs)
#' print(p)
#' }
#' @export
plotDMRsManhattan <- function(dmrs,
                              genome = "hg38",
                              promoter_col = "in_promoter_of",
                              gene_body_col = "in_gene_body_of",
                              point_size = 1.1,
                              point_alpha = 0.75,
                              block_col = "block_id",
                              show_blocks = TRUE,
                              block_alpha = 0.12,
                              block_linewidth = 0.25,
                              output_file = NULL,
                              width = 12,
                              height = 6) {
    dmrs <- convertToGRanges(dmrs, genome = genome)
    if (!( "score" %in% colnames(S4Vectors::mcols(dmrs)))) {
        stop("Column 'score' not found in DMR metadata.")
    }
    if (!is.numeric(point_size) || length(point_size) != 1 || is.na(point_size) || point_size <= 0) {
        stop("point_size must be a single positive numeric value.")
    }
    if (!is.numeric(point_alpha) || length(point_alpha) != 1 || is.na(point_alpha) || point_alpha < 0 || point_alpha > 1) {
        stop("point_alpha must be a single numeric value in [0, 1].")
    }
    if (!is.logical(show_blocks) || length(show_blocks) != 1 || is.na(show_blocks)) {
        stop("show_blocks must be TRUE or FALSE.")
    }
    if (!is.numeric(block_alpha) || length(block_alpha) != 1 || is.na(block_alpha) || block_alpha < 0 || block_alpha > 1) {
        stop("block_alpha must be a single numeric value in [0, 1].")
    }
    if (!is.numeric(block_linewidth) || length(block_linewidth) != 1 || is.na(block_linewidth) || block_linewidth < 0) {
        stop("block_linewidth must be a single non-negative numeric value.")
    }

    scores <- suppressWarnings(as.numeric(S4Vectors::mcols(dmrs)[["score"]]))
    keep <- is.finite(scores)
    if (!any(keep)) {
        stop("No finite score values found in column 'score'.")
    }
    dmrs <- dmrs[keep]
    scores <- scores[keep]

    region_class <- .classifyDMRRegionForManhattan(dmrs, promoter_col = promoter_col, gene_body_col = gene_body_col)
    m <- !is.na(region_class)
    dmrs <- dmrs[m]
    scores <- scores[m]
    region_class <- region_class[m]
    chr <- as.character(GenomicRanges::seqnames(dmrs))
    chr_levels <- .orderChromosomesNaturally(chr)

    dmr_df <- data.frame(
        chr = chr,
        start = GenomicRanges::start(dmrs),
        end = GenomicRanges::end(dmrs),
        score = scores,
        region_class = region_class,
        stringsAsFactors = FALSE
    )
    if (block_col %in% colnames(S4Vectors::mcols(dmrs))) {
        dmr_df$block_id <- as.character(S4Vectors::mcols(dmrs)[[block_col]])
    } else {
        dmr_df$block_id <- NA_character_
    }
    dmr_df$chr <- factor(dmr_df$chr, levels = chr_levels)
    dmr_df$midpoint <- floor((dmr_df$start + dmr_df$end) / 2)

    chr_lengths <- tapply(dmr_df$end, dmr_df$chr, max, na.rm = TRUE)
    chr_lengths <- as.numeric(chr_lengths)
    names(chr_lengths) <- levels(dmr_df$chr)
    chr_lengths[!is.finite(chr_lengths)] <- 0

    chr_offsets <- c(0, cumsum(chr_lengths)[seq_len(max(1L, length(chr_lengths) - 1L))])
    chr_offsets <- chr_offsets[seq_along(chr_lengths)]
    names(chr_offsets) <- names(chr_lengths)

    dmr_df$position <- dmr_df$midpoint + chr_offsets[as.character(dmr_df$chr)]
    axis_df <- stats::aggregate(position ~ chr, data = dmr_df, FUN = function(x) mean(range(x)))
    axis_df <- axis_df[match(chr_levels, as.character(axis_df$chr)), , drop = FALSE]
    axis_df <- axis_df[is.finite(axis_df$position), , drop = FALSE]

    chrom_boundaries <- cumsum(chr_lengths)
    chrom_boundaries <- chrom_boundaries[is.finite(chrom_boundaries)]
    chrom_boundaries <- chrom_boundaries[seq_len(max(0L, length(chrom_boundaries) - 1L))]

    region_colors <- c(
        "Promoter" = "#D73027",
        "Gene Body" = "#3182BD",
        "Intergenic" = "#9E9E9E"
    )
    region_shapes <- c(
        "Promoter" = 16,
        "Gene Body" = 17,
        "Intergenic" = 15
    )

    p <- ggplot2::ggplot(dmr_df, ggplot2::aes(x = position, y = score, color = region_class, shape = region_class))
    block_rects <- data.frame()
    if (show_blocks) {
        block_rects <- .buildManhattanBlockRects(dmr_df, block_col = "block_id")
        if (nrow(block_rects) > 0) {
            p <- p + ggplot2::geom_rect(
                data = block_rects,
                ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                fill = block_rects$fill,
                color = block_rects$fill,
                linewidth = block_linewidth,
                alpha = block_alpha,
                show.legend = FALSE
            )
        }
    }
    subtitle <- NULL
    if (nrow(block_rects) > 0) {
        subtitle <- paste0("Identified blocks: ", nrow(block_rects))
    }
    p <- p +
        ggplot2::geom_point(size = point_size, alpha = point_alpha, stroke = 0) +
        ggplot2::scale_color_manual(values = region_colors, drop = TRUE, name = "Primary Region") +
        ggplot2::scale_shape_manual(values = region_shapes, drop = TRUE, name = "Primary Region") +
        ggplot2::scale_x_continuous(
            breaks = axis_df$position,
            labels = as.character(axis_df$chr),
            expand = ggplot2::expansion(mult = c(0.01, 0.01))
        ) +
        ggplot2::labs(
            x = "Chromosome",
            y = "score",
            title = paste0("DMR Manhattan Plot"),
            subtitle = subtitle
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5),
            legend.title = ggplot2::element_text(face = "bold"),
            legend.position = "top"
        )
    if (length(chrom_boundaries) > 0) {
        p <- p + ggplot2::geom_vline(xintercept = chrom_boundaries, color = "#DADADA", linewidth = 0.3)
    }
    if (!is.null(output_file)) {
        ggplot2::ggsave(filename = output_file, plot = p, width = width, height = height, units = "in", dpi = 300)
    }
    p
}
