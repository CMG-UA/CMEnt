# Avoid NSE warnings from R CMD check for ggplot2 aes()
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("Sample", "Beta", "Position", "x", "xend", "y", "yend", "start", "position", "score", "region_class", "chr"))
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
#' @param genome Character. Genome version (default: "hg19"). Ignored if beta_locs is provided.
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
                              genome = "hg19",
                              beta_locs = NULL,
                              extend_by_dmr_size_ratio = 0.2,
                              min_extension_bp = 50,
                              plot_title = TRUE,
                              .ret_details = FALSE) {
    showtext::showtext_auto()


    # Validate input
    if (dmr_index < 1 || dmr_index > length(dmrs)) {
        stop("dmr_index must be between 1 and ", length(dmrs))
    }

    # Extract the specific DMR
    dmr <- dmrs[dmr_index]
    dmr_data <- S4Vectors::mcols(dmr)

    # Get genomic locations if not provided
    use_abs <- FALSE
    if (is.null(beta_locs)) {
        array <- strex::match_arg(array, ignore_case = TRUE)
        beta_locs <- getSortedGenomicLocs(array = array, genome = genome)
        use_abs <- TRUE
    }

    # Extract DMR information

    supporting_sites <- .getSupportingSitesFromColumns(
        dmr,
        max_sup_cpgs_per_dmr_side = NULL,
        separate_by_section = TRUE,
        use_absolute_indices = use_abs,
        beta_locs = beta_locs
    )
    upstream_sup_cpgs_inds <- supporting_sites[[1]]$upstream
    downstream_sup_cpgs_inds <- supporting_sites[[1]]$downstream
    upstream_sup_cpgs <- beta_locs[upstream_sup_cpgs_inds, , drop = FALSE]
    downstream_sup_cpgs <- beta_locs[downstream_sup_cpgs_inds, , drop = FALSE]
    seeds_inds <- supporting_sites[[1]]$seeds
    start_seed_ind <- min(seeds_inds)
    end_seed_ind <- max(seeds_inds)
    if (length(upstream_sup_cpgs_inds) == 0) {
        start_cpg_ind <- start_seed_ind
    } else {
        start_cpg_ind <- min(upstream_sup_cpgs_inds)
    }
    if (length(downstream_sup_cpgs_inds) == 0) {
        end_cpg_ind <- end_seed_ind
    } else {
        end_cpg_ind <- max(downstream_sup_cpgs_inds)
    }

    chr <- as.character(GenomicRanges::seqnames(dmr))
    dmr_start <- GenomicRanges::start(dmr)
    dmr_end <- GenomicRanges::end(dmr)

    # Extract seed IDs from the comma-separated string


    # Get positions

    start_cpg_pos <- beta_locs[start_cpg_ind, "start"]
    end_cpg_pos <- beta_locs[end_cpg_ind, "start"]
    seed_positions <- beta_locs[seeds_inds, "start"]
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
    nsup_cpgs <- beta_locs[start_seed_ind:end_seed_ind, ]
    nsup_cpgs <- nsup_cpgs[which(nsup_cpgs$start < dmr_end & nsup_cpgs$start > dmr_start & !seq(start_seed_ind, end_seed_ind) %in% seeds_inds), , drop = FALSE]

    downstream_nsup_cpgs <- beta_locs[which(beta_locs$start > dmr_end & beta_locs$start <= plot_end), , drop = FALSE]
    upstream_nsup_cpgs <- beta_locs[which(beta_locs$start < dmr_start & beta_locs$start >= plot_start), , drop = FALSE]
    extended_nsup_cpgs <- rbind(
        nsup_cpgs,
        upstream_nsup_cpgs,
        downstream_nsup_cpgs
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

    extended_sup_cpgs <- rbind(
        upstream_sup_cpgs,
        downstream_sup_cpgs
    )
    if (nrow(extended_sup_cpgs) > 0) {
        extended_sup_cpgs_df <- data.frame(
            start = extended_sup_cpgs$start,
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
    if (nrow(extended_nsup_cpgs) > 0) {
        extended_nsup_cpgs_df <- data.frame(
            start = extended_nsup_cpgs$start,
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
        linewidth = 1.5,
        alpha = 0.8
    )

    # Plot seed stems
    p <- p + ggplot2::geom_segment(
        data = seeds_df,
        ggplot2::aes(x = start, xend = start, y = 0, yend = y),
        color = "#377EB8",
        linewidth = 1.2,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"), type = "closed")
    )

    # Plot supported CpG stems
    if (nrow(extended_sup_cpgs_df) > 0) {
        p <- p + ggplot2::geom_segment(
            data = extended_sup_cpgs_df,
            ggplot2::aes(x = start, xend = start, y = 0, yend = y),
            color = "#377EB8",
            linewidth = 0.8,
            alpha = 0.8
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

    # Plot seed points
    p <- p + ggplot2::geom_point(
        data = seeds_df,
        ggplot2::aes(x = start, y = y),
        color = "#377EB8",
        size = 3,
        shape = 16
    )

    # Plot supporting CpG points
    if (nrow(extended_sup_cpgs_df) > 0) {
        p <- p + ggplot2::geom_point(
            data = extended_sup_cpgs_df,
            ggplot2::aes(x = start, y = y),
            color = "#377EB8",
            size = 2,
            shape = 16,
            alpha = 0.8
        )
    }

    # Plot non-supporting CpG points
    if (nrow(extended_nsup_cpgs_df) > 0) {
        p <- p + ggplot2::geom_point(
            data = extended_nsup_cpgs_df,
            ggplot2::aes(x = start, y = y),
            color = "gray70",
            size = 2,
            shape = 16,
            alpha = 0.5
        )
    }

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
    if (nrow(upstream_sup_cpgs) > 0) {
        p <- p + ggplot2::geom_segment(
            data = dmr_upstream_line,
            ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
            color = "#E41A1C",
            linewidth = 1.5,
            alpha = 0.8
        ) + ggplot2::geom_text(
            data = dmr_upstream_line,
            ggplot2::aes(
                x = (x + xend) / 2, y = yend,
                label = "Upstream Extension"
            ),
            vjust = -0.5,
            color = "#000000",
            size = 3
        )
        p <- p + ggplot2::annotate(
            "polygon",
            x = c(min(upstream_sup_cpgs$start), start_seed_pos, start_seed_pos, min(upstream_sup_cpgs$start)),
            y = c(0, 0, 1, 0.5),
            alpha = 0.1,
            fill = "#E41A1C"
        )
    }
    # if downstream extended CpGs exist add shading in the form of a trapezoid
    if (nrow(downstream_sup_cpgs) > 0) {
        p <- p + ggplot2::geom_segment(
            data = dmr_downstream_line,
            ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
            color = "#E41A1C",
            linewidth = 1.5,
            alpha = 0.8
        ) + ggplot2::geom_text(
            data = dmr_downstream_line,
            ggplot2::aes(x = (x + xend) / 2, y = y, label = "Downstream Extension"),
            vjust = -0.5,
            color = "#000000",
            size = 3
        )
        p <- p + ggplot2::annotate(
            "polygon",
            x = c(end_seed_pos, max(downstream_sup_cpgs$start), max(downstream_sup_cpgs$start), end_seed_pos),
            y = c(0, 0, 0.5, 1),
            alpha = 0.1,
            fill = "#E41A1C"
        )
    }

    # Create title if not provided
    title <- sprintf(
        "DMR #%d: %s:%s-%s\n%d seeds (\u0394\u03b2=%.3f)",
        dmr_index,
        chr,
        format(dmr_start, big.mark = ",", scientific = FALSE),
        format(dmr_end, big.mark = ",", scientific = FALSE),
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
    if (nrow(extended_nsup_cpgs) > 0 || nrow(extended_sup_cpgs) > 0) {
        p <- p +
            ggplot2::scale_y_continuous(
                breaks = c(0.5, 1),
                labels = c("Array CpGs", "seeds"),
                limits = c(-0.1, 1.15)
            )
    } else {
        p <- p +
            ggplot2::scale_y_continuous(
                breaks = c(1),
                labels = c("seeds"),
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

    p <- p + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 0.5)
    )

    # Add ticks for the seeds on the x-axis
    breaks <- c(plot_start, beta_locs[start_cpg_ind:end_cpg_ind, "start"], plot_end)
    cpgs_labs <- paste0(beta_locs[start_cpg_ind:end_cpg_ind, "start"], "\n(", rownames(beta_locs[start_cpg_ind:end_cpg_ind, , drop = FALSE]), ")")
    breaks_labels <- c(
        as.character(plot_start), cpgs_labs,
        as.character(plot_end)
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
        total_shown_positions <- rbind(extended_nsup_cpgs, extended_sup_cpgs, beta_locs[seeds_inds, , drop = FALSE])
        total_shown_positions <- total_shown_positions[order(total_shown_positions$start), ]
        return(invisible(list(structure_plot = p, breaks = breaks, breaks_labels = breaks_labels, chr = chr, total_locs = total_shown_positions)))
    }
    invisible(p)
}


# Create beta heatmap plot
.plotBetaHeatmap <- function(dmr_data, beta_data, total_shown_positions, pheno = NULL, sample_group_col = "Sample_Group") {
    cpg_ids <- rownames(total_shown_positions)
    cpg_locs <- total_shown_positions[, c("chr", "start")]


    # Mark seeds
    seed_ids <- unlist(strsplit(as.character(dmr_data$seeds), ","))
    is_seed <- cpg_ids %in% seed_ids

    # Create heatmap
    # Prepare data
    beta_data <- as.data.frame(beta_data)
    beta_data[, "CpG"] <- rownames(beta_data)
    beta_melted <- suppressWarnings(suppressMessages(reshape2::melt(beta_data, id_vars = "CpG")))
    colnames(beta_melted) <- c("CpG", "Sample", "Beta")
    beta_melted$Position <- cpg_locs[as.character(beta_melted$CpG), "start"]
    beta_melted$is_seed <- is_seed[match(beta_melted$CpG, cpg_ids)]
    if (!is.null(pheno) && !is.null(sample_group_col)) {
        beta_melted$Group <- pheno[as.character(beta_melted$Sample), sample_group_col]
        # Order samples by group
        sample_order <- rownames(pheno)[order(pheno[[sample_group_col]])]
        beta_melted$Sample <- factor(beta_melted$Sample, levels = sample_order)
    }
    heatmap_plot <- ggplot2::ggplot(beta_melted) +
        ggplot2::geom_tile(ggplot2::aes(x = Position, y = Sample, fill = Beta)) +
        ggplot2::scale_fill_gradient(
            low = "white", high = "red",
            limits = c(min(beta_melted$Beta), max(beta_melted$Beta)),
            name = "\u03b2-values"
        ) +
        ggplot2::labs(
            y = "Sample"
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7)
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 0.5, vjust = 0.5)
        )
    # Add vertical lines for seeds
    if (any(is_seed)) {
        seed_positions <- unique(beta_melted$Position[beta_melted$is_seed])
        heatmap_plot <- heatmap_plot +
            ggplot2::geom_vline(
                xintercept = seed_positions,
                color = "yellow",
                linetype = "dashed",
                linewidth = 0.5,
                alpha = 0.7
            )
    }
    heatmap_plot
}


minmaxscale <- function(x) {
    (x - min(x)) / max(max(x) - min(x), 1e-10)
}

.summarizeMotifContext <- function(dmrs, dmr_index, genome, array, beta_locs, motif_flank_size = 5) {
    ret <- list(top_interactions = data.frame(), jaspar = data.frame())
    if (length(dmrs) < 2) {
        return(ret)
    }

    top_n <- as.integer(getOption("DMRsegal.plotDMR_top_motif_interactions", 5L))
    pool_size <- as.integer(getOption("DMRsegal.plotDMR_interaction_pool_size", 300L))
    top_n <- ifelse(is.na(top_n) || top_n < 1, 5L, top_n)
    pool_size <- ifelse(is.na(pool_size) || pool_size < 2, min(300L, length(dmrs)), pool_size)
    pool_size <- min(pool_size, length(dmrs))

    if ("rank" %in% colnames(S4Vectors::mcols(dmrs))) {
        ord <- order(S4Vectors::mcols(dmrs)$rank, na.last = TRUE)
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
            flank_size = motif_flank_size
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

    sim_matrix <- .extractMotifsSimilarity(cand_dmrs, flank_size = motif_flank_size)
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
        if ("rank" %in% colnames(S4Vectors::mcols(dmrs))) {
            top_tbl$rank <- S4Vectors::mcols(dmrs)$rank[top_tbl$dmr_index]
        }
        ret$top_interactions <- top_tbl
    }

    if (isTRUE(getOption("DMRsegal.plotDMR_query_jaspar", TRUE))) {
        ret$jaspar <- tryCatch(
            comparePWMToJaspar(list(target_pwm)),
            error = function(e) data.frame()
        )
    }
    ret
}

.buildMotifContextLines <- function(motif_context, top_n = 5) {
    lines <- c("Motif Context", "")
    top_tbl <- motif_context$top_interactions
    if (nrow(top_tbl) > 0) {
        lines <- c(lines, "Top motif interactions:")
        for (i in seq_len(min(top_n, nrow(top_tbl)))) {
            rank_txt <- if ("rank" %in% colnames(top_tbl) && !is.na(top_tbl$rank[i])) {
                paste0(" rank=", top_tbl$rank[i])
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
                rank_txt
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

.plotPWM <- function(dmr, genome, array, beta_locs, motif_flank_size = 5) {
    # Extract DMR motifs if not already present
    if (!"pwm" %in% colnames(S4Vectors::mcols(dmr))) {
        dmr <- extractDMRMotifs(dmr,
            genome = genome, array = array,
            beta_locs = beta_locs, flank_size = motif_flank_size
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

    position_labels <- c(seq(-motif_flank_size, 0), seq(0, motif_flank_size))

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
                y = "Information Content (bits)",
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
#' @param genome Character. Genome version (default: "hg19").
#' @param array Character. Array platform type (default: "450K"). Ignored if beta_locs is provided.
#' @param beta_locs Data frame. Genomic locations sorted by position (optional). If NULL, will be fetched based on array and genome.
#' @param ncol Integer. Number of columns in the grid (default: 1).
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
                     genome = "hg19",
                     array = c("450K", "27K", "EPIC", "EPICv2"),
                     beta_locs = NULL,
                     ncol = 1,
                     ...) {
    showtext::showtext_auto()
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
#' @param genome Character. Genome version (default: "hg19").
#' @param array Character. Array platform type. Must be NULL if input is not array-based. Ignored if beta_locs is provided. (default: "450K")
#' @param beta_locs Data frame. Genomic locations sorted by position (optional).
#' @param sample_group_col Character. Column in pheno for sample grouping (default: "Sample_Group").
#' @param extend_by_dmr_size_ratio Numeric. Ratio of the DMR width to extend the plot region outside of the DMR in both sides (default: 0.2).
#' @param min_extension_bp Integer. Minimum extension in base pairs (default: 50).
#' @param max_cpgs Integer. Maximum number of CpGs to show in heatmap (default: 100).
#' @param plot_motif Logical. Whether to plot the sequence logo motif (default: TRUE).
#' @param motif_flank_size Integer. Number of base pairs to include as flanking regions around each CpG site for motif extraction (default: 5).
#' @param plot_title Logical. Whether to display the title on the plot. If FALSE, the title is shown in the logs (default: TRUE).
#' @param output_file Character. If provided, saves the plot to the specified file path (PDF format).
#'
#' @return A combined plot object (gridExtra) containing the DMR structure plot, beta values heatmap (if beta is provided),
#'   and sequence logo motif plot (if motif information is available and plot_motif is TRUE).
#'
#' @examples
#' # Using BetaHandler
#' beta_handler <- getBetaHandler(beta = "beta.txt", array = "450K", genome = "hg19")
#' plotDMR(dmrs, 1, beta = beta_handler, pheno = pheno_df)
#'
#' # Using a file path (handler created automatically)
#' plotDMR(dmrs, 1, beta = "beta.txt", pheno = pheno_df)
#'
#' # Using a beta matrix
#' plotDMR(dmrs, 1, beta = beta_matrix, pheno = pheno_df)
#'
#' # With custom flank size for motif extraction
#' plotDMR(dmrs, 1, beta = beta_matrix, pheno = pheno_df, motif_flank_size = 10)
#'
#' # Without motif plot
#' plotDMR(dmrs, 1, beta = beta_matrix, pheno = pheno_df, plot_motif = FALSE)
#'
#' @export
plotDMR <- function(dmrs,
                    dmr_index,
                    beta = NULL,
                    pheno = NULL,
                    genome = "hg19",
                    array = c("450K", "27K", "EPIC", "EPICv2"),
                    beta_locs = NULL,
                    sample_group_col = "Sample_Group",
                    extend_by_dmr_size_ratio = 0.2,
                    min_extension_bp = 50,
                    max_cpgs = 100,
                    plot_motif = TRUE,
                    motif_flank_size = 5,
                    plot_title = TRUE,
                    output_file = NULL,
                    width = 8,
                    height = 12) {

    dmrs <- convertToGRanges(dmrs, genome)
    if (!is.null(array)) {
        if (length(array) > 1) {
            array <- array[[1]]
        }
        array <- strex::match_arg(array, ignore_case = TRUE)
    }

    showtext::showtext_auto(enable = TRUE)
    showtext::showtext_opts(dpi = 300)
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
            sample_group_col = sample_group_col,
            total_shown_positions = total_shown_positions
        )
        heatmap_plot <- heatmap_plot +
            ggplot2::scale_x_continuous(
                breaks = breaks,
                labels = breaks_labels
            ) +
            ggplot2::coord_cartesian(xlim = c(breaks[1], breaks[length(breaks)])) +
            ggplot2::labs(
                x = sprintf("Genomic Position on %s (bp)", ret$chr)
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
        pwm_plot <- .plotPWM(dmr, genome = genome, array = array, beta_locs = beta_locs, motif_flank_size = motif_flank_size)
        if (!is.null(pwm_plot)) {
            motif_context <- .summarizeMotifContext(
                dmrs = dmrs,
                dmr_index = dmr_index,
                genome = genome,
                array = array,
                beta_locs = beta_locs,
                motif_flank_size = motif_flank_size
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
            heights = c(1.1, 1, 0.55)
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

.filterDMRsForCircos <- function(dmrs, max_dmrs_per_chr) {
    dmrs_split <- split(seq_along(dmrs), as.character(GenomicRanges::seqnames(dmrs)))
    filtered_dmrs_list <- lapply(dmrs_split, function(inds) {
        chr_dmrs <- dmrs[inds]
        if (length(chr_dmrs) <= max_dmrs_per_chr) {
            return(chr_dmrs)
        }
        abs_delta_beta <- abs(chr_dmrs$delta_beta)
        score <- minmaxscale(abs_delta_beta)
        ord <- order(score, decreasing = TRUE)
        selected_indices <- ord[1:max_dmrs_per_chr]
        chr_dmrs[selected_indices]
    })
    filtered_dmrs <- unlist(GenomicRanges::GRangesList(filtered_dmrs_list))
    filtered_dmrs
}

.hasNonEmptyString <- function(x) {
    x <- as.character(x)
    !is.na(x) & nzchar(trimws(x))
}

.normalizeCircosRegion <- function(region) {
    if (is.null(region)) {
        return(NULL)
    }

    region_df <- NULL
    if (inherits(region, "GRanges")) {
        if (length(region) == 0) {
            return(NULL)
        }
        region_df <- data.frame(
            chr = as.character(GenomeInfoDb::seqnames(region)),
            start = GenomicRanges::start(region),
            end = GenomicRanges::end(region),
            stringsAsFactors = FALSE
        )
    } else if (is.character(region) && length(region) == 1) {
        clean <- gsub(",", "", trimws(region))
        m <- regexec("^([^:]+):(\\d+)-(\\d+)$", clean)
        parsed <- regmatches(clean, m)[[1]]
        if (length(parsed) != 4) {
            stop("region must be in the format 'chr:start-end' (e.g., 'chr7:100000-200000').")
        }
        region_df <- data.frame(
            chr = parsed[2],
            start = as.numeric(parsed[3]),
            end = as.numeric(parsed[4]),
            stringsAsFactors = FALSE
        )
    } else if (is.data.frame(region)) {
        req_cols <- c("chr", "start", "end")
        if (!all(req_cols %in% colnames(region))) {
            stop("region data frame must contain columns: chr, start, end.")
        }
        region_df <- data.frame(
            chr = as.character(region$chr),
            start = as.numeric(region$start),
            end = as.numeric(region$end),
            stringsAsFactors = FALSE
        )
    } else if (is.list(region) && all(c("chr", "start", "end") %in% names(region))) {
        region_df <- data.frame(
            chr = as.character(region$chr),
            start = as.numeric(region$start),
            end = as.numeric(region$end),
            stringsAsFactors = FALSE
        )
    } else {
        stop("region must be NULL, GRanges, a 'chr:start-end' string, or a data.frame/list with chr/start/end.")
    }

    region_df <- region_df[.hasNonEmptyString(region_df$chr), , drop = FALSE]
    if (nrow(region_df) == 0) {
        return(NULL)
    }
    if (any(!is.finite(region_df$start)) || any(!is.finite(region_df$end))) {
        stop("region start/end must be finite numeric values.")
    }
    region_df$start <- floor(region_df$start)
    region_df$end <- ceiling(region_df$end)
    if (any(region_df$start < 1L) || any(region_df$end < 1L)) {
        stop("region start/end must be >= 1.")
    }
    if (any(region_df$start > region_df$end)) {
        stop("region start must be <= end.")
    }
    unique(region_df)
}

.filterDMRsByScopeForCircos <- function(dmrs, chromosomes = NULL, region_df = NULL) {
    ret <- dmrs
    if (!is.null(chromosomes)) {
        keep_chr <- as.character(GenomicRanges::seqnames(ret)) %in% chromosomes
        ret <- ret[keep_chr]
    }
    if (!is.null(region_df) && nrow(region_df) > 0) {
        region_gr <- GenomicRanges::GRanges(
            seqnames = region_df$chr,
            ranges = IRanges::IRanges(start = region_df$start, end = region_df$end)
        )
        overlaps <- GenomicRanges::findOverlaps(ret, region_gr, ignore.strand = TRUE)
        keep_region <- seq_along(ret) %in% S4Vectors::queryHits(overlaps)
        ret <- ret[keep_region]
    }
    ret
}

.orderChromosomesNaturally <- function(chromosomes) {
    chromosomes <- unique(as.character(chromosomes))
    if (length(chromosomes) == 0) {
        return(chromosomes)
    }
    chr_clean <- gsub("^chr", "", chromosomes, ignore.case = TRUE)
    chr_num <- suppressWarnings(as.numeric(chr_clean))
    chr_special <- match(toupper(chr_clean), c("X", "Y", "M", "MT"))
    chr_special[is.na(chr_special)] <- Inf
    ord <- order(is.na(chr_num), chr_num, chr_special, chr_clean)
    chromosomes[ord]
}

.getCircosSectorBounds <- function(dmrs, unique_chrs, region_df = NULL) {
    if (!is.null(region_df) && nrow(region_df) > 0) {
        bounds <- do.call(rbind, lapply(unique_chrs, function(chr) {
            chr_rows <- region_df[region_df$chr == chr, , drop = FALSE]
            if (nrow(chr_rows) == 0) {
                return(NULL)
            }
            data.frame(
                chr = chr,
                start = min(chr_rows$start),
                end = max(chr_rows$end),
                stringsAsFactors = FALSE
            )
        }))
    } else {
        bounds <- do.call(rbind, lapply(unique_chrs, function(chr) {
            chr_dmrs <- dmrs[as.character(GenomicRanges::seqnames(dmrs)) == chr]
            if (length(chr_dmrs) == 0) {
                return(NULL)
            }
            data.frame(
                chr = chr,
                start = min(GenomicRanges::start(chr_dmrs)),
                end = max(GenomicRanges::end(chr_dmrs)),
                stringsAsFactors = FALSE
            )
        }))
    }
    if (is.null(bounds) || nrow(bounds) == 0) {
        return(NULL)
    }
    bounds$start <- pmax(1L, floor(bounds$start))
    bounds$end <- pmax(bounds$start + 1L, ceiling(bounds$end))
    bounds
}

.clampPointToSectorBounds <- function(point, chr, sector_bounds) {
    if (is.null(sector_bounds) || length(point) == 0) {
        return(point)
    }
    idx <- match(chr, sector_bounds$chr)
    if (is.na(idx)) {
        return(point)
    }
    lower <- sector_bounds$start[idx]
    upper <- sector_bounds$end[idx]
    point <- pmax(lower, pmin(upper, point))
    if (length(point) == 2 && point[1] > point[2]) {
        point <- sort(point)
    }
    point
}

.subsetCytobandForCircos <- function(cytoband, unique_chrs, region_df = NULL) {
    if (is.null(cytoband) || nrow(cytoband) == 0) {
        return(NULL)
    }
    cb <- cytoband[cytoband$V1 %in% unique_chrs, , drop = FALSE]
    if (nrow(cb) == 0) {
        return(NULL)
    }

    if (!is.null(region_df) && nrow(region_df) > 0) {
        trimmed <- do.call(rbind, lapply(seq_len(nrow(region_df)), function(i) {
            chr <- region_df$chr[i]
            chr_cb <- cb[cb$V1 == chr, , drop = FALSE]
            if (nrow(chr_cb) == 0) {
                return(NULL)
            }
            reg_start <- region_df$start[i]
            reg_end <- region_df$end[i]
            overlap <- chr_cb[chr_cb$V3 >= reg_start & chr_cb$V2 <= reg_end, , drop = FALSE]
            if (nrow(overlap) == 0) {
                return(NULL)
            }
            overlap$V2 <- pmax(overlap$V2, reg_start)
            overlap$V3 <- pmin(overlap$V3, reg_end)
            overlap[overlap$V2 < overlap$V3, , drop = FALSE]
        }))
        if (is.null(trimmed) || nrow(trimmed) == 0) {
            return(NULL)
        }
        cb <- unique(trimmed)
    }
    cb <- cb[order(match(cb$V1, unique_chrs), cb$V2, cb$V3), , drop = FALSE]
    cb
}

.orderCircosInteractionRows <- function(link_data) {
    if (is.null(link_data) || nrow(link_data) == 0) {
        return(link_data)
    }
    if ("component_best_rank" %in% colnames(link_data) && any(is.finite(link_data$component_best_rank))) {
        rank_vec <- link_data$component_best_rank
        rank_vec[!is.finite(rank_vec)] <- Inf
        ord <- order(rank_vec, -link_data$sim, link_data$component_id)
    } else {
        ord <- order(-link_data$sim, link_data$component_id)
    }
    link_data[ord, , drop = FALSE]
}

.selectCircosInteractions <- function(link_data, max_components) {
    link_data <- .orderCircosInteractionRows(link_data)
    if (is.null(link_data) || nrow(link_data) == 0) {
        return(link_data)
    }
    if (is.null(max_components) || !is.finite(max_components)) {
        return(link_data)
    }
    max_components <- as.integer(max_components)
    if (max_components < 1L) {
        return(link_data[0, , drop = FALSE])
    }
    if (nrow(link_data) <= max_components) {
        return(link_data)
    }

    matched_mask <- !is.na(link_data$has_jaspar_match) & link_data$has_jaspar_match
    matched <- .orderCircosInteractionRows(link_data[matched_mask, , drop = FALSE])
    unmatched <- .orderCircosInteractionRows(link_data[!matched_mask, , drop = FALSE])
    unique_matched <- unique(matched$component_id) # order is preserved from matched
    n_matched_take <- min(length(unique_matched), max_components)
    matched_comps <- unique_matched[seq_len(n_matched_take)]

    ret <- if (n_matched_take > 0) {
        matched[matched$component_id %in% matched_comps, , drop = FALSE]
    } else {
        matched[0, , drop = FALSE]
    }
    if (length(matched_comps) < max_components && nrow(unmatched) > 0) {
        unique_unmatched <- unique(unmatched$component_id)
        n_unmatched_take <- min(length(unique_unmatched), max_components - length(matched_comps))
        unmatched_comps <- unique_unmatched[seq_len(n_unmatched_take)]
        ret <- rbind(ret, unmatched[unmatched$component_id %in% unmatched_comps, , drop = FALSE])
    }
    ret
}

.wrapCircosLegendLabel <- function(label, width = 88L) {
    paste(strwrap(label, width = width), collapse = "\n")
}

#' Plot Circos Visualization of DMRs
#'
#' @description Creates a circular genome plot using circlize showing DMRs with multiple tracks:
#' ideogram track with chromosome bands, DMR arcs colored by delta beta,
#' beta value heatmaps for each sample, and motif-based interaction links between DMRs.
#'
#' @param dmrs GRanges object or data frame. DMR results from findDMRsFromSeeds.
#' @param beta BetaHandler object, character path to beta file, or beta values matrix.
#' @param pheno Data frame or character path to phenotype file. Sample information with
#'   rownames matching beta column names (required for beta track).
#' @param genome Character. Genome version (e.g., "hg19").
#' @param array Character. Array platform type (default: "450K"). Ignored if sorted_locs is provided.
#' @param sorted_locs Data frame. Genomic locations sorted by position (optional). If NULL, will be fetched based on array and genome.
#' @param sample_group_col Character. Column in pheno for sample grouping (default: "Sample_Group").
#' @param min_similarity Numeric. Minimum motifs PWM similarity threshold for considering DMRs are related (default: 0.7).
#' @param flank_size Integer. Flanking region size for motif extraction in bp (default: 5).
#' @param max_num_samples_per_group Integer. Maximum number of samples to show per group in heatmap (default: 5).
#' @param max_dmrs_per_chr Integer. Maximum number of DMRs to use per chromosome (default: 10). The DMRs with highest absolute delta beta will be selected.
#' @param max_cpgs_per_dmr Integer. Maximum number of CpGs to show per DMR in scatter/heatmap (default: 5).
#' @param max_components Integer. Maximum number of interactions to plot (default: 30).
#' @param chromosomes Character vector. Subset of chromosomes to display (default: NULL, show all available).
#' @param region Genomic region to display. Can be NULL, a GRanges, a string in the form `chr:start-end`,
#'   or a data.frame/list with columns `chr`, `start`, `end`.
#' @param link_geometry Character. Geometry of interaction links: `"midpoint"` (compact, file-size friendly)
#'   or `"range"` (uses full DMR span where possible). Default: `"midpoint"`.
#' @param unmatched_interaction_color Character. Color used for interaction components without JASPAR matches.
#'   These links are shown but omitted from the interaction legend (default: `"#B3B3B3"`).
#' @param legend_width_ratio Numeric. Fraction of horizontal canvas reserved for legends (default: 0.34).
#' @param degenerate_resolution Integer. Resolution in base pairs for showing lines instead of ribbons/lines instead of rectangles (default: 1e6).
#' @param ... Additional arguments (currently unused).
#'
#' @return NULL (creates plot in graphics device).
#'
#' @examples
#' \dontrun{
#' # Load DMR results
#' dmrs <- readRDS("dmrs.rds")
#'
#' # Basic circos plot
#' plotDMRsCircos(dmrs, beta = "beta.txt", pheno = pheno_df)
#'
#' # With custom genome
#' plotDMRsCircos(dmrs, beta = "beta.txt", pheno = pheno_df, genome = "hg38")
#' }
#'
#' @importFrom circlize circos.initializeWithIdeogram circos.trackPlotRegion circos.genomicHeatmap
#' @importFrom circlize circos.genomicLink circos.clear colorRamp2 circos.rect circos.link CELL_META
#' @export
plotDMRsCircos <- function(dmrs,
                           beta,
                           pheno,
                           genome = "hg19",
                           array = c("450K", "27K", "EPIC", "EPICv2"),
                           sorted_locs = NULL,
                           sample_group_col = "Sample_Group",
                           min_similarity = 0.7,
                           flank_size = 5,
                           max_num_samples_per_group = 5,
                           max_dmrs_per_chr = 10,
                           max_cpgs_per_dmr = 5,
                           min_component_size = 2,
                           max_components = 30,
                           chromosomes = NULL,
                           region = NULL,
                           link_geometry = c("midpoint", "range"),
                           unmatched_interaction_color = "#B3B3B3",
                           legend_width_ratio = 0.34,
                           degenerate_resolution = 1e6,
                           output_file = NULL,
                           verbose = NULL,
                           ...) {
    dmrs <- convertToGRanges(dmrs, genome)
    link_geometry <- strex::match_arg(link_geometry, ignore_case = TRUE)
    if (!is.null(max_components)) {
        if (!is.numeric(max_components) || length(max_components) != 1 || is.na(max_components)) {
            stop("max_components must be NULL or a single numeric value.")
        }
        if (max_components < 1) {
            stop("max_components must be >= 1 when provided.")
        }
        max_components <- as.integer(max_components)
    }
    if (!is.character(unmatched_interaction_color) || length(unmatched_interaction_color) != 1 || !nzchar(unmatched_interaction_color)) {
        stop("unmatched_interaction_color must be a single non-empty color string.")
    }
    if (!is.numeric(legend_width_ratio) || length(legend_width_ratio) != 1 || is.na(legend_width_ratio) ||
        legend_width_ratio <= 0.05 || legend_width_ratio >= 0.80) {
        stop("legend_width_ratio must be a single numeric value in (0.05, 0.80).")
    }

    requested_chrs <- NULL
    if (!is.null(chromosomes)) {
        if (!is.character(chromosomes)) {
            stop("chromosomes must be a character vector when provided.")
        }
        requested_chrs <- unique(trimws(chromosomes))
        requested_chrs <- requested_chrs[nzchar(requested_chrs)]
        if (length(requested_chrs) == 0) {
            stop("chromosomes must contain at least one non-empty chromosome name.")
        }
    }
    region_df <- .normalizeCircosRegion(region)
    dmrs <- .filterDMRsByScopeForCircos(dmrs, chromosomes = requested_chrs, region_df = region_df)
    if (length(dmrs) == 0) {
        stop("No DMRs remain after applying chromosome/region filters.")
    }
    present_chrs <- unique(as.character(GenomicRanges::seqnames(dmrs)))
    unique_chrs <- if (is.null(requested_chrs)) {
        .orderChromosomesNaturally(present_chrs)
    } else {
        requested_chrs[requested_chrs %in% present_chrs]
    }
    if (length(unique_chrs) == 0) {
        unique_chrs <- .orderChromosomesNaturally(present_chrs)
    }
    if (!is.null(region_df)) {
        region_df <- region_df[region_df$chr %in% unique_chrs, , drop = FALSE]
        if (nrow(region_df) == 0) {
            stop("No region entries overlap the selected chromosomes/DMRs.")
        }
    }

    if (!is.null(verbose)) {
        options("DMRsegal.verbose" = verbose)
    }
    verbose <- getOption("DMRsegal.verbose", default = 2)
    beta_handler <- getBetaHandler(
        beta = beta,
        array = array,
        genome = genome,
        sorted_locs = sorted_locs
    )

    if (is.character(pheno) && length(pheno) == 1 && file.exists(pheno)) {
        pheno <- read.table(pheno, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
    }
    if (!is.data.frame(pheno)) {
        stop("pheno must be a data frame or a valid file path to a phenotype table")
    }
    if (!(sample_group_col %in% colnames(pheno))) {
        stop(sprintf("sample_group_col '%s' not found in pheno data frame", sample_group_col))
    }
    beta_locs <- beta_handler$getBetaLocs()

    .log_step("Preparing data for Circos plot...")

    .log_step("Filtering DMRs for plotting by maximum absolute delta beta...", level = 2)
    heatmap_dmrs <- .filterDMRsForCircos(dmrs, max_dmrs_per_chr)
    .log_success("DMRs filtered for Circos plot heatmap", level = 2)

    .log_info("Total DMRs to plot on heatmap: ", length(heatmap_dmrs), level = 2)

    cytoband <- .getCytobandData(genome)
    cytoband_subset <- .subsetCytobandForCircos(cytoband, unique_chrs, region_df = region_df)
    sector_bounds <- .getCircosSectorBounds(dmrs, unique_chrs, region_df = region_df)

    .log_step("Preparing DMRs data...", level = 2)
    arc_data <- .prepareCircosArcData(dmrs)
    if (!is.null(arc_data) && length(arc_data$chr) > 0 && !is.null(sector_bounds)) {
        arc_idx <- match(arc_data$chr, sector_bounds$chr)
        keep_arc <- !is.na(arc_idx)
        arc_data$chr <- arc_data$chr[keep_arc]
        arc_data$start <- arc_data$start[keep_arc]
        arc_data$end <- arc_data$end[keep_arc]
        arc_data$delta_beta <- arc_data$delta_beta[keep_arc]
        arc_idx <- arc_idx[keep_arc]
        arc_data$start <- pmax(arc_data$start, sector_bounds$start[arc_idx])
        arc_data$end <- pmin(arc_data$end, sector_bounds$end[arc_idx])
    }
    .log_success("DMR arcs data prepared", level = 2)
    if (getOption("DMRsegal.make_debug_dir", FALSE)) {
        .log_info("Saving DMR arcs data to debug/circos_arc_data.tsv", level = 1)
        dir.create("debug", showWarnings = FALSE)
        write.table(arc_data, file = "debug/circos_arc_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    }

    .log_step("Preparing heatmap data...", level = 2)
    heatmap_data <- .prepareCircosHeatmapData(
        heatmap_dmrs, beta_handler, pheno, sample_group_col,
        max_cpgs_per_dmr, max_num_samples_per_group
    )
    heatmap_df <- heatmap_data$heatmap_df
    reduced_pheno <- heatmap_data$reduced_pheno
    if (!is.null(heatmap_df) && nrow(heatmap_df) > 0 && !is.null(sector_bounds)) {
        hidx <- match(heatmap_df$chr, sector_bounds$chr)
        keep_h <- !is.na(hidx)
        heatmap_df <- heatmap_df[keep_h, , drop = FALSE]
        if (nrow(heatmap_df) > 0) {
            hidx <- hidx[keep_h]
            heatmap_df$start <- pmax(heatmap_df$start, sector_bounds$start[hidx])
            heatmap_df$end <- pmin(heatmap_df$end, sector_bounds$end[hidx])
            heatmap_df <- heatmap_df[heatmap_df$start <= heatmap_df$end, , drop = FALSE]
        }
    }
    .log_success("Heatmap data prepared", level = 2)
    .log_info("Total heatmap entries: ", nrow(heatmap_df), level = 2)
    if (getOption("DMRsegal.make_debug_dir", FALSE)) {
        .log_info("Saving heatmap data to debug/circos_heatmap_data.tsv", level = 1)
        dir.create("debug", showWarnings = FALSE)
        write.table(heatmap_df, file = "debug/circos_heatmap_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    }

    .log_step("Computing motif-based DMR interactions...", level = 2)
    link_data <- .prepareCircosLinkData(
        dmrs, genome, array, beta_locs, min_similarity, flank_size, max_components, min_component_size
    )
    if (!is.null(link_data) && nrow(link_data) > 0 && !is.null(sector_bounds)) {
        lidx1 <- match(link_data$chr1, sector_bounds$chr)
        lidx2 <- match(link_data$chr2, sector_bounds$chr)
        keep_l <- !is.na(lidx1) & !is.na(lidx2)
        link_data <- link_data[keep_l, , drop = FALSE]
        if (nrow(link_data) > 0) {
            lidx1 <- lidx1[keep_l]
            lidx2 <- lidx2[keep_l]
            link_data$start1 <- pmax(link_data$start1, sector_bounds$start[lidx1])
            link_data$end1 <- pmin(link_data$end1, sector_bounds$end[lidx1])
            link_data$start2 <- pmax(link_data$start2, sector_bounds$start[lidx2])
            link_data$end2 <- pmin(link_data$end2, sector_bounds$end[lidx2])
            keep_span <- link_data$start1 <= link_data$end1 & link_data$start2 <= link_data$end2
            link_data <- link_data[keep_span, , drop = FALSE]
        }
    }

    .log_success("DMR interactions data prepared", level = 2)
    if (getOption("DMRsegal.make_debug_dir", FALSE) && !is.null(link_data) && nrow(link_data) > 0) {
        .log_info("Saving DMR interactions data to debug/circos_link_data.tsv", level = 1)
        dir.create("debug", showWarnings = FALSE)
        write.table(link_data, file = "debug/circos_link_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    }
    if (!is.null(output_file)) {
        grDevices::cairo_pdf(output_file, width = 14, height = 10, fallback_resolution = 300)
    }
    if (.Device == "null device") {
        grDevices::cairo_pdf(width = 14, height = 10, fallback_resolution = 300)
    }
    .log_step("Creating Circos plot...")
    on.exit(try(circlize::circos.clear(), silent = TRUE), add = TRUE)
    graphics::plot.new()
    circle_width_ratio <- 1 - legend_width_ratio
    root_layout <- grid::viewport(layout = grid::grid.layout(
        nrow = 1, ncol = 2,
        widths = grid::unit.c(grid::unit(circle_width_ratio, "npc"), grid::unit(legend_width_ratio, "npc"))
    ))
    grid::pushViewport(root_layout)
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::pushViewport(grid::viewport(width = grid::unit(1, "snpc"), height = grid::unit(1, "snpc")))
    graphics::par(omi = gridBase::gridOMI(), new = TRUE)
    circlize::circos.par(gap.after = c(rep(2, length(unique_chrs) - 1), 10))

    heatmap_legend <- NULL
    arc_legend <- NULL
    link_legend <- NULL
    use_manual_init <- FALSE

    if (!is.null(cytoband_subset) && nrow(cytoband_subset) > 0) {
        circlize::circos.initializeWithIdeogram(cytoband = cytoband_subset, plotType = c("labels", "axis"))
    } else if (!is.null(region_df) && !is.null(sector_bounds) && nrow(sector_bounds) > 0) {
        use_manual_init <- TRUE
        circlize::circos.initialize(
            factors = sector_bounds$chr,
            xlim = as.matrix(sector_bounds[, c("start", "end"), drop = FALSE])
        )
        circlize::circos.trackPlotRegion(
            ylim = c(0, 1),
            bg.border = NA,
            track.height = 0.06,
            panel.fun = function(x, y) {
                circlize::circos.axis(h = "top", labels.cex = 0.45, major.tick.percentage = 0.3)
                circlize::circos.text(
                    x = circlize::CELL_META$xcenter,
                    y = circlize::CELL_META$ylim[2] + circlize::convert_y(1.2, "mm"),
                    labels = circlize::CELL_META$sector.index,
                    cex = 0.6,
                    facing = "inside",
                    niceFacing = TRUE,
                    adj = c(0.5, 0)
                )
            }
        )
    } else {
        circlize::circos.initializeWithIdeogram(species = genome, chromosome.index = unique_chrs, plotType = c("labels", "axis"))
    }

    if (!is.null(heatmap_df) && nrow(heatmap_df) > 0) {
        .log_step("Adding heatmap track...", level = 2)
        heatmap_numeric <- as.matrix(heatmap_df[, -(1:3), drop = FALSE])
        storage.mode(heatmap_numeric) <- "numeric"
        valid_vals <- as.numeric(heatmap_numeric)
        valid_vals <- valid_vals[is.finite(valid_vals)]
        if (length(valid_vals) == 0) {
            .log_warn("Heatmap beta values are all missing/non-finite. Skipping heatmap track.")
        } else {
            q <- stats::quantile(valid_vals, probs = c(0.05, 0.5, 0.95), na.rm = TRUE, names = FALSE, type = 8)
            if (length(unique(q)) < 3) {
                q <- c(min(valid_vals), mean(range(valid_vals)), max(valid_vals))
            }
            if (length(unique(q)) < 3) {
                q <- c(q[1] - 1e-6, q[2], q[3] + 1e-6)
            }
            col_fun <- circlize::colorRamp2(q, c("#2b83ba", "#f7f7f7", "#d7191c"))
            heatmap_height <- 0.3
            heatmap_df_plot <- data.frame(
                heatmap_df[, 1:3, drop = FALSE],
                heatmap_numeric,
                check.names = FALSE,
                stringsAsFactors = FALSE
            )
            previous_track_index <- circlize::get.current.track.index()
            circlize::circos.genomicHeatmap(
                bed = heatmap_df_plot,
                col = col_fun,
                border = NA,
                side = "outside",
                border_lwd = 0,
                line_col = "#6E6E6E",
                line_lwd = 0.45,
                heatmap_height = heatmap_height,
            )
            heatmap_track_index <- previous_track_index + 1
        }

        if (length(valid_vals) > 0) {
            suppressMessages(circlize::circos.track(track.index = heatmap_track_index, panel.fun = function(x, y) {
                if (circlize::CELL_META$sector.numeric.index == length(unique_chrs)) { # the last sector
                    groups <- as.character(reduced_pheno[[sample_group_col]])
                    group_runs <- rle(groups)
                    unique_groups <- unique(group_runs$values)

                    group_colors <- colorspace::qualitative_hcl(length(unique_groups), palette = "Pastel 1")
                    names(group_colors) <- unique_groups
                    y_limits <- circlize::CELL_META$cell.ylim
                    n_samples <- sum(group_runs$lengths)
                    if (n_samples <= 0) {
                        return(invisible(NULL))
                    }
                    row_height <- (y_limits[2] - y_limits[1]) / n_samples
                    y_cursor <- y_limits[1]
                    for (i in seq_along(group_runs$values)) {
                        group <- group_runs$values[i]
                        group_size <- group_runs$lengths[i]
                        y_next <- y_cursor + group_size * row_height
                        y_bottom <- min(y_cursor, y_next)
                        y_top <- max(y_cursor, y_next)
                        circlize::circos.rect(
                            circlize::CELL_META$cell.xlim[2] + circlize::convert_x(1, "mm"),
                            y_bottom,
                            circlize::CELL_META$cell.xlim[2] + circlize::convert_x(5, "mm"),
                            y_top,
                            col = group_colors[group],
                            border = NA
                        )
                        circlize::circos.text(
                            circlize::CELL_META$cell.xlim[2] + circlize::convert_x(3, "mm"),
                            y_bottom + (y_top - y_bottom) / 2,
                            group,
                            cex = 0.5,
                            facing = "clockwise"
                        )
                        y_cursor <- y_next
                    }
                }
            }, bg.border = NA))


            # Make legend for heatmap
            legend_at <- signif(q, 2)
            heatmap_legend <- ComplexHeatmap::Legend(
                title = "Samples \u03b2-values",
                at = legend_at,
                col_fun = col_fun,
                title_position = "topleft",
                legend_height = grid::unit(4, "cm"),
                labels_gp = grid::gpar(fontsize = 8),
                title_gp = grid::gpar(fontsize = 10, fontface = "bold")
            )
            .log_success("Heatmap track added", level = 2)
        }
    }

    if (!use_manual_init && !is.null(cytoband_subset) && nrow(cytoband_subset) > 0) {
        circlize::circos.genomicIdeogram(cytoband = cytoband_subset)
    } else if (!use_manual_init) {
        circlize::circos.genomicIdeogram(species = genome)
    }


    if (!is.null(arc_data)) {
        .log_step("Adding DMR arc track...", level = 2)

        delta_beta <- as.numeric(arc_data$delta_beta)
        valid_vals <- delta_beta[is.finite(delta_beta)]
        if (length(valid_vals) == 0) {
            .log_warn("DMR delta_beta values are all missing/non-finite. Skipping arc track.")
            arc_data <- NULL
        } else {
            neg_min <- suppressWarnings(min(valid_vals[valid_vals < 0], na.rm = TRUE))
            pos_max <- suppressWarnings(max(valid_vals[valid_vals > 0], na.rm = TRUE))
            if (!is.finite(neg_min)) {
                neg_min <- -max(1e-6, abs(pos_max) * 0.05)
            }
            if (!is.finite(pos_max)) {
                pos_max <- max(1e-6, abs(neg_min) * 0.05)
            }
            if (neg_min >= 0) {
                neg_min <- -max(1e-6, abs(pos_max) * 0.05)
            }
            if (pos_max <= 0) {
                pos_max <- max(1e-6, abs(neg_min) * 0.05)
            }
            q <- c(neg_min, 0, pos_max)
            col_fun <- circlize::colorRamp2(q, c("#055709", "white", "#801414"))
        }

        if (!is.null(arc_data)) {
            arc_colors <- col_fun(arc_data$delta_beta)

            arc_df <- data.frame(
                chr = arc_data$chr,
                start = arc_data$start,
                end = arc_data$end,
                delta_beta = arc_data$delta_beta,
                color = arc_colors,
                stringsAsFactors = FALSE
            )

            circlize::circos.trackPlotRegion(
                sectors = unique_chrs,
                bg.border = NA,
                bg.col = "white",
                track.height = 0.05,
                ylim = c(0, 1),
                panel.fun = function(x, y) {
                    chr <- circlize::CELL_META$sector.index
                    dmr_chr <- arc_df[arc_df$chr == chr, ]
                    if (nrow(dmr_chr) > 0) {
                        for (i in seq_len(nrow(dmr_chr))) {
                            if (dmr_chr$end[i] - dmr_chr$start[i] < degenerate_resolution) {
                                circlize::circos.lines(
                                    x = c((dmr_chr$start[i] + dmr_chr$end[i]) / 2, (dmr_chr$start[i] + dmr_chr$end[i]) / 2),
                                    y = c(0, 1),
                                    col = dmr_chr$color[i],
                                    border = NA
                                )
                            } else {
                                circlize::circos.rect(
                                    xleft = dmr_chr$start[i],
                                    xright = dmr_chr$end[i],
                                    ybottom = 0,
                                    ytop = 1,
                                    col = dmr_chr$color[i],
                                    border = NA,
                                    density = NULL
                                )
                            }
                        }
                    }
                }
            )
            # Add legend for arc colors
            q <- signif(q, 2)
            arc_legend <- ComplexHeatmap::Legend(
                title = "DMR \u0394\u03b2",
                at = q,
                col_fun = col_fun,
                title_position = "topleft",
                legend_height = grid::unit(4, "cm"),
                labels_gp = grid::gpar(fontsize = 8),
                title_gp = grid::gpar(fontsize = 10, fontface = "bold")
            )

            .log_success("Arc track added", level = 2)
        }
    }


    if (!is.null(link_data) && nrow(link_data) > 0) {
        .log_step("Adding link track...", level = 2)

        link_data$has_jaspar_match <- !is.na(link_data$has_jaspar_match) & link_data$has_jaspar_match
        comp_data <- link_data[
            !duplicated(link_data$component_id),
            c("component_id", "size", "consensus_seq", "jaspar_names", "jaspar_corr", "has_jaspar_match", "component_best_rank"),
            drop = FALSE
        ]

        component_colors <- setNames(
            rep(unmatched_interaction_color, nrow(comp_data)),
            as.character(comp_data$component_id)
        )
        matched_components <- comp_data[comp_data$has_jaspar_match, , drop = FALSE]
        if (nrow(matched_components) > 0) {
            rank_vec <- matched_components$component_best_rank
            rank_vec[!is.finite(rank_vec)] <- Inf
            matched_components <- matched_components[
                order(rank_vec, -matched_components$size, matched_components$component_id),
                ,
                drop = FALSE
            ]
            hues <- seq(15, 375, length.out = nrow(matched_components) + 1)
            matched_cols <- grDevices::hcl(h = hues, l = 62, c = 95)[seq_len(nrow(matched_components))]
            component_colors[as.character(matched_components$component_id)] <- matched_cols
        }

        link_data$colors <- vapply(seq_len(nrow(link_data)), function(i) {
            comp_id <- as.character(link_data$component_id[i])
            if (!isTRUE(link_data$has_jaspar_match[i])) {
                return(grDevices::adjustcolor(unmatched_interaction_color, alpha.f = 0.85))
            }
            base_color <- component_colors[comp_id]
            sim_val <- link_data$sim[i]
            if (!is.finite(sim_val)) {
                sim_val <- 0
            }
            sim_val <- pmax(0, pmin(1, sim_val))
            shade_factor <- 0.45 + 0.55 * sim_val
            base_rgb <- grDevices::col2rgb(base_color) / 255
            grDevices::rgb(
                red = base_rgb[1] * shade_factor,
                green = base_rgb[2] * shade_factor,
                blue = base_rgb[3] * shade_factor,
                alpha = 0.85
            )
        }, character(1))

        legend_components <- comp_data[comp_data$has_jaspar_match, , drop = FALSE]
        if (nrow(legend_components) > 0) {
            rank_vec <- legend_components$component_best_rank
            rank_vec[!is.finite(rank_vec)] <- Inf
            legend_components <- legend_components[order(rank_vec, legend_components$component_id), , drop = FALSE]
            link_legend_colors <- component_colors[as.character(legend_components$component_id)]
            link_legend_labels <- vapply(seq_len(nrow(legend_components)), function(i) {
                rank_prefix <- if (is.finite(legend_components$component_best_rank[i])) {
                    paste0("[rank=", as.integer(legend_components$component_best_rank[i]), "] ")
                } else {
                    ""
                }
                label <- paste0(
                    rank_prefix,
                    "[n=", legend_components$size[i], "] ",
                    legend_components$consensus_seq[i]
                )
                if (.hasNonEmptyString(legend_components$jaspar_names[i])) {
                    jas_names <- trimws(strsplit(legend_components$jaspar_names[i], ",", fixed = TRUE)[[1]])
                    jas_corr <- if (.hasNonEmptyString(legend_components$jaspar_corr[i])) {
                        trimws(strsplit(legend_components$jaspar_corr[i], ",", fixed = TRUE)[[1]])
                    } else {
                        rep("", length(jas_names))
                    }
                    n_show <- min(3L, length(jas_names))
                    for (j in seq_len(n_show)) {
                        corr_val <- suppressWarnings(as.numeric(jas_corr[j]))
                        corr_txt <- if (is.finite(corr_val)) {
                            paste0(" (", signif(corr_val, 3), ")")
                        } else {
                            ""
                        }
                        label <- paste0(label, " | ", jas_names[j], corr_txt)
                    }
                    if (length(jas_names) > n_show) {
                        label <- paste0(label, " ...")
                    }
                }
                .wrapCircosLegendLabel(label)
            }, character(1))
            link_legend <- ComplexHeatmap::Legend(
                labels = link_legend_labels,
                legend_gp = grid::gpar(fill = link_legend_colors),
                title = "DMR Interaction Components\n(JASPAR-annotated)",
                title_position = "topleft",
                title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                labels_gp = grid::gpar(fontsize = 8)
            )
        }

        has_directionality <- "rank" %in% colnames(S4Vectors::mcols(dmrs)) &&
            any(is.finite(suppressWarnings(as.numeric(S4Vectors::mcols(dmrs)$rank))))
        directional_flag <- if (has_directionality) 1 else 0
        for (i in seq_len(nrow(link_data))) {
            point1 <- c(link_data$start1[i], link_data$end1[i])
            point2 <- c(link_data$start2[i], link_data$end2[i])
            if (identical(link_geometry, "midpoint")) {
                point1 <- mean(point1)
                point2 <- mean(point2)
            } else {
                if (point1[2] - point1[1] < degenerate_resolution) {
                    point1 <- mean(point1)
                }
                if (point2[2] - point2[1] < degenerate_resolution) {
                    point2 <- mean(point2)
                }
            }
            point1 <- .clampPointToSectorBounds(point1, link_data$chr1[i], sector_bounds)
            point2 <- .clampPointToSectorBounds(point2, link_data$chr2[i], sector_bounds)
            if (has_directionality) {
                circlize::circos.link(
                    sector.index1 = link_data$chr1[i],
                    point1 = point1,
                    sector.index2 = link_data$chr2[i],
                    point2 = point2,
                    col = link_data$colors[i],
                    directional = directional_flag,
                    arr.type = "triangle",
                    arr.length = 0.3,
                    arr.width = 0.05,
                    border = NA
                )
            } else {
                circlize::circos.link(
                    sector.index1 = link_data$chr1[i],
                    point1 = point1,
                    sector.index2 = link_data$chr2[i],
                    point2 = point2,
                    col = link_data$colors[i],
                    directional = 0,
                    border = NA
                )
            }
        }
        .log_success("Link track added", level = 2)
    }
    circlize::circos.clear()
    grid::upViewport(2)
    top_legends <- Filter(Negate(is.null), list(heatmap_legend, arc_legend))
    has_any_legend <- length(top_legends) > 0 || !is.null(link_legend)
    if (has_any_legend) {
        grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
        top_packed <- NULL
        if (length(top_legends) > 0) {
            top_args <- c(top_legends, list(
                direction = "horizontal",
                column_gap = grid::unit(6, "mm"),
                row_gap = grid::unit(6, "mm"),
                max_width = grid::unit(0.98, "npc")
            ))
            top_packed <- do.call(ComplexHeatmap::packLegend, top_args)
        }
        lgd_list <- top_packed
        if (!is.null(link_legend)) {
            if (is.null(top_packed)) {
                lgd_list <- link_legend
            } else {
                lgd_list <- ComplexHeatmap::packLegend(
                    top_packed,
                    link_legend,
                    direction = "vertical",
                    row_gap = grid::unit(3, "mm"),
                    max_height = grid::unit(0.98, "npc")
                )
            }
        }
        ComplexHeatmap::draw(
            lgd_list,
            x = grid::unit(0, "npc"),
            y = grid::unit(1, "npc"),
            just = c("left", "top")
        )
        grid::upViewport()
    }
    grid::upViewport()
    .log_success("Circos plot created successfully")
    if (!is.null(output_file)) {
        grDevices::dev.off()
    }


    invisible(NULL)
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

#' Plot Manhattan-Style View of DMR Scores
#'
#' @description Creates a Manhattan-style genome-wide scatter plot of DMR scores.
#' DMRs are colored by their dominant genomic region class inferred from DMR annotations.
#'
#' @param dmrs GRanges object or data frame. DMR results from findDMRsFromSeeds.
#' @param score_col Character. Metadata column used for y-axis values (default: `"score"`).
#' @param genome Character. Genome version passed to `convertToGRanges` (default: `"hg19"`).
#' @param promoter_col Character. Metadata column indicating promoter overlap (default: `"in_promoter_of"`).
#' @param gene_body_col Character. Metadata column indicating gene-body overlap (default: `"in_gene_body_of"`).
#' @param point_size Numeric. Point size (default: `1.1`).
#' @param point_alpha Numeric. Point alpha in `[0, 1]` (default: `0.75`).
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
                              score_col = "score",
                              genome = "hg19",
                              promoter_col = "in_promoter_of",
                              gene_body_col = "in_gene_body_of",
                              point_size = 1.1,
                              point_alpha = 0.75) {
    dmrs <- convertToGRanges(dmrs, genome = genome)
    if (!(score_col %in% colnames(S4Vectors::mcols(dmrs)))) {
        stop("score_col '", score_col, "' not found in DMR metadata.")
    }
    if (!is.numeric(point_size) || length(point_size) != 1 || is.na(point_size) || point_size <= 0) {
        stop("point_size must be a single positive numeric value.")
    }
    if (!is.numeric(point_alpha) || length(point_alpha) != 1 || is.na(point_alpha) || point_alpha < 0 || point_alpha > 1) {
        stop("point_alpha must be a single numeric value in [0, 1].")
    }

    scores <- suppressWarnings(as.numeric(S4Vectors::mcols(dmrs)[[score_col]]))
    keep <- is.finite(scores)
    if (!any(keep)) {
        stop("No finite score values found in score_col '", score_col, "'.")
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

    p <- ggplot2::ggplot(dmr_df, ggplot2::aes(x = position, y = score, color = region_class)) +
        ggplot2::geom_point(size = point_size, alpha = point_alpha, stroke = 0) +
        ggplot2::scale_color_manual(values = region_colors, drop = TRUE, name = "Primary Region") +
        ggplot2::scale_x_continuous(
            breaks = axis_df$position,
            labels = as.character(axis_df$chr),
            expand = ggplot2::expansion(mult = c(0.01, 0.01))
        ) +
        ggplot2::labs(
            x = "Chromosome",
            y = score_col,
            title = paste0("DMR Manhattan Plot (", score_col, ")")
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
    p
}


.getCytobandData <- function(genome) {
    cache_dir <- tools::R_user_dir("DMRsegal", which = "cache")
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE)
    }

    cache_file <- file.path(cache_dir, paste0("cytoband_", genome, ".rds"))

    if (file.exists(cache_file)) {
        .log_info("Loading cached cytoband data for ", genome, level = 3)
        return(readRDS(cache_file))
    }

    .log_step("Downloading cytoband data from UCSC for ", genome, "...", level = 3)

    tryCatch(
        {
            url <- paste0("https://hgdownload.cse.ucsc.edu/goldenpath/", genome, "/database/cytoBandIdeo.txt.gz")

            temp_file <- tempfile(fileext = ".txt.gz")
            result <- tryCatch(
                {
                    utils::download.file(url, temp_file, quiet = TRUE, mode = "wb")
                    TRUE
                },
                error = function(e) {
                    .log_warn("cytoBandIdeo not found, trying cytoBand table...")
                    FALSE
                }
            )

            if (!result) {
                url <- paste0("https://hgdownload.cse.ucsc.edu/goldenpath/", genome, "/database/cytoBand.txt.gz")
                utils::download.file(url, temp_file, quiet = TRUE, mode = "wb")
            }

            cytoband <- read.table(temp_file, sep = "\t", stringsAsFactors = FALSE, header = FALSE)

            if (ncol(cytoband) >= 5) {
                colnames(cytoband) <- c("V1", "V2", "V3", "V4", "V5")
            } else if (ncol(cytoband) == 4) {
                colnames(cytoband) <- c("V1", "V2", "V3", "V4")
                cytoband$V5 <- ""
            } else {
                stop("Unexpected cytoband format")
            }

            unlink(temp_file)

            saveRDS(cytoband, cache_file)
            .log_success("Cytoband data downloaded and cached", level = 3)

            return(cytoband)
        },
        error = function(e) {
            .log_warn("Failed to download cytoband data: ", e$message, ". Using default ideogram.")
            NULL
        }
    )
}


.getGenomeLengths <- function(genome, chromosomes) {
    supported_organisms <- Organism.dplyr::supportedOrganisms()
    matched_row <- supported_organisms[grepl(tolower(genome), tolower(supported_organisms$TxDb)), , drop = FALSE]

    if (nrow(matched_row) == 0) {
        .log_warn("Unsupported genome: ", genome, ". Using default genome string.")
        return(genome)
    }

    txdb_pkg <- matched_row$TxDb[1]

    if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
        .log_warn("TxDb package ", txdb_pkg, " not installed. Using default genome string.")
        return(genome)
    }

    txdb <- tryCatch(
        {
            if (!isNamespaceLoaded(txdb_pkg)) {
                loadNamespace(txdb_pkg)
            }
            getExportedValue(txdb_pkg, txdb_pkg)
        },
        error = function(e) {
            .log_warn("Failed to load TxDb: ", e$message)
            NULL
        }
    )

    if (is.null(txdb)) {
        return(genome)
    }

    seqinfo <- GenomeInfoDb::seqinfo(txdb)
    seqlengths <- GenomeInfoDb::seqlengths(seqinfo)

    unique_chrs <- unique(chromosomes)
    unique_chrs <- unique_chrs[unique_chrs %in% names(seqlengths)]

    if (length(unique_chrs) == 0) {
        .log_warn("No valid chromosomes found in DMRs. Using default genome.")
        return(genome)
    }

    genome_lengths <- as.list(seqlengths[unique_chrs])

    genome_lengths
}

.closest_rows_indices_to_centroids <- function(beta_pcs, centers, verbose = FALSE) { # nolint
    stopifnot(nrow(beta_pcs) >= nrow(centers))
    selection <- integer()
    assessed_pcs <- beta_pcs
    assessed_centers <- centers
    assessed_indices <- seq_len(nrow(beta_pcs))
    while (length(selection) < nrow(centers)) {
        nn <- FNN::get.knnx(assessed_pcs, assessed_centers, k = 1)
        center_indices <- as.vector(nn$nn.index)
        duplicated_centers <- duplicated(center_indices)
        chosen_indices <- center_indices[!duplicated_centers]
        selection <- c(selection, assessed_indices[chosen_indices])
        if (any(duplicated_centers)) {
            .log_info(sum(duplicated_centers), " duplicated nearest neighbours, retrying on the duplicates...", level = 3)
            assessed_pcs <- assessed_pcs[-chosen_indices, , drop = FALSE]
            assessed_indices <- assessed_indices[-chosen_indices]
            assessed_centers <- assessed_centers[duplicated_centers, , drop = FALSE]
        }
    }
    sort(selection)
}

.prepareCircosHeatmapData <- function(dmrs, beta_handler, pheno, sample_group_col, max_sup_cpgs_per_dmr_side = 2, max_num_samples_per_group = 10) {
    beta_col_names <- beta_handler$getBetaColNames()
    pheno <- pheno[rownames(pheno) %in% beta_col_names, , drop = FALSE]


    if (nrow(pheno) == 0) {
        .log_warn("No samples in pheno match the samples in beta values. Skipping heatmap track.")
        return(list(heatmap_df = NULL, reduced_pheno = NULL))
    }
    # Order pheno by sample group
    pheno <- pheno[order(pheno[[sample_group_col]]), , drop = FALSE]
    dmrs_cpgs_list <- .getSupportingSitesFromColumns(
        dmrs,
        max_sup_cpgs_per_dmr_side = max_sup_cpgs_per_dmr_side,
        separate_by_section = FALSE,
        use_absolute_indices = FALSE,
        beta_locs = beta_handler$getBetaLocs()
    )
    dmrs_cpgs_inds <- unique(unlist(dmrs_cpgs_list, use.names = FALSE))
    if (length(dmrs_cpgs_inds) == 0) {
        .log_warn("No supporting CpGs available for selected DMRs. Skipping heatmap track.")
        return(list(heatmap_df = NULL, reduced_pheno = NULL))
    }


    shown_locs <- beta_handler$getBetaLocs()[dmrs_cpgs_inds, c("chr", "start", "end"), drop = FALSE]
    shown_locs$chr <- as.character(shown_locs$chr)
    shown_locs$start <- as.numeric(shown_locs$start)
    shown_locs$end <- as.numeric(shown_locs$end)
    beta_data <- beta_handler$getBeta(
        row_names = dmrs_cpgs_inds,
        col_names = rownames(pheno)
    )
    if (max_num_samples_per_group > 0) {
        selected_samples <- c()
        groups <- pheno[[sample_group_col]]
        unique_groups <- unique(groups)
        reduced_pheno <- data.frame()
        for (group in unique_groups) {
            group_samples <- rownames(pheno)[groups == group]
            if (length(group_samples) > max_num_samples_per_group) {
                # Do K-means KNN sampling, selecting max_num_samples_per_group samples that best represent the group
                group_beta <- beta_data[, group_samples, drop = FALSE]
                non_zero_cols <- apply(group_beta, 2, function(v) stats::var(v, na.rm = TRUE) != 0)
                group_beta <- group_beta[, non_zero_cols, drop = FALSE]
                # transpose to have samples as rows
                group_beta <- t(group_beta)
                pcs <- stats::prcomp(group_beta, center = TRUE, scale. = TRUE, rank. = min(10, ncol(group_beta)))$x
                set.seed(getOption("DMRsegal.random_seed", 42))
                kmeans <- stats::kmeans(pcs, centers = max_num_samples_per_group, algorithm = "Lloyd", iter.max = 1000, nstart = 5)
                group_samples <- group_samples[.closest_rows_indices_to_centroids(pcs, kmeans$centers)]
            }
            reduced_pheno <- rbind(reduced_pheno, pheno[group_samples, , drop = FALSE])
            selected_samples <- c(selected_samples, group_samples)
        }
        beta_data <- beta_data[, selected_samples, drop = FALSE]
    } else {
        reduced_pheno <- pheno
    }
    beta_data <- as.matrix(beta_data)
    storage.mode(beta_data) <- "numeric"
    heatmap_df <- data.frame(
        shown_locs,
        beta_data,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    list(heatmap_df = heatmap_df, reduced_pheno = reduced_pheno)
}


.prepareCircosArcData <- function(dmrs) {
    list(
        chr = as.character(GenomicRanges::seqnames(dmrs)),
        start = GenomicRanges::start(dmrs),
        end = GenomicRanges::end(dmrs),
        delta_beta = as.numeric(S4Vectors::mcols(dmrs)$delta_beta)
    )
}


.prepareCircosLinkData <- function(dmrs, genome, array, beta_locs, min_similarity, flank_size, max_components, min_component_size) {
    ret <- tryCatch(
        {
            computeDMRsInteraction(
                dmrs = dmrs,
                genome = genome,
                array = array,
                min_similarity = min_similarity,
                beta_locs = beta_locs,
                flank_size = flank_size,
                min_component_size = min_component_size
            )
        },
        error = function(e) {
            .log_warn("Failed to compute motif-based interactions: ", e$message)
            NULL
        }
    )

    if (is.null(ret) || nrow(ret$interactions) == 0) {
        .log_warn("No significant interactions found at similarity >=", min_similarity, ". Skipping link track.")
        return(NULL)
    }

    .log_success("Found ", nrow(ret$interactions), " motif-based interactions")
    if (!is.null(ret$components) && nrow(ret$components) > 0) {
        largest_component <- max(ret$components$size, na.rm = TRUE)
        if (is.finite(largest_component) && largest_component > 0.8 * length(dmrs)) {
            .log_warn(
                "Largest interaction component contains ", largest_component, " / ", length(dmrs),
                " DMRs. This often indicates dense generic CpG-context motif connectivity."
            )
        }
    }
    link_data <- ret$interactions
    if ("component_id" %in% colnames(link_data)) {
        link_data <- link_data[!is.na(link_data$component_id), , drop = FALSE]
    }
    if (nrow(link_data) == 0) {
        .log_warn("No motif interactions remain after component filtering. Skipping link track.")
        return(NULL)
    }

    if (!is.null(ret$components) && nrow(ret$components) > 0) {
        link_data <- merge(
            link_data,
            ret$components,
            by = "component_id",
            all.x = TRUE,
            suffixes = c("", "_comp"),
            sort = FALSE
        )
    } else {
        link_data$size <- NA_real_
        link_data$consensus_seq <- NA_character_
        link_data$jaspar_names <- NA_character_
        link_data$jaspar_corr <- NA_character_
    }

    link_data$has_jaspar_match <- .hasNonEmptyString(link_data$jaspar_names)
    if ("rank" %in% colnames(S4Vectors::mcols(dmrs))) {
        ranks <- suppressWarnings(as.numeric(S4Vectors::mcols(dmrs)$rank))
        pair_rank <- pmin(ranks[link_data$index1], ranks[link_data$index2], na.rm = TRUE)
        pair_rank[!is.finite(pair_rank)] <- NA_real_
        link_data$interaction_rank <- pair_rank
        comp_rank_df <- stats::aggregate(
            interaction_rank ~ component_id,
            data = link_data,
            FUN = function(x) {
                x <- x[is.finite(x)]
                if (length(x) == 0) {
                    NA_real_
                } else {
                    min(x)
                }
            }
        )
        colnames(comp_rank_df)[colnames(comp_rank_df) == "interaction_rank"] <- "component_best_rank"
        link_data <- merge(link_data, comp_rank_df, by = "component_id", all.x = TRUE, sort = FALSE)
    } else {
        link_data$component_best_rank <- NA_real_
    }

    original_n <- nrow(link_data)
    link_data <- .selectCircosInteractions(link_data, max_components)
    if (!is.null(max_components) && nrow(link_data) < original_n) {
        .log_info("Selected ", nrow(link_data), " interactions (max_components=", max_components, ").", level = 2)
    }
    if (nrow(link_data) > 0) {
        matched_n <- sum(link_data$has_jaspar_match, na.rm = TRUE)
        .log_info(
            "Selected interactions: ",
            matched_n, " JASPAR-matched + ",
            nrow(link_data) - matched_n, " unmatched.",
            level = 2
        )
    }
    link_data <- .orderCircosInteractionRows(link_data)

    if (nrow(link_data) == 0) {
        .log_warn("No motif interactions remain after filtering. Skipping link track.")
        return(NULL)
    }
    rownames(link_data) <- NULL
    link_data
}
