# Avoid NSE warnings from R CMD check for ggplot2 aes()
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("Sample", "Beta", "Position", "x", "xend", "y", "yend", "start"))
}

#' Plot DMR Structure with DMPs and Extended CpGs
#'
#' @description Visualizes the structure of Differentially Methylated Regions (DMRs)
#' identified by findDMRsFromSeeds, showing the underlying DMPs as stem plots connected
#' by horizontal lines to form DMRs, with extended CpG regions shown as vertical lines.
#' The plot distinguishes between DMPs (differentially methylated positions), supporting
#' CpGs that extend the DMR, and non-supporting CpGs in the surrounding region.
#'
#' @param dmrs GRanges object. Output from findDMRsFromSeeds containing DMR information.
#' @param dmr_index Integer. Which DMR to plot (default: 1).
#' @param array Character. Array platform type: "450K", "27K", "EPIC", or "EPICv2" (default: "450K"). Ignored if sorted_locs is provided.
#' @param genome Character. Genome version (default: "hg19"). Ignored if sorted_locs is provided.
#' @param sorted_locs Data frame. Genomic locations sorted by position (optional). If NULL, will be fetched based on array and genome.
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
                              sorted_locs = NULL,
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
    if (is.null(sorted_locs)) {
        array <- strex::match_arg(array, ignore_case = TRUE)
        sorted_locs <- getSortedGenomicLocs(array = array, genome = genome)
    }

    # Extract DMR information
    chr <- as.character(GenomicRanges::seqnames(dmr))
    sorted_locs <- sorted_locs[sorted_locs$chr == chr, , drop = FALSE]
    dmr_start <- GenomicRanges::start(dmr)
    dmr_end <- GenomicRanges::end(dmr)
    start_cpg <- dmr_data$start_cpg
    end_cpg <- dmr_data$end_cpg
    start_cpg_ind <- which(rownames(sorted_locs) == start_cpg)
    start_cpg_pos <- sorted_locs[start_cpg_ind, "start"]
    end_cpg_ind <- which(rownames(sorted_locs) == end_cpg)
    end_cpg_pos <- sorted_locs[end_cpg_ind, "start"]
    start_dmp <- dmr_data$start_dmp
    end_dmp <- dmr_data$end_dmp
    start_dmp_ind <- which(rownames(sorted_locs) == start_dmp)
    end_dmp_ind <- which(rownames(sorted_locs) == end_dmp)

    # Extract DMP IDs from the comma-separated string
    dmp_ids <- unlist(strsplit(as.character(dmr_data$dmps), ","))

    # Get DMP positions
    dmp_positions <- sorted_locs[dmp_ids, "start"]
    start_dmp_pos <- dmr_data$start_dmp_pos
    end_dmp_pos <- dmr_data$end_dmp_pos

    if (start_dmp_ind > start_cpg_ind) {
        upstream_sup_cpgs <- sorted_locs[start_cpg_ind:(start_dmp_ind - 1), ]
    } else {
        upstream_sup_cpgs <- data.frame()
    }
    if (end_dmp_ind < end_cpg_ind) {
        downstream_sup_cpgs <- sorted_locs[(end_dmp_ind + 1):end_cpg_ind, ]
    } else {
        downstream_sup_cpgs <- data.frame()
    }
    if (extend_by_dmr_size_ratio > 0) {
        dmr_size <- dmr_end - dmr_start
        ext <- round(dmr_size * extend_by_dmr_size_ratio)
    } else {
        ext <- 0
    }
    ext <- max(ext, min_extension_bp)
    plot_start <- max(1, start_cpg_pos - ext)
    plot_end <- end_cpg_pos + ext
    nsup_cpgs <- sorted_locs[start_dmp_ind:end_dmp_ind, ]
    nsup_cpgs <- nsup_cpgs[which(nsup_cpgs$start < dmr_end & nsup_cpgs$start > dmr_start & !rownames(nsup_cpgs) %in% dmp_ids), , drop = FALSE]

    downstream_nsup_cpgs <- sorted_locs[which(sorted_locs$start > dmr_end & sorted_locs$start <= plot_end), , drop = FALSE]
    upstream_nsup_cpgs <- sorted_locs[which(sorted_locs$start < dmr_start & sorted_locs$start >= plot_start), , drop = FALSE]
    extended_nsup_cpgs <- rbind(
        nsup_cpgs,
        upstream_nsup_cpgs,
        downstream_nsup_cpgs
    )


    # Create plotting data frame
    # 1. DMPs (stem plots at y=1)
    dmp_df <- data.frame(
        cpg_id = dmp_ids,
        start = dmp_positions,
        y = 1,
        type = "DMP",
        stringsAsFactors = FALSE
    )

    # 2. DMR connection (horizontal line connecting DMPs at y=1)
    dmr_line <- data.frame(
        x = start_dmp_pos,
        xend = end_dmp_pos,
        y = 1,
        yend = 1,
        type = "DMR",
        stringsAsFactors = FALSE
    )

    # 3. Extended supporting CpGs (vertical lines at y=0.5)
    if (start_cpg_pos != start_dmp_pos) {
        dmr_upstream_line <- data.frame(
            x = start_cpg_pos,
            xend = start_dmp_pos,
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
    if (end_cpg_pos != end_dmp_pos) {
        dmr_downstream_line <- data.frame(
            x = end_dmp_pos,
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
            cpg_id = rownames(extended_sup_cpgs),
            start = extended_sup_cpgs$start,
            y = 0.5,
            type = "Extended_CpG",
            stringsAsFactors = FALSE
        )
    } else {
        extended_sup_cpgs_df <- data.frame(
            cpg_id = character(0),
            start = numeric(0),
            y = numeric(0),
            type = character(0),
            stringsAsFactors = FALSE
        )
    }
    # 4. Non-supporting CpGs in extended region
    if (nrow(extended_nsup_cpgs) > 0) {
        extended_nsup_cpgs_df <- data.frame(
            cpg_id = rownames(extended_nsup_cpgs),
            start = extended_nsup_cpgs$start,
            y = 0.5,
            type = "Extended_CpG",
            stringsAsFactors = FALSE
        )
    } else {
        extended_nsup_cpgs_df <- data.frame(
            cpg_id = character(0),
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

    # Plot DMP stems
    p <- p + ggplot2::geom_segment(
        data = dmp_df,
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

    # Plot DMP points
    p <- p + ggplot2::geom_point(
        data = dmp_df,
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
        xmin = start_dmp_pos,
        xmax = end_dmp_pos,
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
            x = c(min(upstream_sup_cpgs$start), start_dmp_pos, start_dmp_pos, min(upstream_sup_cpgs$start)),
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
            x = c(end_dmp_pos, max(downstream_sup_cpgs$start), max(downstream_sup_cpgs$start), end_dmp_pos),
            y = c(0, 0, 0.5, 1),
            alpha = 0.1,
            fill = "#E41A1C"
        )
    }

    # Create title if not provided
    title <- sprintf(
        "DMR #%d: %s:%s-%s\n%d DMPs, %d Sequence CpGs (\u0394\u03b2=%.3f)",
        dmr_index,
        chr,
        format(dmr_start, big.mark = ",", scientific = FALSE),
        format(dmr_end, big.mark = ",", scientific = FALSE),
        dmr_data$dmps_num,
        dmr_data$cpgs_num,
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
                labels = c("Array CpGs", "DMPs"),
                limits = c(-0.1, 1.15)
            )
    } else {
        p <- p +
            ggplot2::scale_y_continuous(
                breaks = c(1),
                labels = c("DMPs"),
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

    # Add ticks for the DMPs on the x-axis
    breaks <- c(plot_start, sorted_locs[start_cpg_ind:end_cpg_ind, "start"], plot_end)
    cpgs_labs <- paste0(sorted_locs[start_cpg_ind:end_cpg_ind, "start"], "\n(", rownames(sorted_locs[start_cpg_ind:end_cpg_ind, , drop = FALSE]), ")")
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
        total_shown_positions <- rbind(extended_nsup_cpgs, extended_sup_cpgs, sorted_locs[dmp_ids, , drop = FALSE])
        total_shown_positions <- total_shown_positions[order(total_shown_positions$start), ]
        return(invisible(list(structure_plot = p, breaks = breaks, breaks_labels = breaks_labels, chr = chr, total_locs = total_shown_positions)))
    }
    invisible(p)
}


# Create beta heatmap plot
.plotBetaHeatmap <- function(dmr_data, beta_data, total_shown_positions, pheno = NULL, sample_group_col = "Sample_Group") {
    cpg_ids <- rownames(total_shown_positions)
    cpg_locs <- total_shown_positions[, c("chr", "start")]


    # Mark DMPs
    dmp_ids <- unlist(strsplit(as.character(dmr_data$dmps), ","))
    is_dmp <- cpg_ids %in% dmp_ids

    # Create heatmap
    # Prepare data
    beta_data <- as.data.frame(beta_data)
    beta_data[, "CpG"] <- rownames(beta_data)
    beta_melted <- suppressWarnings(suppressMessages(reshape2::melt(beta_data, id_vars = "CpG")))
    colnames(beta_melted) <- c("CpG", "Sample", "Beta")
    beta_melted$Position <- cpg_locs[as.character(beta_melted$CpG), "start"]
    beta_melted$is_DMP <- is_dmp[match(beta_melted$CpG, cpg_ids)]
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
    # Add vertical lines for DMPs
    if (any(is_dmp)) {
        dmp_positions <- unique(beta_melted$Position[beta_melted$is_DMP])
        heatmap_plot <- heatmap_plot +
            ggplot2::geom_vline(
                xintercept = dmp_positions,
                color = "yellow",
                linetype = "dashed",
                linewidth = 0.5,
                alpha = 0.7
            )
    }
    return(heatmap_plot)
}


minmaxscale <- function(x) {
    (x - min(x)) / max(max(x) - min(x), 1e-10)
}

.plotPWM <- function(dmr, genome, array, sorted_locs, motif_flank_size = 5) {
    # Extract DMR motifs if not already present
    if (!"pwm" %in% colnames(S4Vectors::mcols(dmr))) {
        dmr <- extractDMRMotifs(dmr,
            genome = genome, array = array,
            genomic_locs = sorted_locs, flank_size = motif_flank_size
        )
    }
    pwm <- mcols(dmr)$pwm[[1]]
    
    if (is.null(pwm) || !is.matrix(pwm)) {
        return(NULL)
    }

    n_positions <- ncol(pwm)
    if (n_positions == 0) {
        return(NULL)
    }

    position_labels <- c(seq(-motif_flank_size, 0), seq(0, motif_flank_size))

    rownames(pwm) <- BASE_LEVELS

    base_colors <- c("A" = "#109648", "C" = "#255C99", "G" = "#F7B32B", "T" = "#D62839")

    pwm_plot <- ggseqlogo::ggseqlogo(pwm, method = "custom", seq_type = "dna") +
        ggplot2::scale_fill_manual(values = base_colors) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = 10),
            axis.text.y = ggplot2::element_text(size = 10),
            axis.title = ggplot2::element_text(size = 11, face = "bold"),
            plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5)
        ) +
        ggplot2::labs(
            x = "Position Relative to CpG",
            y = "Information Content (bits)",
            title = "Supporting CpG Motif PWM"
        ) +
        ggplot2::scale_x_continuous(breaks = 1:n_positions, labels = as.character(position_labels))

    return(pwm_plot)
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
#' @param array Character. Array platform type (default: "450K").
#' @param ncol Integer. Number of columns in the grid (default: 1).
#' @param ... Additional arguments passed to .plotDMRStructure or plotDMRWithBeta.
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
                     genome = c("hg19", "hg38", "mm10", "mm39"),
                     array = c("450K", "27K", "EPIC", "EPICv2"),
                     ncol = 1,
                     ...) {
    showtext::showtext_auto()
    array <- strex::match_arg(array, ignore_case = TRUE)
    genome <- strex::match_arg(genome, ignore_case = TRUE)
    if (is.null(dmr_indices)) {
        score <- minmaxscale(abs(dmrs$delta_beta))
        ord <- order(score, decreasing = TRUE)
        dmr_indices <- ord[seq_len(min(top_n, length(dmrs)))]
    }

    # Create individual plots
    if (!is.null(beta)) {
        if (is.character(beta)) {
            beta <- getBetaHandler(
                beta = beta,
                array = array,
                genome = genome
            )
        }
        plot_list <- lapply(seq_along(dmr_indices), function(idx) {
            i <- dmr_indices[idx]
            plotDMR(
                dmrs = dmrs,
                dmr_index = i,
                beta = beta,
                pheno = pheno,
                sample_group_col = sample_group_col,
                ...
            )
        })
        invisible(plot_list)
    }
}


#' Plot DMR
#'
#' @description Creates a detailed DMR plot with an integrated heatmap showing
#' beta values across samples for DMPs and surrounding CpGs. The plot consists of
#' two panels: the top panel shows the DMR structure with DMPs and extended CpGs,
#' and the bottom panel displays a heatmap of beta values for all samples, if beta values are provided.
#' Additionally, if motif information is available or can be extracted, a sequence logo
#' plot is added showing the nucleotide composition and information content around CpG sites in the DMR.
#'
#' @param dmrs GRanges object. Output from findDMRsFromSeeds.
#' @param dmr_index Integer. Which DMR to plot.
#' @param beta BetaHandler object, character path to beta file, or beta values matrix.
#'   If a character path or matrix is provided, a BetaHandler will be created automatically.
#' @param pheno Data frame or character path to phenotype file. Sample information with rownames matching beta column names (required).
#' @param sorted_locs Data frame. Genomic locations sorted by position (optional).
#' @param array Character. Array platform type (default: "450K").
#' @param genome Character. Genome version (default: "hg19").
#' @param sample_group_col Character. Column in pheno for sample grouping (default: "Sample_Group").
#' @param extend_by_dmr_size_ratio Numeric. Ratio of the DMR width to extend the plot region outside of the DMR in both sides (default: 0.2).
#' @param min_extension_bp Integer. Minimum extension in base pairs (default: 50).
#' @param max_cpgs Integer. Maximum number of CpGs to show in heatmap (default: 100).
#' @param plot_motif Logical. Whether to plot the sequence logo motif (default: TRUE).
#' @param motif_flank_size Integer. Number of base pairs to include as flanking regions around each CpG site for motif extraction (default: 5).
#' @param plot_title Logical. Whether to display the title on the plot. If FALSE, the title is shown in the logs (default: TRUE).
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
                    sorted_locs = NULL,
                    array = c("450K", "27K", "EPIC", "EPICv2"),
                    genome = c("hg19", "hg38", "mm10", "mm39"),
                    sample_group_col = "Sample_Group",
                    extend_by_dmr_size_ratio = 0.2,
                    min_extension_bp = 50,
                    max_cpgs = 100,
                    plot_motif = TRUE,
                    motif_flank_size = 5,
                    plot_title = TRUE) {
    showtext::showtext_auto()
    array <- strex::match_arg(array, ignore_case = TRUE)
    genome <- strex::match_arg(genome, ignore_case = TRUE)

    # Create BetaHandler if a file path or matrix was provided
    if (!is.null(beta)) {
        if ((is.character(beta) && length(beta) == 1 && file.exists(beta)) || is.matrix(beta) || is.data.frame(beta)) {
            beta_handler <- getBetaHandler(
                beta = beta,
                array = array,
                genome = genome
            )
        } else if (!"BetaHandler" %in% class(beta)) {
            stop("beta_handler must be either a file path (character) or a BetaHandler object")
        } else {
            beta_handler <- beta
        }
    }

    if (!is.null(pheno) && !(sample_group_col %in% colnames(pheno))) {
        stop(sprintf("sample_group_col '%s' not found in pheno data frame", sample_group_col))
    }
    dmr <- dmrs[dmr_index]
    dmr_data <- S4Vectors::mcols(dmr)

    if (!is.null(pheno)) {
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


    # Create structure plot
    ret <- .plotDMRStructure(
        dmrs = dmrs,
        dmr_index = dmr_index,
        array = array,
        genome = genome,
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
        beta_data <- beta_handler$getBeta(
            row_names = rownames(total_shown_positions),
            col_names = if (is.null(pheno)) NULL else rownames(pheno)
        )
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
        warning("Beta values not provided. Only the DMR structure plot will be returned.")
        grobs <- list(ggplot2::ggplotGrob(structure_plot))
    }

    if (plot_motif) {
        pwm_plot <- .plotPWM(dmr, genome = genome, array = array, sorted_locs = sorted_locs, motif_flank_size = motif_flank_size)
        if (!is.null(pwm_plot)) {
            grobs <- c(grobs, list(ggplot2::ggplotGrob(pwm_plot)))
        }
    }

    max_width <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
    combined <- do.call(gridExtra::gtable_rbind, grobs)
    combined$widths <- max_width
    grid::grid.draw(combined)
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
#' @param array Character. Array platform type (default: "450K").
#' @param genome Character. Genome version (e.g., "hg19", "hg38", "mm10", "mm39").
#' @param sample_group_col Character. Column in pheno for sample grouping (default: "Sample_Group").
#' @param min_sim Numeric. Minimum motifs PWM similarity threshold for considering DMRs are related (default: 0.7).
#' @param flank_size Integer. Flanking region size for motif extraction in bp (default: 5).
#' @param max_cpgs_per_dmr Integer. Maximum number of CpGs to show per DMR in scatter/heatmap (default: 100).
#' @param max_interactions Integer. Maximum number of interactions to plot (default: 30).
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
                           array = "450K",
                           sorted_locs = NULL,
                           sample_group_col = "Sample_Group",
                           min_sim = 0.7,
                           flank_size = 5,
                           max_num_samples = 20, # @TODO implement that
                           max_dmrs_per_chr = 5,
                           max_cpgs_per_dmr = 5,
                           max_interactions = 30,
                           ...) {
    dmrs <- convertToGRanges(dmrs, genome)

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
    if (is.null(sorted_locs)) {
        sorted_locs <- beta_handler$getGenomicLocs()
    }

    .log_step("Preparing data for Circos plot...")

    .log_step("Filtering DMRs for plotting by maximum absolute delta beta...", level = 2)
    dmrs <- .filterDMRsForCircos(dmrs, max_dmrs_per_chr)
    .log_success("DMRs filtered for Circos plot", level = 2)

    .log_info("Total DMRs to plot: ", length(dmrs), level = 2)

    cytoband <- .getCytobandData(genome)

    .log_step("Preparing DMRs data...", level = 2)
    arc_data <- .prepareCircosArcData(dmrs)
    .log_success("DMR arcs data prepared", level = 2)
    if (verbose >= 3) {
        dir.create("debug", showWarnings = FALSE)
        write.table(arc_data, file = "debug/circos_arc_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    }

    .log_step("Preparing heatmap data...", level = 2)
    heatmap_data <- .prepareCircosHeatmapData(
        dmrs, beta_handler, pheno, sample_group_col,
        sorted_locs, max_cpgs_per_dmr
    )
    .log_success("Heatmap data prepared", level = 2)
    .log_info("Total heatmap entries: ", nrow(heatmap_data), level = 2)
    if (verbose >= 3) {
        dir.create("debug", showWarnings = FALSE)
        write.table(heatmap_data, file = "debug/circos_heatmap_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    }

    .log_step("Computing motif-based DMR interactions...", level = 2)
    link_data <- .prepareCircosLinkData(
        dmrs, genome, array, min_sim, flank_size, sorted_locs
    )
    .log_success("DMR interactions data prepared", level = 2)
    if (verbose >= 3) {
        dir.create("debug", showWarnings = FALSE)
        write.table(link_data, file = "debug/circos_link_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    }

    .log_step("Creating Circos plot...")

    unique_chrs <- unique(as.character(GenomicRanges::seqnames(dmrs)))

    if (!is.null(cytoband)) {
        cytoband_subset <- cytoband[cytoband$V1 %in% unique_chrs, ]
        circlize::circos.initializeWithIdeogram(cytoband = cytoband_subset, plotType = c("ideogram", "labels"))
    } else {
        circlize::circos.initializeWithIdeogram(species = genome, plotType = c("ideogram", "labels"))
    }

    if (!is.null(arc_data)) {
        .log_step("Adding DMR arc track...", level = 2)

        positive_delta <- grDevices::colorRampPalette(c("white", "#801414"))(100)
        negative_delta <- grDevices::colorRampPalette(c("#055709", "white"))(100)

        arc_colors <- ifelse(
            arc_data$delta_beta > 0,
            positive_delta[pmin(100, ceiling(abs(arc_data$delta_beta) * 100))],
            negative_delta[pmin(100, ceiling(abs(arc_data$delta_beta) * 100))]
        )

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
            bg.col = "#ffe7c1",
            track.height = 0.05,
            ylim = c(0, 1),
            panel.fun = function(x, y) {
                chr <- circlize::CELL_META$sector.index
                dmr_chr <- arc_df[arc_df$chr == chr, ]
                if (nrow(dmr_chr) > 0) {
                    for (i in seq_len(nrow(dmr_chr))) {
                        circlize::circos.rect(
                            xleft = dmr_chr$start[i],
                            xright = dmr_chr$end[i],
                            ybottom = 0,
                            ytop = 1,
                            col = dmr_chr$color[i],
                            border = dmr_chr$color[i]
                        )
                    }
                }
            }
        )
        .log_success("Arc track added", level = 2)
    }

    if (!is.null(heatmap_data) && nrow(heatmap_data) > 0) {
        .log_step("Adding heatmap track...", level = 2)
        col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

        circlize::circos.genomicHeatmap(
            bed = heatmap_data,
            col = col_fun,
            border = "white",
            side = "inside",
            border_lwd = 0.2,
            line_lwd = 0.2,
            heatmap_height = 0.2
        )
        .log_success("Heatmap track added", level = 2)
    }

    if (!is.null(link_data) && nrow(link_data) > 0) {
        .log_step("Adding link track...", level = 2)
        link_data <- link_data[order(link_data$sim, decreasing = TRUE), ]
        if (nrow(link_data) > max_interactions) {
            link_data <- link_data[1:max_interactions, ]
            .log_info("Limiting to top ", max_interactions, " interactions based on similarity", level = 2)
        }

        link_colors <- circlize::colorRamp2(
            c(min(link_data$sim), median(link_data$sim), max(link_data$sim)),
            c("#89adc2", "#62b0f0", "#1696e0")
        )

        for (i in seq_len(nrow(link_data))) {
            point1 <- c(link_data$start1[i], link_data$end1[i])
            point2 <- c(link_data$start2[i], link_data$end2[i])
            if (point1[2] - point1[1] < 1e5) {
                point1 <- mean(point1)
            }
            if (point2[2] - point2[1] < 1e5) {
                point2 <- mean(point2)
            }
            circlize::circos.link(
                sector.index1 = link_data$chr1[i],
                point1 = point1,
                sector.index2 = link_data$chr2[i],
                point2 = point2,
                col = link_colors(link_data$sim[i]),
                border = NA
            )
        }
        .log_success("Link track added", level = 2)
    }

    .log_success("Circos plot created successfully")

    circlize::circos.clear()

    invisible(NULL)
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
                    download.file(url, temp_file, quiet = TRUE, mode = "wb")
                    TRUE
                },
                error = function(e) {
                    .log_warn("cytoBandIdeo not found, trying cytoBand table...")
                    FALSE
                }
            )

            if (!result) {
                url <- paste0("https://hgdownload.cse.ucsc.edu/goldenpath/", genome, "/database/cytoBand.txt.gz")
                download.file(url, temp_file, quiet = TRUE, mode = "wb")
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
            return(NULL)
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

.prepareCircosHeatmapData <- function(dmrs, beta_handler, pheno, sample_group_col, sorted_locs, max_sup_cpgs_per_dmr_side = 2) {
    beta_col_names <- beta_handler$getBetaColNames()
    pheno <- pheno[rownames(pheno) %in% beta_col_names, , drop = FALSE]

    if (nrow(pheno) == 0) {
        .log_warn("No samples in pheno match the samples in beta values. Skipping heatmap track.")
        return(NULL)
    }
    available_cpgs <- beta_handler$getBetaRowNames()
    dmrs_cpgs_list <- getSupportingSites(
        dmrs, 
        available_cpgs,
        max_sup_cpgs_per_dmr_side = max_sup_cpgs_per_dmr_side,
        ret_index = TRUE,
        separate_by_section = FALSE
    )
    dmrs_cpgs <- unlist(dmrs_cpgs_list)


    shown_locs <- sorted_locs[dmrs_cpgs, c("chr", "start", "end"), drop = FALSE]
    beta_data <- beta_handler$getBeta(
        row_names = dmrs_cpgs,
        col_names = rownames(pheno)
    )
    cbind(shown_locs, beta_data)
}


.prepareCircosArcData <- function(dmrs) {
    list(
        chr = as.character(GenomicRanges::seqnames(dmrs)),
        start = GenomicRanges::start(dmrs),
        end = GenomicRanges::end(dmrs),
        delta_beta = abs(S4Vectors::mcols(dmrs)$delta_beta)
    )
}


.prepareCircosLinkData <- function(dmrs, genome, array, min_sim, flank_size, sorted_locs) {
    ret <- tryCatch(
        {
            computeDMRsInteraction(
                dmrs = dmrs,
                genome = genome,
                array = array,
                min_sim = min_sim,
                genomic_locs = sorted_locs,
                flank_size = flank_size
            )
        },
        error = function(e) {
            .log_warn("Failed to compute motif-based interactions: ", e$message)
            NULL
        }
    )

    if (is.null(ret) || nrow(ret$interactions) == 0) {
        .log_warn("No significant interactions found at similarity >=", min_sim, ". Skipping link track.")
        return(NULL)
    }

    .log_success("Found ", nrow(ret$interactions), " motif-based interactions")
    ret$interactions
}
