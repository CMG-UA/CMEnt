# Avoid NSE warnings from R CMD check for ggplot2 aes()
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("Sample", "Beta", "Position", "x", "xend", "y", "yend", "pos"))
}

#' Plot DMR Structure with DMPs and Extended CpGs
#'
#' @description Visualizes the structure of Differentially Methylated Regions (DMRs)
#' identified by findDMRsFromSeeds, showing the underlying DMPs as stem plots connected
#' by horizontal lines to form DMRs, with extended CpG regions shown as vertical lines.
#'
#' @param dmrs GRanges object. Output from findDMRsFromSeeds containing DMR information.
#' @param dmr_index Integer. Which DMR to plot (default: 1).
#' @param beta_file Character. Path to the methylation beta values file (optional).
#' @param sorted_locs Data frame. Genomic locations sorted by position (optional).
#' @param array Character. Array platform type (default: "450K").
#' @param genome Character. Genome version (default: "hg19").
#' @param show_pvalues Logical. Whether to show p-values on DMP stems (default: TRUE).
#' @param title Character. Plot title (default: NULL, auto-generated).
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' # Load DMR results
#' dmrs <- readRDS("dmrs.rds")
#'
#' # Plot first DMR
#' plotDMR(dmrs, dmr_index = 1)
#'
#' # Plot with custom options
#' plotDMR(dmrs, dmr_index = 2, show_pvalues = TRUE)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_text theme_minimal
#' @importFrom ggplot2 labs scale_color_gradient2 theme element_text scale_y_continuous
#' @importFrom GenomicRanges mcols start end seqnames
#' @export
plotDMR <- function(dmrs,
                    dmr_index = 1,
                    beta_file = NULL,
                    sorted_locs = NULL,
                    array = c("450K", "27K", "EPIC", "EPICv2"),
                    genome = c("hg19", "hg38", "mm10", "mm39"),
                    show_pvalues = TRUE,
                    title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting. Please install it.")
    }

    array <- match.arg(array)
    genome <- match.arg(genome)

    # Validate input
    if (dmr_index < 1 || dmr_index > length(dmrs)) {
        stop("dmr_index must be between 1 and ", length(dmrs))
    }

    # Extract the specific DMR
    dmr <- dmrs[dmr_index]
    dmr_data <- S4Vectors::mcols(dmr)

    # Get genomic locations if not provided
    if (is.null(sorted_locs)) {
        sorted_locs <- getSortedGenomicLocs(array = array, genome = genome)
    }

    # Extract DMR information
    chr <- as.character(GenomicRanges::seqnames(dmr))
    dmr_start <- GenomicRanges::start(dmr)
    dmr_end <- GenomicRanges::end(dmr)
    start_cpg <- dmr_data$start_cpg
    end_cpg <- dmr_data$end_cpg
    start_dmp <- dmr_data$start_dmp
    end_dmp <- dmr_data$end_dmp
    if (start_dmp > start_cpg) {
        upstream_sup_cpgs <- sorted_locs[start_cpg - 1:start_dmp, ]
    } else {
        upstream_sup_cpgs <- data.frame()
    }
    if (end_dmp < end_cpg) {
        downstream_sup_cpgs <- sorted_locs[end_dmp + 1:end_cpg, ]
    } else {
        downstream_sup_cpgs <- data.frame()
    }

    # Extract DMP IDs from the comma-separated string
    dmp_ids <- unlist(strsplit(as.character(dmr_data$dmps), ","))

    # Get all CpGs in the extended region
    start_ind <- dmr_data$start_ind
    end_ind <- dmr_data$end_ind

    # Get DMP positions
    dmp_positions <- sorted_locs[dmp_ids, "pos"]
    start_dmp_pos <- dmr_data$start_dmp_pos
    end_dmp_pos <- dmr_data$end_dmp_pos

    # Create plotting data frame
    # 1. DMPs (stem plots at y=1)
    dmp_df <- data.frame(
        cpg_id = dmp_ids,
        pos = dmp_positions,
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

    # 3. Extended CpGs (vertical lines at y=0.5)
    extended_cpgs <- rbind(
        upstream_sup_cpgs,
        downstream_sup_cpgs
    )
    if (nrow(extended_cpgs) > 0) {
        extended_cpgs_df <- data.frame(
            cpg_id = rownames(extended_cpgs),
            pos = extended_cpgs$pos,
            y = 0.5,
            type = "Extended_CpG",
            stringsAsFactors = FALSE
        )
    } else {
        extended_cpgs_df <- data.frame(
            cpg_id = character(0),
            pos = numeric(0),
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
        ggplot2::aes(x = pos, xend = pos, y = 0, yend = y),
        color = "#377EB8",
        linewidth = 1.2,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"), type = "closed")
    )

    # Plot CpG stems
    if (nrow(extended_cpgs_df) > 0) {
        p <- p + ggplot2::geom_segment(
            data = extended_cpgs_df,
            ggplot2::aes(x = pos, xend = pos, y = 0, yend = y),
            color = "gray50",
            linewidth = 0.8,
            alpha = 0.7
        )
    }

    # Plot DMP points
    p <- p + ggplot2::geom_point(
        data = dmp_df,
        ggplot2::aes(x = pos, y = y),
        color = "#377EB8",
        size = 3,
        shape = 16
    )

    # Plot CpG points
    if (nrow(extended_cpgs_df) > 0) {
        p <- p + ggplot2::geom_point(
            data = extended_cpgs_df,
            ggplot2::aes(x = pos, y = y),
            color = "gray50",
            size = 2,
            shape = 16,
            alpha = 0.7
        )
    }   

    # Add DMR region shading
    p <- p + ggplot2::annotate(
        "rect",
        xmin = start_dmp_pos,
        xmax = end_dmp_pos,
        ymin = -0.05,
        ymax = 1.1,
        alpha = 0.1,
        fill = "#E41A1C"
    )
    # if upstream extended CpGs exist add shading in the form of a trapezoid
    if (nrow(upstream_sup_cpgs) > 0) {
        p <- p + ggplot2::annotate(
            "polygon",
            x = c(min(upstream_sup_cpgs$pos), start_dmp_pos, start_dmp_pos, min(upstream_sup_cpgs$pos)),
            y = c(0, 0, 1, 0.5),
            alpha = 0.1,
            fill = "gray50"
        )
    }
    # if downstream extended CpGs exist add shading in the form of a trapezoid
    if (nrow(downstream_sup_cpgs) > 0) {
        p <- p + ggplot2::annotate(
            "polygon",
            x = c(end_dmp_pos, max(downstream_sup_cpgs$pos), max(downstream_sup_cpgs$pos), end_dmp_pos),
            y = c(0, 0, 0.5, 1),
            alpha = 0.1,
            fill = "gray50"
        )
    }
    # Add extended region shading
    p <- p + ggplot2::annotate(
        "rect",
        xmin = dmr_start,
        xmax = dmr_end,
        ymin = -0.05,
        ymax = 0.9,
        alpha = 0.05,
        fill = "gray50"
    )

    # Create title if not provided
    if (is.null(title)) {
        title <- sprintf(
            "DMR #%d: %s:%s-%s\n%d DMPs, %d CpGs (delta_beta=%.3f, p=%.2e)",
            dmr_index,
            chr,
            format(dmr_start, big.mark = ",", scientific = FALSE),
            format(dmr_end, big.mark = ",", scientific = FALSE),
            dmr_data$dmps_num,
            dmr_data$cpgs_num,
            dmr_data$delta_beta,
            dmr_data$pval_adj
        )
    }

    # Styling
    p <- p +
        ggplot2::scale_y_continuous(
            breaks = c(0, 0.5, 1),
            labels = c("Extended\nCpGs", "", "DMPs/DMR"),
            limits = c(-0.1, 1.15)
        ) +
        ggplot2::labs(
            title = title,
            x = sprintf("Genomic Position on %s (bp)", chr),
            y = ""
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 11),
            axis.text.y = ggplot2::element_text(hjust = 0.5),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank()
        )

    # Add legend manually using annotation
    legend_x <- dmr_start + (dmr_end - dmr_start) * 0.02
    legend_y_start <- 0.5

    p <- p +
        ggplot2::annotate("text", x = legend_x, y = legend_y_start + 0.45,
            label = "Legend:", hjust = 0, fontface = "bold", size = 3.5) +
        ggplot2::annotate("segment", x = legend_x, xend = legend_x + (dmr_end - dmr_start) * 0.03,
            y = legend_y_start + 0.35, yend = legend_y_start + 0.35,
            color = "#377EB8", linewidth = 1.2,
            arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm"), type = "closed")) +
        ggplot2::annotate("text", x = legend_x + (dmr_end - dmr_start) * 0.04,
            y = legend_y_start + 0.35, label = "DMP (seed site)", hjust = 0, size = 3) +
        ggplot2::annotate("segment", x = legend_x, xend = legend_x + (dmr_end - dmr_start) * 0.03,
            y = legend_y_start + 0.25, yend = legend_y_start + 0.25,
            color = "#E41A1C", linewidth = 1.5) +
        ggplot2::annotate("text", x = legend_x + (dmr_end - dmr_start) * 0.04,
            y = legend_y_start + 0.25, label = "DMR (connected region)", hjust = 0, size = 3) +
        ggplot2::annotate("segment", x = legend_x + (dmr_end - dmr_start) * 0.015,
            xend = legend_x + (dmr_end - dmr_start) * 0.015,
            y = legend_y_start + 0.12, yend = legend_y_start + 0.18,
            color = "gray70", linewidth = 0.5, alpha = 0.7) +
        ggplot2::annotate("text", x = legend_x + (dmr_end - dmr_start) * 0.04,
            y = legend_y_start + 0.15, label = "Extended CpG", hjust = 0, size = 3)

    return(p)
}


#' Plot Multiple DMRs in a Grid
#'
#' @description Creates a grid of DMR plots for multiple regions.
#'
#' @param dmrs GRanges object. Output from findDMRsFromSeeds.
#' @param dmr_indices Integer vector. Which DMRs to plot (default: 1:min(4, length(dmrs))).
#' @param ncol Integer. Number of columns in the grid (default: 2).
#' @param ... Additional arguments passed to plotDMR.
#'
#' @return A combined ggplot2 object (requires patchwork or gridExtra).
#'
#' @examples
#' \dontrun{
#' dmrs <- readRDS("dmrs.rds")
#' plotDMRs(dmrs, dmr_indices = 1:6, ncol = 3)
#' }
#'
#' @export
plotDMRs <- function(dmrs,
                     dmr_indices = NULL,
                     ncol = 2,
                     ...) {
    if (is.null(dmr_indices)) {
        dmr_indices <- 1:min(4, length(dmrs))
    }

    # Create individual plots
    plot_list <- lapply(dmr_indices, function(i) {
        plotDMR(dmrs, dmr_index = i, ...)
    })

    # Try to use patchwork if available, otherwise gridExtra, otherwise return list
    if (requireNamespace("patchwork", quietly = TRUE)) {
        return(patchwork::wrap_plots(plot_list, ncol = ncol))
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
        return(gridExtra::grid.arrange(grobs = plot_list, ncol = ncol))
    } else {
        message("Install 'patchwork' or 'gridExtra' package for combined plots. Returning list of plots.")
        return(plot_list)
    }
}


#' Plot DMR with Beta Values Heatmap
#'
#' @description Creates a detailed DMR plot with an integrated heatmap showing
#' beta values across samples for DMPs and surrounding CpGs.
#'
#' @param dmrs GRanges object. Output from findDMRsFromSeeds.
#' @param dmr_index Integer. Which DMR to plot.
#' @param beta_handler BetaFileHandler object OR character path to beta file.
#'   If a character path is provided, a BetaFileHandler will be created automatically.
#' @param pheno Data frame. Phenotype data with sample information (required).
#' @param sorted_locs Data frame. Genomic locations (optional).
#' @param array Character. Array platform type (default: "450K").
#' @param genome Character. Genome version (default: "hg19").
#' @param sample_group_col Character. Column in pheno for sample grouping.
#' @param max_cpgs Integer. Maximum number of CpGs to show in heatmap (default: 100).
#'
#' @return A combined plot object.
#'
#' @examples
#' \dontrun{
#' # Using BetaFileHandler
#' beta_handler <- BetaFileHandler$new(beta_file = "beta.txt", array = "450K", genome = "hg19")
#' plotDMRWithBeta(dmrs, 1, beta_handler = beta_handler, pheno = pheno_df)
#'
#' # Or using a file path (handler created automatically)
#' plotDMRWithBeta(dmrs, 1, beta_handler = "beta.txt", pheno = pheno_df)
#' }
#'
#' @export
plotDMRWithBeta <- function(dmrs,
                            dmr_index,
                            beta_handler,
                            pheno,
                            sorted_locs = NULL,
                            array = c("450K", "27K", "EPIC", "EPICv2"),
                            genome = c("hg19", "hg38", "mm10", "mm39"),
                            sample_group_col = "Sample_Group",
                            max_cpgs = 100) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required. Please install it.")
    }
    if (!requireNamespace("reshape2", quietly = TRUE)) {
        stop("Package 'reshape2' is required. Please install it.")
    }

    array <- match.arg(array)
    genome <- match.arg(genome)

    # Create BetaFileHandler if a file path was provided
    if (is.character(beta_handler)) {
        beta_handler <- BetaFileHandler$new(
            beta_file = beta_handler,
            array = array,
            genome = genome,
            verbose = 0
        )
    } else if (!inherits(beta_handler, "BetaFileHandler")) {
        stop("beta_handler must be either a file path (character) or a BetaFileHandler object")
    }

    # Get genomic locations
    if (is.null(sorted_locs)) {
        sorted_locs <- getSortedGenomicLocs(array = array, genome = genome)
    }

    # Update the handler's sorted_locs if not already set
    if (is.null(beta_handler$sorted_locs)) {
        beta_handler$sorted_locs <- sorted_locs
    }

    # Extract DMR
    dmr <- dmrs[dmr_index]
    dmr_data <- S4Vectors::mcols(dmr)

    # Get CpG IDs in the region
    start_ind <- dmr_data$start_ind
    end_ind <- dmr_data$end_ind

    # Subsample if too many CpGs
    if ((end_ind - start_ind + 1) > max_cpgs) {
        step <- ceiling((end_ind - start_ind + 1) / max_cpgs)
        cpg_indices <- seq(start_ind, end_ind, by = step)
    } else {
        cpg_indices <- start_ind:end_ind
    }

    cpg_ids <- rownames(sorted_locs)[cpg_indices]
    cpg_locs <- sorted_locs[cpg_ids, ]

    # Read beta values using BetaFileHandler
    beta_data <- beta_handler$getBeta(
        row_names = cpg_ids,
        col_names = rownames(pheno)
    )

    # Mark DMPs
    dmp_ids <- unlist(strsplit(as.character(dmr_data$dmps), ","))
    is_dmp <- cpg_ids %in% dmp_ids

    # Create structure plot
    structure_plot <- plotDMR(
        dmrs = dmrs,
        dmr_index = dmr_index,
        sorted_locs = sorted_locs,
        array = array,
        genome = genome,
        show_pvalues = FALSE,
        title = NULL
    )

    # Create heatmap
    # Prepare data
    beta_melted <- reshape2::melt(beta_data)
    colnames(beta_melted) <- c("CpG", "Sample", "Beta")
    beta_melted$Position <- cpg_locs[as.character(beta_melted$CpG), "pos"]
    beta_melted$is_DMP <- is_dmp[match(beta_melted$CpG, cpg_ids)]
    beta_melted$Group <- pheno[as.character(beta_melted$Sample), sample_group_col]

    # Order samples by group
    sample_order <- rownames(pheno)[order(pheno[[sample_group_col]])]
    beta_melted$Sample <- factor(beta_melted$Sample, levels = sample_order)

    heatmap_plot <- ggplot2::ggplot(beta_melted) +
        ggplot2::geom_tile(ggplot2::aes(x = Position, y = Sample, fill = Beta)) +
        ggplot2::scale_fill_gradient2(
            low = "blue", mid = "white", high = "red",
            midpoint = 0.5,
            limits = c(0, 1),
            name = "Beta Value"
        ) +
        ggplot2::labs(
            x = "Genomic Position (bp)",
            y = "Sample",
            title = sprintf("Beta Values for DMR #%d", dmr_index)
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
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

    # Combine plots
    if (requireNamespace("patchwork", quietly = TRUE)) {
        combined <- structure_plot / heatmap_plot +
            patchwork::plot_layout(heights = c(1, 2))
        return(combined)
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
        return(gridExtra::grid.arrange(structure_plot, heatmap_plot, nrow = 2, heights = c(1, 2)))
    } else {
        message("Install 'patchwork' or 'gridExtra' for combined plots. Returning list.")
        return(list(structure = structure_plot, heatmap = heatmap_plot))
    }
}
