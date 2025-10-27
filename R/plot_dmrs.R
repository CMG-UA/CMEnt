# Avoid NSE warnings from R CMD check for ggplot2 aes()
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("Sample", "Beta", "Position", "x", "xend", "y", "yend", "pos"))
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
#' @param sorted_locs Data frame. Genomic locations sorted by position (optional). If NULL, will be fetched based on array and genome.
#' @param pval_col Character. The name of the column containing the aggregated p-value used in the title (default: "pval_adj").
#' @param array Character. Array platform type: "450K", "27K", "EPIC", or "EPICv2" (default: "450K").
#' @param genome Character. Genome version: "hg19", "hg38", "mm10", or "mm39" (default: "hg19").
#' @param extend_by_dmr_size_ratio Numeric. Ratio of the DMR width to extend the plot region outside of the DMR on both sides (default: 0.2).
#' @param min_extension_bp Integer. Minimum extension in base pairs for the plot region (default: 50).
#' @param plot_title Logical. Whether to display the title on the plot. If FALSE, the title is logged instead (default: TRUE).
#' @param .ret_details Logical. Internal parameter to return additional details (breaks, labels, chromosome, locations) for use by plotDMRWithBeta (default: FALSE).
#'
#' @return A ggplot2 object showing the DMR structure. If .ret_details is TRUE, returns a list containing the plot and additional information.
#'
#' @examples
#' \dontrun{
#' # Load DMR results
#' dmrs <- readRDS("dmrs.rds")
#'
#' # Plot first DMR with title
#' plotDMR(dmrs, dmr_index = 1)
#'
#' # Plot without title (title will be logged)
#' plotDMR(dmrs, dmr_index = 1, plot_title = FALSE)
#'
#' # Plot with custom p-value column
#' plotDMR(dmrs, dmr_index = 1, pval_col = "pval_raw")
#'
#' # Plot with extended region
#' plotDMR(dmrs, dmr_index = 1, extend_by_dmr_size_ratio = 0.5, min_extension_bp = 100)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_text theme_minimal
#' @importFrom ggplot2 labs scale_color_gradient2 theme element_text scale_y_continuous
#' @importFrom GenomicRanges mcols start end seqnames
#' @export
plotDMR <- function(dmrs,
                    dmr_index = 1,
                    sorted_locs = NULL,
                    pval_col = "pval_adj",
                    array = c("450K", "27K", "EPIC", "EPICv2"),
                    genome = c("hg19", "hg38", "mm10", "mm39"),
                    extend_by_dmr_size_ratio = 0.2,
                    min_extension_bp = 50,
                    plot_title = TRUE,
                    .ret_details = FALSE) {
    showtext::showtext_auto()
    array <- strex::match_arg(array, ignore_case = TRUE)
    genome <- strex::match_arg(genome, ignore_case = TRUE)

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
    sorted_locs <- sorted_locs[sorted_locs$chr == chr, , drop = FALSE]
    dmr_start <- GenomicRanges::start(dmr)
    dmr_end <- GenomicRanges::end(dmr)
    start_cpg <- dmr_data$start_cpg
    end_cpg <- dmr_data$end_cpg
    start_cpg_ind <- which(rownames(sorted_locs) == start_cpg)
    start_cpg_pos <- sorted_locs[start_cpg_ind, "pos"]
    end_cpg_ind <- which(rownames(sorted_locs) == end_cpg)
    end_cpg_pos <- sorted_locs[end_cpg_ind, "pos"]
    start_dmp <- dmr_data$start_dmp
    end_dmp <- dmr_data$end_dmp
    start_dmp_ind <- which(rownames(sorted_locs) == start_dmp)
    end_dmp_ind <- which(rownames(sorted_locs) == end_dmp)

    # Extract DMP IDs from the comma-separated string
    dmp_ids <- unlist(strsplit(as.character(dmr_data$dmps), ","))

    # Get DMP positions
    dmp_positions <- sorted_locs[dmp_ids, "pos"]
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
    nsup_cpgs <- nsup_cpgs[which(nsup_cpgs$pos < dmr_end & nsup_cpgs$pos > dmr_start & !rownames(nsup_cpgs) %in% dmp_ids), , drop = FALSE]

    downstream_nsup_cpgs <- sorted_locs[which(sorted_locs$pos > dmr_end & sorted_locs$pos <= plot_end), , drop = FALSE]
    upstream_nsup_cpgs <- sorted_locs[which(sorted_locs$pos < dmr_start & sorted_locs$pos >= plot_start), , drop = FALSE]
    extended_nsup_cpgs <- rbind(
        nsup_cpgs,
        upstream_nsup_cpgs,
        downstream_nsup_cpgs
    )


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
            pos = extended_sup_cpgs$pos,
            y = 0.5,
            type = "Extended_CpG",
            stringsAsFactors = FALSE
        )
    } else {
        extended_sup_cpgs_df <- data.frame(
            cpg_id = character(0),
            pos = numeric(0),
            y = numeric(0),
            type = character(0),
            stringsAsFactors = FALSE
        )
    }
    # 4. Non-supporting CpGs in extended region
    if (nrow(extended_nsup_cpgs) > 0) {
        extended_nsup_cpgs_df <- data.frame(
            cpg_id = rownames(extended_nsup_cpgs),
            pos = extended_nsup_cpgs$pos,
            y = 0.5,
            type = "Extended_CpG",
            stringsAsFactors = FALSE
        )
    } else {
        extended_nsup_cpgs_df <- data.frame(
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

    # Plot supported CpG stems
    if (nrow(extended_sup_cpgs_df) > 0) {
        p <- p + ggplot2::geom_segment(
            data = extended_sup_cpgs_df,
            ggplot2::aes(x = pos, xend = pos, y = 0, yend = y),
            color = "#377EB8",
            linewidth = 0.8,
            alpha = 0.8
        )
    }

    # Plot non-supported CpG stems
    if (nrow(extended_nsup_cpgs_df) > 0) {
        p <- p + ggplot2::geom_segment(
            data = extended_nsup_cpgs_df,
            ggplot2::aes(x = pos, xend = pos, y = 0, yend = y),
            color = "gray50",
            linewidth = 0.8,
            alpha = 0.5
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

    # Plot supporting CpG points
    if (nrow(extended_sup_cpgs_df) > 0) {
        p <- p + ggplot2::geom_point(
            data = extended_sup_cpgs_df,
            ggplot2::aes(x = pos, y = y),
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
            ggplot2::aes(x = pos, y = y),
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
            x = c(min(upstream_sup_cpgs$pos), start_dmp_pos, start_dmp_pos, min(upstream_sup_cpgs$pos)),
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
            x = c(end_dmp_pos, max(downstream_sup_cpgs$pos), max(downstream_sup_cpgs$pos), end_dmp_pos),
            y = c(0, 0, 0.5, 1),
            alpha = 0.1,
            fill = "#E41A1C"
        )
    }

    # Create title if not provided
    title <- sprintf(
        "DMR #%d: %s:%s-%s\n%d DMPs, %d Sequence CpGs (\u0394\u03b2=%.3f, p=%.2e)",
        dmr_index,
        chr,
        format(dmr_start, big.mark = ",", scientific = FALSE),
        format(dmr_end, big.mark = ",", scientific = FALSE),
        dmr_data$dmps_num,
        dmr_data$cpgs_num,
        dmr_data$delta_beta,
        dmr_data[[pval_col]]
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
    breaks <- c(plot_start, sorted_locs[start_cpg_ind:end_cpg_ind, "pos"], plot_end)
    cpgs_labs <- paste0(sorted_locs[start_cpg_ind:end_cpg_ind, "pos"], "\n(", rownames(sorted_locs[start_cpg_ind:end_cpg_ind, , drop = FALSE]), ")")
    breaks_labels <- c(
        as.character(plot_start), cpgs_labs,
        as.character(plot_end)
    )
    p <- p + ggplot2::scale_x_continuous(
        breaks = breaks,
        labels = breaks_labels
    )
    p <- p + ggplot2::coord_cartesian(xlim = c(breaks[1], breaks[length(breaks)]))
    p <- p + ggplot2::labs(
        x = sprintf("Genomic Position on %s (bp)", chr)
    )
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

minmaxscale <- function(x) {
    (x - min(x)) / max(max(x) - min(x), 1e-10)
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
#' @param pval_col Character. The name of the column containing the aggregated p-value (default: "pval_adj").
#' @param sample_group_col Character. Column in pheno for sample grouping (default: "Sample_Group").
#' @param genome Character. Genome version (default: "hg19").
#' @param array Character. Array platform type (default: "450K").
#' @param ncol Integer. Number of columns in the grid (default: 1).
#' @param ... Additional arguments passed to plotDMR or plotDMRWithBeta.
#'
#' @return If beta is NULL: A combined ggplot2 object (requires patchwork or gridExtra), or a list of ggplot objects.
#'   If beta is provided: A list of combined plot objects with structure and heatmap.
#'
#' @examples
#' \dontrun{
#' # Plot structure only
#' dmrs <- readRDS("dmrs.rds")
#' plotDMRs(dmrs, dmr_indices = 1:6, ncol = 3)
#'
#' # Plot with beta values heatmap
#' plotDMRs(dmrs, top_n = 4, beta = "beta.txt", pheno = pheno_df)
#' }
#'
#' @export
plotDMRs <- function(dmrs,
                     dmr_indices = NULL,
                     top_n = 4,
                     beta = NULL,
                     pheno = NULL,
                     pval_col = "pval_adj",
                     sample_group_col = "Sample_Group",
                     genome = c("hg19", "hg38", "mm10", "mm39"),
                     array = c("450K", "27K", "EPIC", "EPICv2"),
                     ncol = 1,
                     ...) {
    showtext::showtext_auto()
    array <- strex::match_arg(array, ignore_case = TRUE)
    genome <- strex::match_arg(genome, ignore_case = TRUE)
    if (is.null(dmr_indices)) {
        score <- minmaxscale(abs(dmrs$delta_beta)) * minmaxscale(-log10(mcols(dmrs)[, pval_col] + 1e-10))
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
            plotDMRWithBeta(
                dmrs = dmrs,
                dmr_index = i,
                pval_col = pval_col,
                beta = beta,
                pheno = pheno,
                sample_group_col = sample_group_col,
                ...
            )
        })
        invisible(plot_list)
    } else {
        plot_list <- lapply(dmr_indices, function(i) {
            plotDMR(dmrs, dmr_index = i, pval_col = pval_col, ...)
        })
        invisible(gridExtra::grid.arrange(grobs = plot_list, ncol = ncol))
    }
}


#' Plot DMR with Beta Values Heatmap
#'
#' @description Creates a detailed DMR plot with an integrated heatmap showing
#' beta values across samples for DMPs and surrounding CpGs. The plot consists of
#' two panels: the top panel shows the DMR structure with DMPs and extended CpGs,
#' and the bottom panel displays a heatmap of beta values for all samples.
#'
#' @param dmrs GRanges object. Output from findDMRsFromSeeds.
#' @param dmr_index Integer. Which DMR to plot.
#' @param beta BetaHandler object, character path to beta file, or beta values matrix.
#'   If a character path or matrix is provided, a BetaHandler will be created automatically.
#' @param pheno Data frame or character path to phenotype file. Sample information with rownames matching beta column names (required).
#' @param sorted_locs Data frame. Genomic locations sorted by position (optional).
#' @param array Character. Array platform type (default: "450K").
#' @param genome Character. Genome version (default: "hg19").
#' @param pval_col Character. The name of the column containing the aggregated p-value (default: "pval_adj").
#' @param sample_group_col Character. Column in pheno for sample grouping (default: "Sample_Group").
#' @param extend_by_dmr_size_ratio Numeric. Ratio of the DMR width to extend the plot region outside of the DMR in both sides (default: 0.2).
#' @param min_extension_bp Integer. Minimum extension in base pairs (default: 50).
#' @param max_cpgs Integer. Maximum number of CpGs to show in heatmap (default: 100).
#' @param plot_title Logical. Whether to display the title on the plot. If FALSE, the title is shown in the logs (default: TRUE).
#'
#' @return A combined plot object (gtable) containing the DMR structure plot and beta values heatmap,
#'   or a list of plots if required packages are not available.
#'
#' @examples
#' \dontrun{
#' # Using BetaHandler
#' beta_handler <- getBetaHandler(beta = "beta.txt", array = "450K", genome = "hg19")
#' plotDMRWithBeta(dmrs, 1, beta = beta_handler, pheno = pheno_df)
#'
#' # Using a file path (handler created automatically)
#' plotDMRWithBeta(dmrs, 1, beta = "beta.txt", pheno = pheno_df)
#'
#' # Using a beta matrix
#' plotDMRWithBeta(dmrs, 1, beta = beta_matrix, pheno = pheno_df)
#' }
#' @export
plotDMRWithBeta <- function(dmrs,
                            dmr_index,
                            beta,
                            pheno,
                            sorted_locs = NULL,
                            array = c("450K", "27K", "EPIC", "EPICv2"),
                            genome = c("hg19", "hg38", "mm10", "mm39"),
                            pval_col = "pval_adj",
                            sample_group_col = "Sample_Group",
                            extend_by_dmr_size_ratio = 0.2,
                            min_extension_bp = 50,
                            max_cpgs = 100,
                            plot_title = TRUE) {
    showtext::showtext_auto()
    array <- strex::match_arg(array, ignore_case = TRUE)
    genome <- strex::match_arg(genome, ignore_case = TRUE)

    # Create BetaHandler if a file path or matrix was provided
    if (is.character(beta) && length(beta) == 1 && file.exists(beta) || is.matrix(beta) || is.data.frame(beta)) {
        beta_handler <- getBetaHandler(
            beta = beta,
            array = array,
            genome = genome
        )
    } else if (!"BetaHandler" %in% class(beta)) {
        stop("beta_handler must be either a file path (character) or a BetaHandler object")
    }
    if (is.character(pheno) && length(pheno) == 1 && file.exists(pheno)) {
        pheno <- read.table(pheno, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
    }
    if (!is.data.frame(pheno)) {
        stop("pheno must be a data frame or a valid file path to a phenotype table")
    }
    if (!(sample_group_col %in% colnames(pheno))) {
        stop(sprintf("sample_group_col '%s' not found in pheno data frame", sample_group_col))
    }

    # Extract DMR
    dmr <- dmrs[dmr_index]
    dmr_data <- S4Vectors::mcols(dmr)

    pheno <- pheno[rownames(pheno) %in% beta_handler$getBetaColNames(), , drop = FALSE]
    if (nrow(pheno) == 0) {
        stop("No samples in pheno match the samples in beta values")
    }


    # Create structure plot
    ret <- plotDMR(
        dmrs = dmrs,
        dmr_index = dmr_index,
        array = array,
        genome = genome,
        pval_col = pval_col,
        plot_title = plot_title,
        extend_by_dmr_size_ratio = extend_by_dmr_size_ratio,
        min_extension_bp = min_extension_bp,
        .ret_details = TRUE
    )
    structure_plot <- ret$structure_plot
    breaks <- ret$breaks
    breaks_labels <- ret$breaks_labels
    total_shown_positions <- ret$total_locs
    cpg_ids <- rownames(total_shown_positions)
    cpg_locs <- total_shown_positions[, c("chr", "pos")]

    # Read beta values using BetaHandler
    beta_data <- beta_handler$getBeta(
        row_names = rownames(total_shown_positions),
        col_names = rownames(pheno)
    )

    # Mark DMPs
    dmp_ids <- unlist(strsplit(as.character(dmr_data$dmps), ","))
    is_dmp <- cpg_ids %in% dmp_ids


    chr <- ret$chr
    # Remove x-axis title for combined plot
    structure_plot <- structure_plot +
        ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank()
        )

    # Create heatmap
    # Prepare data
    beta_data <- as.data.frame(beta_data)
    beta_data[, "CpG"] <- rownames(beta_data)
    beta_melted <- suppressWarnings(suppressMessages(reshape2::melt(beta_data, id_vars = "CpG")))
    colnames(beta_melted) <- c("CpG", "Sample", "Beta")
    beta_melted$Position <- cpg_locs[as.character(beta_melted$CpG), "pos"]
    beta_melted$is_DMP <- is_dmp[match(beta_melted$CpG, cpg_ids)]
    beta_melted$Group <- pheno[as.character(beta_melted$Sample), sample_group_col]

    # Order samples by group
    sample_order <- rownames(pheno)[order(pheno[[sample_group_col]])]
    beta_melted$Sample <- factor(beta_melted$Sample, levels = sample_order)

    heatmap_plot <- ggplot2::ggplot(beta_melted) +
        ggplot2::geom_tile(ggplot2::aes(x = Position, y = Sample, fill = Beta)) +
        ggplot2::scale_fill_gradient(
            low = "white", high = "red",
            limits = c(min(beta_melted$Beta), max(beta_melted$Beta)),
            name = "\u03b2-values"
        ) +
        ggplot2::labs(
            y = "Sample",
            x = sprintf("Genomic Position on %s (bp)", chr),
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7),
        ) +
        ggplot2::scale_x_continuous(
            breaks = breaks,
            labels = breaks_labels
        ) +
        ggplot2::coord_cartesian(xlim = c(breaks[1], breaks[length(breaks)])) +
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
    g1 <- ggplot2::ggplotGrob(structure_plot)
    g2 <- ggplot2::ggplotGrob(heatmap_plot)
    max_width <- grid::unit.pmax(g1$widths, g2$widths)
    combined <- gridExtra::gtable_rbind(g1, g2)
    combined$widths <- max_width
    grid::grid.draw(combined)
    invisible(combined)
}
