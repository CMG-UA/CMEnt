#' Find Differentially Methylated Regions from DMPs
#'
#' Identifies Differentially Methylated Regions (DMRs) from pre-computed 
#' Differentially Methylated Positions (DMPs) using correlation analysis and 
#' region expansion.
#'
#' @param beta.file Path to beta value file, or NULL if using tabix
#' @param tabix.file Path to tabix file, or NULL if using beta file
#' @param dmps.tsv.file Path to DMPs TSV file
#' @param pheno Data frame with sample phenotype information
#' @param pval.col Column name for p-values in DMPs file
#' @param sample_group.col Column name for sample groups in phenotype
#' @param min.cpg.delta_beta Minimum delta beta for CpG connectivity
#' @param expansion.step Number of CpGs to check in each expansion step
#' @param max.pval Maximum p-value for correlation
#' @param max.lookup.dist Maximum distance for DMR lookup
#' @param min.dmps Minimum number of DMPs per region
#' @param min.cpgs Minimum number of CpGs per region
#' @param ignored.sample.groups Sample groups to ignore
#' @param output.id Identifier for output files
#' @param njobs Number of parallel jobs
#' @param verbose Whether to print progress messages
#' @param beta.row.names.file Optional file with beta matrix row names
#'
#' @return A \code{GRanges} object containing the identified DMRs with metadata
#'   including number of CpGs, number of DMPs, mean p-value, and expansion
#'   stop reasons.
#'
#' @examples
#' \dontrun{
#' # Example will be added with example data
#' }
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
#' @importFrom stats p.adjust
#' @importFrom data.table fread
#' @importFrom parallel mclapply
#'
#' @export
findDMRsFromDMPs <- function(beta.file = NULL,
                            tabix.file = NULL,
                            dmps.tsv.file,
                            pheno,
                            pval.col = "pval_adj",
                            sample_group.col,
                            min.cpg.delta_beta = 0,
                            expansion.step = 5,
                            max.pval = 0.05,
                            max.lookup.dist = 10000,
                            min.dmps = 1,
                            min.cpgs = 50,
                            ignored.sample.groups = NULL,
                            output.id = NULL,
                            njobs = 1,
                            verbose = FALSE,
                            beta.row.names.file = NULL) {
    
    # Input validation
    if (is.null(beta.file) && is.null(tabix.file)) {
        stop("Either beta.file or tabix.file must be provided")
    }
    if (!is.null(beta.file) && !is.null(tabix.file)) {
        stop("Only one of beta.file or tabix.file should be provided")
    }
    
    # Process phenotype data
    if (!sample_group.col %in% colnames(pheno)) {
        stop("sample_group.col not found in phenotype data")
    }
    sample.groups <- as.factor(pheno[[sample_group.col]])
    
    # Read DMPs
    dmps <- data.table::fread(dmps.tsv.file)
    if (!pval.col %in% colnames(dmps)) {
        stop("pval.col not found in DMPs file")
    }
    
    # Initialize beta data handling
    if (!is.null(beta.file)) {
        if (!file.exists(beta.file)) {
            stop("Beta file not found: ", beta.file)
        }
        beta.info <- .get.beta.col.names.and.inds(beta.file, colnames(pheno))
    } else {
        if (!file.exists(tabix.file)) {
            stop("Tabix file not found: ", tabix.file)
        }
        beta.info <- .get.beta.col.names.and.inds(tabix.file, colnames(pheno),
                                                 is.tabix = TRUE)
    }
    
    # Process DMPs and find initial regions
    dmps <- dmps[order(dmps$chr, dmps$pos), ]
    initial.regions <- .findInitialRegions(
        dmps = dmps,
        max.lookup.dist = max.lookup.dist,
        min.dmps = min.dmps
    )
    
    if (verbose) {
        message("Found ", nrow(initial.regions), " initial regions")
    }
    
    # Expand regions in parallel
    expanded.regions <- parallel::mclapply(
        seq_len(nrow(initial.regions)),
        function(i) {
            region <- initial.regions[i, ]
            .expandDMRs(
                dmr = region,
                beta.file = beta.file,
                beta.row.names = beta.info$row.names,
                beta.col.names = beta.info$col.names,
                sample.groups = sample.groups,
                sorted.locs = dmps,
                max.pval = max.pval,
                min.cpg.delta_beta = min.cpg.delta_beta,
                expansion.step = expansion.step
            )
        },
        mc.cores = njobs
    )
    
    # Filter regions
    valid.regions <- do.call(rbind, expanded.regions)
    valid.regions <- valid.regions[valid.regions$n_cpgs >= min.cpgs, ]
    
    # Convert to GRanges with all metadata from original implementation
    gr <- GenomicRanges::GRanges(
        seqnames = dmrs$chr,
        ranges = IRanges::IRanges(
            start = dmrs$start,
            end = dmrs$end
        ),
        mcols = S4Vectors::DataFrame(
            # DMP information
            start_dmp = dmrs$start_dmp,
            end_dmp = dmrs$end_dmp,
            start_dmp_pos = dmrs$start_dmp_pos,
            end_dmp_pos = dmrs$end_dmp_pos,
            dmps_num = dmrs$dmps_num,
            
            # Delta beta statistics
            delta_beta = dmrs$delta_beta,
            delta_beta_sd = dmrs$delta_beta_sd,
            delta_beta_se = dmrs$delta_beta_se,
            delta_beta_min = dmrs$delta_beta_min,
            delta_beta_max = dmrs$delta_beta_max,
            delta_beta_start = dmrs$delta_beta_start,
            delta_beta_mid = dmrs$delta_beta_mid,
            delta_beta_end = dmrs$delta_beta_end,
            
            # Cases beta statistics
            cases_beta = dmrs$cases_beta,
            cases_beta_max = dmrs$cases_beta_max,
            cases_beta_min = dmrs$cases_beta_min,
            cases_beta_sd = dmrs$cases_beta_sd,
            cases_beta_se = dmrs$cases_beta_se,
            cases_beta_start = dmrs$cases_beta_start,
            cases_beta_mid = dmrs$cases_beta_mid,
            cases_beta_end = dmrs$cases_beta_end,
            
            # Cases beta SD statistics
            cases_beta_dmps_sd = dmrs$cases_beta_dmps_sd,
            cases_beta_dmps_sd_max = dmrs$cases_beta_dmps_sd_max,
            cases_beta_dmps_sd_min = dmrs$cases_beta_dmps_sd_min,
            cases_beta_dmps_sd_start = dmrs$cases_beta_dmps_sd_start,
            cases_beta_dmps_sd_mid = dmrs$cases_beta_dmps_sd_mid,
            cases_beta_dmps_sd_end = dmrs$cases_beta_dmps_sd_end,
            
            # Controls beta statistics
            controls_beta = dmrs$controls_beta,
            controls_beta_max = dmrs$controls_beta_max,
            controls_beta_min = dmrs$controls_beta_min,
            controls_beta_sd = dmrs$controls_beta_sd,
            controls_beta_se = dmrs$controls_beta_se,
            controls_beta_start = dmrs$controls_beta_start,
            controls_beta_mid = dmrs$controls_beta_mid,
            controls_beta_end = dmrs$controls_beta_end,
            
            # Controls beta SD statistics
            controls_beta_dmps_sd = dmrs$controls_beta_dmps_sd,
            controls_beta_dmps_sd_max = dmrs$controls_beta_dmps_sd_max,
            controls_beta_dmps_sd_min = dmrs$controls_beta_dmps_sd_min,
            controls_beta_dmps_sd_start = dmrs$controls_beta_dmps_sd_start,
            controls_beta_dmps_sd_mid = dmrs$controls_beta_dmps_sd_mid,
            controls_beta_dmps_sd_end = dmrs$controls_beta_dmps_sd_end,
            
            # Sample counts and correlation
            cases_num = dmrs$cases_num,
            controls_num = dmrs$controls_num,
            corr_pval = dmrs$corr_pval,
            
            # DMP statistics
            dmps_pval_adj = dmrs$dmps_pval_adj,
            dmps_pval_adj_min = dmrs$dmps_pval_adj_min,
            dmps_pval_adj_max = dmrs$dmps_pval_adj_max,
            dmps_pval = dmrs$dmps_pval,
            dmps_pval_min = dmrs$dmps_pval_min,
            dmps_pval_max = dmrs$dmps_pval_max,
            dmps_qval = dmrs$dmps_qval,
            dmps_qval_min = dmrs$dmps_qval_min,
            dmps_qval_max = dmrs$dmps_qval_max,
            
            # Region information
            stop_connection_reason = dmrs$stop_connection_reason,
            dmps = dmrs$dmps
        )
    )
    
    # Save results if output.id provided
    if (!is.null(output.id)) {
        saveRDS(gr, file = paste0(output.id, "_dmrs.rds"))
    }
    
    gr
}