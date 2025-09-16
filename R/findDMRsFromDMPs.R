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
                            expansion.step = 500,
                            max.pval = 0.05,
                            max.lookup.dist = 10000,
                            min.dmps = 1,
                            min.cpgs = 50,
                            ignored.sample.groups = NULL,
                            output.id = NULL,
                            njobs = 1,
                            verbose = FALSE,
                            extreme.verbosity = FALSE,
                            beta.row.names.file = NULL) {
    verbose <- verbose || extreme.verbosity
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
    
    # Create case-control vector from sample groups
    casecontrol <- NULL
    if (length(unique(sample.groups)) == 2 && 
        !any(levels(sample.groups) %in% ignored.sample.groups)) {
        if (verbose) message("Creating case-control vector from binary sample groups")
        casecontrol <- as.numeric(sample.groups == levels(sample.groups)[2])
    }
    
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
        
        # Read beta row names if not provided
        if (is.null(beta.row.names.file)) {
            beta.row.names <- unlist(data.table::fread(
                file = beta.file,
                select = 1,
                header = TRUE
            ))
        } else {
            beta.row.names <- readLines(beta.row.names.file)
        }
        
        beta.info <- .get.beta.col.names.and.inds(beta.file, rownames(pheno))
        beta.info$row.names <- beta.row.names
    } else {
        if (!file.exists(tabix.file)) {
            stop("Tabix file not found: ", tabix.file)
        }
        beta.info <- .get.beta.col.names.and.inds(tabix.file, rownames(pheno),
                                                 is.tabix = TRUE)
        
        # For tabix files, read row names from first column
        if (is.null(beta.row.names.file)) {
            beta.row.names <- unlist(data.table::fread(
                file = tabix.file,
                select = 1,
                header = TRUE
            ))
        } else {
            beta.row.names <- readLines(beta.row.names.file)
        }
        beta.info$row.names <- beta.row.names
    }
    if (is.null(beta.row.names.file) && !is.null(output.id)) {
        output.beta.row.names.file <- paste0(output.id, "_row_names.txt")
        beta.row.names.file <- output.beta.row.names.file # load from previous run
        if (verbose) message("Saving beta file row names to beta row names file: ", beta.row.names.file)
        writeLines(paste(beta.row.names, collapse='\n'), beta.row.names.file)        
    }
    # Process DMPs and create sorted locations
    # First check if all DMPs are in beta file
    unique_dmps <- unique(dmps$dmp)
    if (!all(unique_dmps %in% beta.info$row.names)) {
        stop("Some of the DMPs are not present in the beta file. DMPs: ", 
             paste(unique_dmps[!(unique_dmps %in% beta.info$row.names)], collapse=','))
    }
    
    # Create full CpG locations dataframe from methylation data and DMPs
    # First get genomic locations for all CpGs in beta matrix
    if (!is.null(beta.file)) {
        # Try loading EPIC annotations first
        cpg_locations <- try({
            if (requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
                IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
            } else {
                NULL
            }
        })
        
        # If EPIC fails or missing CpGs, try 450k annotations
        if (inherits(cpg_locations, "try-error") || is.null(cpg_locations) || 
            !all(beta.info$row.names %in% rownames(cpg_locations))) {
            cpg_locations_450k <- try({
                if (requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
                    IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
                } else {
                    NULL
                }
            })
            
            if (!inherits(cpg_locations_450k, "try-error") && !is.null(cpg_locations_450k)) {
                if (is.null(cpg_locations)) {
                    cpg_locations <- cpg_locations_450k
                } else {
                    # Merge EPIC and 450k annotations, keeping EPIC versions where duplicated
                    missing_cpgs <- setdiff(rownames(cpg_locations_450k), rownames(cpg_locations))
                    if (length(missing_cpgs) > 0) {
                        cpg_locations <- rbind(cpg_locations,
                                            cpg_locations_450k[missing_cpgs,])
                    }
                }
            }
        }
        
        if (is.null(cpg_locations)) {
            stop("Could not load either EPIC or 450k annotations. Please install ",
                 "IlluminaHumanMethylationEPICanno.ilm10b4.hg19 and/or ",
                 "IlluminaHumanMethylation450kanno.ilmn12.hg19")
        }
        
        # Check if we have all CpG locations
        missing_cpgs <- setdiff(beta.info$row.names, rownames(cpg_locations))
        if (length(missing_cpgs) > 0) {
            warning("Could not find genomic locations for ", length(missing_cpgs), 
                   " CpGs in annotation packages. These will be excluded unless ",
                   "they are DMPs with positions provided.")
            if (verbose) {
                message("First few missing CpGs: ", 
                       paste(head(missing_cpgs, 5), collapse=", "))
            }
        }
        
        # Create initial sorted locations from available annotations
        available_cpgs <- intersect(beta.info$row.names, rownames(cpg_locations))
        sorted.locs <- data.frame(
            chr = cpg_locations[available_cpgs, "chr"],
            pos = cpg_locations[available_cpgs, "pos"],
            cpg = available_cpgs,
            stringsAsFactors = FALSE
        )
    } else {
        # For tabix files, positions should be in file
        stop("Tabix file support for full CpG locations not yet implemented")
    }
    
    # Ensure all DMPs are in the sorted locations
    unique_dmps <- unique(dmps$dmp)
    missing_dmps <- unique_dmps[!unique_dmps %in% sorted.locs$cpg]
    if (length(missing_dmps) > 0) {
        # Add missing DMP locations from DMP file
        missing_locs <- as.data.frame(dmps[dmps$dmp %in% missing_dmps, c("chr", "pos", "dmp")])
        names(missing_locs)[names(missing_locs) == "dmp"] <- "cpg"
        sorted.locs <- rbind(sorted.locs, missing_locs)
    }
    
    # Sort by chromosome and position
    sorted.locs <- sorted.locs[order(sorted.locs$chr, sorted.locs$pos), ]
    rownames(sorted.locs) <- sorted.locs$cpg
    
    # Process DMPs and find initial regions
    initial.regions <- .findInitialRegions(
        dmps = dmps,
        max.lookup.dist = max.lookup.dist,
        min.dmps = min.dmps
    )
    
    if (verbose) {
        message("Found ", nrow(initial.regions), " initial regions")
    }
    browser()
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
                sorted.locs = sorted.locs,
                max.pval = max.pval,
                min.cpg.delta_beta = min.cpg.delta_beta,
                casecontrol = casecontrol,  # Add case-control vector
                expansion.step = expansion.step,
                verbose = verbose,
                extreme.verbosity = extreme.verbosity
            )
        },
        mc.cores = njobs
    )
    
    # Filter regions
    valid.regions <- do.call(rbind, expanded.regions)
    valid.regions <- valid.regions[valid.regions$n_cpgs >= min.cpgs, ]
    
    # Convert to GRanges with all metadata from original implementation
    gr <- GenomicRanges::GRanges(
        seqnames = valid.regions$chr,
        ranges = IRanges::IRanges(
            start = valid.regions$start,
            end = valid.regions$end
        ),
        mcols = S4Vectors::DataFrame(
            # CpG and DMP information
            start_cpg = valid.regions$start_cpg,  # Added CpG ID fields
            end_cpg = valid.regions$end_cpg,      # Added CpG ID fields
            start_dmp = valid.regions$start_dmp,  # Keep DMP fields for compatibility
            end_dmp = valid.regions$end_dmp,      # Keep DMP fields for compatibility
            start_dmp_pos = valid.regions$start,  # Use actual genomic positions
            end_dmp_pos = valid.regions$end,      # Use actual genomic positions
            dmps_num = valid.regions$dmps_num,
            
            # Delta beta statistics
            delta_beta = valid.regions$delta_beta,
            delta_beta_sd = valid.regions$delta_beta_sd,
            delta_beta_se = valid.regions$delta_beta_se,
            delta_beta_min = valid.regions$delta_beta_min,
            delta_beta_max = valid.regions$delta_beta_max,
            delta_beta_start = valid.regions$delta_beta_start,
            delta_beta_mid = valid.regions$delta_beta_mid,
            delta_beta_end = valid.regions$delta_beta_end,
            
            # Cases beta statistics
            cases_beta = valid.regions$cases_beta,
            cases_beta_max = valid.regions$cases_beta_max,
            cases_beta_min = valid.regions$cases_beta_min,
            cases_beta_sd = valid.regions$cases_beta_sd,
            cases_beta_se = valid.regions$cases_beta_se,
            cases_beta_start = valid.regions$cases_beta_start,
            cases_beta_mid = valid.regions$cases_beta_mid,
            cases_beta_end = valid.regions$cases_beta_end,
            
            # Cases beta SD statistics
            cases_beta_dmps_sd = valid.regions$cases_beta_dmps_sd,
            cases_beta_dmps_sd_max = valid.regions$cases_beta_dmps_sd_max,
            cases_beta_dmps_sd_min = valid.regions$cases_beta_dmps_sd_min,
            cases_beta_dmps_sd_start = valid.regions$cases_beta_dmps_sd_start,
            cases_beta_dmps_sd_mid = valid.regions$cases_beta_dmps_sd_mid,
            cases_beta_dmps_sd_end = valid.regions$cases_beta_dmps_sd_end,
            
            # Controls beta statistics
            controls_beta = valid.regions$controls_beta,
            controls_beta_max = valid.regions$controls_beta_max,
            controls_beta_min = valid.regions$controls_beta_min,
            controls_beta_sd = valid.regions$controls_beta_sd,
            controls_beta_se = valid.regions$controls_beta_se,
            controls_beta_start = valid.regions$controls_beta_start,
            controls_beta_mid = valid.regions$controls_beta_mid,
            controls_beta_end = valid.regions$controls_beta_end,
            
            # Controls beta SD statistics
            controls_beta_dmps_sd = valid.regions$controls_beta_dmps_sd,
            controls_beta_dmps_sd_max = valid.regions$controls_beta_dmps_sd_max,
            controls_beta_dmps_sd_min = valid.regions$controls_beta_dmps_sd_min,
            controls_beta_dmps_sd_start = valid.regions$controls_beta_dmps_sd_start,
            controls_beta_dmps_sd_mid = valid.regions$controls_beta_dmps_sd_mid,
            controls_beta_dmps_sd_end = valid.regions$controls_beta_dmps_sd_end,
            
            # Sample counts and correlation
            cases_num = valid.regions$cases_num,
            controls_num = valid.regions$controls_num,
            corr_pval = valid.regions$corr_pval,
            
            # DMP statistics
            dmps_pval_adj = valid.regions$dmps_pval_adj,
            dmps_pval_adj_min = valid.regions$dmps_pval_adj_min,
            dmps_pval_adj_max = valid.regions$dmps_pval_adj_max,
            dmps_pval = valid.regions$dmps_pval,
            dmps_pval_min = valid.regions$dmps_pval_min,
            dmps_pval_max = valid.regions$dmps_pval_max,
            dmps_qval = valid.regions$dmps_qval,
            dmps_qval_min = valid.regions$dmps_qval_min,
            dmps_qval_max = valid.regions$dmps_qval_max,
            
            # Region information
            stop_connection_reason = valid.regions$stop_connection_reason,
            dmps = valid.regions$dmps
        )
    )
    
    # Save results if output.id provided
    if (!is.null(output.id)) {
        saveRDS(gr, file = paste0(output.id, "_dmrs.rds"))
    }
    
    gr
}