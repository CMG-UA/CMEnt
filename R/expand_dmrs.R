#' DMR Expansion Functions
#' @keywords internal

#' Expand DMRs to nearby CpG regions
#'
#' @param dmr DMR data frame with initial boundaries
#' @param beta.file Path to beta value file
#' @param beta.row.names Vector of beta matrix row names
#' @param beta.col.names Vector of beta matrix column names
#' @param sample.groups Factor of sample groups
#' @param sorted.locs Data frame of sorted CpG locations
#' @param max.pval Maximum p-value for correlation
#' @param min.cpg.delta_beta Minimum delta beta threshold
#' @param casecontrol Optional case-control vector
#' @param tabix.file Optional path to tabix file
#' @param expansion.step Number of CpGs to check in each step
#' @return Expanded DMR data frame
#' @importFrom data.table fread
#' @keywords internal
.expandDMRs <- function(dmr,
                       beta.file,
                       beta.row.names,
                       beta.col.names,
                       sample.groups,
                       sorted.locs,
                       max.pval,
                       min.cpg.delta_beta = 0,
                       casecontrol = NULL,
                       tabix.file = NULL,
                       expansion.step = 500) {
    
    # Get beta value column information
    if (!is.null(beta.file)) {
        ret <- .get.beta.col.names.and.inds(beta.file, beta.col.names)
    } else {
        ret <- .get.beta.col.names.and.inds(tabix.file, beta.col.names,
                                           is.tabix = TRUE)
    }
    
    cols.inds <- ret$beta.col.inds
    file.beta.col.names <- ret$file.beta.col.names
    beta.col.names <- ret$beta.col.names
    sorted.locs <- sorted.locs[beta.row.names,]
    
    # Initialize DMR boundaries
    dmr.start <- dmr$start_dmp
    dmr.end <- dmr$end_dmp
    dmr.start.ind <- which(beta.row.names == dmr.start)
    
    if (length(dmr.start.ind) == 0) {
        stop("Could not find the start CpG ", dmr.start,
             " in the beta file row names.")
    }
    
    # Initialize statistics tracking
    all_betas <- list()  # Store beta values during expansion
    all_correlations <- numeric()  # Store correlation p-values
    
    # Expand upstream
    end.site.ind <- dmr.start.ind[[1]]
    upstream.exp <- end.site.ind
    upstream.stop.found <- FALSE
    
    while (TRUE) {
        if (end.site.ind < 0) {
            upstream.stop.reason <- "end-of-input"
            upstream.exp <- 1
            break
        }
        
        start.site.ind <- max(0, end.site.ind - expansion.step) + 1
        x <- which(sorted.locs[start.site.ind:end.site.ind, "chr"] == dmr$chr)
        
        if (length(x) == 0) {
            upstream.stop.reason <- "end-of-input"
            break
        }
        
        x <- x[[1]]
        start.site.ind <- start.site.ind + x - 1
        exp.step <- end.site.ind - start.site.ind + 1
        
        # Read beta values
        if (!is.null(beta.file)) {
            upstream.betas <- data.table::fread(
                file = beta.file,
                skip = start.site.ind,
                nrows = exp.step,
                header = FALSE,
                data.table = FALSE,
                colClasses = c("character",
                             rep("numeric", length(file.beta.col.names)))
            )
            upstream.betas <- upstream.betas[, cols.inds]
        } else {
            upstream.region <- paste0(
                sorted.locs[start.site.ind, 'chr'], ':',
                sorted.locs[start.site.ind, 'pos'], '-',
                sorted.locs[end.site.ind, 'pos'] + 1
            )
            upstream.betas <- try(Rsamtools::tabix(
                upstream.region,
                tabix.file,
                check.valid = FALSE,
                verbose = FALSE
            ))
            if (inherits(upstream.betas, "try-error")) {
                warning("Error reading upstream region ", upstream.region)
                upstream.stop.reason <- "error-reading-tabix"
                break
            }
            upstream.betas <- as.data.frame(
                sapply(upstream.betas[, beta.col.names], as.numeric)
            )
        }
        
        upstream.betas <- upstream.betas[rev(seq_len(nrow(upstream.betas))), ,
                                       drop = FALSE]
        
        # Test connectivity and collect statistics
        i <- 1
        while (TRUE) {
            corr.ret <- .test.connectivity(
                site1.beta = unlist(upstream.betas[i, ]),
                site2.beta = unlist(upstream.betas[i + 1, ]),
                sample.groups = sample.groups,
                max.pval = max.pval,
                casecontrol = casecontrol,
                min.delta_beta = min.cpg.delta_beta
            )
            
            # Store statistics
            all_correlations <- c(all_correlations, corr.ret[[2]])  # p-value
            all_betas[[length(all_betas) + 1]] <- upstream.betas[i, ]
            
            if (!corr.ret[[1]]) {
                upstream.exp <- end.site.ind - i + 1
                upstream.stop.found <- TRUE
                upstream.stop.reason <- corr.ret$reason
                break
            }
            
            i <- i + 1
            if (i == nrow(upstream.betas)) {
                # Store last beta values
                all_betas[[length(all_betas) + 1]] <- upstream.betas[i, ]
                break
            }
        }
        
        if (upstream.stop.found)
            break
            
        end.site.ind <- end.site.ind - expansion.step
    }
    
    # Expand downstream
    dmr.end.ind <- which(beta.row.names == dmr.end)
    if (length(dmr.end.ind) == 0) {
        stop("Could not find the end CpG ", dmr.end,
             " in the beta file row names.")
    }
    
    start.site.ind <- dmr.end.ind[[1]]
    downstream.exp <- start.site.ind
    downstream.stop.found <- FALSE
    
    while (TRUE) {
        if (start.site.ind > length(beta.row.names)) {
            downstream.stop.reason <- "end-of-input"
            downstream.exp <- start.site.ind - 1
            break
        }
        
        end.site.ind <- min(start.site.ind + expansion.step,
                           length(beta.row.names))
        x <- which(rev(sorted.locs[start.site.ind:end.site.ind, "chr"]) == dmr$chr)
        
        if (length(x) == 0) {
            downstream.stop.reason <- "end-of-input"
            downstream.exp <- start.site.ind - 1
            break
        }
        
        x <- x[[1]]
        end.site.ind <- end.site.ind - x + 1
        
        if (end.site.ind <= start.site.ind + 1) {
            downstream.stop.reason <- "end-of-input"
            downstream.exp <- start.site.ind - 1
            break
        }
        
        # Read beta values
        if (!is.null(beta.file)) {
            downstream.betas <- data.table::fread(
                file = beta.file,
                skip = start.site.ind,
                nrows = expansion.step - x + 1,
                header = FALSE,
                data.table = FALSE,
                colClasses = c("character",
                             rep("numeric", length(file.beta.col.names)))
            )
            downstream.betas <- downstream.betas[, cols.inds]
        } else {
            downstream.region <- paste0(
                sorted.locs[start.site.ind, "chr"], ":",
                sorted.locs[start.site.ind, "pos"], "-",
                sorted.locs[end.site.ind, "pos"]
            )
            downstream.betas <- try(Rsamtools::tabix(
                downstream.region,
                tabix.file,
                check.valid = FALSE,
                verbose = FALSE
            ))
            if (inherits(downstream.betas, "try-error")) {
                warning("Error reading downstream region ", downstream.region,
                       " from tabix file, stopping extension.")
                downstream.stop.reason <- "error-reading-tabix"
                break
            }
            downstream.betas <- as.data.frame(
                sapply(downstream.betas[, beta.col.names], as.numeric)
            )
        }
        
        # Test connectivity and collect statistics
        i <- 1
        while (TRUE) {
            corr.ret <- .test.connectivity(
                site1.beta = unlist(downstream.betas[i, ]),
                site2.beta = unlist(downstream.betas[i + 1, ]),
                sample.groups = sample.groups,
                max.pval = max.pval,
                casecontrol = casecontrol,
                min.delta_beta = min.cpg.delta_beta
            )
            
            # Store statistics
            all_correlations <- c(all_correlations, corr.ret[[2]])  # p-value
            all_betas[[length(all_betas) + 1]] <- downstream.betas[i, ]
            
            if (!corr.ret[[1]]) {
                downstream.exp <- start.site.ind + i - 1
                downstream.stop.reason <- corr.ret$reason
                downstream.stop.found <- TRUE
                break
            }
            
            i <- i + 1
            if (i == nrow(downstream.betas)) {
                # Store last beta values
                all_betas[[length(all_betas) + 1]] <- downstream.betas[i, ]
                break
            }
        }
        
        if (downstream.stop.found)
            break
            
        prev.start.site.ind <- start.site.ind
        start.site.ind <- start.site.ind + expansion.step - x + 1
        
        if (start.site.ind < prev.start.site.ind) {
            stop("BUG: start.site.ind < prev.start.site.ind during downstream ",
                 "expansion. start.site.ind: ", start.site.ind,
                 " prev.start.site.ind: ", prev.start.site.ind,
                 " expansion.step: ", expansion.step,
                 " x: ", x)
        }
    }
    
    # Calculate final statistics and update DMR fields to match original implementation
    region_cpgs <- which(beta.row.names >= beta.row.names[upstream.exp] & 
                        beta.row.names <= beta.row.names[downstream.exp])
    
    # Calculate beta values for DMPs in the region
    beta_matrix <- do.call(rbind, all_betas)
    region_dmps <- which(beta.row.names >= dmr$start_dmp & 
                        beta.row.names <= dmr$end_dmp)
    
    # Get the actual DMPs beta values
    dmp_indices <- match(beta.row.names[region_dmps], rownames(beta_matrix))
    dmp_indices <- dmp_indices[!is.na(dmp_indices)]
    dmp_betas <- beta_matrix[dmp_indices, , drop=FALSE]
    
    if (!is.null(casecontrol) && length(dmp_indices) > 0) {
        # Calculate case and control betas
        cases_betas <- rowMeans(dmp_betas[, casecontrol, drop=FALSE])
        controls_betas <- rowMeans(dmp_betas[, !casecontrol, drop=FALSE])
        delta_betas <- cases_betas - controls_betas
        
        # Calculate case and control SDs
        cases_sds <- apply(dmp_betas[, casecontrol, drop=FALSE], 1, sd)
        controls_sds <- apply(dmp_betas[, !casecontrol, drop=FALSE], 1, sd)
        
        # Update all DMR fields
        dmr$start_cpg <- beta.row.names[upstream.exp]
        dmr$end_cpg <- beta.row.names[downstream.exp]
        dmr$start <- sorted.locs[dmr$start_cpg, "pos"]
        dmr$end <- sorted.locs[dmr$end_cpg, "pos"]
        dmr$n_cpgs <- length(region_cpgs)
        dmr$dmps_num <- length(dmp_indices)
        dmr$stop_connection_reason <- downstream.stop.reason  # We use the last stop reason
        
        # Delta beta statistics
        dmr$delta_beta <- mean(abs(delta_betas)) * sign(sum(sign(delta_betas)))
        dmr$delta_beta_sd <- sd(delta_betas)
        dmr$delta_beta_se <- sd(delta_betas)/sqrt(length(delta_betas))
        dmr$delta_beta_min <- min(delta_betas)
        dmr$delta_beta_max <- max(delta_betas)
        dmr$delta_beta_start <- delta_betas[1]
        dmr$delta_beta_mid <- delta_betas[ceiling(length(delta_betas)/2)]
        dmr$delta_beta_end <- delta_betas[length(delta_betas)]
        
        # Cases beta statistics
        dmr$cases_beta <- mean(abs(cases_betas)) * sign(sum(sign(cases_betas)))
        dmr$cases_beta_max <- max(cases_betas)
        dmr$cases_beta_min <- min(cases_betas)
        dmr$cases_beta_sd <- sd(cases_betas)
        dmr$cases_beta_se <- sd(cases_betas)/sqrt(length(cases_betas))
        dmr$cases_beta_start <- cases_betas[1]
        dmr$cases_beta_mid <- cases_betas[ceiling(length(cases_betas)/2)]
        dmr$cases_beta_end <- cases_betas[length(cases_betas)]
        
        # Cases beta SD statistics
        dmr$cases_beta_dmps_sd <- mean(cases_sds)
        dmr$cases_beta_dmps_sd_max <- max(cases_sds)
        dmr$cases_beta_dmps_sd_min <- min(cases_sds)
        dmr$cases_beta_dmps_sd_start <- cases_sds[1]
        dmr$cases_beta_dmps_sd_mid <- cases_sds[ceiling(length(cases_sds)/2)]
        dmr$cases_beta_dmps_sd_end <- cases_sds[length(cases_sds)]
        
        # Controls beta statistics
        dmr$controls_beta <- mean(abs(controls_betas)) * sign(sum(sign(controls_betas)))
        dmr$controls_beta_max <- max(controls_betas)
        dmr$controls_beta_min <- min(controls_betas)
        dmr$controls_beta_sd <- sd(controls_betas)
        dmr$controls_beta_se <- sd(controls_betas)/sqrt(length(controls_betas))
        dmr$controls_beta_start <- controls_betas[1]
        dmr$controls_beta_mid <- controls_betas[ceiling(length(controls_betas)/2)]
        dmr$controls_beta_end <- controls_betas[length(controls_betas)]
        
        # Controls beta SD statistics
        dmr$controls_beta_dmps_sd <- mean(controls_sds)
        dmr$controls_beta_dmps_sd_max <- max(controls_sds)
        dmr$controls_beta_dmps_sd_min <- min(controls_sds)
        dmr$controls_beta_dmps_sd_start <- controls_sds[1]
        dmr$controls_beta_dmps_sd_mid <- controls_sds[ceiling(length(controls_sds)/2)]
        dmr$controls_beta_dmps_sd_end <- controls_sds[length(controls_sds)]
        
        # Sample counts
        dmr$cases_num <- sum(casecontrol)
        dmr$controls_num <- sum(!casecontrol)
        
        # Correlation statistics
        dmr$corr_pval <- max(all_correlations, na.rm=TRUE)  # Same as original max_pval_in_region
    } else {
        # Handle the case where we don't have case/control info or no DMPs
        dmr$start_cpg <- beta.row.names[upstream.exp]
        dmr$end_cpg <- beta.row.names[downstream.exp]
        dmr$start <- sorted.locs[dmr$start_cpg, "pos"]
        dmr$end <- sorted.locs[dmr$end_cpg, "pos"]
        dmr$n_cpgs <- length(region_cpgs)
        dmr$dmps_num <- length(dmp_indices)
        dmr$stop_connection_reason <- downstream.stop.reason
    }
    
    dmr
}