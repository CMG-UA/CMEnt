#' Test connectivity between CpG sites
#'
#' Tests the correlation between two CpG sites within each sample group and
#' optionally checks delta beta thresholds.
#'
#' @param site1.beta Beta values for first site
#' @param site2.beta Beta values for second site
#' @param sample.groups Factor of sample groups
#' @param max.pval Maximum p-value threshold
#' @param casecontrol Optional case-control vector
#' @param min.delta_beta Minimum delta beta threshold
#' @param extreme.verbosity Whether to print detailed messages
#' @return List with connectivity test results
#' @keywords internal
.test.connectivity <- function(site1.beta, site2.beta, sample.groups, max.pval,
                             casecontrol = NULL, min.delta_beta = 0,
                             extreme.verbosity = FALSE) {
    pval <- 0
    delta_beta <- NULL

    # Apply Bonferroni correction
    max.pval.corrected <- max.pval / length(unique(sample.groups))
    
    # Test correlation within each sample group
    for (g in levels(sample.groups)) {
        # Skip groups with too few samples
        if (sum(sample.groups == g) < 3)
            next
        
        # Try correlation test with error handling
        op <- options(warn = 2)$warn
        group_mask <- sample.groups == g
        corr.ret <- try(stats::cor.test(
            site1.beta[group_mask],
            site2.beta[group_mask],
            method = "pearson"
        ))
        options(warn = op)
        
        # Handle correlation test errors
        if (inherits(corr.ret, "try-error")) {
            if (extreme.verbosity) {
                message(".test.connectivity: Error occurred in cor.test while processing the following:")
                message("casecontrol:", paste(casecontrol, collapse=","))
                message("site2.beta:", paste(site2.beta, collapse=","))
                message("sample.groups:", paste(sample.groups, collapse=","))
                message("max.pval.corrected:", max.pval.corrected)
                message("Error message:", corr.ret)
                browser()
            }
            return(list(FALSE, pval, delta_beta, failing = g, reason = "error occurred"))
        }
        
        # Check p-value
        r <- max(pval, corr.ret$p.value)
        if (is.null(r) || is.na(r)) {
            return(list(FALSE, pval, delta_beta, failing = g, reason = "na pval"))
        }
        pval <- r
        if (pval > max.pval.corrected) {
            return(list(FALSE, pval, delta_beta, failing = g, 
                       reason = "pval>max.pval (corrected)"))
        }
    }
    
    # Check delta beta if case-control info provided
    if (!is.null(casecontrol) && (min.delta_beta > 0)) {
        if (length(casecontrol) != length(site2.beta)) {
            if (extreme.verbosity) {
                message(".test.connectivity: Error occurred while computing delta beta for the following:")
                message("casecontrol:", paste(casecontrol, collapse=","))
                message("site2.beta:", paste(site2.beta, collapse=","))
                message("sample.groups:", paste(sample.groups, collapse=","))
                message("max.pval.corrected:", max.pval.corrected)
            }
            stop(paste0(
                "The provided casecontrol vector has length ", length(casecontrol),
                " while site2.beta has length ", length(site2.beta)
            ))
        }
        
        # Calculate delta beta with NA handling for both sites
        delta_beta1 <- mean(site1.beta[casecontrol == 1], na.rm = TRUE) -
                      mean(site1.beta[casecontrol == 0], na.rm = TRUE)
        delta_beta2 <- mean(site2.beta[casecontrol == 1], na.rm = TRUE) -
                      mean(site2.beta[casecontrol == 0], na.rm = TRUE)
        
        # Check both delta betas have same direction and meet minimum threshold
        if (is.null(delta_beta1) || is.na(delta_beta1) || 
            is.null(delta_beta2) || is.na(delta_beta2) ||
            (abs(delta_beta1) < min.delta_beta) ||
            (abs(delta_beta2) < min.delta_beta) ||
            (sign(delta_beta1) != sign(delta_beta2))) {
            delta_beta <- min(abs(delta_beta1), abs(delta_beta2))
            return(list(FALSE, pval, delta_beta, reason = "delta_beta<min.delta_beta or inconsistent direction"))
        }
        
        # Use the smaller delta beta as the connectivity measure
        delta_beta <- min(abs(delta_beta1), abs(delta_beta2)) * sign(delta_beta1)
    }
    
    list(TRUE, pval, delta_beta)
}