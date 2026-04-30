#' @title Helper functions for CMEnt
#' @description A collection of helper functions for the CMEnt package.


#' Find DMPs using DSS on BSseq objects
#'
#' This helper function identifies differentially methylated positions (DMPs) from a BSseq object using the DSS package. It allows for flexible specification of sample groups, covariates, and chromosome filtering.
#'
#' @param bsseq A BSseq object or a file path to a saved BSseq object (RDS format).
#' @param samplesheet A data frame or a file path to a tab-delimited text file containing sample metadata. Must include columns for sample IDs and group labels.
#' @param samplesheet_sep The separator used in the samplesheet file if a file path is provided. Default is tab ("\\t").
#' @param group_col The name of the column in the samplesheet that contains the group labels for comparison. Default is "Sample_Group".
#' @param id_col The name of the column in the samplesheet that contains the sample IDs. Default is "Sample_ID".
#' @param chr A character vector of chromosome names to include in the analysis, or "auto" to automatically include chr1-chr22, or "all" to include chr1-chr22 plus chrX and chrY. Default is "auto".
#' @param case_group The specific group label in the group_col to treat as the "case" group for comparison. If NULL, the first unique group in group_col will be used as the case group. Default is NULL.
#' @param covariates A character vector of additional covariate column names from the samplesheet to include in the DSS model, or a comma-separated string of covariate names. Default is NULL (no additional covariates).
#' @param fdr_thres The false discovery rate threshold for calling DMPs. Default is 0.05.
#' @param output_file An optional file path to save the DMP results as a tab-delimited text file. If the file name ends with ".gz", the output will be gzipped. Default is NULL (no file output).
#' @param njobs The number of parallel jobs to use for chromosome-level analysis. Default is the number of available CPU cores minus one.
#'
#' @return A data frame of identified DMPs with columns for chromosome, position, site ID, p-value, q-value, delta beta, and DMP score.
#' 
#' @examples
#' \dontrun{
#' # Load example BSseq data
#' data("BSobj", package = "bsseq")
#' # Create a sample metadata data frame
#' samplesheet <- data.frame(
#'    Sample_ID = colnames(BSobj),
#'   Sample_Group = c(rep("Condition1", 3), rep("Condition2", 3)),
#'   Age = c(30, 32, 31, 28, 29, 27)
#' )
#' # Find DMPs with DSS
#' dmps <- findDMPsBSSeq(
#'    bsseq = BSobj,
#'    samplesheet = samplesheet,
#'    group_col = "Sample_Group",
#'    id_col = "Sample_ID",
#'    case_group = "Condition2",
#'    covariates = "Age",
#'    fdr_thres = 0.05,
#'    output_file = "dmp_results.tsv.gz",
#'    njobs = 4
#' )
#' }
#' @importFrom stats as.formula
#' @export
findDMPsBSSeq <- function(
    bsseq,
    samplesheet,
    samplesheet_sep = "\t",
    group_col = "Sample_Group",
    id_col = "Sample_ID",
    chr = "auto",
    case_group = NULL,
    covariates = NULL,
    fdr_thres = 0.05,
    output_file = NULL,
    njobs = max(1L, BiocParallel::bpnworkers(BiocParallel::bpparam()) - 1L)
) {
    if (!requireNamespace("bsseq", quietly = TRUE)) {
        BiocManager::install("bsseq", update = FALSE, ask = FALSE)
    }
    if (!requireNamespace("DSS", quietly = TRUE)) {
        BiocManager::install("DSS", update = FALSE, ask = FALSE)
    }
    require(DSS)
    if (chr == "auto") {
        chr <- paste0("chr", c(1:22))
    } else if (chr == "all") {
        chr <- c(paste0("chr", c(1:22)), "chrX", "chrY")
    } else {
        chr <- trimws(unlist(strsplit(chr, ",")))
    }
    if (is.character(bsseq)) {
        bsseq_obj <- readRDS(bsseq)
    } else if (inherits(bsseq, "BSseq")) {
        bsseq_obj <- bsseq
    } else {
        stop("bsseq argument must be a file path or a BSseq object.")
    }
    seqn <- as.character(seqnames(bsseq_obj))
    bsseq_obj <- bsseq_obj[seqn %in% chr, ]
    seqn <- as.character(seqnames(bsseq_obj))
    if (length(seqn) == 0L) {
        stop("No loci remain after chromosome filtering.")
    }
    if (is.character(samplesheet)) {
        pheno <- read.table(samplesheet, header = TRUE, stringsAsFactors = FALSE, sep = samplesheet_sep, check.names = FALSE)
    } else if (is.data.frame(samplesheet)) {
        pheno <- samplesheet
    } else {
        stop("samplesheet argument must be a file path or a data frame.")
    }
    colnames(pheno) <- trimws(colnames(pheno))

    if (!all(c(group_col, id_col) %in% colnames(pheno))) {
        stop(paste("Group column", group_col, "or ID column", id_col, "not found in samplesheet."))
    }

    pheno[[id_col]] <- as.character(pheno[[id_col]])
    if (anyDuplicated(colnames(bsseq_obj)) > 0) {
        stop("Duplicate sample IDs found in BSseq object column names.")
    }
    if (anyDuplicated(pheno[[id_col]]) > 0) {
        stop("Duplicate sample IDs found in the samplesheet ID column.")
    }

    if (!is.null(case_group) && !case_group %in% pheno[[group_col]]) {
        stop(paste("Specified case group", case_group, "not found in group column", group_col))
    }
    if (!all(colnames(bsseq_obj) %in% pheno[[id_col]])) {
        stop("Not all samples in BSseq object are present in the samplesheet.")
    }
    if (is.null(case_group)) {
        case_group <- unique(pheno[[group_col]])[[1]]
    }

    njobs <- as.integer(njobs)
    if (is.na(njobs) || njobs < 1L) {
        stop("njobs must be a positive integer.")
    }

    # Ensure model rows are in the same order as BSseq sample columns.
    pheno <- pheno[match(colnames(bsseq_obj), pheno[[id_col]]), , drop = FALSE]
    rownames(pheno) <- pheno[[id_col]]

    if (!is.null(covariates)) {
        if (length(covariates) == 1 && is.character(covariates)) {
            covariates <- unlist(strsplit(covariates, ","))
        }
        covariates <- trimws(as.character(covariates))
        covariates <- covariates[nzchar(covariates)]
        if (length(covariates) == 0) {
            covariates <- NULL
        }
    }
    if (!is.null(covariates)) {
        missing_covariates <- setdiff(covariates, colnames(pheno))
        if (length(missing_covariates) > 0) {
            stop(paste("The following covariates are missing from the provided samplesheet:", paste(missing_covariates, collapse = ",")))
        }
    }

    pheno$condition <- factor(
        pheno[[group_col]] == case_group,
        levels = c(FALSE, TRUE),
        labels = c("control", "case")
    )
    if (length(unique(pheno$condition)) < 2) {
        stop("Condition has fewer than two levels after applying case_group; cannot fit contrast.")
    }
    covariates_formula <- covariates
    if (!is.null(covariates_formula)) {
        covariates_formula <- ifelse(
            make.names(covariates_formula) == covariates_formula,
            covariates_formula,
            paste0("`", covariates_formula, "`")
        )
    }
    formula <- as.formula(paste("~ condition", if (!is.null(covariates_formula)) paste("+", paste(covariates_formula, collapse = " + "))))

    chr_in_bsseq <- unique(seqn)
    nworkers <- min(njobs, length(chr_in_bsseq))

    run_dss_for_chr <- function(chr_name) {
        chr_idx <- seqn == chr_name
        bsseq_chr <- bsseq_obj[chr_idx, ]
        if (nrow(bsseq_chr) == 0L) {
            return(NULL)
        }
        suppressMessages({
            fit <- DSS::DMLfit.multiFactor(
                BSobj = bsseq_chr,
                design = pheno,
                formula = formula
            )
            dml_chr <- DSS::DMLtest.multiFactor(fit, term = "condition")
        })
        dml_chr <- as.data.frame(dml_chr)

        case_idx <- pheno$condition == "case"
        control_idx <- pheno$condition == "control"
        meth_chr <- bsseq::getMeth(bsseq_chr, type = "raw")
        case_mean <- rowMeans(meth_chr[, case_idx, drop = FALSE], na.rm = TRUE)
        control_mean <- rowMeans(meth_chr[, control_idx, drop = FALSE], na.rm = TRUE)
        dml_chr$delta_beta <- case_mean - control_mean
        dml_chr$delta_beta[!is.finite(dml_chr$delta_beta)] <- NA_real_

        dml_chr
    }

    bp_param <- .makeBiocParallelParam(nworkers, n_tasks = length(chr_in_bsseq))
    dml_by_chr <- BiocParallel::bplapply(
        chr_in_bsseq,
        run_dss_for_chr,
        BPPARAM = bp_param
    )
    dml_by_chr <- Filter(Negate(is.null), dml_by_chr)
    if (length(dml_by_chr) == 0L) {
        stop("DSS::DMLtest returned no chromosome-level results.")
    }
    dml_test <- do.call(rbind, dml_by_chr)

    dml_test <- dml_test[
        is.finite(dml_test$pvals) &
            is.finite(dml_test$fdrs) &
            is.finite(dml_test$delta_beta),
        ,
        drop = FALSE
    ]
    if (nrow(dml_test) == 0) {
        stop("DSS::DMLtest returned no usable sites for DMP generation.")
    }

    dml_test$site_id <- paste0(dml_test$chr, ":", dml_test$pos)
    dml_test$score <- -log10(dml_test$pvals + .Machine$double.xmin) * abs(dml_test$delta_beta)
    dml_test$score <- (dml_test$score - min(dml_test$score)) / (max(dml_test$score) - min(dml_test$score) + .Machine$double.xmin)

    dml_sig <- dml_test[dml_test$fdrs <= fdr_thres, , drop = FALSE]
    if (nrow(dml_sig) == 0) {
        dmps <- data.frame(
            chr = character(),
            start = integer(),
            end = integer(),
            site_id = character(),
            pval = numeric(),
            qval = numeric(),
            delta_beta = numeric(),
            score = numeric(),
            stringsAsFactors = FALSE
        )
    } else {
        ord_score <- order(-dml_sig$score, dml_sig$pvals, na.last = NA)
        dml_sig <- dml_sig[ord_score, , drop = FALSE]

        dmps <- data.frame(
            chr = as.character(dml_sig$chr),
            start = as.integer(dml_sig$pos),
            end = as.integer(dml_sig$pos),
            site_id = dml_sig$site_id,
            pval = dml_sig$pvals,
            qval = dml_sig$fdrs,
            delta_beta = dml_sig$delta_beta,
            score = dml_sig$score,
            stringsAsFactors = FALSE
        )
    }
    if (nrow(dmps) > 0) {
        dmps <- dmps[order(dmps$chr, dmps$start, na.last = TRUE), , drop = FALSE]
    }
    if (!is.null(output_file)) {
        dmps_out <- dmps
        dmps_out$delta_beta <- formatC(dmps_out$delta_beta, format = "f", digits = 2)
        dmps_out$score <- formatC(dmps_out$score, format = "f", digits = 2)
        dmps_out$pval <- formatC(dmps_out$pval, format = "e", digits = 2)
        dmps_out$qval <- formatC(dmps_out$qval, format = "e", digits = 2)

        if (grepl("\\.gz$", output_file, ignore.case = TRUE)) {
            con <- gzfile(output_file, open = "wt")
            tryCatch(
                {
                    write.table(dmps_out, con, sep = "\t", quote = FALSE, row.names = FALSE)
                },
                finally = {
                    close(con)
                }
            )
        } else {
            write.table(dmps_out, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
        }
    }
    dmps
}
