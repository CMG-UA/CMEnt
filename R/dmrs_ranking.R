#' @keywords internal
#' @noRd
.performCrossPrediction <- function(beta_mat, groups, nfold = getOption("DMRsegal.ranking_nfold", 5)) {
    set.seed(getOption("DMRsegal.random_seed", 42))
    groups <- as.factor(groups)
    if (ncol(beta_mat) != length(groups)) {
        stop(
            "Mismatch between beta matrix columns (", ncol(beta_mat),
            ") and group labels (", length(groups), ")."
        )
    }
    if (nlevels(groups) < 2) {
        stop("Ranking requires at least two classes in '__casecontrol__'.")
    }
    group_folds <- split(seq_len(ncol(beta_mat)), groups)
    for (g in names(group_folds)) {
        if (length(group_folds[[g]]) < nfold) {
            gsize <- length(group_folds[[g]])
            stop(paste0(
                "Number of samples in group (", gsize, ") '", g, "' is less than nfold = ", nfold,
                ". Cannot perform stratified cross-prediction. Reduce nfold using options(DMRsegal.ranking_nfold=", gsize,
                ") or increase number of samples in this group."
            ))
        }
    }
    folds <- integer(ncol(beta_mat))
    for (g in names(group_folds)) {
        idx <- group_folds[[g]]
        folds[idx] <- sample(rep(seq_len(nfold), length.out = length(idx)))
    }
    predictions <- vector("character", ncol(beta_mat))
    decision_values <- rep(NA_real_, ncol(beta_mat))
    for (fold in 1:nfold) {
        train_indices <- which(folds != fold)
        test_indices <- which(folds == fold)
        train_groups <- as.factor(groups[train_indices])
        model <- e1071::svm(
            t(beta_mat[, train_indices]),
            train_groups,
            kernel = "radial",
            scale = TRUE
        )
        fold_pred <- predict(model, t(beta_mat[, test_indices]), decision.values = TRUE)
        fold_pred_chr <- as.character(fold_pred)
        predictions[test_indices] <- fold_pred_chr

        fold_decision <- as.numeric(attr(fold_pred, "decision.values"))
        if (length(fold_decision) != length(test_indices)) {
            fold_decision <- rep(NA_real_, length(test_indices))
        } else {
            train_levels <- levels(train_groups)
            pred_from_sign <- ifelse(fold_decision >= 0, train_levels[1], train_levels[2])
            if (mean(pred_from_sign == fold_pred_chr, na.rm = TRUE) < 0.5) {
                fold_decision <- -fold_decision
            }
        }
        decision_values[test_indices] <- fold_decision
    }
    confusion_matrix <- table(Predicted = predictions, Actual = as.character(groups))
    cv_accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)

    # Margin-sensitive score (in (0, 1]) based on logistic loss of SVM decision values.
    # This breaks ties among perfect-accuracy DMRs by rewarding larger separating margins.
    y <- ifelse(as.character(groups) == levels(groups)[1], 1, -1)
    finite_mask <- is.finite(decision_values)
    if (!any(finite_mask)) {
        margin_score <- cv_accuracy
    } else {
        logistic_loss <- log1p(exp(-y[finite_mask] * decision_values[finite_mask]))
        margin_score <- exp(-mean(logistic_loss))
    }
    c(score = margin_score, cv_accuracy = cv_accuracy)
}

#' Rank DMRs Based on Classification Score
#'
#' @description Ranks Differentially Methylated Regions (DMRs) based on their ability to
#' discriminate between sample groups using cross-validated Support Vector Machine (SVM)
#' classification. For each DMR, this function performs stratified k-fold cross-prediction
#' using an RBF kernel SVM and computes a margin-sensitive classification score based on
#' decision values, which serves as a measure of the DMR's discriminative power.
#'
#' @param dmrs Data frame or GRanges object containing DMR coordinates and metadata
#' @param beta Character. Path to beta value file, tabix file, beta matrix, BetaHandler object, or bed file
#' @param pheno Data frame. Phenotype data containing sample group information
#' @param genome Character. Genome version (e.g., "hg19", "hg38", "mm10"). Default is "hg19"
#' @param array Character. Array platform type (e.g., "450K", "EPIC", "EPICv2"). Default is "450K"
#' @param sorted_locs Data frame. Optional pre-computed sorted genomic locations. Default is NULL
#' @param sample_group_col Character. Column name in pheno containing sample group information. Default is "Sample_Group"
#'
#' @return GRanges object with DMRs ordered by score and additional metadata columns:
#' \itemize{
#'   \item score: Margin-sensitive cross-validated classification score for the DMR
#'   \item cv_accuracy: Raw cross-validated classification accuracy for the DMR
#' }
#'
#' @details
#' The function uses stratified k-fold cross-prediction to ensure balanced representation
#' of sample groups in each fold. The number of folds can be controlled using the
#' option "DMRsegal.ranking_nfold" (default is 5). An RBF (Radial Basis Function) kernel
#' SVM is trained on the beta values of CpG sites within each DMR.
#'
#' The `score` combines classification correctness and margin confidence,
#' making it more sensitive than plain cross-validated accuracy when many DMRs
#' classify perfectly. The `cv_accuracy` column stores the raw cross-validated
#' accuracy for reference.
#'
#' @examples
#' # Load example data
#' beta <- loadExampleInputData("beta")
#' pheno <- loadExampleInputData("pheno")
#'
#' # Load pre-computed DMRs
#' dmrs <- readRDS(system.file("extdata", "example_output.rds", package = "DMRsegal"))
#'
#' # Rank DMRs
#' ranked_dmrs <- rankDMRs(
#'     dmrs = dmrs,
#'     beta = beta,
#'     pheno = pheno,
#'     sample_group_col = "Sample_Group"
#' )
#'
#' @export
rankDMRs <- function(
    dmrs, beta, pheno, covariates = NULL,
    genome = "hg19", array = "450K", sorted_locs = NULL,
    sample_group_col = "Sample_Group"
) {
    verbose <- getOption("DMRsegal.verbose", 1)
    dmrs <- convertToGRanges(dmrs, genome = genome)
    beta_handler <- getBetaHandler(beta, array = array, genome = genome, sorted_locs = sorted_locs)
    beta_col_names <- beta_handler$getBetaColNames()
    missing_pheno_samples <- setdiff(beta_col_names, rownames(pheno))
    if (length(missing_pheno_samples) > 0) {
        stop(
            "The following beta samples are missing from pheno row names: ",
            paste(head(missing_pheno_samples, 10), collapse = ","),
            if (length(missing_pheno_samples) > 10) " ..." else ""
        )
    }
    pheno <- pheno[beta_col_names, , drop = FALSE]
    supporting_sites <- .getSupportingSitesFromColumns(
        dmrs,
        use_absolute_indices = FALSE,
        separate_by_section = FALSE,
        beta_locs = beta_handler$getBetaLocs()
    )
    if (! "__casecontrol__" %in% colnames(pheno)) {
        pheno[, "__casecontrol__"] <- pheno[, sample_group_col] != pheno[1, sample_group_col]
    }
    class_values <- unique(pheno[, "__casecontrol__"])
    class_values <- class_values[!is.na(class_values)]
    if (length(class_values) < 2) {
        stop("Ranking requires at least two classes in '__casecontrol__'.")
    }
    .setupParallel()
    p_con <- NULL
    if (verbose > 0) {
        # check if version of progressr is equal or higher than >= 0.17.0-9002, otherwise p_con will not be used

        if (utils::packageVersion("progressr") >= "0.17.0-9002") {
            p_con <- progressr::progressor(steps = length(dmrs), message = "Ranking DMRs..")
        }
    }
    design <- NULL
    if (!is.null(covariates)) {
        design <- as.data.frame(pheno[, covariates, drop = FALSE])
        design <- data.frame(`(Intercept)` = 1, design, check.names = FALSE)
        # convert any string columns to factors
        for (col in colnames(design)) {
            if (is.character(design[[col]])) {
                design[[col]] <- as.numeric(as.factor(design[[col]]))
            }
        }
        design <- as.matrix(design)
        storage.mode(design) <- "double"
    }
    cv_metrics <- future.apply::future_sapply(
        X = seq_along(dmrs),
        FUN = function(i) {
            site_indices <- supporting_sites[[i]]
            beta_mat <- beta_handler$getBeta(row_names = site_indices, col_names = beta_col_names)
            # Convert to M-values
            beta_mat <- log2(beta_mat / (1 - beta_mat + 1e-6) + 1e-6)
            if (!is.null(covariates)) {
                beta_mat <- .remove_confounder_effect(beta_mat, design)
            }
            cv_results <- .performCrossPrediction(beta_mat, pheno[, "__casecontrol__"])
            if (verbose > 0 && !is.null(p_con)) p_con()
            cv_results
        },
        future.seed = TRUE,
        future.globals = c("supporting_sites", "beta_handler", "pheno", "p_con", ".performCrossPrediction", "covariates", "design"),
        future.stdout = NA
    )
    .finalizeParallel()
    mcols(dmrs)$score <- as.numeric(cv_metrics["score", ])
    mcols(dmrs)$cv_accuracy <- as.numeric(cv_metrics["cv_accuracy", ])
    mcols(dmrs)$rank <- as.numeric(as.factor(rank(-mcols(dmrs)$score, ties.method = "first")))
    dmrs[order(mcols(dmrs)$rank)]
}
