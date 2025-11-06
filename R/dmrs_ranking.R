#' @keywords internal
#' @noRd
.performCrossPrediction <- function(beta_mat, groups, nfold = getOption("DMRsegal.ranking_nfold", 5)) {
    set.seed(getOption("DMRsegal.random_seed", 42))
    group_folds <- split(sample(seq_len(ncol(beta_mat))), as.factor(groups))
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
    folds <- unlist(lapply(group_folds, function(x) sample(rep(1:nfold, length.out = length(x)))))
    predictions <- vector("character", ncol(beta_mat))
    for (fold in 1:nfold) {
        train_indices <- which(folds != fold)
        test_indices <- which(folds == fold)
        model <- e1071::svm(t(beta_mat[, train_indices]), as.factor(groups[train_indices]), kernel = "radial")
        predictions[test_indices] <- predict(model, t(beta_mat[, test_indices]))
    }
    confusion_matrix <- table(Predicted = predictions, Actual = groups)
    accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    return(accuracy)
}

#' Rank DMRs Based on Classification Accuracy
#'
#' @description Ranks Differentially Methylated Regions (DMRs) based on their ability to
#' discriminate between sample groups using cross-validated Support Vector Machine (SVM)
#' classification. For each DMR, this function performs stratified k-fold cross-prediction
#' using an RBF kernel SVM to compute classification accuracy, which serves as a measure
#' of the DMR's discriminative power.
#'
#' @param dmrs Data frame or GRanges object containing DMR coordinates and metadata
#' @param beta Character. Path to beta value file, tabix file, beta matrix, BetaHandler object, or bed file
#' @param pheno Data frame. Phenotype data containing sample group information
#' @param genome Character. Genome version (e.g., "hg19", "hg38", "mm10"). Default is "hg19"
#' @param array Character. Array platform type (e.g., "450K", "EPIC", "EPICv2"). Default is "450K"
#' @param sorted_locs Data frame. Optional pre-computed sorted genomic locations. Default is NULL
#' @param sample_group_col Character. Column name in pheno containing sample group information. Default is "Sample_Group"
#' @param casecontrol_col Character. Column name in pheno for case/control status. If NULL, the first level of sample_group_col is assumed to be control. Default is NULL
#'
#' @return GRanges object with DMRs ordered by p-value and an additional metadata column:
#' \itemize{
#'   \item accuracy: Cross-validated classification accuracy for the DMR
#' }
#'
#' @details
#' The function uses stratified k-fold cross-prediction to ensure balanced representation
#' of sample groups in each fold. The number of folds can be controlled using the
#' option "DMRsegal.ranking_nfold" (default is 5). An RBF (Radial Basis Function) kernel
#' SVM is trained on the beta values of CpG sites within each DMR.
#'
#' Classification accuracy represents how well the methylation pattern of a DMR can
#' distinguish between sample groups. Higher accuracy values indicate stronger
#' discriminative power and potentially more biologically relevant regions.
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
rankDMRs <- function(dmrs, beta, pheno, genome = "hg19", array = "450K", sorted_locs = NULL,
                     sample_group_col = "Sample_Group",
                     casecontrol_col = NULL) {
    dmrs <- convertToGRanges(dmrs, genome = genome)
    beta_handler <- getBetaHandler(beta, array = array, genome = genome, sorted_locs = sorted_locs)
    supporting_sites <- getSupportingSites(dmrs, use_absolute_indices = FALSE, separate_by_section = FALSE)
    if (is.null(casecontrol_col)) {
        casecontrol_col <- "__CASE_CONTROL__"
        pheno[, casecontrol_col] <- pheno[, sample_group_col] != pheno[1, sample_group_col]
    }
    accuracies <- sapply(seq_along(dmrs), function(i) {
        site_indices <- supporting_sites[[i]]
        beta_mat <- beta_handler$getBeta(row_names = site_indices)
        cv_results <- .performCrossPrediction(beta_mat, pheno[, casecontrol_col])
        cv_results
    })
    mcols(dmrs)$accuracy <- accuracies
    mcols(dmrs)$rank <- rank(-mcols(dmrs)$accuracy, ties.method = "first")
    return(dmrs[order(mcols(dmrs)$rank)])
}
