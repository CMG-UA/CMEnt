test_that("SVM scoring maps binary groups to minus one and plus one", {
    values <- rep(c(1, 0), 3L)
    folds <- rep(1:3, each = 2L)

    for (groups in list(values, factor(values), factor(values, levels = c(1, 0)), values == 1)) {
        plan <- CMEnt:::.prepareCrossPredictionFolds(groups, folds = folds)
        expect_identical(plan$groups_y, rep(c(1, -1), 3L))
    }
})

test_that("SVM margin scores retain their value with the binary sign mapping", {
    values <- rep(0:1, each = 6L)
    folds <- rep(1:3, 4L)
    beta <- rbind(
        c(0.10, 0.12, 0.14, 0.11, 0.13, 0.15, 0.85, 0.87, 0.89, 0.86, 0.88, 0.90),
        c(0.20, 0.22, 0.24, 0.21, 0.23, 0.25, 0.75, 0.77, 0.79, 0.76, 0.78, 0.80)
    )

    # The score was captured before changing the sign convention. Exercise both
    # class orders because e1071's raw decision orientation can change with them.
    for (groups in list(factor(values), factor(values, levels = c(1, 0)))) {
        for (ord in list(seq_along(groups), rev(seq_along(groups)))) {
            result <- CMEnt:::.performCrossPrediction(
                beta[, ord, drop = FALSE], groups[ord], folds = folds[ord]
            )
            expect_equal(result, c(score = 0.7302394460844595, cv_accuracy = 1), tolerance = 1e-6)
        }
    }
})
