# ============================================================
# 02_cross_validation.R
#
# NMF cross-validation workflow used for rank sweeps.
#
# Manuscript analysis:
# - Ranks evaluated: k = 2:80
# - Maximum iterations: 50
# - Held-out test MSE traced at every iteration
# - Three replicate CV results were used for the final
#   supplementary rank-sweep tables.
#
# Revealable rank (k_star) was selected manually for each
# replicate from the epoch-50 held-out test-MSE curve as the
# largest rank in the near-minimum stable region immediately
# preceding sustained overfitting.
#
# Automatic minimum-error ranks are retained separately for
# downstream analyses.
# ============================================================

library(singlet)
library(dplyr)


cross_val <- function(
    A,
    ranks = 2:80,
    n_replicates = 3
) {

  cv_data <- singlet::cross_validate_nmf(
    A,
    ranks,
    n_replicates = n_replicates,
    trace_test_mse = 1,
    maxit = 50,
    tol = -1,
    tol_overfit = Inf
  )

  cv_data$test_error <- as.numeric(cv_data$test_error)

  cv_data
}


get_best_rank <- function(cv_data) {

  final_iter <- max(cv_data$iter)

  cv_data |>
    filter(iter == final_iter) |>
    group_by(k) |>
    summarise(
      mean_test_err = mean(test_error),
      .groups = "drop"
    ) |>
    arrange(mean_test_err) |>
    slice(1) |>
    pull(k)
}
