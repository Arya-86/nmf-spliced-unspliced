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
# IMPORTANT:
# The final revealable rank (k_star) was NOT selected
# automatically by minimum test error.
#
# k_star was manually selected for each replicate from the
# epoch-50 held-out test-MSE curve as the largest rank in the
# near-minimum stable region immediately preceding sustained
# overfitting.
#
# ============================================================


library(singlet)
library(dplyr)


# ------------------------------------------------------------
# Run NMF cross-validation
# ------------------------------------------------------------

cross_val <- function(
    A,
    ranks = 2:80,
    n_replicates = 3,
    maxit = 50
) {

  cv_data <- singlet::cross_validate_nmf(
    A,
    ranks,
    n_replicates = n_replicates,
    trace_test_mse = 1,
    maxit = maxit,
    tol = -1,
    tol_overfit = Inf
  )

  cv_data$test_error <- as.numeric(cv_data$test_error)

  cv_data
}


# ------------------------------------------------------------
# Automatic minimum-error rank
# ------------------------------------------------------------
#
# This helper is retained because minimum-error ranks were
# used operationally in parts of the downstream workflow.
#
# IMPORTANT:
# This is NOT the manuscript revealable rank (k_star).
# The reported k_star values were selected manually from
# replicate-specific epoch-50 CV curves.
# ------------------------------------------------------------

get_best_rank <- function(cv_data) {

  final_iter <- max(cv_data$iter)

  final_cv <- cv_data |>
    dplyr::filter(iter == final_iter) |>
    dplyr::group_by(k) |>
    dplyr::summarise(
      mean_test_err = mean(test_error),
      .groups = "drop"
    ) |>
    dplyr::arrange(mean_test_err)

  final_cv$k[1]
}


# ------------------------------------------------------------
# Extract epoch-50 held-out test MSE
# ------------------------------------------------------------
#
# This preserves replicate-level values when replicate
# information is present in the output.
# ------------------------------------------------------------

extract_final_epoch <- function(
    cv_data,
    final_iter = 50
) {

  cv_data |>
    dplyr::filter(iter == final_iter)
}


# ------------------------------------------------------------
# Run CV for all four representations
# ------------------------------------------------------------

run_cv_all_representations <- function(
    dataset,
    ranks = 2:80,
    n_replicates = 3,
    maxit = 50
) {

  message("Running concatenated representation...")

  concat_cv <- cross_val(
    dataset$concat,
    ranks = ranks,
    n_replicates = n_replicates,
    maxit = maxit
  )


  message("Running spliced-only representation...")

  spliced_cv <- cross_val(
    dataset$concat[
      dataset$spliced_indices,
      ,
      drop = FALSE
    ],
    ranks = ranks,
    n_replicates = n_replicates,
    maxit = maxit
  )


  message("Running unspliced-only representation...")

  unspliced_cv <- cross_val(
    dataset$concat[
      dataset$unspliced_indices,
      ,
      drop = FALSE
    ],
    ranks = ranks,
    n_replicates = n_replicates,
    maxit = maxit
  )


  message("Running added representation...")

  added_cv <- cross_val(
    dataset$added,
    ranks = ranks,
    n_replicates = n_replicates,
    maxit = maxit
  )


  list(
    concat = concat_cv,
    spliced = spliced_cv,
    unspliced = unspliced_cv,
    added = added_cv
  )
}


# ------------------------------------------------------------
# Automatic minimum-error ranks
# ------------------------------------------------------------

get_all_best_ranks <- function(cv_results) {

  c(
    concat = get_best_rank(cv_results$concat),
    spliced = get_best_rank(cv_results$spliced),
    unspliced = get_best_rank(cv_results$unspliced),
    added = get_best_rank(cv_results$added)
  )
}


# ------------------------------------------------------------
# Extract final-epoch rank-sweep results
# ------------------------------------------------------------

get_final_epoch_results <- function(
    cv_results,
    final_iter = 50
) {

  list(
    concat = extract_final_epoch(
      cv_results$concat,
      final_iter
    ),

    spliced = extract_final_epoch(
      cv_results$spliced,
      final_iter
    ),

    unspliced = extract_final_epoch(
      cv_results$unspliced,
      final_iter
    ),

    added = extract_final_epoch(
      cv_results$added,
      final_iter
    )
  )
}
