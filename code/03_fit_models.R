# ============================================================
# 03_fit_models.R
#
# Fit final NMF models for each expression representation.
# ============================================================

library(singlet)


fit_nmf_models <- function(
    dataset,
    ranks,
    seed = 123
) {

  set.seed(seed)

  spliced_nmf <- singlet::run_nmf(
    dataset$concat[dataset$spliced_indices, , drop = FALSE],
    rank = ranks[["spliced"]],
    L1 = 0,
    tol = 1e-5
  )

  unspliced_nmf <- singlet::run_nmf(
    dataset$concat[dataset$unspliced_indices, , drop = FALSE],
    rank = ranks[["unspliced"]],
    L1 = 0,
    tol = 1e-5
  )

  added_nmf <- singlet::run_nmf(
    dataset$added,
    rank = ranks[["added"]],
    L1 = 0,
    tol = 1e-5
  )

  concat_nmf <- singlet::run_nmf(
    dataset$concat,
    rank = ranks[["concat"]],
    L1 = 0,
    tol = 1e-5
  )

  list(
    spliced = spliced_nmf,
    unspliced = unspliced_nmf,
    added = added_nmf,
    concat = concat_nmf
  )
}
