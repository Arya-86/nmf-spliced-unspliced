# ============================================================
# 03_fit_models.R
#
# Fit final NMF models for the four expression representations.
# Ranks are supplied from the preceding cross-validation step.
# ============================================================

library(singlet)


# ------------------------------------------------------------
# Fit final NMF models
# ------------------------------------------------------------

fit_nmf_models <- function(
    dataset,
    ranks,
    seed = 123
) {

  required_ranks <- c(
    "spliced",
    "unspliced",
    "added",
    "concat"
  )

  stopifnot(
    all(required_ranks %in% names(ranks))
  )

  set.seed(seed)

  spliced_nmf <- singlet::run_nmf(
    dataset$concat[
      dataset$spliced_indices,
      ,
      drop = FALSE
    ],
    rank = ranks[["spliced"]],
    L1 = 0,
    tol = 1e-5
  )

  unspliced_nmf <- singlet::run_nmf(
    dataset$concat[
      dataset$unspliced_indices,
      ,
      drop = FALSE
    ],
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


  # Check that the expected feature labels are preserved
  stopifnot(
    !any(grepl("_s$|_us$", rownames(added_nmf$w))),
    all(grepl("_s$|_us$", rownames(concat_nmf$w)))
  )


  list(
    spliced = spliced_nmf,
    unspliced = unspliced_nmf,
    added = added_nmf,
    concat = concat_nmf
  )
}
