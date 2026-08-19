# ============================================================
# 04_factor_matching.R
#
# Factor-loading comparisons across NMF representations.
# ============================================================


# Remove modality suffixes and keep unique gene rows
standardize_nmf_features <- function(model) {
  rownames(model) <- sub("_(s|us)$", "", rownames(model))
  model[!duplicated(rownames(model)), , drop = FALSE]
}


# Align models to shared genes
align_nmf_rows <- function(...) {
  models <- lapply(list(...), standardize_nmf_features)

  common_features <- Reduce(
    intersect,
    lapply(models, rownames)
  )

  lapply(
    models,
    function(model) {
      model[common_features, , drop = FALSE]
    }
  )
}


# L2-normalize factor-loading vectors
l2_normalize <- function(model) {
  sweep(
    model,
    2,
    sqrt(colSums(model^2)),
    "/"
  )
}


# Compare Added factors with the spliced and unspliced
# components of the Concatenated model
match_added_to_concat <- function(added_w, concat_w) {

  spliced_idx <- grep("_s$", rownames(concat_w))
  unspliced_idx <- grep("_us$", rownames(concat_w))

  concat_s <- concat_w[spliced_idx, , drop = FALSE]
  concat_us <- concat_w[unspliced_idx, , drop = FALSE]

  rownames(concat_s) <- sub("_s$", "", rownames(concat_s))
  rownames(concat_us) <- sub("_us$", "", rownames(concat_us))

  aligned <- align_nmf_rows(
    added_w,
    concat_s,
    concat_us
  )

  added_n <- l2_normalize(aligned[[1]])
  concat_s_n <- l2_normalize(aligned[[2]])
  concat_us_n <- l2_normalize(aligned[[3]])

  sim_to_spliced <- t(added_n) %*% concat_s_n
  sim_to_unspliced <- t(added_n) %*% concat_us_n

  data.frame(
    added_factor = seq_len(ncol(added_w)),
    max_sim_spliced = apply(sim_to_spliced, 1, max),
    max_sim_unspliced = apply(sim_to_unspliced, 1, max),
    better_match = ifelse(
      apply(sim_to_spliced, 1, max) >
        apply(sim_to_unspliced, 1, max),
      "Spliced",
      "Unspliced"
    )
  )
}


get_concat_block <- function(
    concat_w,
    modality = c(
      "spliced",
      "unspliced"
    )
) {

  modality <- match.arg(
    modality
  )

  suffix <- if (
    modality == "spliced"
  ) {
    "_s$"
  } else {
    "_us$"
  }

  row_idx <- grep(
    suffix,
    rownames(concat_w)
  )

  if (length(row_idx) == 0) {
    stop(
      "No ",
      modality,
      " features found in concatenated loading matrix."
    )
  }

  block <- concat_w[
    row_idx,
    ,
    drop = FALSE
  ]

  rownames(block) <- sub(
    suffix,
    "",
    rownames(block)
  )

  block
}
