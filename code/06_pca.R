# ============================================================
# 06_pca.R
#
# PCA scree analysis for each expression representation.
# ============================================================

library(Matrix)
library(irlba)


get_pca_scree <- function(
    expr_matrix,
    label,
    n_pcs = 50
) {

  expr_matrix <- as.matrix(expr_matrix)

  expr_matrix <- expr_matrix[
    Matrix::rowSums(expr_matrix) > 0,
  ]

  expr_matrix <- t(
    scale(
      t(expr_matrix),
      center = TRUE,
      scale = TRUE
    )
  )

  pca_res <- prcomp_irlba(
    expr_matrix,
    n = n_pcs,
    center = FALSE,
    scale. = FALSE
  )

  var_explained <-
    pca_res$sdev^2 /
    sum(pca_res$sdev^2)

  data.frame(
    PC = seq_along(var_explained),
    VarianceExplained = var_explained,
    Modality = label
  )
}
