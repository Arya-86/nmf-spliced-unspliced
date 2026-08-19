# ============================================================
# 08_factor_similarity.R
#
# Factor similarity and bipartite matching analysis.
# ============================================================

library(RcppML)
library(reshape2)
library(ggplot2)
library(cowplot)


get_bipartite_match_costs <- function(cost) {

  if (ncol(cost) < nrow(cost)) {
    cost <- t(cost)
  }

  match_result <- bipartiteMatch(cost)
  assignment <- match_result$pairs

  if (is.null(assignment)) {

    assignment <- match_result$assignment

    if (is.null(assignment)) {
      stop(
        "bipartiteMatch() returned neither 'pairs' nor 'assignment'."
      )
    }

    assignment <- assignment + 1L
  }

  assignment <- unlist(
    assignment,
    use.names = FALSE
  )

  row_idx <- seq_along(assignment)

  valid_idx <-
    !is.na(assignment) &
    assignment >= 1 &
    assignment <= ncol(cost)

  row_idx <- row_idx[valid_idx]
  assignment <- assignment[valid_idx]

  if (length(assignment) == 0) {
    stop(
      "bipartiteMatch() returned no valid assignments."
    )
  }

  cost[cbind(row_idx, assignment)]
}


factor_similarity_analysis <- function(models) {

  sus_nmf <- models$concat$w
  s_plus_us_nmf <- models$added$w
  s_nmf <- models$spliced$w
  us_nmf <- models$unspliced$w

  concat_s_nmf <- get_concat_block(
    sus_nmf,
    "spliced"
  )

  concat_us_nmf <- get_concat_block(
    sus_nmf,
    "unspliced"
  )

  aligned_models <- align_nmf_rows(
    concat_s_nmf,
    concat_us_nmf,
    s_plus_us_nmf,
    s_nmf,
    us_nmf
  )

  concat_s_nmf <- aligned_models[[1]]
  concat_us_nmf <- aligned_models[[2]]
  s_plus_us_nmf <- aligned_models[[3]]
  s_nmf <- aligned_models[[4]]
  us_nmf <- aligned_models[[5]]

  w <- cbind(
    concat_s_nmf,
    s_plus_us_nmf,
    s_nmf,
    us_nmf
  )

  modalities <- c(
    rep("concatenated", ncol(concat_s_nmf)),
    rep("added", ncol(s_plus_us_nmf)),
    rep("spliced", ncol(s_nmf)),
    rep("unspliced", ncol(us_nmf))
  )

  w <- l2_normalize(w)
  G <- crossprod(w)

  G_adj <- RcppML::cosine(
    concat_us_nmf,
    us_nmf
  )

  G[
    1:ncol(concat_s_nmf),
    (ncol(G) - ncol(us_nmf) + 1):ncol(G)
  ] <- G_adj

  G[lower.tri(G)] <- t(G)[lower.tri(G)]

  order <- c(
    hclust(
      dist(
        t(G[, which(modalities == "concatenated")])
      )
    )$order,

    hclust(
      dist(
        t(G[, which(modalities == "added")])
      )
    )$order +
      min(which(modalities == "added")) - 1,

    hclust(
      dist(
        t(G[, which(modalities == "spliced")])
      )
    )$order +
      min(which(modalities == "spliced")) - 1,

    hclust(
      dist(
        t(G[, which(modalities == "unspliced")])
      )
    )$order +
      min(which(modalities == "unspliced")) - 1
  )

  G <- G[order, order]

  G_ <- G

  rownames(G_) <-
    colnames(G_) <-
    seq_len(ncol(G_))

  rownames(G) <-
    colnames(G) <-
    modalities

  df <- reshape2::melt(G)
  df2 <- reshape2::melt(G_)

  plot_data <- data.frame(
    x = df2$Var1,
    y = df2$Var2,
    sim = df2$value,
    modality1 = df$Var1,
    modality2 = df$Var2
  )

  p1 <- ggplot(
    plot_data,
    aes(x, y, fill = sim)
  ) +
    geom_tile() +
    theme_classic() +
    scale_fill_viridis_c(
      option = "B",
      end = 0.9
    ) +
    scale_y_continuous(
      expand = c(0.01, 0.01)
    ) +
    scale_x_reverse(
      expand = c(0.01, 0.01)
    ) +
    facet_grid(
      rows = vars(modality2),
      cols = vars(modality1),
      space = "free",
      scales = "free"
    ) +
    theme(
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank()
    ) +
    labs(
      fill = "Similarity",
      x = "NMF factors",
      y = "NMF factors"
    )

  models_to_compare <- list(
    concatenated = concat_s_nmf,
    added = s_plus_us_nmf,
    spliced = s_nmf,
    unspliced = us_nmf
  )

  results <- list()

  for (i in seq_along(models_to_compare)) {

    for (j in seq_along(models_to_compare)) {

      if (i != j) {

        model1_name <-
          names(models_to_compare)[[i]]

        model2_name <-
          names(models_to_compare)[[j]]

        model1 <- l2_normalize(
          models_to_compare[[i]]
        )

        model2 <- l2_normalize(
          models_to_compare[[j]]
        )

        if (
          model1_name == "concatenated" &&
          model2_name == "unspliced"
        ) {
          model1 <- l2_normalize(
            concat_us_nmf
          )
        }

        if (
          model2_name == "concatenated" &&
          model1_name == "unspliced"
        ) {
          model2 <- l2_normalize(
            concat_us_nmf
          )
        }

        cost <- 1 - crossprod(
          model1,
          model2
        )

        matched_costs <-
          get_bipartite_match_costs(
            cost
          )

        results[[length(results) + 1]] <-
          data.frame(
            model1 = model1_name,
            model2 = model2_name,
            cost = matched_costs
          )
      }
    }
  }

  results <- do.call(
    rbind,
    results
  )

  p2 <- ggplot(
    results,
    aes(
      model1,
      cost,
      color = model2
    )
  ) +
    geom_point(
      position =
        position_jitterdodge(
          jitter.width = 0.1
        ),
      size = 1
    ) +
    geom_boxplot() +
    theme_classic() +
    labs(
      x = "NMF Model Modality",
      y = "Cost of factor bipartite matching",
      color = "Matched Modality"
    ) +
    theme(
      aspect.ratio = 1.5,
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1,
          vjust = 1
        )
    )

  figure <- cowplot::plot_grid(
    p1,
    p2,
    nrow = 1,
    rel_widths = c(1, 0.7),
    labels = c("A", "B")
  )

  list(
    similarity_matrix = G,
    matching_results = results,
    figure = figure
  )
}
