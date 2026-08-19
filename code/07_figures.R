# ============================================================
# 07_figures.R
#
# Plotting functions used for manuscript figures.
# ============================================================

library(ggplot2)
library(viridis)
library(cowplot)


# ------------------------------------------------------------
# PCA scree plot
# ------------------------------------------------------------

plot_pca_scree <- function(scree_df) {

  ggplot(
    scree_df,
    aes(
      x = PC,
      y = VarianceExplained,
      color = Modality
    )
  ) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    geom_point(size = 1.5, alpha = 0.9) +
    scale_x_continuous(
      breaks = seq(0, 50, 10)
    ) +
    scale_y_continuous(
      expand = c(0.01, 0.01)
    ) +
    labs(
      x = "Principal Component",
      y = "Explained Variance",
      color = "Modality"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}


# ------------------------------------------------------------
# NMF cross-validation figure
# ------------------------------------------------------------

plot_cv <- function(df) {

  lower_limit <- min(df$test_err)
  upper_limit <- max(df$test_err[df$iter == min(df$iter)])

  p1 <- ggplot(
    df,
    aes(
      x = k,
      y = test_err,
      color = iter,
      group = iter
    )
  ) +
    geom_line() +
    scale_color_viridis_c(option = "B") +
    facet_grid(cols = vars(group)) +
    scale_y_continuous(
      limits = c(lower_limit, upper_limit)
    ) +
    labs(
      x = "NMF Model Rank",
      y = "Test Set Loss",
      color = "Epoch"
    ) +
    theme_classic() +
    theme(aspect.ratio = 1)

  plot_data <- subset(
    df,
    k %in% c(5, 30, 80)
  )

  plot_data$k <- factor(
    plot_data$k,
    levels = c(5, 30, 80),
    labels = c(
      "Underfit (k = 5)",
      "Better Fit (k = 30)",
      "Overfit (k = 80)"
    )
  )

  p2 <- ggplot(
    plot_data,
    aes(
      x = iter,
      y = test_err,
      color = k,
      group = k
    )
  ) +
    geom_line(linewidth = 0.5) +
    scale_color_viridis_d(
      option = "B",
      end = 0.9,
      begin = 0.1
    ) +
    facet_grid(cols = vars(group)) +
    scale_y_continuous(
      limits = c(lower_limit, upper_limit)
    ) +
    labs(
      x = "NMF Model Epoch",
      y = "Test Set Loss",
      color = "Rank"
    ) +
    theme_classic() +
    theme(aspect.ratio = 1)

  cowplot::plot_grid(
    p1,
    p2,
    labels = c("A", "B"),
    nrow = 2,
    align = "hv"
  )
}


# ------------------------------------------------------------
# GSEA comparison for spliced and unspliced loadings
# ------------------------------------------------------------

plot_gsea_velocity <- function(
    gsea_results,
    factor_pos = 1,
    pathways_per_panel = 5
) {

  gsea_s <- gsea_results[[factor_pos]]$spliced
  gsea_us <- gsea_results[[factor_pos]]$unspliced

  rownames(gsea_s) <- gsea_s$pathway
  rownames(gsea_us) <- gsea_us$pathway

  gsea_s$pathway <- NULL
  gsea_us$pathway <- NULL

  plot_data <- merge(
    gsea_s,
    gsea_us,
    by = "row.names",
    all = TRUE
  )

  plot_data <- data.frame(
    pathway = plot_data$Row.names,
    padj_s = plot_data$padj.x,
    padj_us = plot_data$padj.y,
    NES_s = plot_data$NES.x,
    NES_us = plot_data$NES.y
  )

  plot_data$padj_s[is.na(plot_data$padj_s)] <- 1
  plot_data$padj_us[is.na(plot_data$padj_us)] <- 1
  plot_data$NES_s[is.na(plot_data$NES_s)] <- 0
  plot_data$NES_us[is.na(plot_data$NES_us)] <- 0

  plot_data <- subset(
    plot_data,
    padj_s < 0.05 | padj_us < 0.05
  )

  plot_data <- plot_data[
    grep("GOBP_", plot_data$pathway),
    ,
    drop = FALSE
  ]

  fold_change <- function(x, y) {
    denom <- x + y
    denom[denom == 0] <- NA
    ifelse(
      x > y,
      (x - y) / denom,
      -(y - x) / denom
    )
  }

  plot_data$NES_fc <- fold_change(
    plot_data$NES_us,
    plot_data$NES_s
  )

  plot_data$pathway <- gsub(
    "GOBP_",
    "",
    plot_data$pathway
  )

  n <- nrow(plot_data)
  ord_fc <- order(
    plot_data$NES_fc,
    decreasing = TRUE
  )

  idx_unspliced <- ord_fc[
    seq_len(pathways_per_panel)
  ]

  idx_spliced <- ord_fc[
    (n - pathways_per_panel + 1):n
  ]

  remaining <- setdiff(
    seq_len(n),
    c(idx_unspliced, idx_spliced)
  )

  idx_shared <- remaining[
    order(abs(plot_data$NES_fc[remaining]))
  ][seq_len(pathways_per_panel)]

  plot_data$group <- NA_character_

  plot_data$group[idx_unspliced] <- "most unspliced"
  plot_data$group[idx_spliced] <- "most spliced"
  plot_data$group[idx_shared] <- "most shared"

  plot_data <- plot_data[
    !is.na(plot_data$group),
    ,
    drop = FALSE
  ]

  p1 <- plot_data[
    ,
    c("pathway", "padj_s", "NES_s", "group")
  ]

  p2 <- plot_data[
    ,
    c("pathway", "padj_us", "NES_us", "group")
  ]

  p1$mode <- "spliced"
  p2$mode <- "unspliced"

  colnames(p1) <-
    colnames(p2) <-
    c("pathway", "padj", "NES", "group", "mode")

  plot_data_long <- rbind(
    p1,
    p2
  )

  plot_data_long$group <- factor(
    plot_data_long$group,
    levels = c(
      "most shared",
      "most spliced",
      "most unspliced"
    )
  )

  ggplot(
    plot_data_long,
    aes(
      x = NES,
      y = pathway,
      color = mode,
      size = -log10(padj)
    )
  ) +
    geom_point(alpha = 0.9) +
    facet_grid(
      rows = vars(group),
      scales = "free_y",
      space = "free_y"
    ) +
    scale_x_continuous(
      limits = c(0, 3),
      expand = c(0, 0)
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 7)
    ) +
    labs(
      x = "NES",
      y = NULL,
      color = "mode",
      size = expression(-log[10](padj))
    )
}
