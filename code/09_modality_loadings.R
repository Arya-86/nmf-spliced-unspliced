# ============================================================
# 09_modality_loadings.R
#
# Modality-specific loading analyses.
# ============================================================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)


plot_selected_modality_loadings <- function(
    concat_nmf,
    spliced_dominant,
    unspliced_dominant,
    n_label = 7
) {

  spliced_indices <- grep(
    "_s$",
    rownames(concat_nmf$w)
  )

  unspliced_indices <- grep(
    "_us$",
    rownames(concat_nmf$w)
  )

  focus_factors <- c(
    spliced_dominant,
    unspliced_dominant
  )

  factor_class <- tibble(
    concat_factor = focus_factors,
    class = c(
      rep(
        "Spliced-dominant",
        length(spliced_dominant)
      ),
      rep(
        "Unspliced-dominant",
        length(unspliced_dominant)
      )
    )
  )

  dissociation_df <- bind_rows(
    lapply(
      focus_factors,
      function(k) {

        w_s <- concat_nmf$w[
          spliced_indices,
          k
        ]

        w_us <- concat_nmf$w[
          unspliced_indices,
          k
        ]

        genes <- sub(
          "_(s|us)$",
          "",
          rownames(concat_nmf$w)[
            spliced_indices
          ]
        )

        data.frame(
          gene = genes,
          w_s = w_s,
          w_us = w_us
        ) %>%
          mutate(
            log2_ratio =
              log2(
                (w_s + 1e-9) /
                  (w_us + 1e-9)
              ),
            max_w =
              pmax(
                w_s,
                w_us
              ),
            concat_factor = k
          ) %>%
          filter(
            max_w > 0
          ) %>%
          arrange(
            desc(max_w)
          ) %>%
          slice_head(
            n = 50
          ) %>%
          left_join(
            factor_class,
            by = "concat_factor"
          ) %>%
          mutate(
            factor_label =
              paste0(
                "Factor ",
                k,
                "\n(",
                class,
                ")"
              ),
            label =
              ifelse(
                rank(-max_w) <= n_label,
                gene,
                NA
              )
          )
      }
    )
  )

  dissociation_df <-
    dissociation_df %>%
    mutate(
      factor_label = factor(
        factor_label,
        levels = paste0(
          "Factor ",
          focus_factors,
          "\n(",
          ifelse(
            focus_factors %in%
              spliced_dominant,
            "Spliced-dominant",
            "Unspliced-dominant"
          ),
          ")"
        )
      )
    )

  figure <- ggplot(
    dissociation_df,
    aes(
      x = w_s,
      y = w_us,
      colour = max_w,
      label = label
    )
  ) +
    geom_point(
      size = 2,
      alpha = 0.85
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey40"
    ) +
    geom_text_repel(
      size = 2.8,
      na.rm = TRUE,
      max.overlaps = 15,
      xlim = c(NA, NA),
      clip = "off"
    ) +
    facet_wrap(
      ~ factor_label,
      ncol = 2
    ) +
    scale_colour_viridis_c(
      option = "magma",
      name = "Max\nloading"
    ) +
    labs(
      x = "Spliced feature weight",
      y = "Unspliced feature weight"
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      strip.background =
        element_blank(),
      strip.text =
        element_text(
          face = "bold",
          size = 10
        )
    )

  list(
    data = dissociation_df,
    figure = figure
  )
}


plot_modality_balance_gsea <- function(
    concat_nmf,
    gene_sets
) {

  spliced_indices <- grep(
    "_s$",
    rownames(concat_nmf$w)
  )

  unspliced_indices <- grep(
    "_us$",
    rownames(concat_nmf$w)
  )

  modality_scores <- data.frame(
    nmf_factor =
      seq_len(
        ncol(concat_nmf$w)
      ),
    spliced =
      colSums(
        concat_nmf$w[
          spliced_indices,
        ]
      ),
    unspliced =
      colSums(
        concat_nmf$w[
          unspliced_indices,
        ]
      )
  )

  modality_scores$spliced_fraction <-
    modality_scores$spliced /
    (
      modality_scores$spliced +
        modality_scores$unspliced
    )

  p1 <- ggplot(
    modality_scores,
    aes(
      x = "factors",
      y = spliced_fraction
    )
  ) +
    geom_violin(
      fill = "white",
      color = "black",
      linewidth = 1
    ) +
    geom_jitter(
      width = 0.08,
      size = 2.8,
      alpha = 0.9
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    theme_classic(
      base_size = 18
    ) +
    labs(
      x = "",
      y = "Proportion of Spliced Weights per Factor"
    ) +
    theme(
      axis.text.x =
        element_text(size = 16),
      axis.title.y =
        element_text(size = 22),
      axis.text.y =
        element_text(size = 16),
      aspect.ratio = 2,
      plot.margin =
        margin(10, 10, 10, 10)
    )

  all_gsea <- list()

  cat("Running GSEA:\n")

  for (
    i in seq_len(
      ncol(concat_nmf$w)
    )
  ) {

    cat(i, " ")

    s_weights <- concat_nmf$w[
      spliced_indices,
      i
    ]

    names(s_weights) <- rownames(
      concat_nmf$w[
        spliced_indices,
      ]
    )

    s_res <- run_gsea_fgsea(
      s_weights,
      gene_sets
    )

    us_weights <- concat_nmf$w[
      unspliced_indices,
      i
    ]

    names(us_weights) <- rownames(
      concat_nmf$w[
        unspliced_indices,
      ]
    )

    us_res <- run_gsea_fgsea(
      us_weights,
      gene_sets
    )

    if (
      !is.null(s_res) &&
      !is.null(us_res)
    ) {

      s_df <- s_res %>%
        dplyr::select(
          pathway,
          padj
        ) %>%
        dplyr::rename(
          padj_s = padj
        )

      us_df <- us_res %>%
        dplyr::select(
          pathway,
          padj
        ) %>%
        dplyr::rename(
          padj_us = padj
        )

      merged <- full_join(
        s_df,
        us_df,
        by = "pathway"
      )

      merged$factor <- i

      all_gsea[[i]] <- merged
    }
  }

  cat("\nDone.\n")

  gsea_df <- bind_rows(
    all_gsea
  )

  gsea_df$padj_s[
    is.na(gsea_df$padj_s)
  ] <- 1

  gsea_df$padj_us[
    is.na(gsea_df$padj_us)
  ] <- 1

  gsea_df$padj_s <-
    -log10(
      gsea_df$padj_s +
        1e-300
    )

  gsea_df$padj_us <-
    -log10(
      gsea_df$padj_us +
        1e-300
    )

  p2 <- ggplot(
    gsea_df,
    aes(
      x = padj_s,
      y = padj_us
    )
  ) +
    geom_point(
      size = 0.5,
      alpha = 0.35
    ) +
    geom_vline(
      xintercept =
        -log10(0.05),
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept =
        -log10(0.05),
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.8
    ) +
    coord_cartesian(
      xlim = c(0, 20),
      ylim = c(0, 20)
    ) +
    theme_classic(
      base_size = 18
    ) +
    labs(
      x = expression(
        -log[10](padj) ~
          "(spliced)"
      ),
      y = expression(
        -log[10](padj) ~
          "(unspliced)"
      )
    ) +
    theme(
      axis.title =
        element_text(size = 22),
      axis.text =
        element_text(size = 16),
      aspect.ratio = 1,
      plot.margin =
        margin(10, 10, 10, 10)
    )

  figure <- cowplot::plot_grid(
    p1,
    p2,
    nrow = 1,
    labels = c("A", "B"),
    rel_widths = c(1, 2)
  )

  list(
    modality_scores =
      modality_scores,
    gsea_results =
      all_gsea,
    gsea_plot_data =
      gsea_df,
    figure =
      figure
  )
}
