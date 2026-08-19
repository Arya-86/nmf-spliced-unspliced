# ============================================================
# 09_modality_loadings.R
#
# Modality-specific loading analyses for Figures 4 and 5.
# ============================================================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)


# ============================================================
# Figure 4
# PBMC spliced vs unspliced feature loadings
# ============================================================

pbmc_concat_nmf <- pbmc_models$concat

spliced_indices <- grep(
  "_s$",
  rownames(pbmc_concat_nmf$w)
)

unspliced_indices <- grep(
  "_us$",
  rownames(pbmc_concat_nmf$w)
)


spliced_dominant <- c(
  34,
  35
)

unspliced_dominant <- c(
  23,
  28
)

focus_factors <- c(
  spliced_dominant,
  unspliced_dominant
)


factor_class <- tibble(
  concat_factor = focus_factors,
  class = c(
    "Spliced-dominant",
    "Spliced-dominant",
    "Unspliced-dominant",
    "Unspliced-dominant"
  )
)


n_label <- 7


dissociation_df <- bind_rows(
  lapply(
    focus_factors,
    function(k) {

      w_s <- pbmc_concat_nmf$w[
        spliced_indices,
        k
      ]

      w_us <- pbmc_concat_nmf$w[
        unspliced_indices,
        k
      ]

      genes <- sub(
        "_(s|us)$",
        "",
        rownames(
          pbmc_concat_nmf$w
        )[spliced_indices]
      )

      df <- data.frame(
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
              rank(-max_w) <=
                n_label,
              gene,
              NA
            )
        )

      df
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


figure4 <- ggplot(
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


# ============================================================
# Figure 5
# Heart modality preference and GSEA comparison
# ============================================================

heart_concat_nmf <- heart_models$concat

h_spliced_indices <- grep(
  "_s$",
  rownames(
    heart_concat_nmf$w
  )
)

h_unspliced_indices <- grep(
  "_us$",
  rownames(
    heart_concat_nmf$w
  )
)


# ------------------------------------------------------------
# Figure 5A
# ------------------------------------------------------------

heart_scores <- data.frame(
  nmf_factor =
    1:ncol(
      heart_concat_nmf$w
    ),
  spliced =
    colSums(
      heart_concat_nmf$w[
        h_spliced_indices,
      ]
    ),
  unspliced =
    colSums(
      heart_concat_nmf$w[
        h_unspliced_indices,
      ]
    )
)


heart_scores$spliced_fraction <-
  heart_scores$spliced /
  (
    heart_scores$spliced +
      heart_scores$unspliced
  )


p1 <- ggplot(
  heart_scores,
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
      margin(
        10,
        10,
        10,
        10
      )
  )


# ------------------------------------------------------------
# Figure 5B
# ------------------------------------------------------------

all_gsea <- list()

cat("Running GSEA:\n")


for (
  i in 1:ncol(
    heart_concat_nmf$w
  )
) {

  cat(i, " ")

  # Spliced weights
  s_weights <- heart_concat_nmf$w[
    h_spliced_indices,
    i
  ]

  names(s_weights) <- rownames(
    heart_concat_nmf$w[
      h_spliced_indices,
    ]
  )

  s_res <- run_gsea_fgsea(
    s_weights,
    mouse_gene_sets
  )


  # Unspliced weights
  us_weights <- heart_concat_nmf$w[
    h_unspliced_indices,
    i
  ]

  names(us_weights) <- rownames(
    heart_concat_nmf$w[
      h_unspliced_indices,
    ]
  )

  us_res <- run_gsea_fgsea(
    us_weights,
    mouse_gene_sets
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
  is.na(
    gsea_df$padj_s
  )
] <- 1

gsea_df$padj_us[
  is.na(
    gsea_df$padj_us
  )
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
      margin(
        10,
        10,
        10,
        10
      )
  )


figure5 <- cowplot::plot_grid(
  p1,
  p2,
  nrow = 1,
  labels = c("A", "B"),
  rel_widths = c(1, 2)
)
