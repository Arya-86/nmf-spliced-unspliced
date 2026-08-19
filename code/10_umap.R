# ============================================================
# 10_umap.R
#
# UMAP comparison across NMF representations.
# ============================================================

library(Seurat)
library(dplyr)
library(uwot)
library(ggplot2)


run_umap <- function(
    nmf_model,
    modality
) {

  set.seed(123)

  umap_data <- uwot::umap(
    t(nmf_model$h),
    metric = "cosine",
    min_dist = 0.3
  )

  umap_data <- as.data.frame(
    umap_data
  )

  colnames(umap_data) <- c(
    "UMAP1",
    "UMAP2"
  )

  umap_data$Modality <-
    modality

  umap_data
}


compare_nmf_umaps <- function(
    dataset,
    models
) {

  added_matched <- dataset$added[
    ,
    colnames(models$added$h)
  ]

  stopifnot(
    identical(
      colnames(added_matched),
      colnames(models$added$h)
    )
  )

  rownames(added_matched) <-
    make.unique(
      rownames(added_matched)
    )

  sc_data <- CreateSeuratObject(
    counts = added_matched
  )

  sc_data[["nmf"]] <-
    CreateDimReducObject(
      embeddings =
        t(models$added$h),
      loadings =
        models$added$w,
      key = "NMF_",
      assay =
        DefaultAssay(sc_data)
    )

  sc_data <- FindNeighbors(
    sc_data,
    reduction = "nmf",
    dims =
      seq_len(
        ncol(
          sc_data[
            ["nmf"]
          ]@cell.embeddings
        )
      )
  ) %>%
    FindClusters()

  clusters <- data.frame(
    sample =
      colnames(
        models$added$h
      ),
    cluster =
      as.integer(
        sc_data$seurat_clusters
      )
  )

  plot_data <- rbind(
    cbind(
      clusters,
      run_umap(
        models$spliced,
        "Spliced"
      )
    ),
    cbind(
      clusters,
      run_umap(
        models$unspliced,
        "Unspliced"
      )
    ),
    cbind(
      clusters,
      run_umap(
        models$added,
        "Added"
      )
    ),
    cbind(
      clusters,
      run_umap(
        models$concat,
        "Concatenated"
      )
    )
  )

  plot_data$cluster <-
    factor(
      plot_data$cluster
    )

  plot_data$Modality <-
    factor(
      plot_data$Modality,
      levels = c(
        "Spliced",
        "Unspliced",
        "Added",
        "Concatenated"
      )
    )

  figure <- ggplot(
    plot_data,
    aes(
      x = UMAP1,
      y = UMAP2,
      color = cluster
    )
  ) +
    geom_point(
      size = 0.25,
      alpha = 0.8
    ) +
    facet_wrap(
      ~ Modality,
      ncol = 2,
      scales = "free"
    ) +
    scale_color_brewer(
      palette = "Set1",
      name = "Cluster"
    ) +
    theme_void(
      base_size = 12
    ) +
    theme(
      strip.text =
        element_text(
          face = "bold",
          size = 12
        ),
      legend.position =
        "bottom",
      legend.box =
        "horizontal"
    )

  list(
    clusters = clusters,
    plot_data = plot_data,
    figure = figure
  )
}
