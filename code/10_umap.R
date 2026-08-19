# ============================================================
# 10_umap.R
#
# UMAP comparison across NMF representations for Figure S4.
# ============================================================

library(Seurat)
library(dplyr)
library(uwot)
library(ggplot2)


# Final Heart models
heart_spliced_nmf <- heart_models$spliced
heart_unspliced_nmf <- heart_models$unspliced
heart_added_nmf <- heart_models$added
heart_concat_nmf <- heart_models$concat


# ------------------------------------------------------------
# Cluster cells using the Added-model NMF representation
# ------------------------------------------------------------

heart_added_matched <- heart$added[
  ,
  colnames(heart_added_nmf$h)
]

stopifnot(
  identical(
    colnames(heart_added_matched),
    colnames(heart_added_nmf$h)
  )
)

rownames(heart_added_matched) <-
  make.unique(
    rownames(heart_added_matched)
  )

sc_data <- CreateSeuratObject(
  counts = heart_added_matched
)

sc_data[["nmf"]] <- CreateDimReducObject(
  embeddings = t(heart_added_nmf$h),
  loadings = heart_added_nmf$w,
  key = "NMF_",
  assay = DefaultAssay(sc_data)
)

sc_data <- FindNeighbors(
  sc_data,
  reduction = "nmf",
  dims = 1:ncol(
    sc_data[["nmf"]]@cell.embeddings
  )
) %>%
  FindClusters()

clusters <- data.frame(
  sample = colnames(
    heart_added_nmf$h
  ),
  cluster = as.integer(
    sc_data$seurat_clusters
  )
)


# ------------------------------------------------------------
# Generate UMAP coordinates
# ------------------------------------------------------------

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

  umap_data$Modality <- modality

  umap_data
}


# Transfer the Added-model cluster labels
# to all four embeddings
plot_data <- rbind(
  cbind(
    clusters,
    run_umap(
      heart_spliced_nmf,
      "Spliced"
    )
  ),
  cbind(
    clusters,
    run_umap(
      heart_unspliced_nmf,
      "Unspliced"
    )
  ),
  cbind(
    clusters,
    run_umap(
      heart_added_nmf,
      "Added"
    )
  ),
  cbind(
    clusters,
    run_umap(
      heart_concat_nmf,
      "Concatenated"
    )
  )
)

plot_data$cluster <- factor(
  plot_data$cluster
)

plot_data$Modality <- factor(
  plot_data$Modality,
  levels = c(
    "Spliced",
    "Unspliced",
    "Added",
    "Concatenated"
  )
)


# ------------------------------------------------------------
# Figure S4
# ------------------------------------------------------------

figure_s4 <- ggplot(
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
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    legend.position = "bottom",
    legend.box = "horizontal"
  )
