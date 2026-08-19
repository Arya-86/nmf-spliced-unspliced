# ============================================================
# 05_gsea.R
#
# Gene set enrichment analysis of NMF factor loadings.
# ============================================================

library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)


# Build Gene Ontology (C5) gene sets
mouse_gene_sets <- msigdbr(
  species = "Mus musculus",
  collection = "C5"
) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::distinct() %>%
  dplyr::group_by(gs_name) %>%
  dplyr::summarise(
    genes = list(gene_symbol),
    .groups = "drop"
  ) %>%
  tibble::deframe()


# Run GSEA for one factor
run_gsea_fgsea <- function(
    gene_weights,
    gene_sets,
    top_fraction = 0.2
) {

  names(gene_weights) <- sub(
    "_(s|us)$",
    "",
    names(gene_weights)
  )

  gene_weights <- gene_weights[
    is.finite(gene_weights) &
      gene_weights > 0
  ]

  if (length(gene_weights) < 5) {
    return(NULL)
  }

  # Collapse duplicated gene symbols by maximum loading
  gene_weights <- tapply(
    gene_weights,
    names(gene_weights),
    max
  )

  gene_weights <- setNames(
    as.numeric(gene_weights),
    names(gene_weights)
  )

  gene_weights <- sort(
    gene_weights,
    decreasing = TRUE
  )

  n_keep <- max(
    1L,
    ceiling(length(gene_weights) * top_fraction)
  )

  gene_weights <- gene_weights[
    seq_len(n_keep)
  ]

  if (length(gene_weights) < 20) {
    return(NULL)
  }

  fgsea_res <- fgsea::fgseaMultilevel(
    pathways = gene_sets,
    stats = gene_weights,
    minSize = 5,
    maxSize = 250,
    eps = 0,
    scoreType = "pos"
  )

  if (nrow(fgsea_res) == 0) {
    return(NULL)
  }

  data.frame(
    pathway = fgsea_res$pathway,
    padj = fgsea_res$padj,
    NES = fgsea_res$NES,
    size = fgsea_res$size
  )
}
