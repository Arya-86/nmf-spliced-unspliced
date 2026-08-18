prepare_dataset <- function(
    data_dir,
    dataset_prefix,
    min_features = 200,
    max_features = 15000,
    min_umis = 500,
    max_umis = 75000
) {

  # Load matrices and metadata
  spliced_counts <- read_npz_into_R(
    file.path(data_dir, paste0(dataset_prefix, "_spliced.npz"))
  )

  unspliced_counts <- read_npz_into_R(
    file.path(data_dir, paste0(dataset_prefix, "_unspliced.npz"))
  )

  gene_data <- read.csv(
    file.path(data_dir, paste0(dataset_prefix, "_genes.tsv")),
    sep = "\t",
    header = FALSE
  )

  cell_data <- read.csv(
    file.path(data_dir, paste0(dataset_prefix, "_barcodes.tsv")),
    sep = "\t",
    header = FALSE
  )

  gene_names <- as.character(gene_data[, 1])

  colnames(spliced_counts) <-
    colnames(unspliced_counts) <-
    as.character(cell_data[, 1])

  rownames(spliced_counts) <- gene_names
  rownames(unspliced_counts) <- gene_names


  # ----------------------------------------------------------
  # Added representation
  # ----------------------------------------------------------

  combined_counts <- spliced_counts + unspliced_counts

  keep_cells <- which(
    colSums(combined_counts != 0) > min_features &
    colSums(combined_counts != 0) < max_features &
    colSums(combined_counts) > min_umis &
    colSums(combined_counts) < max_umis
  )

  combined_counts <- combined_counts[, keep_cells]

  # Normalize exactly as in the original pipeline
  added <- Seurat::LogNormalize(combined_counts)

  # Remove duplicated gene names AFTER normalization,
  # matching the original analysis.
  num_dups <- sum(duplicated(rownames(added)))

  added <- added[
    !duplicated(rownames(added)),
  ]

  stopifnot(
    length(unique(rownames(added))) == nrow(added)
  )

  message(
    "Removed ",
    num_dups,
    " duplicated gene names from Added representation"
  )


  # ----------------------------------------------------------
  # Concatenated representation
  # ----------------------------------------------------------

  # Preserve the original feature rows.
  # Do NOT deduplicate before constructing the concat matrix.
  rownames(spliced_counts) <-
    paste0(gene_names, "_s")

  rownames(unspliced_counts) <-
    paste0(gene_names, "_us")

  concat <- rbind(
    spliced_counts[, keep_cells],
    unspliced_counts[, keep_cells]
  )

  spliced_indices <- seq_len(nrow(spliced_counts))

  unspliced_indices <-
    (nrow(spliced_counts) + 1):nrow(concat)

  # Normalize exactly as in the original pipeline
  concat <- Seurat::LogNormalize(concat)


  # ----------------------------------------------------------
  # Return processed objects
  # ----------------------------------------------------------

  list(
    added = added,
    concat = concat,
    spliced_indices = spliced_indices,
    unspliced_indices = unspliced_indices,
    keep_cells = keep_cells
  )
}
