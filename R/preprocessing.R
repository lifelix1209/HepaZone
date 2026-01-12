# ==============================================================================
# Data Preprocessing Functions
# ==============================================================================

#' Preprocess scRNA-seq Data for Zonation Analysis
#'
#' Performs preprocessing steps specific for liver zonation analysis:
#' 1. Removes mitochondrial genes (Mt-*)
#' 2. Removes Mup genes (major urinary proteins)
#' 3. Cell-level normalization (counts per cell)
#' 4. Gene-level normalization (relative to max expression)
#'
#' @param seurat_obj A Seurat object
#' @param mt_pattern Pattern to identify mitochondrial genes (default: "^Mt-")
#' @param remove_mup Remove Mup genes (default: TRUE)
#' @param mup_pattern Pattern to identify Mup genes (default: "^Mup")
#' @param do_normalize Perform normalization (default: TRUE)
#' @return Preprocessed Seurat object with normalized data in 'data' and 'scale.data' slots
#' @examples
#' \dontrun{
#' seurat_obj <- load_10x_data("./data/ZT00")
#' seurat_obj <- preprocess_zonation(seurat_obj)
#' }
#' @export
preprocess_zonation <- function(seurat_obj,
                                 mt_pattern = "^Mt-",
                                 remove_mup = TRUE,
                                 mup_pattern = "^Mup",
                                 do_normalize = TRUE) {

  # Get the counts matrix
  counts <- seurat_obj@assays$RNA@counts

  # Step 1: Identify and remove mitochondrial genes
  mt_genes <- grep(mt_pattern, rownames(counts), value = TRUE)
  if (length(mt_genes) > 0) {
    message(sprintf("Removing %d mitochondrial genes", length(mt_genes)))
    counts <- counts[!rownames(counts) %in% mt_genes, ]
  }

  # Step 2: Remove Mup genes if requested
  if (remove_mup) {
    mup_genes <- grep(mup_pattern, rownames(counts), value = TRUE)
    if (length(mup_genes) > 0) {
      message(sprintf("Removing %d Mup genes", length(mup_genes)))
      counts <- counts[!rownames(counts) %in% mup_genes, ]
    }
  }

  # Update the Seurat object with filtered counts
  seurat_obj@assays$RNA@counts <- counts

  # Step 3: Cell-level normalization
  # Each cell's expression divided by total expression in that cell
  if (do_normalize) {
    message("Performing cell-level normalization...")

    # Calculate cell totals (colSums)
    cell_totals <- Matrix::colSums(counts)

    # Avoid division by zero
    cell_totals[cell_totals == 0] <- 1

    # Normalize: counts / cell_total
    mat_norm <- t(t(counts) / cell_totals)

    # Store in 'data' slot (normalized but not scaled)
    seurat_obj@assays$RNA@data <- as(mat_norm, "dgCMatrix")

    # Step 4: Gene-level normalization
    # Each gene divided by its maximum expression
    message("Performing gene-level normalization...")

    gene_max <- Matrix::rowMax(mat_norm)
    gene_max[gene_max == 0] <- 1  # Avoid division by zero

    # Normalize: gene_expression / gene_max
    mat_norm_genes <- mat_norm / gene_max

    # Update 'data' slot with double-normalized data
    seurat_obj@assays$RNA@data <- as(mat_norm_genes, "dgCMatrix")

    # Store the original normalized matrix for later use
    if (!"mat_norm" %in% names(seurat_obj@misc)) {
      seurat_obj@misc$mat_norm <- mat_norm_genes
    }
  }

  # Update metadata
  n_genes_after <- nrow(counts)
  n_cells_after <- ncol(counts)
  message(sprintf("Preprocessing complete: %d genes x %d cells",
                  n_genes_after, n_cells_after))

  return(seurat_obj)
}


#' Filter Cells by Mitochondrial Content
#'
#' Identifies and optionally removes cells with high mitochondrial gene expression,
#' which often indicates low-quality or dying cells.
#'
#' @param seurat_obj A Seurat object
#' @param mt_pattern Pattern for mitochondrial genes (default: "^Mt-")
#' @param max_mt_percent Maximum allowed mitochondrial percentage (default: 20)
#' @param remove If TRUE, remove high-mt cells; if FALSE, just add metadata (default: TRUE)
#' @return Seurat object with mitochondrial percentage in meta.data
#' @export
filter_by_mitochondrial <- function(seurat_obj,
                                     mt_pattern = "^Mt-",
                                     max_mt_percent = 20,
                                     remove = TRUE) {
  # Calculate mitochondrial percentage
  counts <- seurat_obj@assays$RNA@counts
  mt_genes <- grep(mt_pattern, rownames(counts), value = TRUE)

  if (length(mt_genes) == 0) {
    message("No mitochondrial genes found with pattern: ", mt_pattern)
    return(seurat_obj)
  }

  mt_counts <- Matrix::colSums(counts[mt_genes, , drop = FALSE])
  total_counts <- Matrix::colSums(counts)

  mt_percent <- (mt_counts / total_counts) * 100

  # Add to metadata
  seurat_obj$mt_percent <- mt_percent

  # Filter if requested
  if (remove) {
    cells_to_keep <- which(mt_percent < max_mt_percent)
    message(sprintf("Removing %d cells with mt > %d%% (keeping %d of %d)",
                    ncol(seurat_obj) - length(cells_to_keep),
                    max_mt_percent,
                    length(cells_to_keep),
                    ncol(seurat_obj)))

    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
  }

  return(seurat_obj)
}


#' Filter Low-Complexity Genes
#'
#' Removes genes with very low expression across cells.
#'
#' @param seurat_obj A Seurat object
#' @param min_cells Express in at least this many cells (default: 3)
#' @param min_counts Minimum total counts across all cells (default: 100)
#' @return Filtered Seurat object
#' @export
filter_low_complexity <- function(seurat_obj,
                                   min_cells = 3,
                                   min_counts = 100) {
  counts <- seurat_obj@assays$RNA@counts

  # Gene expressed in at least min_cells
  expr_per_gene <- Matrix::rowSums(counts > 0)
  genes_by_cells <- expr_per_gene >= min_cells

  # Gene has at least min_counts total
  total_counts <- Matrix::rowSums(counts)
  genes_by_counts <- total_counts >= min_counts

  # Keep genes that pass both filters
  genes_to_keep <- genes_by_cells & genes_by_counts

  n_removed <- sum(!genes_to_keep)
  message(sprintf("Removing %d low-complexity genes", n_removed))

  seurat_obj@assays$RNA@counts <- counts[genes_to_keep, ]

  return(seurat_obj)
}


#' Get Preprocessing Summary
#'
#' Generates a summary of the preprocessing steps applied to a Seurat object.
#'
#' @param seurat_obj A Seurat object after preprocessing
#' @return A data frame with preprocessing statistics
#' @export
get_preprocessing_summary <- function(seurat_obj) {
  summary_df <- data.frame(
    Metric = c(
      "Number of genes",
      "Number of cells",
      "Total counts",
      "Median counts per cell",
      "Mean counts per cell",
      "Median genes per cell",
      "Mean genes per cell"
    ),
    Value = c(
      nrow(seurat_obj),
      ncol(seurat_obj),
      sum(seurat_obj@assays$RNA@counts),
      median(Matrix::colSums(seurat_obj@assays$RNA@counts)),
      mean(Matrix::colSums(seurat_obj@assays$RNA@counts)),
      median(Matrix::colSums(seurat_obj@assays$RNA@counts > 0)),
      mean(Matrix::colSums(seurat_obj@assays$RNA@counts > 0))
    )
  )

  return(summary_df)
}
