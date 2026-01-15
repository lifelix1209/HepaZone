# ==============================================================================
# Spatial Zonation Core Functions
# ==============================================================================

#' Halpern 2017 Central Vein Marker Genes
#'
#' A dataset containing marker genes for central vein hepatocytes
#' from Bahar Halpern et al., Nature 2017.
#'
#' @format A character vector of gene symbols
#' @source Halpern & Shenhav et al., Nature 2017
"cv_markers"


#' Halpern 2017 Portal Marker Genes
#'
#' A dataset containing marker genes for portal hepatocytes
#' from Bahar Halpern et al., Nature 2017.
#'
#' @format A character vector of gene symbols
#' @source Halpern & Shenhav et al., Nature 2017
"pn_markers"


#' Default Central Vein Marker Genes - Mouse (Title Case)
#'
#' Default marker genes for central vein (CV) zone based on Halpern 2017.
#' These genes are highly expressed near the central vein (mouse format: Title Case).
#'
#' @keywords internal
.default_cv_markers_mouse <- function() {
  c("Cyp2e1", "Cyp1a2", "Cyp2f2", "Cyp2a5", "Hamp", "Hamp2",
    "Acsl1", "Cyp8b1", "Slc27a5", "Ugt2b36", "Cyp7a1", "Cyp3a11",
    "Gsta1", "Gsta2", "Gsta3", "Gsta4", "Gstm1", "Gstm2", "Gstm3",
    "Mup3", "Mup5", "Mup20")
}


#' Default Portal Marker Genes - Mouse (Title Case)
#'
#' Default marker genes for portal vein (PV) zone based on Halpern 2017.
#' These genes are highly expressed near the portal vein (mouse format: Title Case).
#'
#' @keywords internal
.default_pn_markers_mouse <- function() {
  c("Alb", "Apob", "Tf", "Fgg", "Fgb", "Fga", "Serpina1a", "Serpina1b",
    "Serpina1c", "Serpina1d", "Serpina1e", "Serpina6", "C3", "Ahsg",
    "Fetub", "Lum", "Spp2", "Igfbp1", "Igfbp2", "Rbp4", "Ttr")
}


#' Default Central Vein Marker Genes - Human (UPPERCASE)
#'
#' Default marker genes for central vein (CV) zone based on Halpern 2017.
#' These genes are highly expressed near the central vein (human format: UPPERCASE).
#'
#' @keywords internal
.default_cv_markers_human <- function() {
  c("CYP2E1", "CYP1A2", "CYP2F2", "CYP2A5", "HAMP", "HAMP2",
    "ACSL1", "CYP8B1", "SLC27A5", "UGT2B36", "CYP7A1", "CYP3A11",
    "GSTA1", "GSTA2", "GSTA3", "GSTA4", "GSTM1", "GSTM2", "GSTM3",
    "MUP3", "MUP5", "MUP20")
}


#' Default Portal Marker Genes - Human (UPPERCASE)
#'
#' Default marker genes for portal vein (PV) zone based on Halpern 2017.
#' These genes are highly expressed near the portal vein (human format: UPPERCASE).
#'
#' @keywords internal
.default_pn_markers_human <- function() {
  c("ALB", "APOB", "TF", "FGG", "FGB", "FGA", "SERPINA1", "SERPINA6",
    "C3", "AHSG", "FETUB", "LUM", "SPP2", "IGFBP1", "IGFBP2", "RBP4", "TTR")
}


#' Default Central Vein Marker Genes - Lowercase (lowercase format)
#'
#' Default marker genes for central vein (CV) zone based on Halpern 2017.
#' These genes are highly expressed near the central vein (lowercase format: all lowercase).
#'
#' @keywords internal
.default_cv_markers_lower <- function() {
  c("cyp2e1", "cyp1a2", "cyp2f2", "cyp2a5", "hamp", "hamp2",
    "acsl1", "cyp8b1", "slc27a5", "ugt2b36", "cyp7a1", "cyp3a11",
    "gsta1", "gsta2", "gsta3", "gsta4", "gstm1", "gstm2", "gstm3",
    "mup3", "mup5", "mup20")
}


#' Default Portal Marker Genes - Lowercase (lowercase format)
#'
#' Default marker genes for portal vein (PV) zone based on Halpern 2017.
#' These genes are highly expressed near the portal vein (lowercase format: all lowercase).
#'
#' @keywords internal
.default_pn_markers_lower <- function() {
  c("alb", "apob", "tf", "fgg", "fgb", "fga", "serpina1a", "serpina1b",
    "serpina1c", "serpina1d", "serpina1e", "serpina6", "c3", "ahsg",
    "fetub", "lum", "spp2", "igfbp1", "igfbp2", "rbp4", "ttr")
}


#' Get Default Marker Genes
#'
#' Returns appropriate default marker genes based on detected gene name format.
#' Automatically detects mouse (Title Case), human (UPPERCASE), or lowercase format.
#'
#' @param gene_names Character vector of gene names from the dataset
#' @return List with cv_markers and pn_markers
#' @keywords internal
.get_default_markers <- function(gene_names) {
  # Detect gene name format
  format_detected <- .detect_gene_name_format(gene_names)

  if (format_detected == "mouse") {
    list(
      cv_markers = .default_cv_markers_mouse(),
      pn_markers = .default_pn_markers_mouse(),
      species = "mouse"
    )
  } else if (format_detected == "human") {
    list(
      cv_markers = .default_cv_markers_human(),
      pn_markers = .default_pn_markers_human(),
      species = "human"
    )
  } else if (format_detected == "lowercase") {
    list(
      cv_markers = .default_cv_markers_lower(),
      pn_markers = .default_pn_markers_lower(),
      species = "lowercase"
    )
  } else {
    # Fallback: try lowercase format (common in some datasets)
    warning("Could not detect gene name format. Trying lowercase markers as fallback.")
    list(
      cv_markers = .default_cv_markers_lower(),
      pn_markers = .default_pn_markers_lower(),
      species = "unknown"
    )
  }
}


#' Detect Gene Name Format
#'
#' Detects whether gene names are in mouse (Title Case), human (UPPERCASE),
#' or lowercase format. Genes starting with numbers are excluded from detection.
#'
#' @param gene_names Character vector of gene names
#' @return "mouse", "human", "lowercase", or "unknown"
#' @keywords internal
.detect_gene_name_format <- function(gene_names) {
  if (length(gene_names) == 0) {
    return("unknown")
  }

  # Filter out genes starting with numbers (like 0610007p14rik)
  valid_genes <- gene_names[!grepl("^[0-9]", gene_names)]

  if (length(valid_genes) == 0) {
    return("unknown")
  }

  # Check a sample of genes
  sample_size <- min(100, length(valid_genes))
  sample_genes <- sample(valid_genes, sample_size)

  # Count patterns
  # Mouse: Title Case (first letter uppercase, rest lowercase, may have numbers)
  #         Example: Cyp2e1, Alb, Serpina1a
  # Human: All uppercase
  #         Example: CYP2E1, ALB, SERPINA1
  # Lowercase: All lowercase
  #         Example: cyp2e1, alb, serpina1a

  n_upper <- sum(grepl("^[A-Z]+$", sample_genes))
  n_title <- sum(grepl("^[A-Z][a-z]", sample_genes))
  n_lower <- sum(grepl("^[a-z]+$", sample_genes))

  n_total <- length(sample_genes)

  # Calculate ratios
  upper_ratio <- n_upper / n_total
  title_ratio <- n_title / n_total
  lower_ratio <- n_lower / n_total

  message(sprintf("Gene format detection: %d uppercase (%.1f%%), %d Title Case (%.1f%%), %d lowercase (%.1f%%)",
                  n_upper, 100 * upper_ratio, n_title, 100 * title_ratio, n_lower, 100 * lower_ratio))

  # Decision based on majority
  if (upper_ratio > 0.5) {
    return("human")
  } else if (title_ratio > 0.5) {
    return("mouse")
  } else if (lower_ratio > 0.5) {
    return("lowercase")
  } else {
    # Try to find known markers in different cases
    mouse_markers <- c("Cyp2e1", "Alb", "Gsta1")
    human_markers <- c("CYP2E1", "ALB", "GSTA1")
    lower_markers <- c("cyp2e1", "alb", "gsta1")

    mouse_found <- sum(mouse_markers %in% gene_names)
    human_found <- sum(human_markers %in% gene_names)
    lower_found <- sum(lower_markers %in% gene_names)

    if (human_found > mouse_found && human_found > lower_found) {
      return("human")
    } else if (mouse_found > human_found && mouse_found > lower_found) {
      return("mouse")
    } else if (lower_found > mouse_found && lower_found > human_found) {
      return("lowercase")
    } else {
      # Return lowercase if no clear match but some lowercase markers found
      if (lower_found > 0) {
        return("lowercase")
      }
      return("unknown")
    }
  }
}


#' Calculate Spatial Position Scores Using PCA-Based Multi-Gene Weighting
#'
#' Computes spatial position scores for each cell using CV and PN marker genes.
#' Uses Z-score normalization followed by PCA to extract the first principal
#' component as the spatial gradient, providing optimal multi-gene weighting.
#'
#' @param seurat_obj Preprocessed Seurat object with normalized data
#' @param cv_markers Character vector of central vein marker genes (optional)
#' @param pn_markers Character vector of portal marker genes (optional)
#' @param use_default_markers Use default marker sets if not provided (default: TRUE)
#' @param use_pca Use PCA for multi-gene weighting (default: TRUE)
#' @return Seurat object with 'pca_gradient', 'cv_score', 'pn_score', and 'CL_score' in meta.data
#' @examples
#' \dontrun{
#' seurat_obj <- preprocess_zonation(seurat_obj)
#' seurat_obj <- calculate_spatial_position(seurat_obj)
#' }
#' @export
calculate_spatial_position <- function(seurat_obj,
                                        cv_markers = NULL,
                                        pn_markers = NULL,
                                        use_default_markers = TRUE,
                                        use_pca = TRUE) {

  # Get normalized data matrix using helper
  mat_norm <- .get_data(seurat_obj)

  # Try to get gene names from misc slot (set during preprocessing)
  gene_names <- seurat_obj@misc$gene_names

  # If we have gene names from misc, ensure the data layer has them
  if (!is.null(gene_names) && length(gene_names) > 0) {
    current_names <- rownames(mat_norm)
    if (is.null(current_names) || length(current_names) == 0 ||
        all(nchar(current_names) == 0)) {
      # Set dimnames on the matrix
      dimnames(mat_norm) <- list(gene_names, colnames(mat_norm))
    }
  }

  # Use default markers if not provided (with auto-detection)
  if (is.null(cv_markers) || is.null(pn_markers)) {
    if (use_default_markers) {
      message("Detecting gene name format and selecting appropriate markers...")

      # Get gene names from multiple sources
      all_gene_names <- rownames(mat_norm)

      # Try to get gene names from misc slot if not available from matrix
      if (is.null(all_gene_names) || length(all_gene_names) == 0) {
        all_gene_names <- seurat_obj@misc$gene_names
      }

      # If still empty, try to get from Seurat object features
      if (is.null(all_gene_names) || length(all_gene_names) == 0) {
        assay_obj <- seurat_obj@assays[["RNA"]]
        if ("features" %in% slotNames(assay_obj)) {
          all_gene_names <- assay_obj@features
        }
      }

      # If we have gene names, try auto-detection
      if (length(all_gene_names) > 0) {
        # Auto-detect and get appropriate markers
        marker_info <- .get_default_markers(all_gene_names)

        if (is.null(cv_markers)) {
          cv_markers <- marker_info$cv_markers
        }
        if (is.null(pn_markers)) {
          pn_markers <- marker_info$pn_markers
        }

        message(sprintf("Detected %s gene name format, using %s markers",
                        marker_info$species, toupper(marker_info$species)))
      } else {
        # Fallback: try both mouse and human markers
        message("Could not detect gene names. Trying both mouse and human markers...")

        # Try mouse markers first
        mouse_cv <- .default_cv_markers_mouse()
        mouse_pn <- .default_pn_markers_mouse()
        mouse_cv_found <- sum(mouse_cv %in% rownames(mat_norm))
        mouse_pn_found <- sum(mouse_pn %in% rownames(mat_norm))

        # Try human markers
        human_cv <- .default_cv_markers_human()
        human_pn <- .default_pn_markers_human()
        human_cv_found <- sum(human_cv %in% rownames(mat_norm))
        human_pn_found <- sum(human_pn %in% rownames(mat_norm))

        if (mouse_cv_found + mouse_pn_found >= human_cv_found + human_pn_found) {
          message("Mouse markers matched better, using mouse format")
          if (is.null(cv_markers)) cv_markers <- mouse_cv
          if (is.null(pn_markers)) pn_markers <- mouse_pn
        } else {
          message("Human markers matched better, using human format")
          if (is.null(cv_markers)) cv_markers <- human_cv
          if (is.null(pn_markers)) pn_markers <- human_pn
        }
      }
    } else {
      stop("cv_markers and pn_markers must be provided if use_default_markers = FALSE")
    }
  }

  # Filter markers that exist in the dataset
  cv_markers <- cv_markers[cv_markers %in% rownames(mat_norm)]
  pn_markers <- pn_markers[pn_markers %in% rownames(mat_norm)]

  message(sprintf("Using %d CV markers and %d PN markers",
                  length(cv_markers), length(pn_markers)))

  if (length(cv_markers) == 0) {
    stop("No CV markers found in dataset")
  }
  if (length(pn_markers) == 0) {
    stop("No PN markers found in dataset")
  }

  # Convert to matrix if needed for subsetting
  mat_norm_mat <- as.matrix(mat_norm)

  # Calculate CV score (sum of CV marker expressions)
  cv_score <- numeric(ncol(mat_norm_mat))
  for (m in cv_markers) {
    if (m %in% rownames(mat_norm_mat)) {
      cv_score <- cv_score + mat_norm_mat[m, ]
    }
  }

  # Calculate PN score (sum of PN marker expressions)
  pn_score <- numeric(ncol(mat_norm_mat))
  for (m in pn_markers) {
    if (m %in% rownames(mat_norm_mat)) {
      pn_score <- pn_score + mat_norm_mat[m, ]
    }
  }

  if (use_pca) {
    # === PCA-BASED MULTI-GENE WEIGHTING ===
    message("Computing PCA-based spatial gradient from marker genes...")

    # Combine CV and PN markers
    all_markers <- c(cv_markers, pn_markers)
    marker_mat <- mat_norm_mat[all_markers, , drop = FALSE]  # genes x cells
    n_markers <- nrow(marker_mat)
    n_cells <- ncol(marker_mat)

    # Z-score normalize each marker gene across cells
    # (column-wise standardization for each gene)
    marker_zscore <- matrix(0, nrow = n_markers, ncol = n_cells)
    for (i in 1:n_markers) {
      x <- marker_mat[i, ]
      x_mean <- mean(x, na.rm = TRUE)
      x_sd <- sd(x, na.rm = TRUE)
      if (x_sd > 0) {
        marker_zscore[i, ] <- (x - x_mean) / x_sd
      } else {
        marker_zscore[i, ] <- x - x_mean
      }
    }
    rownames(marker_zscore) <- all_markers
    colnames(marker_zscore) <- colnames(mat_norm_mat)

    # Transpose: cells x genes for PCA
    marker_cell_mat <- t(marker_zscore)  # cells x genes

    # Run PCA
    pca_result <- stats::prcomp(marker_cell_mat, center = FALSE, scale. = FALSE)

    # Extract first principal component (PC1 captures the CV-PN gradient)
    pc1 <- pca_result$x[, 1]

    # Determine direction: PC1 should increase from CV to PN
    # Check correlation with PN score (sum of PN markers)
    corr_with_pn <- cor(pc1, pn_score, use = "complete.obs")
    if (corr_with_pn < 0) {
      # PC1 is negatively correlated with PN, flip it
      pc1 <- -pc1
      message("Flipping PC1 direction to align with PN markers")
    }

    # Map PC1 to [0, 1] range using rank-based transformation
    # This is more robust than min-max scaling
    pc1_ranks <- rank(pc1, ties.method = "average")
    pca_gradient <- (pc1_ranks - 1) / (length(pc1_ranks) - 1)

    # Store in metadata
    seurat_obj$pca_gradient <- pca_gradient
    seurat_obj$pc1 <- pc1
    seurat_obj$cv_score <- cv_score
    seurat_obj$pn_score <- pn_score
    seurat_obj$CL_score <- pca_gradient

    message(sprintf("PCA gradient range: %.3f - %.3f",
                    min(pca_gradient, na.rm = TRUE),
                    max(pca_gradient, na.rm = TRUE)))
    message(sprintf("PC1 explains %.1f%% of marker gene variance",
                    summary(pca_result)$importance[2, 1] * 100))

  } else {
    # === ORIGINAL RATIO-BASED METHOD ===
    message("Using ratio-based CL score (legacy method)...")

    # Calculate CL score: CL = PN / (CV + PN)
    # Range: 0 (CV) to 1 (PV)
    cl_score <- pn_score / (cv_score + pn_score)

    # Store in metadata
    seurat_obj$cv_score <- cv_score
    seurat_obj$pn_score <- pn_score
    seurat_obj$CL_score <- cl_score

    message(sprintf("CL score range: %.3f - %.3f",
                    min(cl_score, na.rm = TRUE),
                    max(cl_score, na.rm = TRUE)))
  }

  return(seurat_obj)
}


#' KNN-Based Spatial Smoothing of Gradient Scores
#'
#' Smooths spatial gradient scores using transcriptome similarity between cells.
#' Each cell's score is averaged with its k nearest neighbors weighted by similarity.
#'
#' @param seurat_obj Seurat object with gradient scores in meta.data
#' @param score_col Name of the score column to smooth (default: "CL_score")
#' @param k Number of nearest neighbors (default: 20)
#' @param use_pca Use PCA for dimensionality reduction before finding neighbors (default: TRUE)
#' @param n_pcs Number of PCs to use if use_pca = TRUE (default: 30)
#' @param weight_by_distance Weight neighbors by similarity distance (default: TRUE)
#' @return Seurat object with smoothed scores added to meta.data
#' @examples
#' \dontrun{
#' seurat_obj <- calculate_spatial_position(seurat_obj)
#' seurat_obj <- knn_smooth_scores(seurat_obj, k = 20)
#' }
#' @export
knn_smooth_scores <- function(seurat_obj,
                               score_col = "CL_score",
                               k = 20,
                               use_pca = TRUE,
                               n_pcs = 30,
                               weight_by_distance = TRUE) {

  message(sprintf("KNN smoothing with k=%d neighbors...", k))

  # Check if score column exists
  if (!score_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Score column '%s' not found in meta.data", score_col))
  }

  scores <- seurat_obj@meta.data[[score_col]]
  n_cells <- length(scores)

  if (n_cells < k + 1) {
    warning(sprintf("Number of cells (%d) < k+1 (%d), using k=%d",
                    n_cells, k + 1, n_cells - 1))
    k <- max(1, n_cells - 1)
  }

  # Get normalized expression matrix
  mat_norm <- .get_data(seurat_obj)

  # Determine number of PCs to use
  max_pcs <- min(n_cells - 1, nrow(mat_norm) - 1)
  n_pcs_use <- min(n_pcs, max_pcs)
  if (n_pcs_use < 10) {
    warning(sprintf("Low number of PCs (%d), using all available", n_pcs_use))
  }

  if (use_pca) {
    # Run PCA on the expression matrix
    message(sprintf("Running PCA with %d components...", n_pcs_use))

    # Use Seurat's PCA if available
    if ("pca" %in% Seurat::Reductions(seurat_obj)) {
      pca_coords <- Seurat::Embeddings(seurat_obj, reduction = "pca")
      # Select top n_pcs_use dimensions
      if (ncol(pca_coords) >= n_pcs_use) {
        pca_coords <- pca_coords[, 1:n_pcs_use, drop = FALSE]
      } else {
        n_pcs_use <- ncol(pca_coords)
      }
    } else {
      # Compute PCA manually
      mat_t <- t(as.matrix(mat_norm))
      pca_result <- stats::prcomp(mat_t, rank = n_pcs_use, center = TRUE, scale. = TRUE)
      pca_coords <- pca_result$x
    }
  } else {
    # Use raw expression for similarity
    message("Using raw expression for neighbor finding...")
    mat_t <- t(as.matrix(mat_norm))
    pca_coords <- mat_t
  }

  # Calculate cell-cell distance matrix (Euclidean in PCA space)
  message("Computing cell-cell distances...")
  cell_dist <- as.matrix(dist(pca_coords))

  # Find k nearest neighbors for each cell (excluding self)
  message("Finding k-nearest neighbors...")
  smoothed_scores <- numeric(n_cells)

  for (i in 1:n_cells) {
    # Get distances to all other cells
    dists <- cell_dist[i, ]

    # Find k nearest neighbors (excluding self)
    neighbor_order <- order(dists)[2:(k + 1)]  # Skip self (index 1)
    neighbor_dists <- dists[neighbor_order]

    if (weight_by_distance) {
      # Convert distances to weights using Gaussian kernel
      # Smaller distance = higher weight
      sigma <- median(neighbor_dists)  # Use median as bandwidth
      if (sigma == 0) sigma <- 1

      weights <- exp(-neighbor_dists^2 / (2 * sigma^2))
      weights <- weights / sum(weights)  # Normalize

      # Weighted average of self and neighbors
      self_weight <- 0.5  # Weight for original cell
      neighbor_weight <- (1 - self_weight) / k

      smoothed_scores[i] <- self_weight * scores[i] +
                           neighbor_weight * sum(scores[neighbor_order])
    } else {
      # Simple average of self and neighbors
      smoothed_scores[i] <- mean(c(scores[i], scores[neighbor_order]))
    }
  }

  # Map to [0, 1] range using rank-based transformation (more robust)
  smoothed_ranks <- rank(smoothed_scores, ties.method = "average")
  smoothed_scores_scaled <- (smoothed_ranks - 1) / (length(smoothed_ranks) - 1)

  # Add smoothed scores to metadata
  new_col_name <- paste0("smoothed_", score_col)
  seurat_obj@meta.data[[new_col_name]] <- smoothed_scores_scaled

  # Update CL_score to use smoothed version if it was the original score
  if (score_col == "CL_score") {
    seurat_obj$CL_score <- smoothed_scores_scaled
    message("Updated CL_score with smoothed values")
  }

  message(sprintf("Smoothed score range: %.3f - %.3f",
                  min(smoothed_scores_scaled, na.rm = TRUE),
                  max(smoothed_scores_scaled, na.rm = TRUE)))

  return(seurat_obj)
}


#' Generate Gamma Distribution Parameters for Spatial Layers
#'
#' Creates Gamma distribution parameters for each spatial layer.
#' The parameters are designed so that each layer has a distinct mode
#' along the CL score axis (0 to 1).
#'
#' @param n_zones Number of spatial zones (default: 10)
#' @param gamma_shape Shape parameter for Gamma distribution (default: 5)
#' @return A data frame with zone, shape, rate, and mode for each zone
#' @keywords internal
.generate_gamma_params <- function(n_zones = 10, gamma_shape = 5) {
  # Generate modes uniformly distributed from 0.05 to 0.95
  modes <- seq(0.05, 0.95, length.out = n_zones)

  # For Gamma distribution: mode = (shape - 1) / rate, for shape > 1
  # Set rate to get desired mode
  params <- data.frame(
    zone = 1:n_zones,
    shape = rep(gamma_shape, n_zones),
    stringsAsFactors = FALSE
  )
  params$mode <- modes
  params$rate <- (params$shape - 1) / modes
  params$rate[params$shape <= 1] <- 1  # Handle edge case

  return(params)
}


#' Map Cells to Spatial Layers Using Gamma Distribution
#'
#' Assigns each cell to spatial layers based on their CL score using
#' a mixture of Gamma distributions. Each zone has a pre-defined Gamma
#' distribution, and cells are probabilistically assigned to zones.
#'
#' @param seurat_obj Seurat object with CL_score in meta.data
#' @param n_zones Number of spatial zones (default: 10)
#' @param gamma_shape Shape parameter for Gamma distributions (default: 5)
#' @param return_matrix If TRUE, return probability matrix; if FALSE, add to Seurat (default: TRUE)
#' @return If return_matrix = TRUE, returns probability matrix (cells x zones).
#'         Otherwise, returns modified Seurat object.
#' @examples
#' \dontrun{
#' prob_matrix <- map_cells_to_layers(seurat_obj, n_zones = 10)
#' }
#' @export
map_cells_to_layers <- function(seurat_obj,
                                 n_zones = 10,
                                 gamma_shape = 5,
                                 return_matrix = TRUE) {

  # Get CL scores
  cl_scores <- seurat_obj$CL_score
  n_cells <- length(cl_scores)

  if (is.null(cl_scores)) {
    stop("CL_score not found. Run calculate_spatial_position() first.")
  }

  # Generate Gamma distribution parameters
  gamma_params <- .generate_gamma_params(n_zones, gamma_shape)

  message(sprintf("Gamma distribution modes: %s",
                  paste(round(gamma_params$mode, 2), collapse = ", ")))

  # Calculate probability density for each cell at each zone's Gamma mode
  prob_matrix <- matrix(0, nrow = n_cells, ncol = n_zones)

  # Scale CL scores to 0-1 range if needed (Gamma needs positive values)
  cl_scaled <- pmax(pmin(cl_scores, 0.999), 0.001)

  for (z in 1:n_zones) {
    shape_z <- gamma_params$shape[z]
    rate_z <- gamma_params$rate[z]

    # Gamma distribution density
    prob_matrix[, z] <- stats::dgamma(cl_scaled, shape = shape_z, rate = rate_z)
  }

  # Normalize to get probabilities (row-wise softmax)
  row_sums <- Matrix::rowSums(prob_matrix)
  row_sums[row_sums == 0] <- 1  # Avoid division by zero
  prob_matrix <- prob_matrix / row_sums

  # Add zone names
  colnames(prob_matrix) <- paste0("Zone_", 1:n_zones)
  rownames(prob_matrix) <- colnames(seurat_obj)

  if (return_matrix) {
    return(prob_matrix)
  } else {
    # Add to Seurat object
    prob_df <- as.data.frame(prob_matrix)
    seurat_obj <- AddMetaData(seurat_obj, prob_df)

    # Also add dominant zone assignment
    dominant_zone <- apply(prob_matrix, 1, which.max)
    seurat_obj$dominant_zone <- dominant_zone

    return(seurat_obj)
  }
}


#' Reconstruct Spatial Gene Expression Profiles
#'
#' Computes the average expression of each gene at each spatial zone
#' using the cell-to-zone probability matrix.
#'
#' @param seurat_obj Preprocessed Seurat object
#' @param prob_matrix Cell-to-zone probability matrix (cells x zones)
#' @return A list containing:
#'   - mean_expression: Gene x Zone matrix of mean expression
#'   - n_cells: Number of cells contributing to each zone
#' @examples
#' \dontrun{
#' prob_matrix <- map_cells_to_layers(seurat_obj)
#' result <- reconstruct_spatial_expression(seurat_obj, prob_matrix)
#' }
#' @export
reconstruct_spatial_expression <- function(seurat_obj, prob_matrix) {

  # Get normalized data matrix using helper
  mat_norm <- .get_data(seurat_obj)

  # Transpose to cells x genes for matrix multiplication
  mat_t <- t(mat_norm)  # cells x genes

  # Calculate mean expression per zone
  # MeanGeneExp = Pmat * mat_norm (genes x zones)
  # Where Pmat is the transpose of prob_matrix (zones x cells)
  mean_expr <- mat_norm %*% prob_matrix  # genes x zones

  # Number of cells contributing to each zone
  zone_contributions <- Matrix::colSums(prob_matrix)

  # Store results
  result <- list(
    mean_expression = as.matrix(mean_expr),
    prob_matrix = prob_matrix,
    n_cells = zone_contributions,
    n_zones = ncol(prob_matrix)
  )

  return(result)
}


#' Create HepaZone Result Object
#'
#' Internal function to compile full results from zonation analysis.
#'
#' @param seurat_obj Original Seurat object
#' @param spatial_result Result from reconstruct_spatial_expression()
#' @param cv_markers CV marker genes used
#' @param pn_markers PN marker genes used
#' @return HepaZone result S3 object
#' @keywords internal
.create_hepa_zone_result <- function(seurat_obj, spatial_result,
                                      cv_markers, pn_markers) {
  result <- list(
    mean_expression = spatial_result$mean_expression,
    prob_matrix = spatial_result$prob_matrix,
    n_cells = spatial_result$n_cells,
    cv_markers = cv_markers,
    pn_markers = pn_markers,
    cl_scores = seurat_obj$CL_score,
    cells = colnames(seurat_obj),
    genes = rownames(seurat_obj),
    n_zones = spatial_result$n_zones,
    preprocessing = seurat_obj@misc
  )

  class(result) <- "HepaZoneResult"
  return(result)
}


#' Print HepaZoneResult Object
#'
#' @param x HepaZoneResult object
#' @param ... Additional arguments (unused)
#' @method print HepaZoneResult
#' @export
print.HepaZoneResult <- function(x, ...) {
  cat("HepaZone Result Object\n")
  cat("======================\n")
  cat(sprintf("Genes: %d\n", nrow(x$mean_expression)))
  cat(sprintf("Cells: %d\n", ncol(x$prob_matrix)))
  cat(sprintf("Zones: %d\n", x$n_zones))

  # Show species if available
  if ("species" %in% names(x)) {
    cat(sprintf("Species: %s\n", x$species))
  }

  cat(sprintf("CV markers: %d\n", length(x$cv_markers)))
  cat(sprintf("PN markers: %d\n", length(x$pn_markers)))

  # Show method parameters
  if ("method_params" %in% names(x)) {
    mp <- x$method_params
    cat(sprintf("\nMethod: PCA=%s, KNN_smooth=%s (k=%d)\n",
                mp$use_pca, mp$knn_smooth, mp$k_neighbors))
  }

  cat("\nExpression matrix dimensions: ",
      paste(dim(x$mean_expression), collapse = " x "), "\n")
}
