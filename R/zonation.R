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


#' Default Central Vein Marker Genes
#'
#' Default marker genes for central vein (CV) zone based on Halpern 2017.
#' These genes are highly expressed near the central vein.
#'
#' @keywords internal
.default_cv_markers <- function() {
  c("Cyp2e1", "Cyp1a2", "Cyp2f2", "Cyp2a5", "Hamp", "Hamp2",
    "Acsl1", "Cyp8b1", "Slc27a5", "Ugt2b36", "Cyp7a1", "Cyp3a11",
    "Gsta1", "Gsta2", "Gsta3", "Gsta4", "Gstm1", "Gstm2", "Gstm3",
    "Mup3", "Mup5", "Mup20")
}


#' Default Portal Marker Genes
#'
#' Default marker genes for portal vein (PV) zone based on Halpern 2017.
#' These genes are highly expressed near the portal vein.
#'
#' @keywords internal
.default_pn_markers <- function() {
  c("Alb", "Apob", "Tf", "Fgg", "Fgb", "Fga", "Serpina1a", "Serpina1b",
    "Serpina1c", "Serpina1d", "Serpina1e", "Serpina6", "C3", "Ahsg",
    "Fetub", "Lum", "Spp2", "Igfbp1", "Igfbp2", "Rbp4", "Ttr")
}


#' Calculate Spatial Position Scores
#'
#' Computes spatial position scores for each cell using CV and PN marker genes.
#' The CL score (central-portal ratio) indicates the position along the
#' liver lobule axis: low values = near central vein, high values = near portal vein.
#'
#' @param seurat_obj Preprocessed Seurat object with normalized data
#' @param cv_markers Character vector of central vein marker genes (optional)
#' @param pn_markers Character vector of portal marker genes (optional)
#' @param use_default_markers Use default marker sets if not provided (default: TRUE)
#' @return Seurat object with 'cv_score', 'pn_score', and 'CL_score' in meta.data
#' @examples
#' \dontrun{
#' seurat_obj <- preprocess_zonation(seurat_obj)
#' seurat_obj <- calculate_spatial_position(seurat_obj)
#' }
#' @export
calculate_spatial_position <- function(seurat_obj,
                                        cv_markers = NULL,
                                        pn_markers = NULL,
                                        use_default_markers = TRUE) {

  # Get normalized data matrix
  mat_norm <- seurat_obj@assays$RNA@data

  # Use default markers if not provided
  if (is.null(cv_markers)) {
    if (use_default_markers) {
      cv_markers <- .default_cv_markers()
    } else {
      stop("cv_markers must be provided if use_default_markers = FALSE")
    }
  }

  if (is.null(pn_markers)) {
    if (use_default_markers) {
      pn_markers <- .default_pn_markers()
    } else {
      stop("pn_markers must be provided if use_default_markers = FALSE")
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

  # Calculate CV score (sum of CV marker expressions)
  cv_expr <- mat_norm[cv_markers, , drop = FALSE]
  cv_score <- Matrix::colSums(cv_expr)

  # Calculate PN score (sum of PN marker expressions)
  pn_expr <- mat_norm[pn_markers, , drop = FALSE]
  pn_score <- Matrix::colSums(pn_expr)

  # Calculate CL score: CL = PN / (CV + PN)
  # Range: 0 (CV) to 1 (PV)
  cl_score <- pn_score / (cv_score + pn_score)

  # Store in metadata
  seurat_obj$cv_score <- as.numeric(cv_score)
  seurat_obj$pn_score <- as.numeric(pn_score)
  seurat_obj$CL_score <- as.numeric(cl_score)

  message(sprintf("CL score range: %.3f - %.3f",
                  min(cl_score, na.rm = TRUE),
                  max(cl_score, na.rm = TRUE)))

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

  # Get normalized data matrix (genes x cells)
  mat_norm <- seurat_obj@assays$RNA@data

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
  cat(sprintf("CV markers: %d\n", length(x$cv_markers)))
  cat(sprintf("PN markers: %d\n", length(x$pn_markers)))
  cat("\nExpression matrix dimensions: ",
      paste(dim(x$mean_expression), collapse = " x "), "\n")
}
