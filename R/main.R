# ==============================================================================
# Main Wrapper Function
# ==============================================================================

#' HepaZone: Liver Spatial Transcriptomics Reconstruction
#'
#' The main function that orchestrates the complete spatial transcriptomics
#' reconstruction pipeline from preprocessing to statistical analysis.
#'
#' @param seurat_obj Input Seurat object (or path to 10x data if from_10x = TRUE)
#' @param cv_markers Character vector of central vein marker genes (optional)
#' @param pn_markers Character vector of portal marker genes (optional)
#' @param n_zones Number of spatial zones (default: 10)
#' @param n_bootstrap Number of bootstrap iterations for SE estimation (default: 500)
#' @param n_permutations Number of permutations for significance testing (default: 1000)
#' @param from_10x If TRUE, seurat_obj is interpreted as path to 10x data (default: FALSE)
#' @param sample.name Sample name for Seurat project (default: "sample")
#' @param mt_pattern Pattern for mitochondrial genes (default: "^Mt-")
#' @param remove_mup Remove Mup genes (default: TRUE)
#' @param do_filter Filter low-quality cells by mitochondrial content (default: TRUE)
#' @param max_mt_percent Maximum mitochondrial percentage for filtering (default: 20)
#' @param gamma_shape Shape parameter for Gamma distribution (default: 5)
#' @param use_pca Use PCA for multi-gene weighting (default: TRUE)
#' @param knn_smooth Apply KNN smoothing after gradient calculation (default: TRUE)
#' @param k_neighbors Number of nearest neighbors for KNN smoothing (default: 20)
#' @param seed Random seed for reproducibility (default: 42)
#' @param verbose Print progress messages (default: TRUE)
#' @return A HepaZone result object containing:
#'   \item{mean_expression}{Gene x Zone matrix of mean expression}
#'   \item{standard_error}{Gene x Zone matrix of standard errors (from bootstrap)}
#'   \item{qvalues}{Gene x Zone matrix of q-values (from permutation test)}
#'   \item{prob_matrix}{Cell x Zone probability matrix}
#'   \item{cv_markers}{Central vein marker genes used}
#'   \item{pn_markers}{Portal marker genes used}
#'   \item{cl_scores}{CL score for each cell}
#'   \item{n_zones}{Number of spatial zones}
#'   \item{svg_results}{Data frame of spatially variable genes}
#'
#' @examples
#' \dontrun{
#' # From 10x data
#' result <- hepa_zone_reconstruct("./data/ZT00", from_10x = TRUE,
#'                                  sample.name = "ZT00")
#'
#' # From existing Seurat object
#' result <- hepa_zone_reconstruct(seurat_obj, n_zones = 10,
#'                                  n_bootstrap = 500, n_permutations = 1000)
#'
#' # With KNN smoothing
#' result <- hepa_zone_reconstruct(seurat_obj, knn_smooth = TRUE, k_neighbors = 30)
#' }
#' @export
hepa_zone_reconstruct <- function(seurat_obj,
                                    cv_markers = NULL,
                                    pn_markers = NULL,
                                    n_zones = 10,
                                    n_bootstrap = 500,
                                    n_permutations = 1000,
                                    from_10x = FALSE,
                                    sample.name = "sample",
                                    mt_pattern = "^Mt-",
                                    remove_mup = TRUE,
                                    do_filter = TRUE,
                                    max_mt_percent = 20,
                                    gamma_shape = 5,
                                    use_pca = TRUE,
                                    knn_smooth = TRUE,
                                    k_neighbors = 20,
                                    seed = 42,
                                    verbose = TRUE) {

  set.seed(seed)

  # Helper function for verbose output
  msg <- function(...) {
    if (verbose) message(sprintf(...))
  }

  # ===========================================================================
  # Step 1: Load/Create Seurat Object
  # ===========================================================================
  msg("Step 1: Loading data...")

  if (from_10x) {
    seurat_obj <- load_10x_data(path = seurat_obj, sample.name = sample.name)
  }

  if (!is(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object or a valid 10x path")
  }

  msg("Loaded %d cells and %d genes", ncol(seurat_obj), nrow(seurat_obj))

  # ===========================================================================
  # Step 2: Quality Filtering
  # ===========================================================================
  msg("Step 2: Quality filtering...")

  if (do_filter) {
    seurat_obj <- filter_by_mitochondrial(
      seurat_obj,
      mt_pattern = mt_pattern,
      max_mt_percent = max_mt_percent,
      remove = TRUE
    )
  }

  # ===========================================================================
  # Step 3: Preprocessing for Zonation
  # ===========================================================================
  msg("Step 3: Preprocessing for zonation analysis...")

  seurat_obj <- preprocess_zonation(
    seurat_obj,
    mt_pattern = mt_pattern,
    remove_mup = remove_mup,
    do_normalize = TRUE
  )

  # ===========================================================================
  # Step 4: Calculate Spatial Position Scores (PCA-based)
  # ===========================================================================
  msg("Step 4: Calculating spatial position scores...")

  seurat_obj <- calculate_spatial_position(
    seurat_obj,
    cv_markers = cv_markers,
    pn_markers = pn_markers,
    use_default_markers = is.null(cv_markers) || is.null(pn_markers),
    use_pca = use_pca
  )

  # Store markers used
  if (is.null(cv_markers)) cv_markers <- .default_cv_markers()
  if (is.null(pn_markers)) pn_markers <- .default_pn_markers()

  # ===========================================================================
  # Step 4b: KNN Smoothing (Optional)
  # ===========================================================================
  if (knn_smooth) {
    msg("Step 4b: Applying KNN smoothing (k=%d)...", k_neighbors)

    seurat_obj <- knn_smooth_scores(
      seurat_obj,
      score_col = "CL_score",
      k = k_neighbors,
      use_pca = TRUE,
      n_pcs = 30,
      weight_by_distance = TRUE
    )
  }

  # ===========================================================================
  # Step 5: Map Cells to Spatial Layers
  # ===========================================================================
  msg("Step 5: Mapping cells to %d spatial layers...", n_zones)

  prob_matrix <- map_cells_to_layers(
    seurat_obj,
    n_zones = n_zones,
    gamma_shape = gamma_shape,
    return_matrix = TRUE
  )

  # ===========================================================================
  # Step 6: Reconstruct Spatial Expression Profiles
  # ===========================================================================
  msg("Step 6: Reconstructing spatial expression profiles...")

  spatial_result <- reconstruct_spatial_expression(seurat_obj, prob_matrix)

  # ===========================================================================
  # Step 7: Bootstrap Standard Error Estimation
  # ===========================================================================
  msg("Step 7: Bootstrap standard error estimation (%d iterations)...", n_bootstrap)

  # Get normalized data using helper
  mat_norm <- .get_data(seurat_obj)
  boot_result <- bootstrap_se(
    expression_matrix = mat_norm,
    prob_matrix = prob_matrix,
    n_bootstrap = n_bootstrap,
    seed = seed
  )

  # ===========================================================================
  # Step 8: Permutation Test
  # ===========================================================================
  msg("Step 8: Permutation test (%d permutations)...", n_permutations)

  perm_result <- permutation_test(
    expression_matrix = mat_norm,
    prob_matrix = prob_matrix,
    n_permutations = n_permutations,
    seed = seed + 1
  )

  # ===========================================================================
  # Step 9: Calculate Q-values
  # ===========================================================================
  msg("Step 9: Multiple testing correction...")

  qvals <- calculate_qvalues(perm_result$p_values)

  # ===========================================================================
  # Step 10: Identify Spatially Variable Genes
  # ===========================================================================
  msg("Step 10: Identifying spatially variable genes...")

  svg_results <- identify_svg(perm_result, qvals, q_threshold = 0.05)
  n_sig <- sum(svg_results$significant)
  msg("Found %d significant spatially variable genes (q < 0.05)", n_sig)

  # ===========================================================================
  # Compile Results
  # ===========================================================================

  # Check if smoothed scores are available
  smoothed_col <- paste0("smoothed_", "CL_score")
  has_smoothed <- smoothed_col %in% colnames(seurat_obj@meta.data)

  result <- list(
    mean_expression = spatial_result$mean_expression,
    standard_error = boot_result$se,
    qvalues = qvals,
    prob_matrix = prob_matrix,
    n_cells = spatial_result$n_cells,
    cv_markers = cv_markers,
    pn_markers = pn_markers,
    cl_scores = seurat_obj$CL_score,
    pca_gradient = if ("pca_gradient" %in% colnames(seurat_obj@meta.data))
                     seurat_obj$pca_gradient else NULL,
    smoothed_scores = if (has_smoothed)
                        seurat_obj@meta.data[[smoothed_col]] else NULL,
    cells = colnames(seurat_obj),
    genes = rownames(seurat_obj),
    n_zones = n_zones,
    svg_results = svg_results,
    permutation_stats = perm_result,
    preprocessing = seurat_obj@misc,
    method_params = list(
      use_pca = use_pca,
      knn_smooth = knn_smooth,
      k_neighbors = k_neighbors
    )
  )

  class(result) <- "HepaZoneResult"

  msg("\n===== HepaZone Analysis Complete =====")
  msg("Output dimensions: %d genes x %d zones",
      nrow(result$mean_expression), result$n_zones)
  msg("Significant SVGs: %d", n_sig)

  return(result)
}


#' Run HepaZone for Multiple Time Points
#'
#' Convenience function to run HepaZone analysis on multiple samples
#' (e.g., different ZT time points) and return combined results.
#'
#' @param sample.paths Named list of paths to 10x directories
#' @param n_zones Number of spatial zones (default: 10)
#' @param n_bootstrap Bootstrap iterations (default: 500)
#' @param n_permutations Permutation iterations (default: 1000)
#' @param ... Additional arguments passed to hepa_zone_reconstruct()
#' @return A list containing:
#'   - results: Named list of HepaZone results per time point
#'   - interaction: Results from time × space interaction analysis
#' @examples
#' \dontrun{
#' sample_paths <- list(
#'   ZT00 = "./data/ZT00",
#'   ZT06 = "./data/ZT06",
#'   ZT12 = "./data/ZT12",
#'   ZT18 = "./data/ZT18"
#' )
#'
#' multi_result <- run_hepa_zone_multi(sample_paths, n_zones = 10)
#' }
#' @export
run_hepa_zone_multi <- function(sample.paths,
                                 n_zones = 10,
                                 n_bootstrap = 500,
                                 n_permutations = 1000,
                                 ...) {
  message("Running HepaZone for multiple time points...")

  results <- list()

  for (sample_name in names(sample.paths)) {
    message("\n===== Processing %s =====", sample_name)

    results[[sample_name]] <- hepa_zone_reconstruct(
      seurat_obj = sample.paths[[sample_name]],
      from_10x = TRUE,
      sample.name = sample_name,
      n_zones = n_zones,
      n_bootstrap = n_bootstrap,
      n_permutations = n_permutations,
      ...
    )
  }

  # Run interaction analysis
  message("\n===== Time x Space Interaction Analysis =====")
  interaction <- analyze_time_space_interaction(results)

  return(list(
    results = results,
    interaction = interaction
  ))
}
