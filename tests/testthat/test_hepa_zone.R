# Testthat tests for HepaZone package

test_that("preprocess_zonation removes mitochondrial genes", {
  # Create a small test Seurat object
  set.seed(42)
  counts <- matrix(rpois(100, 5), nrow = 10)
  rownames(counts) <- c("Mt-Atp6", "Alb", "Cyp2e1", "Cyp1a2", "Mt-Co1",
                        "Hamp", "Actb", "Gapdh", "Mup3", "Mup5")
  colnames(counts) <- paste0("cell_", 1:10)

  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, min.cells = 1, min.features = 1)
  seurat_obj <- preprocess_zonation(seurat_obj, remove_mup = TRUE)

  # Check that Mt- genes are removed
  remaining_genes <- rownames(seurat_obj)
  expect_false(any(grepl("^Mt-", remaining_genes)))
  expect_false(any(grepl("^Mup", remaining_genes)))
})

test_that("calculate_spatial_position computes CL scores", {
  set.seed(42)
  counts <- matrix(rpois(200, 3), nrow = 20)
  rownames(counts) <- paste0("gene_", 1:20)
  colnames(counts) <- paste0("cell_", 1:10)

  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, min.cells = 1, min.features = 1)
  seurat_obj <- preprocess_zonation(seurat_obj, remove_mup = FALSE, do_normalize = TRUE)

  # Add fake marker genes
  cv_markers <- c("gene_1", "gene_2")
  pn_markers <- c("gene_3", "gene_4")

  seurat_obj <- calculate_spatial_position(seurat_obj, cv_markers, pn_markers)

  # Check that CL_score exists
  expect_true("CL_score" %in% colnames(seurat_obj@meta.data))
  expect_true(all(seurat_obj$CL_score >= 0 & seurat_obj$CL_score <= 1))
})

test_that("map_cells_to_layers returns valid probability matrix", {
  set.seed(42)
  counts <- matrix(rpois(100, 3), nrow = 10)
  rownames(counts) <- paste0("gene_", 1:10)
  colnames(counts) <- paste0("cell_", 1:10)

  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, min.cells = 1, min.features = 1)
  seurat_obj <- preprocess_zonation(seurat_obj, remove_mup = FALSE, do_normalize = TRUE)

  # Add CL scores manually
  seurat_obj$CL_score <- runif(ncol(seurat_obj), 0.1, 0.9)

  prob_matrix <- map_cells_to_layers(seurat_obj, n_zones = 5)

  # Check dimensions
  expect_equal(nrow(prob_matrix), ncol(seurat_obj))
  expect_equal(ncol(prob_matrix), 5)

  # Check probabilities sum to 1
  row_sums <- Matrix::rowSums(prob_matrix)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("reconstruct_spatial_expression returns correct dimensions", {
  set.seed(42)
  counts <- matrix(rpois(200, 5), nrow = 20)
  rownames(counts) <- paste0("gene_", 1:20)
  colnames(counts) <- paste0("cell_", 1:10)

  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, min.cells = 1, min.features = 1)
  seurat_obj <- preprocess_zonation(seurat_obj, remove_mup = FALSE, do_normalize = TRUE)
  seurat_obj$CL_score <- runif(ncol(seurat_obj), 0.1, 0.9)

  prob_matrix <- map_cells_to_layers(seurat_obj, n_zones = 5)
  result <- reconstruct_spatial_expression(seurat_obj, prob_matrix)

  # Check dimensions
  expect_equal(nrow(result$mean_expression), 20)  # genes
  expect_equal(ncol(result$mean_expression), 5)   # zones
})

test_that("bootstrap_se returns valid standard errors", {
  set.seed(42)
  counts <- matrix(rpois(50, 3), nrow = 5)
  rownames(counts) <- paste0("gene_", 1:5)
  colnames(counts) <- paste0("cell_", 1:10)

  prob_matrix <- matrix(runif(100), nrow = 10, ncol = 5)
  prob_matrix <- prob_matrix / rowSums(prob_matrix)

  result <- bootstrap_se(counts, prob_matrix, n_bootstrap = 10, seed = 42)

  expect_equal(nrow(result$se), 5)
  expect_equal(ncol(result$se), 5)
  expect_true(all(result$se >= 0))
})

test_that("permutation_test returns valid p-values", {
  set.seed(42)
  counts <- matrix(rpois(50, 3), nrow = 5)
  rownames(counts) <- paste0("gene_", 1:5)
  colnames(counts) <- paste0("cell_", 1:10)

  prob_matrix <- matrix(runif(100), nrow = 10, ncol = 5)
  prob_matrix <- prob_matrix / rowSums(prob_matrix)

  result <- permutation_test(counts, prob_matrix, n_permutations = 10, seed = 123)

  expect_equal(length(result$p_values), 5)
  expect_true(all(result$p_values >= 0 & result$p_values <= 1))
})

test_that("hepa_zone_reconstruct returns HepaZoneResult object", {
  set.seed(42)
  counts <- matrix(rpois(100, 3), nrow = 10)
  rownames(counts) <- c("Alb", "Cyp2e1", "Actb", "Gapdh", "Hamp",
                        "Tf", "Cyp1a2", "Mup3", "Mup5", "Mt-Co1")
  colnames(counts) <- paste0("cell_", 1:10)

  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, min.cells = 1, min.features = 1)

  result <- hepa_zone_reconstruct(
    seurat_obj,
    n_zones = 3,
    n_bootstrap = 10,
    n_permutations = 10,
    seed = 42,
    verbose = FALSE
  )

  expect_s3_class(result, "HepaZoneResult")
  expect_equal(nrow(result$mean_expression), 10)
  expect_equal(result$n_zones, 3)
})

test_that("calculate_qvalues handles edge cases", {
  p_vals <- c(0.01, 0.05, 0.1, NA, 0.001)
  q_vals <- calculate_qvalues(p_vals)

  expect_equal(length(q_vals), length(p_vals))
  expect_true(is.na(q_vals[4]))  # NA should remain NA
})
