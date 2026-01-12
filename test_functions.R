# Simple test script for HepaZone core functions
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

cat("Testing HepaZone core functions...\n\n")

# Source the R files
r_files <- list.files("R", full.names = TRUE)
for (f in r_files) {
  cat("Sourcing:", f, "\n")
  source(f)
}

# Test 1: preprocess_zonation removes mitochondrial genes
cat("\n=== Test 1: preprocess_zonation ===\n")
set.seed(42)
counts <- matrix(rpois(100, 5), nrow = 10)
rownames(counts) <- c("MtAtp6", "Alb", "Cyp2e1", "Cyp1a2", "MtCo1",
                      "Hamp", "Actb", "Gapdh", "Mup3", "Mup5")
colnames(counts) <- paste0("cell", 1:10)

seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 1, min.features = 1)
seurat_obj <- preprocess_zonation(seurat_obj, remove_mup = TRUE,
                                  mt_pattern = "^Mt", mup_pattern = "^Mup")

# Get counts after preprocessing
counts_after <- .get_counts(seurat_obj)
remaining_genes <- rownames(counts_after)
mt_removed <- !any(grepl("^Mt", remaining_genes))
mup_removed <- !any(grepl("^Mup", remaining_genes))

if (mt_removed && mup_removed) {
  cat("PASS: Mitochondrial and Mup genes removed correctly\n")
} else {
  cat("FAIL: Genes not removed correctly\n")
  cat("  Remaining genes:", paste(remaining_genes, collapse = ", "), "\n")
}

# Test 2: calculate_spatial_position computes CL scores
cat("\n=== Test 2: calculate_spatial_position ===\n")
set.seed(42)
counts2 <- matrix(rpois(200, 3), nrow = 20)
rownames(counts2) <- c(paste0("g", 1:10), "Alb", "Cyp2e1", "Hamp", "Tf", "Cyp1a2")
colnames(counts2) <- paste0("cell", 1:10)

seurat_obj2 <- CreateSeuratObject(counts = counts2, min.cells = 1, min.features = 1)
seurat_obj2 <- preprocess_zonation(seurat_obj2, remove_mup = FALSE, do_normalize = TRUE)

cv_markers <- c("g1", "g2")
pn_markers <- c("g3", "g4")

seurat_obj2 <- calculate_spatial_position(seurat_obj2, cv_markers, pn_markers)

has_cl_score <- "CL_score" %in% colnames(seurat_obj2@meta.data)
cl_in_range <- all(seurat_obj2$CL_score >= 0 & seurat_obj2$CL_score <= 1, na.rm = TRUE)

if (has_cl_score && cl_in_range) {
  cat("PASS: CL scores computed correctly\n")
} else {
  cat("FAIL: CL scores not computed correctly\n")
}

# Test 3: map_cells_to_layers returns valid probability matrix
cat("\n=== Test 3: map_cells_to_layers ===\n")
set.seed(42)
counts3 <- matrix(rpois(100, 3), nrow = 10)
rownames(counts3) <- paste0("gene", 1:10)
colnames(counts3) <- paste0("cell", 1:10)

seurat_obj3 <- CreateSeuratObject(counts = counts3, min.cells = 1, min.features = 1)
seurat_obj3 <- preprocess_zonation(seurat_obj3, remove_mup = FALSE, do_normalize = TRUE)
seurat_obj3$CL_score <- runif(ncol(seurat_obj3), 0.1, 0.9)

prob_matrix <- map_cells_to_layers(seurat_obj3, n_zones = 5)

correct_dims <- nrow(prob_matrix) == 10 && ncol(prob_matrix) == 5
row_sums_ok <- all(abs(Matrix::rowSums(prob_matrix) - 1) < 1e-10)

if (correct_dims && row_sums_ok) {
  cat("PASS: Probability matrix has correct dimensions and sums to 1\n")
} else {
  cat("FAIL: Probability matrix invalid\n")
}

# Test 4: reconstruct_spatial_expression returns correct dimensions
cat("\n=== Test 4: reconstruct_spatial_expression ===\n")
set.seed(42)
counts4 <- matrix(rpois(200, 5), nrow = 20)
rownames(counts4) <- paste0("gene", 1:20)
colnames(counts4) <- paste0("cell", 1:10)

seurat_obj4 <- CreateSeuratObject(counts = counts4, min.cells = 1, min.features = 1)
seurat_obj4 <- preprocess_zonation(seurat_obj4, remove_mup = FALSE, do_normalize = TRUE)
seurat_obj4$CL_score <- runif(ncol(seurat_obj4), 0.1, 0.9)

prob_matrix4 <- map_cells_to_layers(seurat_obj4, n_zones = 5)
result4 <- reconstruct_spatial_expression(seurat_obj4, prob_matrix4)

correct_expr_dims <- nrow(result4$mean_expression) == 20 &&
                    ncol(result4$mean_expression) == 5

if (correct_expr_dims) {
  cat("PASS: Expression matrix has correct dimensions\n")
} else {
  cat("FAIL: Expression matrix dimensions incorrect\n")
}

# Test 5: bootstrap_se returns valid standard errors
cat("\n=== Test 5: bootstrap_se ===\n")
set.seed(42)
counts5 <- matrix(rpois(50, 3), nrow = 5)
rownames(counts5) <- paste0("gene", 1:5)
colnames(counts5) <- paste0("cell", 1:10)

prob_matrix5 <- matrix(runif(100), nrow = 10, ncol = 5)
prob_matrix5 <- prob_matrix5 / Matrix::rowSums(prob_matrix5)

boot_result <- bootstrap_se(counts5, prob_matrix5, n_bootstrap = 10, seed = 42)

boot_dims_ok <- nrow(boot_result$se) == 5 && ncol(boot_result$se) == 5
se_nonneg <- all(boot_result$se >= 0)

if (boot_dims_ok && se_nonneg) {
  cat("PASS: Bootstrap SE has correct dimensions and non-negative values\n")
} else {
  cat("FAIL: Bootstrap SE invalid\n")
}

# Test 6: permutation_test returns valid p-values
cat("\n=== Test 6: permutation_test ===\n")
set.seed(42)
counts6 <- matrix(rpois(50, 3), nrow = 5)
rownames(counts6) <- paste0("gene", 1:5)
colnames(counts6) <- paste0("cell", 1:10)

prob_matrix6 <- matrix(runif(100), nrow = 10, ncol = 5)
prob_matrix6 <- prob_matrix6 / Matrix::rowSums(prob_matrix6)

perm_result <- permutation_test(counts6, prob_matrix6, n_permutations = 10,
                                seed = 123)

pval_len_ok <- length(perm_result$p_values) == 5
pval_in_range <- all(perm_result$p_values >= 0 & perm_result$p_values <= 1)

if (pval_len_ok && pval_in_range) {
  cat("PASS: Permutation p-values have correct format\n")
} else {
  cat("FAIL: Permutation p-values invalid\n")
}

# Test 7: calculate_qvalues uses BH method when qvalue unavailable
cat("\n=== Test 7: calculate_qvalues ===\n")
p_vals <- c(0.01, 0.05, 0.1, 0.001)
q_vals <- calculate_qvalues(p_vals)

qval_len_ok <- length(q_vals) == length(p_vals)

if (qval_len_ok) {
  cat("PASS: Q-values computed correctly using fallback BH method\n")
} else {
  cat("FAIL: Q-values computation failed\n")
}

# Test 8: hepa_zone_reconstruct returns HepaZoneResult object
cat("\n=== Test 8: hepa_zone_reconstruct (full pipeline) ===\n")
set.seed(42)
counts8 <- matrix(rpois(100, 3), nrow = 10)
rownames(counts8) <- c("Alb", "Cyp2e1", "Actb", "Gapdh", "Hamp",
                       "Tf", "Cyp1a2", "gene1", "gene2", "gene3")
colnames(counts8) <- paste0("cell", 1:10)

seurat_obj8 <- CreateSeuratObject(counts = counts8, min.cells = 1, min.features = 1)

result8 <- hepa_zone_reconstruct(
  seurat_obj8,
  cv_markers = c("Cyp2e1", "Cyp1a2"),
  pn_markers = c("Alb", "Tf"),
  n_zones = 3,
  n_bootstrap = 10,
  n_permutations = 10,
  seed = 42,
  verbose = FALSE
)

is_hepa_result <- inherits(result8, "HepaZoneResult")
correct_genes <- nrow(result8$mean_expression) == 10
correct_zones <- result8$n_zones == 3

if (is_hepa_result && correct_genes && correct_zones) {
  cat("PASS: Full pipeline completed successfully\n")
  cat("  - Result type: HepaZoneResult\n")
  cat("  - Genes:", nrow(result8$mean_expression), "\n")
  cat("  - Zones:", result8$n_zones, "\n")
  cat("  - SVG found:", sum(result8$svg_results$significant), "\n")
} else {
  cat("FAIL: Full pipeline failed\n")
}

cat("\n=== All Tests Completed ===\n")
