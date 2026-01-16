# ==============================================================================
# HepaZone Test with Optimized Algorithm (PCA + KNN Smoothing)
# ==============================================================================
# This script tests the complete pipeline using hepa_zone_reconstruct()
# which handles all steps automatically.
# ==============================================================================

library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

cat("=============================================================\n")
cat("  HepaZone Test - Optimized Algorithm\n")
cat("  (PCA-based weighting + KNN smoothing)\n")
cat("=============================================================\n\n")

# Load HepaZone functions
source("R/data_io.R")
source("R/preprocessing.R")
source("R/zonation.R")
source("R/statistics.R")
source("R/visualization.R")
source("R/main.R")

# ==============================================================================
# Run Complete Analysis Pipeline
# ==============================================================================
cat(">>> Loading data and running analysis...\n\n")

data_path <- "GSE145197_data/10x_format/ZT00/ZT00A"
n_zones <- 10

seurat_obj <- load_10x_data(
  path = data_path,
  sample.name = "ZT00A",
  min.cells = 10,
  min.features = 200
)

cat(sprintf("Loaded %d cells x %d genes\n\n", nrow(seurat_obj), ncol(seurat_obj)))

# Run complete pipeline in one call
result <- hepa_zone_reconstruct(
  seurat_obj,
  n_zones = n_zones,
  n_bootstrap = 100,
  n_permutations = 200,
  use_pca = TRUE,
  knn_smooth = TRUE,
  k_neighbors = 20,
  seed = 42,
  verbose = TRUE
)

# ==============================================================================
# Results Summary
# ==============================================================================
cat("\n=============================================================\n")
cat("  Analysis Results\n")
cat("=============================================================\n\n")

cat("Method used:\n")
cat(sprintf("  - PCA-based weighting: %s\n", ifelse(result$method_params$use_pca, "Yes", "No")))
cat(sprintf("  - KNN smoothing: %s (k=%d)\n",
            ifelse(result$method_params$knn_smooth, "Yes", "No"),
            result$method_params$k_neighbors))

cat("\nSignificant Spatial Variable Genes:\n")
sig_005 <- sum(result$svg_results$significant, na.rm = TRUE)
sig_010 <- sum(result$svg_results$q_value < 0.10, na.rm = TRUE)
cat(sprintf("  - q < 0.05: %d genes\n", sig_005))
cat(sprintf("  - q < 0.10: %d genes\n", sig_010))

# Classify SVGs by direction
sig_results <- result$svg_results[result$svg_results$significant, ]
n_z <- result$n_zones
cv_idx <- 1:floor(n_z / 3)
pn_idx <- (ceiling(2 * n_z / 3) + 1):n_z

cv_expr <- apply(result$mean_expression[sig_results$gene, cv_idx, drop = FALSE], 1, mean)
pn_expr <- apply(result$mean_expression[sig_results$gene, pn_idx, drop = FALSE], 1, mean)

n_cv_enriched <- sum(cv_expr > pn_expr)
n_pn_enriched <- sum(pn_expr > cv_expr)

cat(sprintf("  - CV-enriched (central): %d\n", n_cv_enriched))
cat(sprintf("  - PN-enriched (portal): %d\n\n", n_pn_enriched))

# Top SVGs
cat("Top 10 Spatial Variable Genes:\n")
sig_results <- sig_results[order(sig_results$q_value), ]
for (i in 1:min(10, nrow(sig_results))) {
  g <- sig_results$gene[i]
  qv <- sig_results$q_value[i]
  direction <- ifelse(cv_expr[g] > pn_expr[g], "CV", "PN")
  cat(sprintf("  %2d. %-15s q=%.4f [%s]\n", i, g, qv, direction))
}

# ==============================================================================
# Visualization
# ==============================================================================
cat("\n>>> Generating visualizations...\n\n")

output_dir <- "GSE145197_results/optimized_test"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get prob_matrix from result
prob_matrix <- result$prob_matrix

# 1. CL Score distribution
cat("  - CL score distribution...\n")
p1 <- plot_cl_distribution(seurat_obj, prob_matrix = prob_matrix,
                          title = "Optimized Algorithm: Cell Distribution")
ggsave(file.path(output_dir, "cl_score_distribution.png"), p1, width = 8, height = 6)

# 2. Top SVG heatmap
cat("  - Top SVG heatmap...\n")
top_svg_genes <- head(sig_results$gene, 25)
if (length(top_svg_genes) >= 1) {
  p2 <- plot_spatial_heatmap(result, genes = top_svg_genes, n_genes = min(25, length(top_svg_genes)),
                             title = "Top 25 SVGs (Optimized Algorithm)")
  ggsave(file.path(output_dir, "top_svg_heatmap.png"), p2, width = 12, height = 10)
}

# 3. Expression gradients
cat("  - Expression gradients...\n")
cv_genes_plot <- head(sig_results$gene[sapply(sig_results$gene, function(g) cv_expr[g] > pn_expr[g])], 3)
pn_genes_plot <- head(sig_results$gene[sapply(sig_results$gene, function(g) pn_expr[g] > cv_expr[g])], 3)
plot_genes <- c(cv_genes_plot, pn_genes_plot)

if (length(plot_genes) >= 1) {
  n_cv <- sum(sapply(plot_genes, function(g) cv_expr[g] > pn_expr[g]))
  n_pn <- length(plot_genes) - n_cv
  colors_vec <- c(rep("#E64B35", n_cv), rep("#4DBBD5", n_pn))
  p3 <- plot_gradient(result, genes = plot_genes,
                     title = "Expression Gradients of SVGs",
                     colors = colors_vec)
  ggsave(file.path(output_dir, "svg_expression_gradients.png"), p3, width = 10, height = 7)
}

# ==============================================================================
# Export Results
# ==============================================================================
cat("\n>>> Exporting results...\n\n")

export_hepa_zone_results(result, output.dir = output_dir)

# Save comparison data
comparison_df <- data.frame(
  gene = rownames(result$mean_expression),
  cv_mean = rowMeans(result$mean_expression[, cv_idx]),
  pn_mean = rowMeans(result$mean_expression[, pn_idx]),
  q_value = result$svg_results$q_value[match(rownames(result$mean_expression),
                                              result$svg_results$gene)]
)
write.csv(comparison_df, file.path(output_dir, "expression_comparison.csv"))

# ==============================================================================
# Summary
# ==============================================================================
cat("\n=============================================================\n")
cat("  Test Complete!\n")
cat("=============================================================\n\n")

cat("Output directory:", output_dir, "\n\n")
cat("Files generated:\n")
cat("  1. cl_score_distribution.png\n")
cat("  2. top_svg_heatmap.png\n")
cat("  3. svg_expression_gradients.png\n")
cat("  4. mean_expression.csv\n")
cat("  5. qvalues.csv\n")
cat("  6. expression_comparison.csv\n")
cat("  7. hepa_zone_full_results.rds\n\n")

cat("=============================================================\n")
