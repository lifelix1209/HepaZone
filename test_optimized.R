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

# Get the processed seurat_obj with all metadata (cv_score, CL_strength, etc.)
# This is necessary because reference semantics may be broken in some R environments
seurat_obj <- result$seurat_obj

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

# Get normalized data from result
normalized_data <- result$preprocessing$mat_norm

# seurat_obj now has all metadata from result$seurat_obj including:
# - CL_rank, CL_strength, CL_score
# - cv_score, pn_score
# - pca_gradient, pc1
# - smoothed_CL_rank, smoothed_CL_strength

# Verify metadata is present
cat(sprintf("seurat_obj has %d metadata columns\n", ncol(seurat_obj@meta.data)))

# ==============================================================================
# QC Plots - Using new visualization functions
# ==============================================================================

# 1. CV vs PN expression colored by CL score
cat("  - CV vs PN by CL score...\n")
p_qc1 <- plot_cv_pn_by_cl(seurat_obj, normalized_data = normalized_data,
                          use_strength = TRUE,
                          title = "CV vs PN Expression (colored by CL_strength)")
ggsave(file.path(output_dir, "01_cv_pn_by_cl.png"), p_qc1, width = 8, height = 7)

# 2. Raw PCA Gradient vs Final CL Score
cat("  - PCA gradient vs CL score...\n")
p_qc2 <- plot_pca_vs_cl(seurat_obj, use_strength = TRUE,
                        title = "Raw PCA Gradient vs Final CL_strength")
ggsave(file.path(output_dir, "02_pca_vs_cl.png"), p_qc2, width = 8, height = 7)

# 3. Distribution of CL Scores
cat("  - CL score distributions...\n")
p_qc3 <- plot_cl_scores_dist(seurat_obj, bins = 50,
                            title = "Distribution of CL Scores")
ggsave(file.path(output_dir, "03_cl_scores_dist.png"), p_qc3, width = 10, height = 5)

# 4. Marker expression vs CL score (using strength)
cat("  - Marker expression vs CL score...\n")
p_qc4 <- plot_marker_vs_cl(seurat_obj, normalized_data = normalized_data,
                           marker_type = "both",
                           use_strength = TRUE, title = "Marker Expression vs CL_strength")
ggsave(file.path(output_dir, "04_marker_vs_cl.png"), p_qc4, width = 8, height = 7)

# 5. Report marker correlations
cat("  - Computing marker correlations...\n")
corr_results <- report_marker_correlations(seurat_obj, normalized_data = normalized_data,
                                           verbose = TRUE)

# ==============================================================================
# Main Analysis Plots
# ==============================================================================

# 6. CL Score distribution
cat("  - CL score distribution...\n")
p6 <- plot_cl_distribution(seurat_obj, prob_matrix = prob_matrix,
                          title = "Cell Distribution Along Spatial Axis")
ggsave(file.path(output_dir, "05_cl_distribution.png"), p6, width = 8, height = 6)

# 7. Top SVG heatmap using ComplexHeatmap (publication quality)
cat("  - Top SVG heatmap (ComplexHeatmap)...\n")
top_svg_genes <- head(sig_results$gene, 25)
if (length(top_svg_genes) >= 1) {
  p7 <- plot_spatial_heatmap(result, genes = top_svg_genes,
                             title = "Top 25 SVGs",
                             use_ggplot = FALSE)
  png(file.path(output_dir, "06_top_svg_heatmap.png"), width = 14, height = 10,
      units = "in", res = 300)
  print(p7)
  dev.off()
  message("  ComplexHeatmap saved")
}

# 8. Expression gradients with individual gene colors
cat("  - Expression gradients...\n")
cv_genes_plot <- head(sig_results$gene[sapply(sig_results$gene, function(g) cv_expr[g] > pn_expr[g])], 3)
pn_genes_plot <- head(sig_results$gene[sapply(sig_results$gene, function(g) pn_expr[g] > cv_expr[g])], 3)
plot_genes <- c(cv_genes_plot, pn_genes_plot)

if (length(plot_genes) >= 1) {
  n_zones <- result$n_zones
  zone_x <- 1:n_zones
  expr_data <- result$mean_expression[plot_genes, , drop = FALSE]

  # Define distinct colors
  n_genes <- length(plot_genes)
  cv_colors <- c("#E64B35", "#F39B7F", "#FF7F5C")
  pn_colors <- c("#4DBBD5", "#00A087", "#3B4992")

  gene_colors <- character(n_genes)
  for (i in 1:n_genes) {
    g <- plot_genes[i]
    if (cv_expr[g] > pn_expr[g]) {
      gene_colors[i] <- cv_colors[(i - 1) %% length(cv_colors) + 1]
    } else {
      gene_colors[i] <- pn_colors[(i - 1) %% length(pn_colors) + 1]
    }
  }

  # Create data frame
  plot_df <- data.frame()
  for (i in 1:n_genes) {
    g <- plot_genes[i]
    gene_df <- data.frame(
      zone = zone_x,
      expression = as.numeric(expr_data[g, ]),
      gene = g,
      color = gene_colors[i]
    )
    plot_df <- rbind(plot_df, gene_df)
  }

  label_df <- plot_df[plot_df$zone == n_zones, ]
  label_df$zone <- n_zones + 0.3

  p8 <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = 0.5, xmax = n_zones/2 + 0.5, ymin = -Inf, ymax = Inf,
             fill = "#FFF5F5", alpha = 0.4) +
    ggplot2::annotate("rect", xmin = n_zones/2 + 0.5, xmax = n_zones + 0.5, ymin = -Inf, ymax = Inf,
             fill = "#F0F7FF", alpha = 0.4) +
    ggplot2::annotate("text", x = n_zones/4, y = Inf, label = "Central Vein", vjust = 2,
             fontface = "bold", color = "#8B0000", size = 4) +
    ggplot2::annotate("text", x = 3*n_zones/4, y = Inf, label = "Portal Vein", vjust = 2,
             fontface = "bold", color = "#00008B", size = 4) +
    ggplot2::geom_line(data = plot_df, ggplot2::aes(x = zone, y = expression, group = gene),
              color = plot_df$color, linewidth = 1.2, alpha = 0.9) +
    ggplot2::geom_text(data = label_df, ggplot2::aes(x = zone, y = expression, label = gene),
              color = label_df$color, size = 3.5, hjust = 0, fontface = "italic") +
    ggplot2::scale_x_continuous(breaks = 1:n_zones, labels = paste0("Z", 1:n_zones),
                       expand = ggplot2::expansion(mult = c(0.02, 0.15))) +
    ggplot2::labs(
      title = "Expression Gradients of Spatially Variable Genes",
      subtitle = sprintf("%d CV genes (red) | %d PN genes (blue)",
                        sum(sapply(plot_genes, function(g) cv_expr[g] > pn_expr[g])),
                        sum(sapply(plot_genes, function(g) pn_expr[g] > cv_expr[g]))),
      x = "Spatial Zone (Central Vein -> Portal Vein)",
      y = "Normalized Expression"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5, color = "gray40"),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text = ggplot2::element_text(size = 10),
      legend.position = "none",
      panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_line(color = "gray95", linewidth = 0.3),
      plot.margin = ggplot2::margin(r = 80)
    )

  ggsave(file.path(output_dir, "07_svg_gradients.png"), p8, width = 12, height = 7)
  message("  Gradient plot saved")
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
cat("QC Plots (new visualization functions):\n")
cat("  01. 01_cv_pn_by_cl.png              - CV vs PN by CL score\n")
cat("  02. 02_pca_vs_cl.png                - PCA gradient vs CL score\n")
cat("  03. 03_cl_scores_dist.png           - CL rank/strength distributions\n")
cat("  04. 04_marker_vs_cl.png             - Marker expression vs CL\n")
cat("\nMain Analysis Plots:\n")
cat("  05. 05_cl_distribution.png          - Cell distribution\n")
cat("  06. 06_top_svg_heatmap.png          - Top 25 SVGs (ComplexHeatmap)\n")
cat("  07. 07_svg_gradients.png            - Expression gradients\n")
cat("\nData Files:\n")
cat("  08. mean_expression.csv\n")
cat("  09. qvalues.csv\n")
cat("  10. expression_comparison.csv\n")
cat("  11. hepa_zone_full_results.rds\n\n")

cat("=============================================================\n")
