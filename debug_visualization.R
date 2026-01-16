# ==============================================================================
# HepaZone Debug Visualization Script
# ==============================================================================
# Plots:
# 1. Raw spatial positioning vs CL_score
# 2. CV marker expression vs CL_score
# 3. PV marker expression vs CL_score
# 4. Zone mapping with soft probabilities
# 5. Entropy and confidence distributions
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)

cat("=============================================================\n")
cat("  HepaZone Debug Visualization\n")
cat("=============================================================\n\n")

# Load HepaZone functions
source("R/data_io.R")
source("R/preprocessing.R")
source("R/zonation.R")
source("R/statistics.R")
source("R/visualization.R")
source("R/main.R")

# ==============================================================================
# Load Data
# ==============================================================================
cat(">>> Loading data (ZT00A sample)...\n\n")

data_path <- "GSE145197_data/10x_format/ZT00/ZT00A"
seurat_obj <- load_10x_data(
  path = data_path,
  sample.name = "ZT00A",
  min.cells = 10,
  min.features = 200
)

cat(sprintf("Loaded: %d genes x %d cells\n\n", nrow(seurat_obj), ncol(seurat_obj)))

# ==============================================================================
# Preprocessing
# ==============================================================================
cat(">>> Preprocessing...\n\n")

seurat_obj <- preprocess_zonation(
  seurat_obj,
  mt_pattern = "^Mt-",
  remove_mup = TRUE
)

# ==============================================================================
# Calculate Spatial Position
# ==============================================================================
cat(">>> Calculating spatial position (PCA-based)...\n\n")

seurat_obj <- calculate_spatial_position(
  seurat_obj,
  use_default_markers = TRUE,
  use_pca = TRUE
)

# Check what columns are in metadata
cat("\nMetadata columns available:\n")
print(colnames(seurat_obj@meta.data))

# ==============================================================================
# KNN Smoothing
# ==============================================================================
cat("\n>>> Applying KNN smoothing...\n\n")

seurat_obj <- knn_smooth_scores(
  seurat_obj,
  score_col = "CL_score",
  k = 20,
  use_pca = TRUE,
  n_pcs = 30,
  weight_by_distance = TRUE
)

# ==============================================================================
# Debug: Check available columns after KNN smoothing
# ==============================================================================
cat("\nMetadata columns after KNN smoothing:\n")
print(colnames(seurat_obj@meta.data))

# Create output directory
output_dir <- "GSE145197_results/debug_plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Plot 1: Raw spatial positioning vs CL_score
# ==============================================================================
cat("\n>>> Plot 1: Raw spatial positioning vs CL_score\n")

# Get the raw PCA gradient and smoothed CL_score
if ("pca_gradient" %in% colnames(seurat_obj@meta.data)) {
  raw_gradient <- seurat_obj$pca_gradient
  final_cl <- seurat_obj$CL_score

  # Create comparison data frame
  compare_df <- data.frame(
    raw_pca_gradient = raw_gradient,
    final_cl_score = final_cl
  )

  p1 <- ggplot(compare_df, aes(x = raw_pca_gradient, y = final_cl_score)) +
    geom_point(alpha = 0.3, size = 0.5, color = "#3B4992") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "Raw PCA Gradient vs Final CL Score",
      subtitle = "After KNN Smoothing",
      x = "Raw PCA Gradient (Before KNN Smoothing)",
      y = "Final CL Score (After KNN Smoothing)"
    ) +
    .theme_publication() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(output_dir, "01_raw_vs_final_cl_score.png"), p1, width = 8, height = 6)

  # Calculate correlation
  corr <- cor(raw_gradient, final_cl, use = "complete.obs")
  cat(sprintf("  Correlation (raw vs final): %.4f\n", corr))
}

# ==============================================================================
# Plot 2: CV Marker Expression vs CL_score
# ==============================================================================
cat("\n>>> Plot 2: CV Marker Expression vs CL_score\n")

# Get CV markers used
cv_markers <- .default_cv_markers_lower()
pn_markers <- .default_pn_markers_lower()

# Filter to available markers
mat_norm <- .get_data(seurat_obj)
available_cv <- cv_markers[cv_markers %in% rownames(mat_norm)]
available_pn <- pn_markers[pn_markers %in% rownames(mat_norm)]

cat(sprintf("  Available CV markers: %d/%d\n", length(available_cv), length(cv_markers)))
cat(sprintf("  Available PN markers: %d/%d\n", length(available_pn), length(pn_markers)))

# Get CL score
cl_score <- seurat_obj$CL_score

# Calculate mean CV marker expression per cell (ensure proper ordering)
cv_expr <- colMeans(mat_norm[available_cv, , drop = FALSE], na.rm = TRUE)

# Ensure cl_score and cv_expr are in the same order
common_cells <- intersect(names(cl_score), names(cv_expr))
cl_score_aligned <- cl_score[common_cells]
cv_expr_aligned <- cv_expr[common_cells]

# Create data frame
cv_df <- data.frame(
  cell = common_cells,
  cl_score = cl_score_aligned,
  cv_expr = cv_expr_aligned
)

# Plot CV expression vs CL score
p2 <- ggplot(cv_df, aes(x = cv_expr, y = cl_score)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#E64B35") +
  geom_smooth(method = "loess", color = "#D55E00", linewidth = 1.2) +
  labs(
    title = "CV Marker Expression vs CL Score",
    subtitle = sprintf("CV markers: %s", paste(available_cv[1:min(5, length(available_cv))], collapse = ", ")),
    x = "Mean CV Marker Expression (Normalized)",
    y = "CL Score (0=CV, 1=PV)"
  ) +
  .theme_publication() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "02_cv_expression_vs_cl_score.png"), p2, width = 8, height = 6)

# Calculate correlation
corr_cv <- cor(cv_expr_aligned, cl_score_aligned, use = "complete.obs")
cat(sprintf("  Correlation (CV expr vs CL_score): %.4f (expected: negative)\n", corr_cv))

# ==============================================================================
# Plot 3: PV Marker Expression vs CL_score
# ==============================================================================
cat("\n>>> Plot 3: PV Marker Expression vs CL_score\n")

# Calculate mean PN marker expression per cell
pn_expr <- colMeans(mat_norm[available_pn, , drop = FALSE], na.rm = TRUE)

# Ensure same ordering
pn_expr_aligned <- pn_expr[common_cells]

# Create data frame
pn_df <- data.frame(
  cell = common_cells,
  cl_score = cl_score_aligned,
  pn_expr = pn_expr_aligned
)

# Plot PN expression vs CL score
p3 <- ggplot(pn_df, aes(x = pn_expr, y = cl_score)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#4DBBD5") +
  geom_smooth(method = "loess", color = "#00A087", linewidth = 1.2) +
  labs(
    title = "PV Marker Expression vs CL Score",
    subtitle = sprintf("PV markers: %s", paste(available_pn[1:min(5, length(available_pn))], collapse = ", ")),
    x = "Mean PV Marker Expression (Normalized)",
    y = "CL Score (0=CV, 1=PV)"
  ) +
  .theme_publication() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "03_pv_expression_vs_cl_score.png"), p3, width = 8, height = 6)

# Calculate correlation
corr_pn <- cor(pn_expr_aligned, cl_score_aligned, use = "complete.obs")
cat(sprintf("  Correlation (PN expr vs CL_score): %.4f (expected: positive)\n", corr_pn))

# ==============================================================================
# Plot 4: Combined CV and PN expression with CL score coloring
# ==============================================================================
cat("\n>>> Plot 4: CV vs PN expression (colored by CL_score)\n")

combined_df <- data.frame(
  cell = common_cells,
  cv_expr = cv_expr_aligned,
  pn_expr = pn_expr_aligned,
  cl_score = cl_score_aligned
)

p4 <- ggplot(combined_df, aes(x = cv_expr, y = pn_expr, color = cl_score)) +
  geom_point(alpha = 0.4, size = 0.5) +
  scale_color_viridis_c(name = "CL Score", option = "viridis") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "CV vs PV Marker Expression",
    subtitle = "Colored by CL Score",
    x = "CV Marker Expression",
    y = "PV Marker Expression"
  ) +
  .theme_publication() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "04_cv_vs_pn_expression.png"), p4, width = 8, height = 6)

# ==============================================================================
# Plot 5: Distribution of CL scores
# ==============================================================================
cat("\n>>> Plot 5: CL Score Distribution\n")

score_df <- data.frame(
  raw_cl = seurat_obj$pca_gradient,
  cl_rank = seurat_obj$CL_rank,
  cl_strength = seurat_obj$CL_strength
)

p5a <- ggplot(score_df, aes(x = cl_rank)) +
  geom_histogram(bins = 50, fill = "#3B4992", color = "white", alpha = 0.7) +
  labs(
    title = "CL_rank Distribution",
    subtitle = sprintf("Rank-based [0,1]: Mean=%.3f, SD=%.3f", mean(score_df$cl_rank, na.rm = TRUE), sd(score_df$cl_rank, na.rm = TRUE)),
    x = "CL_rank",
    y = "Number of Cells"
  ) +
  .theme_publication()

p5b <- ggplot(score_df, aes(x = cl_strength)) +
  geom_histogram(bins = 50, fill = "#E64B35", color = "white", alpha = 0.7) +
  labs(
    title = "CL_strength Distribution",
    subtitle = sprintf("Sigmoid [0,1]: Mean=%.3f, SD=%.3f", mean(score_df$cl_strength, na.rm = TRUE), sd(score_df$cl_strength, na.rm = TRUE)),
    x = "CL_strength",
    y = "Number of Cells"
  ) +
  .theme_publication()

ggsave(file.path(output_dir, "05a_cl_rank_distribution.png"), p5a, width = 8, height = 6)
ggsave(file.path(output_dir, "05b_cl_strength_distribution.png"), p5b, width = 8, height = 6)

# ==============================================================================
# Plot 6: CL_rank vs CL_strength comparison
# ==============================================================================
cat("\n>>> Plot 6: CL_rank vs CL_strength comparison\n")

p6 <- ggplot(score_df, aes(x = cl_rank, y = cl_strength)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#3B4992") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "CL_rank vs CL_strength",
    subtitle = "Comparison of two transformation methods",
    x = "CL_rank (Rank-based)",
    y = "CL_strength (Robust Sigmoid)"
  ) +
  .theme_publication() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "06_cl_rank_vs_strength.png"), p6, width = 8, height = 6)

# Calculate correlation
corr_rank_strength <- cor(score_df$cl_rank, score_df$cl_strength, use = "complete.obs")
cat(sprintf("  Correlation (CL_rank vs CL_strength): %.4f\n", corr_rank_strength))

# ==============================================================================
# Plot 7: CL_strength vs CV/PV marker expression
# ==============================================================================
cat("\n>>> Plot 7: CL_strength vs CV/PV marker expression\n")

combined_viz_df <- data.frame(
  cell = common_cells,
  cv_expr = cv_expr_aligned,
  pn_expr = pn_expr_aligned,
  cl_rank = seurat_obj$CL_rank[common_cells],
  cl_strength = seurat_obj$CL_strength[common_cells]
)

p7a <- ggplot(combined_viz_df, aes(x = cv_expr, y = cl_strength)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#E64B35") +
  geom_smooth(method = "loess", color = "#D55E00", linewidth = 1.2) +
  labs(
    title = "CV Marker Expression vs CL_strength",
    subtitle = "Robust sigmoid preserves magnitude information",
    x = "Mean CV Marker Expression",
    y = "CL_strength (Sigmoid)"
  ) +
  .theme_publication() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p7b <- ggplot(combined_viz_df, aes(x = pn_expr, y = cl_strength)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#4DBBD5") +
  geom_smooth(method = "loess", color = "#00A087", linewidth = 1.2) +
  labs(
    title = "PV Marker Expression vs CL_strength",
    subtitle = "Shows stronger signal for high-confidence cells",
    x = "Mean PV Marker Expression",
    y = "CL_strength (Sigmoid)"
  ) +
  .theme_publication() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "07a_cv_vs_cl_strength.png"), p7a, width = 8, height = 6)
ggsave(file.path(output_dir, "07b_pv_vs_cl_strength.png"), p7b, width = 8, height = 6)

# ==============================================================================
# NEW: Plot 8: Beta Distribution Zone Mapping
# ==============================================================================
cat("\n>>> Plot 8: Beta Distribution Zone Mapping\n")

# Map cells to zones using Beta distribution
zone_result <- map_cells_to_layers(seurat_obj, n_zones = 10, kappa = 30, return_matrix = TRUE)

# Get probability matrix
prob_matrix <- zone_result$prob_matrix
zone_hard <- zone_result$zone_hard
entropy <- zone_result$entropy
confidence <- zone_result$confidence
max_prob <- zone_result$max_prob

cat(sprintf("  Zone mapping complete: %d zones\n", ncol(prob_matrix)))
cat(sprintf("  Entropy: mean=%.3f, median=%.3f\n", mean(entropy), median(entropy)))
cat(sprintf("  Confidence: mean=%.3f, median=%.3f\n", mean(confidence), median(confidence)))

# Plot entropy distribution
p8a <- ggplot(data.frame(entropy = entropy), aes(x = entropy)) +
  geom_histogram(bins = 50, fill = "#8491B4", color = "white", alpha = 0.7) +
  labs(
    title = "Zone Probability Entropy Distribution",
    subtitle = "Higher entropy = more uncertain zone assignment",
    x = "Shannon Entropy",
    y = "Number of Cells"
  ) +
  .theme_publication()

ggsave(file.path(output_dir, "08a_entropy_distribution.png"), p8a, width = 8, height = 6)

# Plot confidence distribution
p8b <- ggplot(data.frame(confidence = confidence), aes(x = confidence)) +
  geom_histogram(bins = 50, fill = "#00A087", color = "white", alpha = 0.7) +
  labs(
    title = "Zone Confidence Distribution",
    subtitle = "Higher confidence = more certain zone assignment",
    x = "Confidence (1 - normalized entropy)",
    y = "Number of Cells"
  ) +
  .theme_publication()

ggsave(file.path(output_dir, "08b_confidence_distribution.png"), p8b, width = 8, height = 6)

# Plot max probability distribution
p8c <- ggplot(data.frame(max_prob = max_prob), aes(x = max_prob)) +
  geom_histogram(bins = 50, fill = "#3B4992", color = "white", alpha = 0.7) +
  labs(
    title = "Maximum Zone Probability Distribution",
    subtitle = "Higher max_prob = more concentrated probability mass",
    x = "Maximum Zone Probability",
    y = "Number of Cells"
  ) +
  .theme_publication()

ggsave(file.path(output_dir, "08c_max_prob_distribution.png"), p8c, width = 8, height = 6)

# ==============================================================================
# NEW: Plot 9: Zone Assignment Comparison
# ==============================================================================
cat("\n>>> Plot 9: Zone Assignment Comparison (Hard vs Soft)\n")

# Create comparison data frame
zone_compare_df <- data.frame(
  cell = rownames(prob_matrix),
  zone_hard = zone_hard,
  zone_soft = apply(prob_matrix, 1, which.max),
  confidence = confidence,
  cl_rank = seurat_obj$CL_rank[rownames(prob_matrix)]
)

# Calculate agreement
agreement <- mean(zone_compare_df$zone_hard == zone_compare_df$zone_soft)
cat(sprintf("  Hard vs Soft zone agreement: %.1f%%\n", 100 * agreement))

# Plot hard vs soft zone assignment
p9 <- ggplot(zone_compare_df, aes(x = zone_hard, y = zone_soft)) +
  geom_bin2d(bins = 10) +
  scale_fill_viridis_c(name = "Count") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Hard vs Soft Zone Assignment",
    subtitle = sprintf("Agreement: %.1f%%", 100 * agreement),
    x = "Hard Zone (CL_rank-based)",
    y = "Soft Zone (Beta-based)"
  ) +
  .theme_publication() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "09_hard_vs_soft_zone.png"), p9, width = 8, height = 6)

# ==============================================================================
# NEW: Plot 10: Confidence vs CL_strength Relationship
# ==============================================================================
cat("\n>>> Plot 10: Confidence vs CL_strength Relationship\n")

conf_strength_df <- data.frame(
  cl_strength = seurat_obj$CL_strength[rownames(prob_matrix)],
  confidence = confidence,
  max_prob = max_prob
)

p10a <- ggplot(conf_strength_df, aes(x = cl_strength, y = confidence)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#00A087") +
  geom_smooth(method = "loess", color = "#3B4992", linewidth = 1.2) +
  labs(
    title = "CL_strength vs Zone Confidence",
    subtitle = "Higher CL_strength generally leads to higher confidence",
    x = "CL_strength",
    y = "Zone Confidence"
  ) +
  .theme_publication() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "10a_strength_vs_confidence.png"), p10a, width = 8, height = 6)

p10b <- ggplot(conf_strength_df, aes(x = cl_strength, y = max_prob)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#E64B35") +
  geom_smooth(method = "loess", color = "#3B4992", linewidth = 1.2) +
  labs(
    title = "CL_strength vs Max Zone Probability",
    subtitle = "Higher CL_strength leads to sharper probability distributions",
    x = "CL_strength",
    y = "Maximum Zone Probability"
  ) +
  .theme_publication() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "10b_strength_vs_max_prob.png"), p10b, width = 8, height = 6)

# ==============================================================================
# NEW: Plot 11: Probability Heatmap (sample)
# ==============================================================================
cat("\n>>> Plot 11: Probability Heatmap (sample cells)\n")

# Sample 200 cells for visualization
set.seed(42)
sample_cells <- sample(rownames(prob_matrix), min(200, nrow(prob_matrix)))

# Order by CL_rank
sample_order <- order(seurat_obj$CL_rank[sample_cells])
sample_cells <- sample_cells[sample_order]

# Create heatmap data
prob_sample <- prob_matrix[sample_cells, ]

# Reshape for ggplot
prob_long <- as.data.frame(prob_sample)
prob_long$cell <- rownames(prob_sample)
prob_long <- prob_long %>%
  tidyr::pivot_longer(cols = -cell, names_to = "zone", values_to = "probability")

# Add CL_rank for ordering
prob_long$cl_rank <- seurat_obj$CL_rank[prob_long$cell]

# Plot heatmap
p11 <- ggplot(prob_long, aes(x = zone, y = cell, fill = probability)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Probability", limits = c(0, 1)) +
  labs(
    title = "Cell-to-Zone Probability Matrix (Sample)",
    subtitle = "Cells ordered by CL_rank",
    x = "Zone",
    y = "Cell"
  ) +
  .theme_publication() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(file.path(output_dir, "11_probability_heatmap.png"), p11, width = 10, height = 8)

# ==============================================================================
# Summary Statistics
# ==============================================================================
cat("\n=============================================================\n")
cat("  Summary Statistics\n")
cat("=============================================================\n\n")

cat("CL Score Statistics:\n")
cat(sprintf("  CL_rank    - Mean: %.3f, SD: %.3f, Range: [%.3f, %.3f]\n",
            mean(seurat_obj$CL_rank, na.rm = TRUE),
            sd(seurat_obj$CL_rank, na.rm = TRUE),
            min(seurat_obj$CL_rank, na.rm = TRUE),
            max(seurat_obj$CL_rank, na.rm = TRUE)))
cat(sprintf("  CL_strength- Mean: %.3f, SD: %.3f, Range: [%.3f, %.3f]\n",
            mean(seurat_obj$CL_strength, na.rm = TRUE),
            sd(seurat_obj$CL_strength, na.rm = TRUE),
            min(seurat_obj$CL_strength, na.rm = TRUE),
            max(seurat_obj$CL_strength, na.rm = TRUE)))

cat("\nCorrelations:\n")
cat(sprintf("  CL_rank vs CL_strength:            %.4f (should be high but not 1.0)\n", corr_rank_strength))
cat(sprintf("  CV marker expression vs CL_rank:   %.4f (expected: negative)\n", corr_cv))
cat(sprintf("  PV marker expression vs CL_rank:   %.4f (expected: positive)\n", corr_pn))

cat("\nZone Mapping Statistics:\n")
cat(sprintf("  Zone entropy:    Mean=%.3f, Median=%.3f, Range=[%.3f, %.3f]\n",
            mean(entropy), median(entropy), min(entropy), max(entropy)))
cat(sprintf("  Zone confidence: Mean=%.3f, Median=%.3f, Range=[%.3f, %.3f]\n",
            mean(confidence), median(confidence), min(confidence), max(confidence)))
cat(sprintf("  Max probability: Mean=%.3f, Median=%.3f, Range=[%.3f, %.3f]\n",
            mean(max_prob), median(max_prob), min(max_prob), max(max_prob)))
cat(sprintf("  Hard vs Soft agreement: %.1f%%\n", 100 * agreement))

cat("\nMarker Genes Used:\n")
cat(sprintf("  CV markers: %s\n", paste(available_cv, collapse = ", ")))
cat(sprintf("  PN markers: %s\n", paste(available_pn, collapse = ", ")))

cat("\n=============================================================\n")
cat(sprintf("  Debug plots saved to: %s\n", output_dir))
cat("=============================================================\n")
cat("\nKey Differences:\n")
cat("  - CL_rank: Uniform distribution, for zoning/bin comparison\n")
cat("  - CL_strength: Normal-like distribution, shows confidence\n")
cat("  - Zone confidence: Derived from probability distribution entropy\n")
cat("=============================================================\n")
