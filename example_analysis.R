# ==============================================================================
# HepaZone Analysis Example
# Using GSE145197 real liver single-cell RNA-seq data
# ==============================================================================
# Data source: GEO accession GSE145197
# Reference: Droin et al. Nat Metab 2021
# Description: Molecular characterization of liver hepatocytes at different
#              circadian time points, enabling spatial transcriptomics reconstruction
#
# IMPORTANT: This script requires GSE145197 data. You can either:
#   1. Download automatically (requires GEOquery package)
#   2. Download manually and load from local files
# ==============================================================================

library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

cat("=============================================================\n")
cat("  HepaZone Analysis - GSE145197 Real Data\n")
cat("=============================================================\n\n")

# ==============================================================================
# Data Loading Options
# ==============================================================================
cat(">>> Loading data...\n\n")

# Try to load from local 10x files first (recommended if data is already downloaded)
# Or use automatic download

# Option A: Load from local directory (if you already have the data)
# Uncomment and modify the following line:
# seurat_obj <- load_10x_data(path = "/path/to/GSE145197_raw/ZT00", sample.name = "ZT00")

# Option B: Automatic download from GEO (if GEOquery is available)
USE_AUTO_DOWNLOAD <- FALSE

if (USE_AUTO_DOWNLOAD) {
  # Check if GEOquery is available
  if (requireNamespace("GEOquery", quietly = TRUE)) {
    # Load HepaZone functions
    source("R/data_io.R")
    source("R/preprocessing.R")
    source("R/zonation.R")
    source("R/statistics.R")
    source("R/visualization.R")
    source("R/main.R")

    # Download and load data
    seurat_obj <- load_gse145197(
      destdir = "GSE145197_data",
      time_point = "ZT00",
      min.cells = 10,
      min.features = 200
    )
    cat("\nData loaded from GEO\n")
  } else {
    cat("GEOquery not available. Please install it or use manual data loading.\n")
    cat("To install: install.packages('GEOquery')\n")
    quit(status = 1)
  }
} else {
  # Option C: Manual data loading
  # You must have downloaded GSE145197 from GEO first

  # Load HepaZone functions
  source("R/data_io.R")
  source("R/preprocessing.R")
  source("R/zonation.R")
  source("R/statistics.R")
  source("R/visualization.R")
  source("R/main.R")

  # Check for local data directory (use converted 10x format)
  data_path <- "GSE145197_data/10x_format/ZT00/ZT00A"

  if (dir.exists(data_path)) {
    cat("Found local data at:", data_path, "\n")
    seurat_obj <- load_10x_data(
      path = data_path,
      sample.name = "ZT00",
      min.cells = 10,
      min.features = 200
    )
  } else {
    cat("ERROR: Local data not found!\n\n")
    cat("=============================================================\n")
    cat("  DATA NOT FOUND\n")
    cat("=============================================================\n\n")
    cat("Please download GSE145197 data manually:\n\n")
    cat("1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145197\n\n")
    cat("2. Download the supplementary file:\n")
    cat("   - Click on 'GSE145197_RAW.tar' (19.8 MB)\n\n")
    cat("3. Extract the archive:\n")
    cat("   tar -xvf GSE145197_RAW.tar\n\n")
    cat("4. The extracted folders contain 10x Genomics data:\n")
    cat("   - ZT00/, ZT06/, ZT12/, ZT18/ (time points)\n\n")
    cat("5. Create the directory structure:\n")
    cat("   mkdir -p GSE145197_data/GSE145197/suppl/ZT00\n")
    cat("   cp -r ZT00/* GSE145197_data/GSE145197/suppl/ZT00/\n\n")
    cat("6. Re-run this script\n\n")
    cat("Alternatively, set USE_AUTO_DOWNLOAD <- TRUE at line 49\n")
    cat("to attempt automatic download (requires GEOquery package).\n")
    cat("=============================================================\n")
    quit(status = 1)
  }
}

cat("\nSeurat object info:\n")
cat(sprintf("  - Genes: %d\n", nrow(seurat_obj)))
cat(sprintf("  - Cells: %d\n", ncol(seurat_obj)))
cat(sprintf("  - Sample: %s\n", seurat_obj@project.name))

# ==============================================================================
# 2. Preprocess Data
# ==============================================================================
cat("\n>>> Step 2: Preprocessing data\n\n")

seurat_obj <- preprocess_zonation(
  seurat_obj,
  mt_pattern = "^Mt-",
  remove_mup = TRUE,
  mup_pattern = "^Mup"
)

# ==============================================================================
# 3. Define Spatial Marker Genes
# ==============================================================================
cat("\n>>> Step 3: Defining spatial marker genes\n\n")

# Central vein (CV) marker genes - Zone 1 (peri-central)
cv_markers <- c(
  "Cyp2e1", "Cyp1a2", "Cyp2f2", "Cyp2a5", "Hamp",
  "Gulo", "Hyou1", "Acsl1", "Mfsd2a", "Sds"
)

# Portal vein (PN) marker genes - Zone 10 (peri-portal)
pn_markers <- c(
  "Alb", "Ass1", "Fgb", "Fgg", "Apoa1",
  "Mup3", "Mup20", "C3", "Serpina1a", "F10"
)

# Check which markers are in the dataset
cv_found <- cv_markers[cv_markers %in% rownames(seurat_obj)]
pn_found <- pn_markers[pn_markers %in% rownames(seurat_obj)]

cat("CV markers found:", length(cv_found), "/", length(cv_markers), "\n")
cat("  ", paste(cv_found, collapse = ", "), "\n")
cat("PN markers found:", length(pn_found), "/", length(pn_markers), "\n")
cat("  ", paste(pn_found, collapse = ", "), "\n\n")

# If few markers found, use data-driven marker selection
if (length(cv_found) < 3 || length(pn_found) < 3) {
  cat("Insufficient known markers. Using data-driven marker selection...\n")
  mat_norm <- .get_data(seurat_obj)
  gene_means <- rowMeans(mat_norm)
  top_genes <- names(sort(gene_means, decreasing = TRUE))[1:100]

  cv_found <- top_genes[1:10]
  pn_found <- top_genes[11:20]
  cat("Using top expressed genes as markers:\n")
  cat("  CV markers:", paste(cv_found, collapse = ", "), "\n")
  cat("  PN markers:", paste(pn_found, collapse = ", "), "\n\n")
}

# ==============================================================================
# 4. Calculate Spatial Position (CL Score)
# ==============================================================================
cat("\n>>> Step 4: Calculating CL Score (spatial position)\n\n")

seurat_obj <- calculate_spatial_position(
  seurat_obj,
  cv_markers = cv_found,
  pn_markers = pn_found
)

cat("CL Score statistics:\n")
cat(sprintf("  - Mean: %.3f\n", mean(seurat_obj$CL_score, na.rm = TRUE)))
cat(sprintf("  - SD: %.3f\n", sd(seurat_obj$CL_score, na.rm = TRUE)))
cat(sprintf("  - Range: %.3f - %.3f\n",
            min(seurat_obj$CL_score, na.rm = TRUE),
            max(seurat_obj$CL_score, na.rm = TRUE)))

# ==============================================================================
# 5. Map Cells to Spatial Layers
# ==============================================================================
cat("\n>>> Step 5: Mapping cells to spatial layers\n\n")

n_zones <- 10
prob_matrix <- map_cells_to_layers(seurat_obj, n_zones = n_zones)

cat("Probability matrix:", nrow(prob_matrix), "x", ncol(prob_matrix), "\n")
cat("Cells per zone:\n")
zone_counts <- colSums(prob_matrix)
for (z in 1:n_zones) {
  cat(sprintf("  Zone %d: %.1f cells (%.1f%%)\n",
              z, zone_counts[z], 100 * zone_counts[z] / sum(zone_counts)))
}

# ==============================================================================
# 6. Run Complete Analysis Pipeline
# ==============================================================================
cat("\n>>> Step 6: Running HepaZone reconstruction\n\n")

result <- hepa_zone_reconstruct(
  seurat_obj,
  cv_markers = cv_found,
  pn_markers = pn_found,
  n_zones = n_zones,
  n_bootstrap = 200,
  n_permutations = 500,
  seed = 42,
  verbose = FALSE
)

cat("\nAnalysis results:\n")
cat(sprintf("  - Genes analyzed: %d\n", nrow(result$mean_expression)))
cat(sprintf("  - Spatial layers: %d\n", result$n_zones))
cat(sprintf("  - Significant SVGs (q<0.05): %d\n",
            sum(result$svg_results$significant, na.rm = TRUE)))
cat(sprintf("  - Significant SVGs (q<0.10): %d\n",
            sum(result$svg_results$q_value < 0.10, na.rm = TRUE)))
cat(sprintf("  - Significant SVGs (q<0.20): %d\n",
            sum(result$svg_results$q_value < 0.20, na.rm = TRUE)))

# ==============================================================================
# 7. Analyze Discovered SVGs
# ==============================================================================
cat("\n>>> Step 7: Analyzing discovered spatial variable genes\n\n")

# Get significant SVGs
sig_results <- result$svg_results[result$svg_results$significant, ]
sig_results <- sig_results[order(sig_results$q_value), ]

# Classify SVGs by expression pattern
n_z <- result$n_zones
cv_idx <- 1:floor(n_z / 3)
pn_idx <- (ceiling(2 * n_z / 3) + 1):n_z

# Calculate mean expression in CV and PN zones
cv_expr <- apply(result$mean_expression[sig_results$gene, cv_idx, drop = FALSE], 1, mean)
pn_expr <- apply(result$mean_expression[sig_results$gene, pn_idx, drop = FALSE], 1, mean)

sig_results$cv_mean <- cv_expr
sig_results$pn_mean <- pn_expr
sig_results$direction <- ifelse(cv_expr > pn_expr, "CV-enriched", "PN-enriched")

# Count by direction
n_cv_enriched <- sum(sig_results$direction == "CV-enriched")
n_pn_enriched <- sum(sig_results$direction == "PN-enriched")

cat("SVG classification:\n")
cat(sprintf("  - CV-enriched (central vein): %d genes\n", n_cv_enriched))
cat(sprintf("  - PN-enriched (portal vein): %d genes\n", n_pn_enriched))

# Show top SVGs
cat("\nTop 20 Spatial Variable Genes:\n")
for (i in 1:min(20, nrow(sig_results))) {
  g <- sig_results$gene[i]
  qv <- sig_results$q_value[i]
  dir <- sig_results$direction[i]
  cv_m <- sig_results$cv_mean[i]
  pn_m <- sig_results$pn_mean[i]
  cat(sprintf("  %2d. %-15s q=%.4f %-12s (CV:%.2f, PN:%.2f)\n",
              i, g, qv, dir, cv_m, pn_m))
}

# ==============================================================================
# 8. Visualization
# ==============================================================================
cat("\n>>> Step 8: Generating visualizations\n\n")

output_dir <- "GSE145197_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 8.1 CL Score distribution
cat("  - CL Score distribution...\n")
p1 <- plot_cl_distribution(seurat_obj, prob_matrix = prob_matrix,
                          title = "GSE145197: Cell Distribution Along Spatial Axis")
ggsave(file.path(output_dir, "cl_score_distribution.png"), p1, width = 8, height = 6)

# 8.2 Top SVG heatmap
cat("  - Top SVG heatmap...\n")
top_svg_genes <- head(sig_results$gene, 25)
p2 <- plot_spatial_heatmap(result, genes = top_svg_genes, n_genes = 25,
                           title = "Top 25 Spatially Variable Genes (GSE145197)")
ggsave(file.path(output_dir, "top_svg_heatmap.png"), p2, width = 12, height = 10)

# 8.3 Expression gradients
cat("  - Expression gradients...\n")
# Select diverse genes: some CV-enriched, some PN-enriched
cv_genes_plot <- head(sig_results$gene[sig_results$direction == "CV-enriched"], 3)
pn_genes_plot <- head(sig_results$gene[sig_results$direction == "PN-enriched"], 3)
plot_genes <- c(cv_genes_plot, pn_genes_plot)

if (length(plot_genes) >= 1) {
  n_cv <- length(cv_genes_plot)
  n_pn <- length(pn_genes_plot)
  colors_vec <- c(rep("#E64B35", n_cv), rep("#4DBBD5", n_pn))
  p3 <- plot_gradient(result, genes = plot_genes,
                     title = "Expression Gradients of SVGs",
                     colors = colors_vec)
  ggsave(file.path(output_dir, "svg_expression_gradients.png"), p3, width = 10, height = 7)
}

# 8.4 SVG classification barplot
cat("  - SVG classification...\n")
svg_class_df <- data.frame(
  Gene = sig_results$gene,
  Direction = sig_results$direction
)

p4 <- ggplot2::ggplot(svg_class_df, ggplot2::aes(x = Direction, fill = Direction)) +
  ggplot2::geom_bar(color = "white") +
  ggplot2::scale_fill_manual(values = c("CV-enriched" = "#E64B35", "PN-enriched" = "#4DBBD5")) +
  ggplot2::labs(title = "Classification of Significant SVGs",
                x = "Expression Pattern", y = "Number of Genes") +
  .theme_publication() +
  ggplot2::theme(legend.position = "none")

ggsave(file.path(output_dir, "svg_classification.png"), p4, width = 6, height = 5)

# 8.5 Marker gene validation
cat("  - Marker gene validation...\n")
# Plot known CV and PN markers if available
known_cv <- cv_found[1:min(5, length(cv_found))]
known_pn <- pn_found[1:min(5, length(pn_found))]
marker_genes <- c(known_cv, known_pn)
marker_genes <- marker_genes[marker_genes %in% rownames(result$mean_expression)]

if (length(marker_genes) >= 1) {
  n_cv_m <- sum(known_cv %in% rownames(result$mean_expression))
  n_pn_m <- sum(known_pn %in% rownames(result$mean_expression))
  colors_markers <- c(rep("#E64B35", n_cv_m), rep("#4DBBD5", n_pn_m))
  p5 <- plot_gradient(result, genes = marker_genes,
                     title = "Known Marker Gene Validation",
                     colors = colors_markers)
  ggsave(file.path(output_dir, "marker_validation.png"), p5, width = 10, height = 7)
}

# ==============================================================================
# 9. Export Results
# ==============================================================================
cat("\n>>> Step 9: Exporting results\n\n")

# Export complete results
export_hepa_zone_results(result, output.dir = output_dir)

# Export detailed SVG list
cat("  - Exporting detailed SVG list...\n")
svg_export <- sig_results
svg_export$cv_mean <- NULL
svg_export$pn_mean <- NULL

for (z in 1:result$n_zones) {
  svg_export[[paste0("Zone_", z)]] <- result$mean_expression[svg_export$gene, z]
}

write.csv(svg_export, file.path(output_dir, "significant_svg_list.csv"))
cat("    Saved to:", file.path(output_dir, "significant_svg_list.csv"), "\n")

# ==============================================================================
# Summary
# ==============================================================================
cat("\n=============================================================\n")
cat("  Analysis Complete!\n")
cat("=============================================================\n\n")

cat("Data source: GEO GSE145197 (Droin et al. Nat Metab 2021)\n")
cat(sprintf("Cells analyzed: %d\n", ncol(seurat_obj)))
cat(sprintf("Genes analyzed: %d\n", nrow(seurat_obj)))
cat(sprintf("Spatial layers: %d\n\n", result$n_zones))

cat("Significant Spatial Variable Genes:\n")
cat(sprintf("  - Total (q<0.05): %d\n", sum(result$svg_results$significant, na.rm = TRUE)))
cat(sprintf("  - CV-enriched: %d\n", n_cv_enriched))
cat(sprintf("  - PN-enriched: %d\n\n", n_pn_enriched))

cat("Output directory:", output_dir, "\n\n")
cat("Generated files:\n")
cat("  1. cl_score_distribution.png   - Cell distribution histogram\n")
cat("  2. top_svg_heatmap.png         - Top SVG heatmap\n")
cat("  3. svg_expression_gradients.png - SVG expression gradients\n")
cat("  4. svg_classification.png      - SVG classification barplot\n")
cat("  5. marker_validation.png       - Known marker validation\n")
cat("  6. mean_expression.csv         - Spatial expression matrix\n")
cat("  7. qvalues.csv                 - Q-value matrix\n")
cat("  8. standard_error.csv          - Standard error matrix\n")
cat("  9. significant_svg_list.csv    - Detailed SVG list\n")
cat(" 10. hepa_zone_full_results.rds  - Complete results object\n\n")

cat("=============================================================\n")
