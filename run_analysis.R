# ==============================================================================
# HepaZone - Quick Run Script
# ==============================================================================
# This script provides a simple interface to run the complete HepaZone analysis.
# Usage: Rscript run_analysis.R [data_path] [n_zones]
# ==============================================================================

# Default settings
data_path <- if (interactive()) "GSE145197_data/10x_format/ZT00/ZT00A" else commandArgs(TRUE)[1]
n_zones <- if (interactive()) 10 else as.integer(commandArgs(TRUE)[2])

if (is.na(n_zones)) n_zones <- 10

# Load required packages
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

# Source all HepaZone functions
source("R/data_io.R")
source("R/preprocessing.R")
source("R/zonation.R")
source("R/statistics.R")
source("R/visualization.R")
source("R/main.R")

cat("=============================================================\n")
cat("  HepaZone Analysis - Quick Run\n")
cat("=============================================================\n")
cat(sprintf("  Data: %s\n", data_path))
cat(sprintf("  Zones: %d\n\n", n_zones))

# Check data exists
if (!dir.exists(data_path)) {
  cat("ERROR: Data directory not found!\n")
  quit(status = 1)
}

# Load data
cat("Loading data...\n")
seurat_obj <- load_10x_data(
  path = data_path,
  sample.name = "ZT00",
  min.cells = 10,
  min.features = 200
)
cat(sprintf("  Loaded %d cells x %d genes\n\n", ncol(seurat_obj), nrow(seurat_obj)))

# Run complete analysis
cat("Running HepaZone analysis...\n")
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

# Quick summary
cat("\n=============================================================\n")
cat("  Results Summary\n")
cat("=============================================================\n")
cat(sprintf("  Genes: %d\n", nrow(result$mean_expression)))
cat(sprintf("  Zones: %d\n", result$n_zones))
cat(sprintf("  Significant SVGs (q<0.05): %d\n\n", sum(result$svg_results$significant, na.rm = TRUE)))

# Export results
output_dir <- "GSE145197_results/quick_run"
dir.create(output_dir, recursive = TRUE)
export_hepa_zone_results(result, output.dir = output_dir)

cat(sprintf("Results saved to: %s\n", output_dir))
cat("=============================================================\n")
