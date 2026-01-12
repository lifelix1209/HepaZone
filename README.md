# HepaZone

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R 4.0+](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)

**HepaZone** is an R package for reconstructing liver spatial transcriptomics profiles from scRNA-seq data. It identifies regionalized marker genes and reconstructs gene expression distribution along the liver lobule spatial axis (from central vein to portal vein), with support for circadian rhythm analysis.

## Key Features

- **Spatial Zonation Reconstruction**: Map hepatocytes to spatial layers along the porto-central axis
- **Circadian Rhythm Analysis**: Analyze time × space interactions across different ZT time points
- **Statistical Testing**: Bootstrap standard errors and permutation tests for significance
- **Seurat Integration**: Built on the Seurat ecosystem for seamless single-cell analysis
- **Flexible Export**: Export results to CSV, TSV, or RDS formats

## Background

Liver zonation refers to the spatial organization of hepatocytes along the porto-central axis, where cells exhibit distinct metabolic functions based on their position. HepaZone enables researchers to:

- Map zonation patterns in healthy and diseased liver tissue
- Track zonation changes during development or pathology
- Identify zone-specific cell populations
- Study circadian rhythm effects on spatial gene expression

## Installation

```r
# Install development version from GitHub
devtools::install_github("yourusername/HepaZone")

# Install from source
install.packages("path/to/HepaZone_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Dependencies

- R >= 4.0.0
- Seurat >= 4.0.0
- Matrix >= 1.3.0
- dplyr >= 1.0.0
- ggplot2 >= 3.3.0
- GEOquery >= 2.60.0
- qvalue >= 2.24.0

## Quick Start

### Single Sample Analysis

```r
library(HepaZone)

# Run complete zonation analysis
result <- hepa_zone_reconstruct(
  seurat_obj = "./data/ZT00",  # Path to 10x output
  from_10x = TRUE,
  sample.name = "ZT00",
  n_zones = 10,
  n_bootstrap = 500,
  n_permutations = 1000
)

# View summary
print(result)

# Plot spatial expression heatmap
plot_spatial_heatmap(result, genes = c("Alb", "Cyp2e1", "Hamp"))

# Plot gene expression gradient
plot_gradient(result, genes = c("Alb", "Cyp2e1", "Hamp"))

# Export results
export_hepa_zone_results(result, output.dir = "./results")
```

### Multi-Time Point Analysis (Circadian Rhythm)

```r
# Define paths for each time point
sample_paths <- list(
  ZT00 = "./data/ZT00",
  ZT06 = "./data/ZT06",
  ZT12 = "./data/ZT12",
  ZT18 = "./data/ZT18"
)

# Run analysis for all time points
multi_result <- run_hepa_zone_multi(
  sample.paths = sample_paths,
  n_zones = 10
)

# Compare expression across time points
plot_time_comparison(
  multi_result$results,
  gene = "Alb",
  time_points = c("ZT00", "ZT06", "ZT12", "ZT18")
)

# View genes with significant time × space interaction
interaction_genes <- multi_result$interaction$high_interaction_genes
```

### Step-by-Step Analysis

```r
# Load and preprocess data
seurat_obj <- load_10x_data("./data/ZT00", sample.name = "ZT00")
seurat_obj <- preprocess_zonation(seurat_obj)

# Calculate spatial position scores
seurat_obj <- calculate_spatial_position(seurat_obj)

# Map cells to spatial layers
prob_matrix <- map_cells_to_layers(seurat_obj, n_zones = 10)

# Reconstruct spatial expression profiles
spatial_result <- reconstruct_spatial_expression(seurat_obj, prob_matrix)

# Bootstrap for standard errors
boot_result <- bootstrap_se(
  expression_matrix = seurat_obj@assays$RNA@data,
  prob_matrix = prob_matrix,
  n_bootstrap = 500
)

# Permutation test
perm_result <- permutation_test(
  expression_matrix = seurat_obj@assays$RNA@data,
  prob_matrix = prob_matrix,
  n_permutations = 1000
)

# Calculate q-values
qvals <- calculate_qvalues(perm_result$p_values)
```

## Output

The main output `HepaZoneResult` object contains:

| Component | Description |
|-----------|-------------|
| `mean_expression` | Gene × Zone matrix of mean expression |
| `standard_error` | Gene × Zone matrix of standard errors |
| `qvalues` | Gene × Zone matrix of q-values |
| `prob_matrix` | Cell × Zone probability matrix |
| `cl_scores` | CL score for each cell |
| `svg_results` | Data frame of spatially variable genes |

## Reference Data

This package uses marker genes from:
- Bahar Halpern & Shenhav et al., Nature 2017
- GEO accession: GSE145197 (mouse liver scRNA-seq)

## Citation

If you use HepaZone in your research, please cite:

```
HepaZone: Liver Spatial Transcriptomics Reconstruction R Package
```

## License

MIT License
