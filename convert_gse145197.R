#!/usr/bin/env Rscript
# Convert GSE145197 TXT format to 10x MTX format
# This script converts the downloaded GEO data to Seurat-compatible format

library(Matrix)

cat("=============================================================\n")
cat("  Converting GSE145197 to 10x Format\n")
cat("=============================================================\n\n")

# Set paths
raw_dir <- "GSE145197_data/GSE145197/suppl"
output_base <- "GSE145197_data/10x_format"

# Create output directory
if (!dir.exists(output_base)) {
  dir.create(output_base, recursive = TRUE)
}

# Sample mapping (time point -> files)
samples <- list(
  ZT00 = c("GSM4308343_UMI_tab_ZT00A.txt.gz", "GSM4308344_UMI_tab_ZT00B.txt.gz", "GSM4308345_UMI_tab_ZT00C.txt.gz"),
  ZT06 = c("GSM4308346_UMI_tab_ZT06A.txt.gz", "GSM4308347_UMI_tab_ZT06B.txt.gz"),
  ZT12 = c("GSM4308348_UMI_tab_ZT12A.txt.gz", "GSM4308349_UMI_tab_ZT12B.txt.gz", "GSM4308350_UMI_tab_ZT12C.txt.gz"),
  ZT18 = c("GSM4308351_UMI_tab_ZT18A.txt.gz", "GSM4308352_UMI_tab_ZT18B.txt.gz")
)

convert_to_10x <- function(input_file, output_dir, sample_name) {
  cat(sprintf("Converting %s -> %s...\n", input_file, sample_name))

  # Read header line to get cell barcodes
  con <- gzfile(input_file, "rt")
  header <- readLines(con, n = 1)
  close(con)

  # Parse header - first element is "gene_names", rest are cell barcodes
  header_parts <- strsplit(header, "\t")[[1]]
  cell_barcodes <- header_parts[-1]  # Remove "gene_names"
  n_cells <- length(cell_barcodes)

  # Read gene names from first column (skip header)
  gene_names <- read.table(gzfile(input_file),
                           sep = "\t",
                           header = FALSE,
                           skip = 1,
                           colClasses = "character",
                           check.names = FALSE,
                           nrows = -1)[, 1]

  # Make gene names unique
  gene_names <- make.unique(gene_names, sep = "_")

  # Read the data (genes x cells matrix) - skip header
  data <- read.table(gzfile(input_file),
                     sep = "\t",
                     header = FALSE,
                     skip = 1,
                     check.names = FALSE,
                     nrows = -1,
                     colClasses = c("character", rep("integer", n_cells)))

  # Remove first column (gene names) and convert to matrix
  # This is genes (rows) x cells (cols) - standard 10x format
  data_matrix <- as.matrix(data[, -1, drop = FALSE])
  mode(data_matrix) <- "integer"
  rownames(data_matrix) <- gene_names

  # Create sparse matrix (genes x cells)
  sparse_mat <- Matrix(data_matrix, sparse = TRUE)

  # Get gene names and barcodes
  genes <- gene_names
  barcodes <- cell_barcodes

  # Create output directory
  sample_dir <- file.path(output_dir, sample_name)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE)
  }

  # Write MTX file
  Matrix::writeMM(sparse_mat, file.path(sample_dir, "matrix.mtx"))

  # Write genes file
  gene_df <- data.frame(gene_id = genes, gene_name = genes)
  write.table(gene_df, file.path(sample_dir, "genes.tsv"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Write barcodes file
  writeLines(barcodes, file.path(sample_dir, "barcodes.tsv"))

  cat(sprintf("  Created: %s (%d cells x %d genes)\n",
              sample_name, nrow(data), ncol(data)))

  invisible(nrow(data))
}

# Convert each time point
for (time_point in names(samples)) {
  cat(sprintf("\n>>> Processing %s:\n", time_point))

  output_dir <- file.path(output_base, time_point)
  total_cells <- 0

  for (file in samples[[time_point]]) {
    input_path <- file.path(raw_dir, file)
    if (file.exists(input_path)) {
      # Extract sample name from filename
      sample_name <- gsub("GSM430834[0-9]+_UMI_tab_", "", file)
      sample_name <- gsub(".txt.gz", "", sample_name)
      cells <- convert_to_10x(input_path, output_dir, sample_name)
      total_cells <- total_cells + cells
    } else {
      cat(sprintf("  Warning: %s not found\n", file))
    }
  }

  cat(sprintf("  Total cells for %s: %d\n", time_point, total_cells))
}

cat("\n=============================================================\n")
cat("  Conversion Complete!\n")
cat("=============================================================\n\n")
cat(sprintf("Output directory: %s\n", output_base))
cat("Files can now be loaded with load_10x_data()\n\n")
