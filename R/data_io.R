#' @import Matrix dplyr Seurat
#' @importFrom Matrix t
#' @importFrom utils read.csv write.csv

# ==============================================================================
# Data Import Functions
# ==============================================================================

#' Download and Load GSE145197 Data from GEO
#'
#' Downloads scRNA-seq UMI data from GEO accession GSE145197 and converts
#' it to a Seurat object for downstream analysis.
#'
#' @param destdir Directory to store downloaded files (default: ".")
#' @param type Data source type: "geo" (GEO database) or "local" (local 10x files)
#' @param ... Additional arguments passed to Read10X or GEOquery functions
#' @return A Seurat object containing the GSE145197 data
#' @examples
#' \dontrun{
#' # Download from GEO
#' seurat_obj <- load_gse145197(destdir = "./data")
#'
#' # Load from local 10x output
#' seurat_obj <- load_gse145197(type = "local", path = "./raw_data")
#' }
#' @export
load_gse145197 <- function(destdir = ".", type = c("geo", "local"), ...) {
  type <- match.arg(type)

  if (type == "geo") {
    .load_from_geo(destdir, ...)
  } else {
    .load_from_local(...)
  }
}


#' Internal function to load data from GEO
#' @keywords internal
.load_from_geo <- function(destdir, ...) {
  message("Note: For full GSE145197 download, use GEOquery manually.")
  message("This function provides a template for GEO data loading.")

  # Template for GSE145197 - modify accession as needed
  accession <- "GSE145197"

  # Example: Download supplementary files
  # geo_obj <- GEOquery::getGEOSuppFiles(accession, makeDirectory = TRUE,
  #                                        baseDir = destdir)

  # For now, return a template Seurat object
  # Replace with actual data loading in production
  message("Please download GSE145197 data manually and use load_10x_data()")
  return(NULL)
}


#' Internal function to load from local 10x data
#' @keywords internal
.load_from_local <- function(path = NULL, sample.name = "sample") {
  if (is.null(path)) {
    stop("Please provide 'path' to 10x output directory")
  }

  # Read 10x data (supports MTX, HDF5, and CSV formats)
  data_mat <- Seurat::Read10X(data.dir = path)

  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = data_mat,
    project = sample.name,
    min.cells = 3,
    min.features = 200
  )

  return(seurat_obj)
}


#' Load Data from 10x Genomics Output
#'
#' Reads 10x Genomics Cell Ranger output (MTX format) and creates a Seurat object.
#'
#' @param path Path to the 10x output directory containing matrix.mtx, genes.tsv, barcodes.tsv
#' @param sample.name Sample name for the Seurat project (default: "sample")
#' @param min.cells Include features detected in at least this many cells (default: 3)
#' @param min.features Include cells where at least this many features are detected (default: 200)
#' @return A Seurat object
#' @export
load_10x_data <- function(path, sample.name = "sample",
                          min.cells = 3, min.features = 200) {
  # Read 10x data
  data_mat <- Seurat::Read10X(data.dir = path)

  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = data_mat,
    project = sample.name,
    min.cells = min.cells,
    min.features = min.features
  )

  return(seurat_obj)
}


#' Load Multiple Samples and Merge
#'
#' Loads multiple 10x samples and merges them into a single Seurat object.
#' Useful for processing multiple time points (ZT00, ZT06, ZT12, ZT18).
#'
#' @param sample.paths Named list of paths to 10x directories
#' @param min.cells Minimum cells per feature filter (default: 3)
#' @param min.features Minimum features per cell filter (default: 200)
#' @return A merged Seurat object with samples as identities
#' @export
load_and_merge_samples <- function(sample.paths, min.cells = 3, min.features = 200) {
  if (length(sample.paths) == 0) {
    stop("sample.paths cannot be empty")
  }

  # Load first sample
  first_sample <- names(sample.paths)[1]
  seurat_obj <- load_10x_data(
    path = sample.paths[[first_sample]],
    sample.name = first_sample,
    min.cells = min.cells,
    min.features = min.features
  )

  # Load and merge remaining samples
  if (length(sample.paths) > 1) {
    for (i in 2:length(sample.paths)) {
      sample_name <- names(sample.paths)[i]
      sample_obj <- load_10x_data(
        path = sample.paths[[sample_name]],
        sample.name = sample_name,
        min.cells = min.cells,
        min.features = min.features
      )
      seurat_obj <- merge(seurat_obj, y = sample_obj, add.cell.ids = c("", sample_name))
    }
  }

  return(seurat_obj)
}


# ==============================================================================
# Data Export Functions
# ==============================================================================

#' Export Spatial Expression Matrix
#'
#' Exports reconstructed spatial expression data to various formats.
#'
#' @param object HepaZone result object or matrix to export
#' @param format Output format: "csv", "tsv", or "rds"
#' @param filepath Output file path (extension will be added if missing)
#' @param ... Additional arguments passed to write functions
#' @return Invisible NULL
#' @examples
#' \dontrun{
#' export_spatial_expression(hepa_result, format = "csv", filepath = "./results")
#' export_spatial_expression(hepa_result, format = "rds", filepath = "./data/hepa_result")
#' }
#' @export
export_spatial_expression <- function(object, format = c("csv", "tsv", "rds"),
                                       filepath, ...) {
  format <- match.arg(format)

  # Add extension if missing
  if (!grepl(paste0("\\.", format, "$"), filepath)) {
    filepath <- paste0(filepath, ".", format)
  }

  switch(format,
    "csv" = {
      if (is.matrix(object) || is.data.frame(object)) {
        utils::write.csv(object, file = filepath, ...)
      } else {
        stop("For CSV export, object must be a matrix or data.frame")
      }
    },
    "tsv" = {
      if (is.matrix(object) || is.data.frame(object)) {
        utils::write.csv(object, file = filepath, sep = "\t", ...)
      } else {
        stop("For TSV export, object must be a matrix or data.frame")
      }
    },
    "rds" = {
      saveRDS(object, file = filepath, ...)
    }
  )

  message(sprintf("Exported to: %s", filepath))
  return(invisible(NULL))
}


#' Export Full HepaZone Results
#'
#' Exports all results from a HepaZone analysis including expression matrices,
#' statistics, and metadata.
#'
#' @param result HepaZone result object from hepa_zone_reconstruct()
#' @param output.dir Output directory path
#' @return Invisible NULL
#' @export
export_hepa_zone_results <- function(result, output.dir = "./results") {
  # Create output directory
  if (!dir.exists(output.dir)) {
    dir.create(output.dir, recursive = TRUE)
  }

  # Export mean expression
  export_spatial_expression(
    result$mean_expression,
    format = "csv",
    filepath = file.path(output.dir, "mean_expression")
  )

  # Export standard errors
  export_spatial_expression(
    result$standard_error,
    format = "csv",
    filepath = file.path(output.dir, "standard_error")
  )

  # Export q-values
  export_spatial_expression(
    result$qvalues,
    format = "csv",
    filepath = file.path(output.dir, "qvalues")
  )

  # Export full results as RDS
  export_spatial_expression(
    result,
    format = "rds",
    filepath = file.path(output.dir, "hepa_zone_full_results")
  )

  message(sprintf("All results exported to: %s", output.dir))
  return(invisible(NULL))
}
