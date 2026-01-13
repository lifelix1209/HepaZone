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
#' @param destdir Directory to store downloaded files (default: "data")
#' @param time_point Time point to load: "ZT00", "ZT06", "ZT12", "ZT18", or "all"
#'                    (default: "ZT00")
#' @param min.cells Include features detected in at least this many cells (default: 10)
#' @param min.features Include cells where at least this many features are detected (default: 200)
#' @param redownload Force re-download of data (default: FALSE)
#' @return A Seurat object (or list of Seurat objects if time_point="all")
#' @examples
#' \dontrun{
#' # Load ZT00 sample
#' seurat_obj <- load_gse145197(destdir = "./data", time_point = "ZT00")
#'
#' # Load all time points
#' all_samples <- load_gse145197(destdir = "./data", time_point = "all")
#' }
#' @export
load_gse145197 <- function(destdir = "data",
                           time_point = c("ZT00", "ZT06", "ZT12", "ZT18", "all"),
                           min.cells = 10,
                           min.features = 200,
                           redownload = FALSE) {

  time_point <- match.arg(time_point)
  geo_info <- NULL

  # Check if data already exists
  raw_dir <- file.path(destdir, "GSE145197", "suppl")

  if (!redownload && dir.exists(raw_dir)) {
    message("Using existing GSE145197 data")
    sample_dirs <- list.dirs(raw_dir, full.names = TRUE, recursive = FALSE)
    geo_info <- list(
      accession = "GSE145197",
      sample_dirs = sample_dirs,
      raw_dir = raw_dir
    )
  } else {
    # Download from GEO
    geo_info <- .load_from_geo(destdir)
  }

  if (is.null(geo_info) || length(geo_info$sample_dirs) == 0) {
    stop("No sample data available")
  }

  # Select time point
  if (time_point == "all") {
    selected_dirs <- geo_info$sample_dirs
  } else {
    # Find directories matching the time point pattern
    zt_pattern <- paste0(".*", time_point, ".*")
    selected_dirs <- grep(zt_pattern, geo_info$sample_dirs, value = TRUE, ignore.case = TRUE)

    if (length(selected_dirs) == 0) {
      # Try exact match
      selected_dirs <- grep(time_point, geo_info$sample_dirs, value = TRUE, ignore.case = TRUE)
    }

    if (length(selected_dirs) == 0) {
      # Use first available sample
      message(sprintf("Time point '%s' not found. Using first available sample.", time_point))
      selected_dirs <- geo_info$sample_dirs[1]
    }
  }

  # Load samples
  if (length(selected_dirs) == 1) {
    # Single sample - return Seurat object
    sample_dir <- selected_dirs[1]
    sample_name <- basename(sample_dir)

    message(sprintf("Loading sample: %s", sample_name))
    seurat_obj <- load_10x_data(
      path = sample_dir,
      sample.name = sample_name,
      min.cells = min.cells,
      min.features = min.features
    )
    return(seurat_obj)

  } else {
    # Multiple samples - return list of Seurat objects
    message(sprintf("Loading %d time points...", length(selected_dirs)))
    result_list <- list()

    for (i in seq_along(selected_dirs)) {
      sample_dir <- selected_dirs[i]
      sample_name <- basename(sample_dir)

      message(sprintf("  Loading %s...", sample_name))
      result_list[[sample_name]] <- load_10x_data(
        path = sample_dir,
        sample.name = sample_name,
        min.cells = min.cells,
        min.features = min.features
      )
    }

    return(result_list)
  }
}


#' Internal function to load data from GEO
#' @keywords internal
.load_from_geo <- function(destdir = ".", ...) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Package 'GEOquery' is required. Install with: install.packages('GEOquery')")
  }

  accession <- "GSE145197"
  message(sprintf("Downloading GSE145197 data from GEO to: %s", destdir))

  # Download supplementary files
  geo_obj <- GEOquery::getGEOSuppFiles(accession,
                                        makeDirectory = TRUE,
                                        baseDir = destdir,
                                        filter_regex = NULL)

  # Find and extract the RAW.tar file
  raw_dir <- file.path(destdir, accession, "suppl")
  raw_tar <- list.files(raw_dir, pattern = "RAW.tar$", full.names = TRUE)

  if (length(raw_tar) == 0) {
    stop("RAW.tar file not found in downloaded files")
  }

  message(sprintf("Extracting %s", basename(raw_tar[1])))
  untar(raw_tar[1], exdir = raw_dir)

  # Find sample directories (ZT00, ZT06, etc.)
  sample_dirs <- list.dirs(raw_dir, full.names = TRUE, recursive = FALSE)
  message(sprintf("Found %d sample directories", length(sample_dirs)))

  if (length(sample_dirs) == 0) {
    stop("No sample directories found in RAW data")
  }

  # Return list of sample paths for further processing
  return(list(
    accession = accession,
    sample_dirs = sample_dirs,
    raw_dir = raw_dir
  ))
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
