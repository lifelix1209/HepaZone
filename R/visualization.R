#==============================================================================
# Visualization Functions
# ==============================================================================

#' Plot Spatial Expression Heatmap
#'
#' Creates a heatmap of gene expression across spatial zones.
#'
#' @param hepa_result HepaZone result object
#' @param genes Character vector of genes to plot (default: top 20 by variance)
#' @param n_genes Number of top genes to show if genes is NULL (default: 20)
#' @param cluster_genes Whether to cluster genes by expression pattern (default: TRUE)
#' @param show_zone_labels Whether to show zone labels on x-axis (default: TRUE)
#' @param title Plot title (default: "Spatial Gene Expression")
#' @return A ggplot object
#' @examples
#' \dontrun{
#' plot_heatmap(hepa_result, genes = c("Alb", "Cyp2e1", "Hamp"))
#' }
#' @export
plot_spatial_heatmap <- function(hepa_result,
                                   genes = NULL,
                                   n_genes = 20,
                                   cluster_genes = TRUE,
                                   show_zone_labels = TRUE,
                                   title = "Spatial Gene Expression") {

  # Get expression matrix
  expr_mat <- hepa_result$mean_expression

  # Select genes if not provided
  if (is.null(genes)) {
    # Select top variable genes
    gene_vars <- apply(expr_mat, 1, var)
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(n_genes, length(gene_vars))]
    genes <- top_genes
  }

  # Subset expression matrix
  expr_subset <- expr_mat[genes, , drop = FALSE]

  # Convert to long format for ggplot
  expr_df <- as.data.frame(expr_subset)
  expr_df$gene <- rownames(expr_subset)
  expr_long <- tidyr::pivot_longer(
    expr_df,
    cols = -gene,
    names_to = "zone",
    values_to = "expression"
  )

  # Create heatmap
  p <- ggplot2::ggplot(expr_long, ggplot2::aes(x = zone, y = gene, fill = expression)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "Expression") +
    ggplot2::labs(title = title, x = "Spatial Zone", y = "Gene") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = if (!show_zone_labels) ggplot2::element_blank() else ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  return(p)
}


#' Plot Gene Expression Gradient
#'
#' Creates a line plot showing gene expression across spatial zones.
#'
#' @param hepa_result HepaZone result object
#' @param genes Character vector of genes to plot
#' @param se_matrix Standard error matrix (from bootstrap) (optional)
#' @param add_ci Add confidence intervals if se_matrix provided (default: TRUE)
#' @param colors Color palette for genes (default: ggplot2 hue)
#' @param title Plot title (default: "Gene Expression Along Spatial Axis")
#' @return A ggplot object
#' @examples
#' \dontrun{
#' plot_gradient(hepa_result, genes = c("Alb", "Cyp2e1", "Hamp"))
#' }
#' @export
plot_gradient <- function(hepa_result,
                           genes,
                           se_matrix = NULL,
                           add_ci = TRUE,
                           colors = NULL,
                           title = "Gene Expression Along Spatial Axis") {

  expr_mat <- hepa_result$mean_expression

  # Validate genes
  genes <- genes[genes %in% rownames(expr_mat)]
  if (length(genes) == 0) {
    stop("None of the specified genes found in the dataset")
  }

  # Prepare data frame
  zones <- 1:hepa_result$n_zones
  expr_df <- data.frame()

  for (g in genes) {
    gene_expr <- expr_mat[g, ]
    gene_df <- data.frame(
      gene = g,
      zone = zones,
      expression = as.numeric(gene_expr)
    )
    expr_df <- rbind(expr_df, gene_df)
  }

  # Add standard errors if provided
  if (!is.null(se_matrix) && add_ci) {
    se_df <- data.frame()
    for (g in genes) {
      if (g %in% rownames(se_matrix)) {
        gene_se <- se_matrix[g, ]
        gene_se_df <- data.frame(
          gene = g,
          zone = zones,
          se = as.numeric(gene_se)
        )
        se_df <- rbind(se_df, gene_se_df)
      }
    }
    expr_df <- merge(expr_df, se_df, by = c("gene", "zone"))
  }

  # Create plot
  p <- ggplot2::ggplot(expr_df, ggplot2::aes(x = zone, y = expression, color = gene)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = zones, labels = paste0("Z", zones)) +
    ggplot2::labs(title = title, x = "Spatial Zone", y = "Normalized Expression") +
    ggplot2::theme_minimal()

  # Add confidence intervals if available
  if ("se" %in% colnames(expr_df)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = expression - 1.96 * se, ymax = expression + 1.96 * se,
                   fill = gene),
      alpha = 0.2,
      color = NA
    )
  }

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }

  return(p)
}


#' Plot CL Score Distribution
#'
#' Visualizes the distribution of cells across spatial zones based on CL scores.
#'
#' @param seurat_obj Seurat object with CL_score
#' @param prob_matrix Cell-to-zone probability matrix (optional)
#' @param fill_color Fill color for histogram (default: "steelblue")
#' @return A ggplot object
#' @export
plot_cl_distribution <- function(seurat_obj, prob_matrix = NULL,
                                  fill_color = "steelblue") {

  cl_scores <- seurat_obj$CL_score

  p <- ggplot2::ggplot(data.frame(CL = cl_scores), ggplot2::aes(x = CL)) +
    ggplot2::geom_histogram(bins = 50, fill = fill_color, color = "white", alpha = 0.8) +
    ggplot2::labs(title = "Distribution of CL Scores (Spatial Position)",
                  x = "CL Score (0=CV, 1=PV)",
                  y = "Number of Cells") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # Add zone boundaries if probability matrix provided
  if (!is.null(prob_matrix)) {
    n_zones <- ncol(prob_matrix)
    zone_centers <- seq(0.05, 0.95, length.out = n_zones)

    for (z in zone_centers) {
      p <- p + ggplot2::geom_vline(xintercept = z, linetype = "dashed", color = "red", alpha = 0.5)
    }
  }

  return(p)
}


#' Plot Spatial Marker Expression
#'
#' Shows the expression of CV and PN marker genes across cells,
#' colored by their CL score.
#'
#' @param seurat_obj Seurat object with CV and PN scores
#' @param cv_genes CV marker genes to highlight
#' @param pn_genes PN marker genes to highlight
#' @param n_cells Number of cells to sample (default: 2000)
#' @return A ggplot object
#' @export
plot_marker_expression <- function(seurat_obj,
                                    cv_genes = NULL,
                                    pn_genes = NULL,
                                    n_cells = 2000) {

  if (is.null(cv_genes)) cv_genes <- .default_cv_markers()
  if (is.null(pn_genes)) pn_genes <- .default_pn_markers()

  # Sample cells if too many
  if (ncol(seurat_obj) > n_cells) {
    set.seed(42)
    cell_idx <- sample(1:ncol(seurat_obj), n_cells)
    plot_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cell_idx])
  } else {
    plot_obj <- seurat_obj
  }

  # Create data frame for plotting
  plot_df <- data.frame(
    cv_score = plot_obj$cv_score,
    pn_score = plot_obj$pn_score,
    CL_score = plot_obj$CL_score
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = cv_score, y = pn_score, color = CL_score)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.5) +
    ggplot2::scale_color_viridis_c(name = "CL Score") +
    ggplot2::labs(title = "CV vs PN Marker Expression",
                  x = "CV Score (Central Vein)",
                  y = "PN Score (Portal Vein)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  return(p)
}


#' Plot Time Point Comparison
#'
#' Compares spatial expression patterns across multiple time points.
#'
#' @param hepa_result_list List of HepaZone results by time point
#' @param gene Gene to compare across time points
#' @param time_points Time point names to include (default: all)
#' @param colors Color palette for time points
#' @return A ggplot object
#' @export
plot_time_comparison <- function(hepa_result_list,
                                   gene,
                                   time_points = NULL,
                                   colors = NULL) {

  if (is.null(time_points)) {
    time_points <- names(hepa_result_list)
  }

  # Collect expression data
  expr_df <- data.frame()
  for (tp in time_points) {
    if (tp %in% names(hepa_result_list)) {
      result <- hepa_result_list[[tp]]
      if (gene %in% rownames(result$mean_expression)) {
        gene_expr <- result$mean_expression[gene, ]
        tp_df <- data.frame(
          time_point = tp,
          zone = 1:result$n_zones,
          expression = as.numeric(gene_expr)
        )
        expr_df <- rbind(expr_df, tp_df)
      }
    }
  }

  p <- ggplot2::ggplot(expr_df, ggplot2::aes(x = zone, y = expression, color = time_point, group = time_point)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = 1:max(expr_df$zone), labels = paste0("Z", 1:max(expr_df$zone))) +
    ggplot2::labs(title = sprintf("Spatial Expression of %s Across Time Points", gene),
                  x = "Spatial Zone",
                  y = "Normalized Expression",
                  color = "Time Point") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }

  return(p)
}


#' Plot Zonation Summary
#'
#' Creates a comprehensive summary figure for a HepaZone analysis.
#'
#' @param seurat_obj Seurat object after zonation analysis
#' @param hepa_result HepaZone result object
#' @param top_n Number of top variable genes to show (default: 10)
#' @return A patchwork or grid plot
#' @export
plot_zonation_summary <- function(seurat_obj, hepa_result, top_n = 10) {

  # CL score distribution
  p1 <- plot_cl_distribution(seurat_obj)

  # Top genes heatmap
  p2 <- plot_spatial_heatmap(hepa_result, n_genes = top_n,
                              show_zone_labels = TRUE,
                              title = paste0("Top ", top_n, " Spatially Variable Genes"))

  # Top genes gradient
  gene_vars <- apply(hepa_result$mean_expression, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(6, length(gene_vars))]
  p3 <- plot_gradient(hepa_result, genes = top_genes,
                      title = "Top Spatially Variable Genes")

  # Combine plots
  # Using base R for simple layout
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

  print(p1)
  print(p2)
  print(p3)

  par(old_par)

  message("Summary plots displayed.")
}
