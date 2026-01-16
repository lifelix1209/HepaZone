#==============================================================================
# Professional Publication-Grade Visualization Functions
#==============================================================================

#' Professional Color Palettes
#' @keywords internal
.hepazone_palettes <- function() {
  list(
    # Nature-style palette
    nature = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", 
               "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"),
    
    # Science-style palette  
    science = c("#3B4992", "#EE0000", "#008B45", "#631879", "#008280",
                "#BB0021", "#5F559B", "#A20056", "#808180", "#1B1919"),
    
    # Cell-style palette
    cell = c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442",
             "#56B4E9", "#E69F00", "#000000", "#999999", "#654321"),
    
    # Zonation-specific gradient
    zonation = c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C")
  )
}


#' Publication Theme
#' @keywords internal
.theme_publication <- function(base_size = 11, base_family = "Arial") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Grid and background
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 1),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.25),
      
      # Axes
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.5),
      axis.text = ggplot2::element_text(colour = "black", size = ggplot2::rel(0.9)),
      axis.title = ggplot2::element_text(colour = "black", size = ggplot2::rel(1.0), 
                                         face = "bold"),
      
      # Legend
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key.size = ggplot2::unit(1.2, "lines"),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.9), face = "bold"),
      
      # Title
      plot.title = ggplot2::element_text(size = ggplot2::rel(1.2), face = "bold",
                                         hjust = 0.5, margin = ggplot2::margin(b = 10)),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1.0), hjust = 0.5,
                                           margin = ggplot2::margin(b = 10)),
      
      # Strip (for facets)
      strip.background = ggplot2::element_rect(fill = "grey90", colour = "black",
                                              linewidth = 0.5),
      strip.text = ggplot2::element_text(size = ggplot2::rel(0.9), face = "bold")
    )
}


#' Plot Spatial Expression Heatmap (Enhanced)
#'
#' Creates a publication-quality heatmap of gene expression across spatial zones.
#'
#' @param hepa_result HepaZone result object
#' @param genes Character vector of genes to plot (default: top 20 by variance)
#' @param n_genes Number of top genes to show if genes is NULL (default: 20)
#' @param cluster_genes Whether to cluster genes hierarchically (default: TRUE)
#' @param cluster_zones Whether to cluster zones hierarchically (default: FALSE)
#' @param scale_method Scaling method: "row", "column", or "none" (default: "row")
#' @param color_palette Color palette: "RdYlBu", "viridis", "magma", or custom (default: "RdYlBu")
#' @param show_values Show expression values in cells (default: FALSE)
#' @param show_zone_labels Whether to show zone labels (default: TRUE)
#' @param title Plot title (default: "Spatial Gene Expression")
#' @param fontsize_row Row label font size (default: 10)
#' @param fontsize_col Column label font size (default: 11)
#' @param use_ggplot Force using ggplot2 instead of ComplexHeatmap (default: TRUE)
#' @return A ggplot object
#' @export
plot_spatial_heatmap <- function(hepa_result,
                                  genes = NULL,
                                  n_genes = 20,
                                  cluster_genes = TRUE,
                                  cluster_zones = FALSE,
                                  scale_method = "row",
                                  color_palette = "RdYlBu",
                                  show_values = FALSE,
                                  show_zone_labels = TRUE,
                                  title = "Spatial Gene Expression",
                                  fontsize_row = 10,
                                  fontsize_col = 11,
                                  use_ggplot = TRUE) {
  
  # Get expression matrix
  expr_mat <- hepa_result$mean_expression
  
  # Select genes
  if (is.null(genes)) {
    gene_vars <- apply(expr_mat, 1, var)
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(n_genes, nrow(expr_mat))]
    genes <- top_genes
  }
  
  # Validate genes
  genes <- genes[genes %in% rownames(expr_mat)]
  if (length(genes) == 0) {
    stop("None of the specified genes found in the dataset")
  }
  
  # Subset and scale expression matrix
  expr_subset <- expr_mat[genes, , drop = FALSE]
  
  if (scale_method == "row") {
    expr_scaled <- t(scale(t(expr_subset)))
  } else if (scale_method == "column") {
    expr_scaled <- scale(expr_subset)
  } else {
    expr_scaled <- expr_subset
  }
  
  # Prepare zone labels - force CV to PV order (Zone 1 to Zone n)
  n_zones <- ncol(expr_scaled)
  zone_labels <- if (show_zone_labels) {
    paste0("Zone ", 1:n_zones)
  } else {
    paste0("Z", 1:n_zones)
  }

  # Ensure columns are in CV->PV order (Zone 1 to n)
  # Set column names and create ordered factor for ggplot version
  colnames(expr_scaled) <- zone_labels
  
  # Use ggplot2 by default (use_ggplot = TRUE for compatibility with ggsave)
  use_complex <- requireNamespace("ComplexHeatmap", quietly = TRUE) && !use_ggplot

  if (use_complex) {

    # Define color scheme
    if (color_palette == "RdYlBu") {
      col_fun <- circlize::colorRamp2(
        c(min(expr_scaled, na.rm = TRUE),
          0,
          max(expr_scaled, na.rm = TRUE)),
        c("#0571B0", "#F7F7F7", "#CA0020")
      )
    } else if (color_palette == "viridis") {
      col_fun <- circlize::colorRamp2(
        seq(min(expr_scaled, na.rm = TRUE), max(expr_scaled, na.rm = TRUE), length.out = 100),
        viridisLite::viridis(100)
      )
    } else if (color_palette == "magma") {
      col_fun <- circlize::colorRamp2(
        seq(min(expr_scaled, na.rm = TRUE), max(expr_scaled, na.rm = TRUE), length.out = 100),
        viridisLite::magma(100)
      )
    } else {
      col_fun <- color_palette
    }

    # Create heatmap annotation - CV to PV gradient
    ha_top <- ComplexHeatmap::HeatmapAnnotation(
      Zone = factor(zone_labels, levels = zone_labels),
      col = list(Zone = circlize::colorRamp2(
        c(1, ncol(expr_scaled)),
        c("#2C7BB6", "#D7191C")  # Blue (CV) to Red (PV)
      )),
      annotation_name_gp = grid::gpar(fontsize = fontsize_col, fontface = "bold"),
      annotation_legend_param = list(
        title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 9)
      )
    )

    # Create heatmap with CV->PV column order (disable column clustering)
    ht <- ComplexHeatmap::Heatmap(
      expr_scaled,
      name = if (scale_method == "none") "Expression" else "Scaled\nExpression",
      col = col_fun,

      # Clustering - DISABLE column clustering to preserve CV->PV order
      cluster_rows = cluster_genes,
      cluster_columns = FALSE,  # Force CV to PV order
      show_row_dend = cluster_genes,
      show_column_dend = FALSE,

      # Row settings
      row_names_gp = grid::gpar(fontsize = fontsize_row),
      row_names_side = "left",
      row_dend_width = grid::unit(15, "mm"),

      # Column settings - explicit column order
      column_order = zone_labels,  # Force CV->PV order
      column_names_gp = grid::gpar(fontsize = fontsize_col, fontface = "bold"),
      column_names_rot = 45,
      column_dend_height = grid::unit(15, "mm"),
      top_annotation = ha_top,

      # Cell settings
      cell_fun = if (show_values) {
        function(j, i, x, y, width, height, fill) {
          grid::grid.text(sprintf("%.2f", expr_scaled[i, j]),
                         x, y, gp = grid::gpar(fontsize = 8))
        }
      } else {
        NULL
      },

      # Border
      border = TRUE,
      rect_gp = grid::gpar(col = "white", lwd = 1),

      # Legend
      heatmap_legend_param = list(
        title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 10),
        legend_height = grid::unit(40, "mm"),
        legend_direction = "vertical"
      ),

      # Title
      column_title = title,
      column_title_gp = grid::gpar(fontsize = 14, fontface = "bold")
    )

    return(ht)
    
  } else {
    # Fallback to ggplot2
    expr_df <- as.data.frame(expr_scaled)
    expr_df$gene <- factor(rownames(expr_scaled), levels = rownames(expr_scaled))

    # Pivot to long format - ensure zone is an ordered factor (CV to PV)
    expr_long <- tidyr::pivot_longer(
      expr_df,
      cols = -gene,
      names_to = "zone",
      values_to = "expression"
    )

    # Convert zone to ordered factor to preserve CV->PV order (Zone 1 to n)
    zone_order <- paste0(if (show_zone_labels) "Zone " else "Z", 1:n_zones)
    expr_long$zone <- factor(expr_long$zone, levels = zone_order)

    # Order genes by clustering if requested
    if (cluster_genes) {
      gene_dist <- dist(expr_scaled)
      gene_hclust <- hclust(gene_dist)
      expr_long$gene <- factor(expr_long$gene, levels = rownames(expr_scaled)[gene_hclust$order])
    }
    
    # Select color palette
    if (color_palette == "RdYlBu") {
      scale_fill <- ggplot2::scale_fill_gradient2(
        low = "#0571B0", mid = "#F7F7F7", high = "#CA0020",
        midpoint = 0, name = if (scale_method == "none") "Expression" else "Scaled\nExpression"
      )
    } else if (color_palette == "viridis") {
      scale_fill <- ggplot2::scale_fill_viridis_c(
        option = "viridis",
        name = if (scale_method == "none") "Expression" else "Scaled\nExpression"
      )
    } else {
      scale_fill <- ggplot2::scale_fill_viridis_c(
        option = "magma",
        name = if (scale_method == "none") "Expression" else "Scaled\nExpression"
      )
    }
    
    p <- ggplot2::ggplot(expr_long, ggplot2::aes(x = zone, y = gene, fill = expression)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      scale_fill +
      ggplot2::labs(title = title, x = "Spatial Zone", y = NULL) +
      .theme_publication(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = ggplot2::element_text(size = fontsize_row),
        legend.position = "right"
      )
    
    if (show_values) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", expression)),
        size = 2.5, color = "black"
      )
    }
    
    return(p)
  }
}


#' Plot Gene Expression Gradient (Enhanced)
#'
#' Creates a publication-quality line plot showing gene expression across spatial zones.
#'
#' @param hepa_result HepaZone result object
#' @param genes Character vector of genes to plot
#' @param se_matrix Standard error matrix (optional)
#' @param add_ci Add 95% confidence intervals (default: TRUE)
#' @param smooth Add smooth spline curve (default: FALSE)
#' @param colors Color palette: "nature", "science", "cell", or custom vector
#' @param alpha_ribbon Transparency for confidence ribbon (default: 0.2)
#' @param line_size Line width (default: 1.2)
#' @param point_size Point size (default: 3)
#' @param title Plot title
#' @param subtitle Plot subtitle (optional)
#' @param show_legend Show legend (default: TRUE)
#' @param legend_position Legend position (default: "right")
#' @return A ggplot object
#' @export
plot_gradient <- function(hepa_result,
                          genes,
                          se_matrix = NULL,
                          add_ci = TRUE,
                          smooth = FALSE,
                          colors = "nature",
                          alpha_ribbon = 0.2,
                          line_size = 1.2,
                          point_size = 3,
                          title = "Gene Expression Along Spatial Axis",
                          subtitle = NULL,
                          show_legend = TRUE,
                          legend_position = "right") {
  
  expr_mat <- hepa_result$mean_expression
  
  # Validate genes
  genes <- genes[genes %in% rownames(expr_mat)]
  if (length(genes) == 0) {
    stop("None of the specified genes found in the dataset")
  }
  
  # Prepare data
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
  
  # Add standard errors
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
    expr_df <- merge(expr_df, se_df, by = c("gene", "zone"), all.x = TRUE)
    expr_df$lower <- expr_df$expression - 1.96 * expr_df$se
    expr_df$upper <- expr_df$expression + 1.96 * expr_df$se
  }
  
  # Set colors
  if (length(colors) == 1 && colors %in% names(.hepazone_palettes())) {
    color_values <- .hepazone_palettes()[[colors]][1:length(genes)]
  } else if (length(colors) >= length(genes)) {
    color_values <- colors[1:length(genes)]
  } else {
    color_values <- .hepazone_palettes()$nature[1:length(genes)]
  }
  
  # Create base plot
  p <- ggplot2::ggplot(expr_df, ggplot2::aes(x = zone, y = expression, 
                                              color = gene, fill = gene, group = gene))
  
  # Add confidence intervals
  if ("se" %in% colnames(expr_df)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = alpha_ribbon,
      color = NA
    )
  }
  
  # Add lines and points
  if (smooth) {
    p <- p + ggplot2::geom_smooth(
      method = "loess", 
      se = FALSE, 
      linewidth = line_size,
      span = 0.5
    )
  } else {
    p <- p + ggplot2::geom_line(linewidth = line_size)
  }
  
  p <- p + ggplot2::geom_point(size = point_size, shape = 21, 
                               color = "white", stroke = 0.8)
  
  # Scales and labels
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = zones, 
      labels = paste0("Z", zones),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_color_manual(
      values = color_values,
      name = "Gene"
    ) +
    ggplot2::scale_fill_manual(
      values = color_values,
      name = "Gene"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Spatial Zone (CV → PV)",
      y = "Normalized Expression (AU)"
    ) +
    .theme_publication(base_size = 11) +
    ggplot2::theme(
      legend.position = if (show_legend) legend_position else "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}


#' Plot CL Score Distribution (Enhanced)
#'
#' Publication-quality visualization of cell distribution across spatial zones.
#'
#' @param seurat_obj Seurat object with CL_score
#' @param prob_matrix Cell-to-zone probability matrix (optional)
#' @param fill_color Fill color (default: "#3B4992")
#' @param add_density Add density curve overlay (default: TRUE)
#' @param add_rug Add rug plot (default: TRUE)
#' @param n_bins Number of histogram bins (default: 50)
#' @param show_zones Show zone boundaries (default: TRUE)
#' @param title Plot title
#' @return A ggplot object
#' @export
plot_cl_distribution <- function(seurat_obj, 
                                  prob_matrix = NULL,
                                  fill_color = "#3B4992",
                                  add_density = TRUE,
                                  add_rug = TRUE,
                                  n_bins = 50,
                                  show_zones = TRUE,
                                  title = "Distribution of Cells Along Spatial Axis") {
  
  cl_scores <- seurat_obj$CL_score
  
  # Remove NA values
  cl_scores <- cl_scores[!is.na(cl_scores)]
  
  plot_df <- data.frame(CL = cl_scores)
  
  # Create histogram with density
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = CL))
  
  # Add histogram
  if (add_density) {
    p <- p + ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = n_bins,
      fill = fill_color,
      color = "white",
      alpha = 0.7,
      linewidth = 0.3
    )
    
    # Add density curve
    p <- p + ggplot2::geom_density(
      color = "#D55E00",
      linewidth = 1.2,
      fill = "#D55E00",
      alpha = 0.1
    )
  } else {
    p <- p + ggplot2::geom_histogram(
      bins = n_bins,
      fill = fill_color,
      color = "white",
      alpha = 0.8,
      linewidth = 0.3
    )
  }
  
  # Add rug plot
  if (add_rug) {
    p <- p + ggplot2::geom_rug(alpha = 0.3, linewidth = 0.5)
  }
  
  # Add zone boundaries
  if (show_zones && !is.null(prob_matrix)) {
    n_zones <- ncol(prob_matrix)
    zone_boundaries <- seq(0, 1, length.out = n_zones + 1)
    
    for (i in 2:n_zones) {
      p <- p + ggplot2::geom_vline(
        xintercept = zone_boundaries[i],
        linetype = "dashed",
        color = "#E64B35",
        linewidth = 0.8,
        alpha = 0.6
      )
    }
    
    # Add zone labels
    zone_centers <- seq(0, 1, length.out = n_zones + 1)[-1] - 1/(2*n_zones)
    for (i in 1:n_zones) {
      p <- p + ggplot2::annotate(
        "text",
        x = zone_centers[i],
        y = Inf,
        label = paste0("Z", i),
        vjust = 1.5,
        size = 3.5,
        fontface = "bold",
        color = "#666666"
      )
    }
  }
  
  # Formatting
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      labels = c("0\n(CV)", "0.2", "0.4", "0.6", "0.8", "1.0\n(PV)"),
      expand = c(0.01, 0.01)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      title = title,
      x = "CL Score (Central Vein ← → Portal Vein)",
      y = if (add_density) "Density" else "Number of Cells"
    ) +
    .theme_publication(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}


#' Plot Spatial Marker Expression (Enhanced)
#'
#' Publication-quality scatter plot of CV vs PN marker expression.
#'
#' @param seurat_obj Seurat object with scores
#' @param cv_genes CV marker genes
#' @param pn_genes PN marker genes
#' @param n_cells Number of cells to plot (default: all if < 5000)
#' @param plot_type Type: "scatter", "hex", or "density" (default: "hex")
#' @param color_palette Color palette for CL score (default: "viridis")
#' @param point_size Point size for scatter plot (default: 0.8)
#' @param alpha Point transparency (default: 0.5)
#' @param add_contours Add density contours (default: TRUE)
#' @param title Plot title
#' @return A ggplot object
#' @export
plot_marker_expression <- function(seurat_obj,
                                    cv_genes = NULL,
                                    pn_genes = NULL,
                                    n_cells = NULL,
                                    plot_type = "hex",
                                    color_palette = "viridis",
                                    point_size = 0.8,
                                    alpha = 0.5,
                                    add_contours = TRUE,
                                    title = "Spatial Marker Expression") {
  
  if (is.null(cv_genes)) cv_genes <- .default_cv_markers()
  if (is.null(pn_genes)) pn_genes <- .default_pn_markers()
  
  # Sample cells
  n_total <- ncol(seurat_obj)
  if (is.null(n_cells)) {
    n_cells <- min(n_total, 5000)
  }
  
  if (n_total > n_cells) {
    set.seed(42)
    cell_idx <- sample(1:n_total, n_cells)
    plot_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cell_idx])
  } else {
    plot_obj <- seurat_obj
  }
  
  # Prepare data
  plot_df <- data.frame(
    cv_score = plot_obj$cv_score,
    pn_score = plot_obj$pn_score,
    CL_score = plot_obj$CL_score
  )
  
  # Remove NA
  plot_df <- plot_df[complete.cases(plot_df), ]
  
  # Create base plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = cv_score, y = pn_score))
  
  # Select plot type
  if (plot_type == "hex") {
    p <- p + ggplot2::geom_hex(ggplot2::aes(fill = ggplot2::after_stat(count)), bins = 50) +
      ggplot2::scale_fill_viridis_c(
        option = "plasma",
        name = "Cell\nDensity",
        trans = "log10"
      )
  } else if (plot_type == "density") {
    p <- p + ggplot2::stat_density_2d(
      ggplot2::aes(fill = ggplot2::after_stat(density)),
      geom = "raster",
      contour = FALSE
    ) +
      ggplot2::scale_fill_viridis_c(
        option = "plasma",
        name = "Density"
      )
  } else {
    # Scatter plot with CL score coloring
    p <- p + ggplot2::geom_point(
      ggplot2::aes(color = CL_score),
      size = point_size,
      alpha = alpha
    )
    
    if (color_palette == "viridis") {
      p <- p + ggplot2::scale_color_viridis_c(
        option = "viridis",
        name = "CL Score"
      )
    } else {
      p <- p + ggplot2::scale_color_gradientn(
        colors = .hepazone_palettes()$zonation,
        name = "CL Score"
      )
    }
  }
  
  # Add contours
  if (add_contours && plot_type %in% c("scatter", "hex")) {
    p <- p + ggplot2::stat_density_2d(
      color = "black",
      linewidth = 0.3,
      alpha = 0.5,
      bins = 6
    )
  }
  
  # Add diagonal reference line
  p <- p + ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.8
  )
  
  # Formatting
  p <- p +
    ggplot2::labs(
      title = title,
      x = "CV Score (Central Vein Markers)",
      y = "PN Score (Portal Vein Markers)"
    ) +
    .theme_publication(base_size = 11) +
    ggplot2::coord_fixed() +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}


#' Plot Time Point Comparison (Enhanced)
#'
#' Publication-quality comparison of spatial patterns across time points.
#'
#' @param hepa_result_list List of HepaZone results by time point
#' @param gene Gene to compare
#' @param time_points Time point names (default: all)
#' @param colors Color palette
#' @param add_ci Add confidence intervals if available (default: TRUE)
#' @param facet_by_time Facet by time point instead of overlay (default: FALSE)
#' @param line_size Line width (default: 1.2)
#' @param point_size Point size (default: 3)
#' @param title Plot title
#' @return A ggplot object
#' @export
plot_time_comparison <- function(hepa_result_list,
                                  gene,
                                  time_points = NULL,
                                  colors = "science",
                                  add_ci = TRUE,
                                  facet_by_time = FALSE,
                                  line_size = 1.2,
                                  point_size = 3,
                                  title = NULL) {
  
  if (is.null(time_points)) {
    time_points <- names(hepa_result_list)
  }
  
  # Collect data
  expr_df <- data.frame()
  
  for (tp in time_points) {
    if (!tp %in% names(hepa_result_list)) next
    
    result <- hepa_result_list[[tp]]
    
    if (!gene %in% rownames(result$mean_expression)) next
    
    gene_expr <- result$mean_expression[gene, ]
    
    tp_df <- data.frame(
      time_point = tp,
      zone = 1:result$n_zones,
      expression = as.numeric(gene_expr)
    )
    
    # Add SE if available
    if (!is.null(result$se_expression) && gene %in% rownames(result$se_expression)) {
      gene_se <- result$se_expression[gene, ]
      tp_df$se <- as.numeric(gene_se)
      tp_df$lower <- tp_df$expression - 1.96 * tp_df$se
      tp_df$upper <- tp_df$expression + 1.96 * tp_df$se
    }
    
    expr_df <- rbind(expr_df, tp_df)
  }
  
  if (nrow(expr_df) == 0) {
    stop("No data found for the specified gene and time points")
  }
  
  # Set colors
  if (length(colors) == 1 && colors %in% names(.hepazone_palettes())) {
    color_values <- .hepazone_palettes()[[colors]][1:length(time_points)]
  } else {
    color_values <- colors
  }
  
  # Order time points
  expr_df$time_point <- factor(expr_df$time_point, levels = time_points)
  
  # Create plot
  p <- ggplot2::ggplot(expr_df, 
                       ggplot2::aes(x = zone, y = expression, 
                                   color = time_point, fill = time_point,
                                   group = time_point))
  
  # Add confidence intervals
  if ("se" %in% colnames(expr_df) && add_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.15,
      color = NA
    )
  }
  
  # Add lines and points
  p <- p +
    ggplot2::geom_line(linewidth = line_size) +
    ggplot2::geom_point(size = point_size, shape = 21, 
                       color = "white", stroke = 0.8)
  
  # Scales
  p <- p +
    ggplot2::scale_x_continuous(
      breaks = unique(expr_df$zone),
      labels = paste0("Z", unique(expr_df$zone))
    ) +
    ggplot2::scale_color_manual(
      values = color_values,
      name = "Time Point"
    ) +
    ggplot2::scale_fill_manual(
      values = color_values,
      name = "Time Point"
    )
  
  # Faceting
  if (facet_by_time) {
    p <- p + ggplot2::facet_wrap(~time_point, nrow = 1)
  }
  
  # Labels
  if (is.null(title)) {
    title <- sprintf("Spatial Expression of %s Across Time Points", gene)
  }
  
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("n = %d time points", length(unique(expr_df$time_point))),
      x = "Spatial Zone",
      y = "Normalized Expression (AU)"
    ) +
    .theme_publication(base_size = 11) +
    ggplot2::theme(
      legend.position = if (facet_by_time) "none" else "right",
      panel.grid.major.x = ggplot2::element_blank()
    )
  
  return(p)
}


#' Plot Zonation Summary (Enhanced)
#'
#' Creates a comprehensive publication-quality summary figure.
#'
#' @param seurat_obj Seurat object after zonation
#' @param hepa_result HepaZone result object
#' @param top_n Number of top genes (default: 10)
#' @param layout Layout: "grid" or "patchwork" (default: "patchwork")
#' @param output_file Output file path (optional, for saving)
#' @param width Figure width in inches (default: 14)
#' @param height Figure height in inches (default: 10)
#' @param dpi Resolution (default: 300)
#' @return Combined plot object
#' @export
plot_zonation_summary <- function(seurat_obj,
                                   hepa_result,
                                   top_n = 10,
                                   layout = "patchwork",
                                   output_file = NULL,
                                   width = 14,
                                   height = 10,
                                   dpi = 300) {
  
  # Panel A: CL score distribution
  p1 <- plot_cl_distribution(
    seurat_obj,
    title = "A. Cell Distribution Along Spatial Axis",
    add_density = TRUE,
    add_rug = FALSE
  )
  
  # Panel B: Marker expression
  p2 <- plot_marker_expression(
    seurat_obj,
    plot_type = "hex",
    title = "B. CV vs PN Marker Expression"
  )
  
  # Panel C: Top genes heatmap
  gene_vars <- apply(hepa_result$mean_expression, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(top_n, length(gene_vars))]
  
  p3 <- plot_spatial_heatmap(
    hepa_result,
    genes = top_genes,
    cluster_genes = TRUE,
    title = sprintf("C. Top %d Spatially Variable Genes", top_n),
    fontsize_row = 9
  )
  
  # Panel D: Top genes gradient
  top_genes_gradient <- names(sort(gene_vars, decreasing = TRUE))[1:min(6, length(gene_vars))]
  
  p4 <- plot_gradient(
    hepa_result,
    genes = top_genes_gradient,
    title = "D. Expression Gradients",
    colors = "nature",
    smooth = FALSE
  )
  
  # Combine plots
  if (layout == "patchwork" && requireNamespace("patchwork", quietly = TRUE)) {
    
    # For ComplexHeatmap object, convert to grob
    if (inherits(p3, "Heatmap")) {
      p3_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(p3))
      p3 <- ggplot2::ggplot() + 
        ggplot2::annotation_custom(p3_grob) +
        ggplot2::theme_void()
    }
    
    combined_plot <- (p1 | p2) / (p3 | p4) +
      patchwork::plot_annotation(
        title = "HepaZone Spatial Analysis Summary",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5)
        )
      )
    
  } else {
    # Fallback to grid layout
    combined_plot <- gridExtra::grid.arrange(
      p1, p2, p3, p4,
      ncol = 2,
      top = grid::textGrob(
        "HepaZone Spatial Analysis Summary",
        gp = grid::gpar(fontsize = 16, fontface = "bold")
      )
    )
  }
  
  # Save if requested
  if (!is.null(output_file)) {
    ggplot2::ggsave(
      filename = output_file,
      plot = combined_plot,
      width = width,
      height = height,
      dpi = dpi,
      units = "in"
    )
    message(sprintf("Summary plot saved to: %s", output_file))
  }
  
  return(combined_plot)
}


#' Save Publication Figure
#'
#' Helper function to save figures in publication-ready formats.
#'
#' @param plot ggplot or ComplexHeatmap object
#' @param filename Output filename (without extension)
#' @param formats Vector of formats: "pdf", "png", "tiff", "eps" (default: c("pdf", "png"))
#' @param width Width in inches (default: 7)
#' @param height Height in inches (default: 5)
#' @param dpi Resolution for raster formats (default: 300)
#' @param output_dir Output directory (default: current directory)
#' @export
save_publication_figure <- function(plot,
                                     filename,
                                     formats = c("pdf", "png"),
                                     width = 7,
                                     height = 5,
                                     dpi = 300,
                                     output_dir = ".") {
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (fmt in formats) {
    output_path <- file.path(output_dir, paste0(filename, ".", fmt))
    
    if (inherits(plot, "Heatmap") || inherits(plot, "HeatmapList")) {
      # ComplexHeatmap object
      if (fmt == "pdf") {
        pdf(output_path, width = width, height = height)
        ComplexHeatmap::draw(plot)
        dev.off()
      } else if (fmt == "png") {
        png(output_path, width = width * dpi, height = height * dpi, res = dpi)
        ComplexHeatmap::draw(plot)
        dev.off()
      } else if (fmt == "tiff") {
        tiff(output_path, width = width * dpi, height = height * dpi, res = dpi,
             compression = "lzw")
        ComplexHeatmap::draw(plot)
        dev.off()
      }
    } else {
      # ggplot object
      ggplot2::ggsave(
        filename = output_path,
        plot = plot,
        width = width,
        height = height,
        dpi = dpi,
        units = "in",
        device = fmt
      )
    }
    
    message(sprintf("Figure saved: %s", output_path))
  }
  
  invisible(NULL)
}
