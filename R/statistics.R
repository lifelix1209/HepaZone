# ==============================================================================
# Statistical Analysis Functions
# ==============================================================================

#' Bootstrap Standard Error Estimation
#'
#' Performs bootstrap resampling to estimate standard errors for
#' spatial expression estimates.
#'
#' @param expression_matrix Gene x Cell normalized expression matrix
#' @param prob_matrix Cell x Zone probability matrix
#' @param n_bootstrap Number of bootstrap iterations (default: 500)
#' @param seed Random seed for reproducibility (default: 42)
#' @param n_cores Number of cores for parallel processing (default: 1)
#' @return A list containing:
#'   - mean: Gene x Zone matrix of mean expression
#'   - se: Gene x Zone matrix of standard errors
#' @examples
#' \dontrun{
#' boot_result <- bootstrap_se(mat_norm, prob_matrix, n_bootstrap = 500)
#' }
#' @export
bootstrap_se <- function(expression_matrix, prob_matrix,
                         n_bootstrap = 500,
                         seed = 42,
                         n_cores = 1) {

  set.seed(seed)

  n_genes <- nrow(expression_matrix)
  n_zones <- ncol(prob_matrix)
  n_cells <- nrow(prob_matrix)

  # Store bootstrap estimates
  boot_means <- array(0, dim = c(n_bootstrap, n_genes, n_zones))

  message(sprintf("Running %d bootstrap iterations...", n_bootstrap))

  pb <- utils::txtProgressBar(min = 0, max = n_bootstrap, style = 3)

  for (b in 1:n_bootstrap) {
    # Resample cells with replacement
    cell_indices <- sample(1:n_cells, size = n_cells, replace = TRUE)

    # Get resampled probability matrix
    prob_boot <- prob_matrix[cell_indices, ]

    # Calculate mean expression for this bootstrap sample
    # MeanGeneExp = mat_norm %*% Pmat
    mat_boot <- expression_matrix[, cell_indices, drop = FALSE]
    boot_means[b, , ] <- as.matrix(mat_boot %*% prob_boot)

    if (b %% 100 == 0) {
      utils::setTxtProgressBar(pb, b)
    }
  }

  close(pb)

  # Calculate mean and standard error across bootstrap iterations
  boot_mean <- apply(boot_means, c(2, 3), mean)
  boot_se <- apply(boot_means, c(2, 3), sd)

  # Add gene and zone names
  rownames(boot_mean) <- rownames(expression_matrix)
  rownames(boot_se) <- rownames(expression_matrix)
  colnames(boot_mean) <- colnames(prob_matrix)
  colnames(boot_se) <- colnames(prob_matrix)

  message("Bootstrap complete.")

  return(list(
    mean = boot_mean,
    se = boot_se
  ))
}


#' Permutation Test for Spatial Expression Significance
#'
#' Performs permutation tests to assess the significance of spatial
#' expression patterns. Shuffles cell identities and recalculates
#' expression profiles to generate null distribution.
#'
#' @param expression_matrix Gene x Cell normalized expression matrix
#' @param prob_matrix Cell x Zone probability matrix
#' @param n_permutations Number of permutations (default: 1000)
#' @param seed Random seed for reproducibility (default: 123)
#' @return A list containing:
#'   - obs_stat: Observed statistic (sum of squared differences across zones)
#'   - perm_stats: Permuted statistics
#'   - p_values: P-value matrix (gene x zone)
#' @examples
#' \dontrun{
#' perm_result <- permutation_test(mat_norm, prob_matrix, n_permutations = 1000)
#' }
#' @export
permutation_test <- function(expression_matrix, prob_matrix,
                              n_permutations = 1000,
                              seed = 123) {

  set.seed(seed)

  n_genes <- nrow(expression_matrix)
  n_zones <- ncol(prob_matrix)
  n_cells <- nrow(prob_matrix)

  # Original expression profile
  obs_expr <- expression_matrix %*% prob_matrix  # genes x zones

  # Calculate observed statistic: sum of squared differences between adjacent zones
  # This captures how much expression changes across the spatial axis
  obs_stat <- matrix(0, nrow = n_genes, ncol = n_zones - 1)
  for (z in 1:(n_zones - 1)) {
    obs_stat[, z] <- (obs_expr[, z+1] - obs_expr[, z])^2
  }
  obs_stat <- rowSums(obs_stat)

  # Store permutation statistics
  perm_stats <- matrix(0, nrow = n_permutations, ncol = n_genes)

  message(sprintf("Running %d permutations...", n_permutations))

  pb <- utils::txtProgressBar(min = 0, max = n_permutations, style = 3)

  for (p in 1:n_permutations) {
    # Shuffle cells
    cell_indices <- sample(1:n_cells)

    # Recalculate expression with shuffled cells
    perm_prob <- prob_matrix[cell_indices, ]
    perm_expr <- expression_matrix %*% perm_prob

    # Calculate statistic for permuted data
    perm_stat <- matrix(0, nrow = n_genes, ncol = n_zones - 1)
    for (z in 1:(n_zones - 1)) {
      perm_stat[, z] <- (perm_expr[, z+1] - perm_expr[, z])^2
    }
    perm_stats[p, ] <- rowSums(perm_stat)

    if (p %% 100 == 0) {
      utils::setTxtProgressBar(pb, p)
    }
  }

  close(pb)

  # Calculate p-values
  p_values <- matrix(0, nrow = n_genes, ncol = 1)
  for (g in 1:n_genes) {
    p_values[g] <- sum(perm_stats[, g] >= obs_stat[g]) / n_permutations
  }

  # Add names
  names(obs_stat) <- rownames(expression_matrix)
  colnames(p_values) <- "p_value"
  rownames(p_values) <- rownames(expression_matrix)

  message("Permutation test complete.")

  return(list(
    obs_stat = obs_stat,
    perm_stats = perm_stats,
    p_values = p_values
  ))
}


#' Calculate Q-values (Multiple Testing Correction)
#'
#' Converts p-values to q-values using the Benjamini-Hochberg procedure.
#' Uses qvalue package if available, otherwise falls back to stats::p.adjust.
#'
#' @param p_values A vector or matrix of p-values
#' @return Q-values with same dimensions as input
#' @examples
#' \dontrun{
#' qvals <- calculate_qvalues(p_values)
#' }
#' @export
calculate_qvalues <- function(p_values) {
  # Check if qvalue is available
  use_qvalue <- requireNamespace("qvalue", quietly = TRUE)

  # Handle vector input
  if (is.vector(p_values)) {
    if (all(is.na(p_values))) {
      return(rep(NA, length(p_values)))
    }
    if (use_qvalue) {
      qvals <- qvalue::qvalue(p_values)$qvalues
    } else {
      qvals <- stats::p.adjust(p_values, method = "BH")
    }
    return(qvals)
  }

  # Handle matrix input
  if (is.matrix(p_values)) {
    qvals <- matrix(NA, nrow = nrow(p_values), ncol = ncol(p_values))

    # Process each column separately
    for (j in 1:ncol(p_values)) {
      col_pvals <- p_values[, j]
      if (!all(is.na(col_pvals))) {
        if (use_qvalue) {
          qvals[, j] <- qvalue::qvalue(col_pvals)$qvalues
        } else {
          qvals[, j] <- stats::p.adjust(col_pvals, method = "BH")
        }
      }
    }

    rownames(qvals) <- rownames(p_values)
    colnames(qvals) <- colnames(p_values)
    return(qvals)
  }

  stop("p_values must be a vector or matrix")
}


#' Identify Spatially Variable Genes
#'
#' Identifies genes with significant spatial expression patterns
#' based on permutation test results.
#'
#' @param perm_result Result from permutation_test()
#' @param qvals Q-values from calculate_qvalues()
#' @param q_threshold Q-value threshold for significance (default: 0.05)
#' @return A data frame of significant spatially variable genes
#' @examples
#' \dontrun{
#' sig_genes <- identify_svg(perm_result, qvals, q_threshold = 0.05)
#' }
#' @export
identify_svg <- function(perm_result, qvals, q_threshold = 0.05) {
  if (is.matrix(qvals)) {
    p_col <- qvals[, 1]
  } else {
    p_col <- qvals
  }

  # Create results data frame
  results_df <- data.frame(
    gene = names(perm_result$obs_stat),
    observed_stat = as.numeric(perm_result$obs_stat),
    p_value = as.numeric(perm_result$p_values),
    q_value = as.numeric(p_col),
    significant = p_col < q_threshold
  )

  # Sort by q-value
  results_df <- results_df[order(results_df$q_value), ]

  return(results_df)
}


#' Time × Space Interaction Analysis
#'
#' Analyzes interaction effects between time points and spatial zones.
#' This is useful for identifying genes whose spatial expression pattern
#' changes over the circadian cycle.
#'
#' @param hepa_result_list List of HepaZone results for different time points
#' @param time_points Character vector of time point names (e.g., c("ZT00", "ZT06"))
#' @param n_cores Number of cores for parallel processing (default: 1)
#' @return List containing interaction statistics and p-values
#' @export
analyze_time_space_interaction <- function(hepa_result_list, time_points = NULL,
                                            n_cores = 1) {
  if (is.null(time_points)) {
    time_points <- names(hepa_result_list)
  }

  if (length(hepa_result_list) < 2) {
    stop("Need at least 2 time points for interaction analysis")
  }

  # Get common genes across all time points
  gene_lists <- lapply(hepa_result_list, function(x) x$genes)
  common_genes <- Reduce(intersect, gene_lists)

  message(sprintf("Found %d genes common across all time points", length(common_genes)))

  # For each gene, test if spatial pattern differs across time points
  n_genes <- length(common_genes)
  n_zones <- hepa_result_list[[1]]$n_zones

  # Calculate variance of expression across time points for each zone
  # High variance suggests time × space interaction
  interaction_stats <- matrix(0, nrow = n_genes, ncol = n_zones)
  rownames(interaction_stats) <- common_genes
  colnames(interaction_stats) <- paste0("Zone_", 1:n_zones)

  for (g in common_genes) {
    for (z in 1:n_zones) {
      # Get expression across all time points for this gene and zone
      expr_vals <- sapply(hepa_result_list, function(res) {
        res$mean_expression[g, z]
      })
      interaction_stats[g, z] <- var(expr_vals)
    }
  }

  # Identify genes with high interaction (top 10% variance)
  mean_interaction <- rowMeans(interaction_stats)
  threshold <- quantile(mean_interaction, 0.9)

  high_interaction_genes <- names(which(mean_interaction > threshold))

  message(sprintf("Found %d genes with high time × space interaction",
                  length(high_interaction_genes)))

  return(list(
    interaction_scores = interaction_stats,
    mean_interaction = mean_interaction,
    high_interaction_genes = high_interaction_genes,
    time_points = time_points
  ))
}
