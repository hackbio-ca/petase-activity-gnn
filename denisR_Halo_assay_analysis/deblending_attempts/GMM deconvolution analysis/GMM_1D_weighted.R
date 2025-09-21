#' 1D Gaussian Mixture Model Colony Quantification
#'
#' Fits a 1D GMM over x-coordinates using intensity weights to quantify
#' colony areas. Uses a wide rolling quantile baseline and responsibility-
#' weighted integration to capture full colony tails.
#'
#' @param washed_path `character(1)`. Path to the washed plate image.
#' @param background_path `character(1)`. Path to the background image.
#' @param gitter_dat_path `character(1)` or `NULL`. Path to gitter .dat file.
#' @param out_dir `character(1)`. Output directory for results.
#' @param row_index `integer(1)`. Which gitter row to analyze (1-based).
#' @param y_row `integer(1)` or `NULL`. Explicit y-coordinate for scanline.
#' @param baseline_window `integer(1)`. Window size for rolling quantile baseline.
#' @param baseline_quantile `numeric(1)`. Quantile for baseline (0.15 = 15th percentile).
#' @param margin `integer(1)`. Extra margin around gitter bounding box for crop.
#' @param min_mass `numeric(1)`. Minimum component mass to retain.
#' @param min_snr `numeric(1)`. Minimum signal-to-noise ratio to retain.
#' @param save_plots `logical(1)`. Whether to save individual colony plots.
#'
#' @return A list containing colony areas, GMM results, and metadata.
#' @export
gmm_1d_weighted <- function(
  washed_path,
  background_path,
  gitter_dat_path = NULL,
  out_dir = "results/gmm_1d",
  row_index = 1,
  y_row = NULL,
  baseline_window = 500,
  baseline_quantile = 0.15,
  margin = 25,
  min_mass = 0.01,
  min_snr = 2.0,
  save_plots = TRUE
) {
  
  # Load required libraries
  if (!require(imager, quietly = TRUE)) stop("imager package required")
  if (!require(mixtools, quietly = TRUE)) stop("mixtools package required")
  if (!require(readr, quietly = TRUE)) stop("readr package required")
  if (!require(dplyr, quietly = TRUE)) stop("dplyr package required")
  if (!require(ggplot2, quietly = TRUE)) stop("ggplot2 package required")
  if (!require(zoo, quietly = TRUE)) stop("zoo package required for rolling statistics")
  
  # Create output directories
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  if (save_plots && !dir.exists(file.path(out_dir, "plots"))) {
    dir.create(file.path(out_dir, "plots"), recursive = TRUE)
  }
  
  # Infer gitter path if not provided
  if (is.null(gitter_dat_path)) {
    gitter_dat_path <- paste0(washed_path, ".dat")
  }
  
  # Read gitter data
  gitter_df <- read_gitter_dat(gitter_dat_path)
  
  # Load and process images
  cat("Loading images...\n")
  washed_img <- imager::load.image(washed_path)
  background_img <- imager::load.image(background_path)
  
  # Convert to grayscale
  washed_gray <- convert_to_grayscale(washed_img)
  background_gray <- convert_to_grayscale(background_img)
  
  # Crop images using gitter bounds
  crop_info <- compute_gitter_crop(gitter_df, dim(washed_gray), margin)
  
  washed_cropped <- crop_image(washed_gray, crop_info)
  background_cropped <- crop_image(background_gray, crop_info)
  
  # Get target row data and adjust coordinates
  row_data <- gitter_df[gitter_df$row == row_index, ]
  if (nrow(row_data) == 0) {
    stop(paste("No colonies found for row_index", row_index))
  }
  
  # Adjust gitter coordinates for crop
  row_data$x_adj <- row_data$x - crop_info$left + 1
  row_data$y_adj <- row_data$y - crop_info$top + 1
  
  # Determine y_row for scanline
  if (is.null(y_row)) {
    y_row <- round(mean(row_data$y_adj))
  } else {
    y_row <- y_row - crop_info$top + 1
  }
  
  cat("Extracting scanline at y =", y_row, "\n")
  
  # Extract scanlines
  washed_line <- extract_scanline(washed_cropped, y_row)
  background_line <- extract_scanline(background_cropped, y_row)
  
  # Background correction
  adjusted <- washed_line - background_line + mean(background_line)
  
  # Compute rolling quantile baseline
  cat("Computing rolling quantile baseline with window =", baseline_window, "\n")
  baseline <- compute_rolling_quantile_baseline(adjusted, baseline_window, baseline_quantile)
  
  # Compute weights: max(adjusted - baseline, 0)
  weights <- pmax(adjusted - baseline, 0)
  
  # Create x-coordinates
  x_coords <- seq_along(adjusted)
  
  # Filter out zero weights for GMM fitting
  nonzero_idx <- weights > 0
  if (sum(nonzero_idx) < 10) {
    stop("Insufficient positive signal for GMM fitting")
  }
  
  x_data <- x_coords[nonzero_idx]
  weight_data <- weights[nonzero_idx]
  
  # Initialize GMM components from gitter centers
  cat("Initializing GMM with", nrow(row_data), "components from gitter centers\n")
  init_means <- row_data$x_adj
  init_sds <- rep(50, length(init_means))  # Initial standard deviation
  init_weights <- rep(1/length(init_means), length(init_means))
  
  # Fit weighted 1D GMM
  cat("Fitting weighted 1D GMM...\n")
  gmm_result <- fit_weighted_gmm_1d(x_data, weight_data, init_means, init_sds, init_weights)
  
  if (is.null(gmm_result)) {
    stop("GMM fitting failed")
  }
  
  # Prune components by mass and SNR
  cat("Pruning components by mass and SNR...\n")
  pruned_gmm <- prune_gmm_components(gmm_result, x_coords, weights, min_mass, min_snr)
  
  # Compute colony areas using responsibility-weighted integration
  cat("Computing colony areas...\n")
  colony_areas <- compute_colony_areas(x_coords, weights, baseline, pruned_gmm)
  
  # Create results data frame
  results_df <- create_results_dataframe(colony_areas, pruned_gmm, row_data, 
                                         crop_info, y_row, row_index)
  
  # Save results
  readr::write_csv(results_df, file.path(out_dir, "gmm_1d_colony_areas.csv"))
  
  # Create plots
  if (save_plots) {
    cat("Creating plots...\n")
    
    # Overview plot
    overview_plot <- plot_gmm_overview(x_coords, adjusted, baseline, weights, 
                                       pruned_gmm, colony_areas)
    ggplot2::ggsave(file.path(out_dir, "plots", "overview.png"), overview_plot, 
                    width = 12, height = 6, dpi = 300)
    
    # Individual colony plots
    for (i in seq_along(colony_areas)) {
      colony_plot <- plot_individual_colony(x_coords, adjusted, baseline, weights,
                                            pruned_gmm, colony_areas, i)
      filename <- paste0("colony_", sprintf("%02d", i), ".png")
      ggplot2::ggsave(file.path(out_dir, "plots", filename), colony_plot,
                      width = 10, height = 8, dpi = 300)
    }
  }
  
  # Return results
  list(
    colony_areas = colony_areas,
    results_df = results_df,
    gmm_result = pruned_gmm,
    scanline_data = list(
      x = x_coords + crop_info$left - 1,  # Convert back to original coordinates
      adjusted = adjusted,
      baseline = baseline,
      weights = weights,
      y_row = y_row + crop_info$top - 1
    ),
    metadata = list(
      baseline_window = baseline_window,
      baseline_quantile = baseline_quantile,
      min_mass = min_mass,
      min_snr = min_snr,
      crop_info = crop_info
    )
  )
}

#' Convert Image to Grayscale
convert_to_grayscale <- function(img) {
  if (dim(img)[4] > 1) {
    return(0.299 * img[,,1,1] + 0.587 * img[,,1,2] + 0.114 * img[,,1,3])
  } else {
    return(img[,,1,1])
  }
}

#' Compute Crop Boundaries from Gitter Data
compute_gitter_crop <- function(gitter_df, img_dim, margin) {
  left <- max(1, min(gitter_df$xl, na.rm = TRUE) - margin)
  right <- min(img_dim[1], max(gitter_df$xr, na.rm = TRUE) + margin)
  top <- max(1, min(gitter_df$yt, na.rm = TRUE) - margin)
  bottom <- min(img_dim[2], max(gitter_df$yb, na.rm = TRUE) + margin)
  
  list(left = as.integer(left), right = as.integer(right), 
       top = as.integer(top), bottom = as.integer(bottom))
}

#' Crop Image Using Boundaries
crop_image <- function(img, crop_info) {
  img[crop_info$left:crop_info$right, crop_info$top:crop_info$bottom]
}

#' Extract Horizontal Scanline
extract_scanline <- function(img, y) {
  as.numeric(img[, y])
}

#' Compute Rolling Quantile Baseline
compute_rolling_quantile_baseline <- function(signal, window_size, quantile) {
  n <- length(signal)
  baseline <- numeric(n)
  
  # Use zoo::rollapply for rolling quantile
  half_window <- floor(window_size / 2)
  
  # Pad signal for edge handling
  padded_signal <- c(rep(signal[1], half_window), signal, rep(signal[n], half_window))
  
  # Compute rolling quantile
  rolling_q <- zoo::rollapply(padded_signal, window_size, 
                              function(x) quantile(x, quantile, na.rm = TRUE),
                              align = "center", fill = NA)
  
  # Remove padding
  baseline <- rolling_q[(half_window + 1):(half_window + n)]
  
  # Handle any remaining NAs
  baseline[is.na(baseline)] <- quantile(signal, quantile, na.rm = TRUE)
  
  return(baseline)
}

#' Fit Weighted 1D GMM
fit_weighted_gmm_1d <- function(x_data, weights, init_means, init_sds, init_props) {
  
  # Expand data points according to weights for mixtools
  # This is a way to handle weighted data in mixtools
  max_weight <- max(weights)
  weight_factor <- 100  # Scale factor for discretization
  
  expanded_x <- c()
  for (i in seq_along(x_data)) {
    n_reps <- max(1, round(weights[i] / max_weight * weight_factor))
    expanded_x <- c(expanded_x, rep(x_data[i], n_reps))
  }
  
  cat("Expanded data from", length(x_data), "to", length(expanded_x), "points\n")
  
  tryCatch({
    # Fit GMM using mixtools
    fit <- mixtools::normalmixEM(expanded_x, k = length(init_means),
                                 lambda = init_props, mu = init_means, 
                                 sigma = init_sds, maxit = 1000, epsilon = 1e-6)
    
    list(
      means = fit$mu,
      sds = fit$sigma,
      weights = fit$lambda,
      loglik = fit$loglik,
      posterior = fit$posterior
    )
  }, error = function(e) {
    cat("mixtools failed:", conditionMessage(e), "\n")
    
    # Fallback: simple initialization
    list(
      means = init_means,
      sds = init_sds,
      weights = init_props,
      loglik = NA,
      posterior = NULL
    )
  })
}

#' Prune GMM Components
prune_gmm_components <- function(gmm_result, x_coords, weights, min_mass, min_snr) {
  
  # Compute component masses and SNRs
  n_components <- length(gmm_result$means)
  keep_components <- logical(n_components)
  
  noise_level <- mad(weights[weights > 0], na.rm = TRUE)
  
  for (i in seq_len(n_components)) {
    # Compute responsibility-weighted mass
    mu <- gmm_result$means[i]
    sigma <- gmm_result$sds[i]
    pi <- gmm_result$weights[i]
    
    # Compute responsibilities for this component
    responsibilities <- pi * dnorm(x_coords, mu, sigma)
    total_resp <- sum(responsibilities, na.rm = TRUE)
    
    # Mass is sum of weight * responsibility
    mass <- sum(weights * responsibilities / total_resp, na.rm = TRUE)
    
    # SNR is peak height / noise
    peak_weight <- max(weights[abs(x_coords - mu) < 2 * sigma], na.rm = TRUE)
    snr <- peak_weight / max(noise_level, 1e-6)
    
    # Keep component if it meets criteria
    keep_components[i] <- (mass >= min_mass) && (snr >= min_snr)
    
    cat(sprintf("Component %d: mu=%.1f, mass=%.3f, SNR=%.2f, keep=%s\n", 
                i, mu, mass, snr, keep_components[i]))
  }
  
  if (!any(keep_components)) {
    cat("Warning: No components meet pruning criteria, keeping all\n")
    keep_components[] <- TRUE
  }
  
  # Return pruned GMM
  list(
    means = gmm_result$means[keep_components],
    sds = gmm_result$sds[keep_components],
    weights = gmm_result$weights[keep_components],
    n_components = sum(keep_components),
    kept_indices = which(keep_components)
  )
}

#' Compute Colony Areas
compute_colony_areas <- function(x_coords, weights, baseline, gmm_result) {
  
  n_components <- length(gmm_result$means)
  areas <- numeric(n_components)
  
  for (i in seq_len(n_components)) {
    mu <- gmm_result$means[i]
    sigma <- gmm_result$sds[i]
    pi <- gmm_result$weights[i]
    
    # Compute responsibilities for this component
    responsibilities <- pi * dnorm(x_coords, mu, sigma)
    
    # Normalize responsibilities
    total_density <- rowSums(outer(x_coords, gmm_result$means, function(x, mu_j) {
      idx <- match(mu_j, gmm_result$means)
      gmm_result$weights[idx] * dnorm(x, mu_j, gmm_result$sds[idx])
    }))
    
    responsibilities <- responsibilities / total_density
    responsibilities[is.na(responsibilities)] <- 0
    
    # Area is responsibility-weighted sum of positive signal
    positive_signal <- pmax(weights, 0)
    areas[i] <- sum(responsibilities * positive_signal, na.rm = TRUE)
    
    cat(sprintf("Colony %d: center=%.1f, area=%.2f\n", i, mu, areas[i]))
  }
  
  areas
}

#' Create Results Data Frame
create_results_dataframe <- function(areas, gmm_result, row_data, crop_info, y_row, row_index) {
  
  n_colonies <- length(areas)
  
  data.frame(
    colony_id = paste0("C", seq_len(n_colonies)),
    center_x = gmm_result$means + crop_info$left - 1,  # Convert back to original coords
    center_y = y_row + crop_info$top - 1,
    sigma_x = gmm_result$sds,
    component_weight = gmm_result$weights,
    area_weighted = areas,
    row_index = row_index,
    gitter_centers = paste(row_data$x, collapse = ";"),
    crop_left = crop_info$left,
    crop_right = crop_info$right,
    crop_top = crop_info$top,
    crop_bottom = crop_info$bottom
  )
}

#' Plot GMM Overview
plot_gmm_overview <- function(x_coords, adjusted, baseline, weights, gmm_result, areas) {
  
  plot_df <- data.frame(
    x = x_coords,
    adjusted = adjusted,
    baseline = baseline,
    weights = weights
  )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = adjusted), color = "blue", alpha = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = baseline), color = "red", alpha = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = weights), color = "green", size = 1) +
    ggplot2::labs(
      title = "1D GMM Colony Quantification Overview",
      subtitle = paste("Found", length(gmm_result$means), "colonies"),
      x = "X coordinate (pixels)",
      y = "Intensity"
    ) +
    ggplot2::theme_minimal()
  
  # Add GMM components
  for (i in seq_along(gmm_result$means)) {
    mu <- gmm_result$means[i]
    sigma <- gmm_result$sds[i]
    
    # Add vertical line at component center
    p <- p + ggplot2::geom_vline(xintercept = mu, color = "purple", 
                                 linetype = "dashed", alpha = 0.8)
    
    # Add component label
    max_y <- max(plot_df$adjusted)
    p <- p + ggplot2::annotate("text", x = mu, y = max_y * 0.9, 
                               label = paste0("C", i, "\nA=", round(areas[i], 1)),
                               size = 3, color = "purple")
  }
  
  p
}

#' Plot Individual Colony
plot_individual_colony <- function(x_coords, adjusted, baseline, weights, 
                                   gmm_result, areas, colony_idx) {
  
  mu <- gmm_result$means[colony_idx]
  sigma <- gmm_result$sds[colony_idx]
  area <- areas[colony_idx]
  
  # Define window around colony
  window_size <- 4 * sigma
  x_min <- max(1, mu - window_size)
  x_max <- min(length(x_coords), mu + window_size)
  
  window_idx <- which(x_coords >= x_min & x_coords <= x_max)
  
  plot_df <- data.frame(
    x = x_coords[window_idx],
    adjusted = adjusted[window_idx],
    baseline = baseline[window_idx],
    weights = weights[window_idx]
  )
  
  # Compute responsibilities for this colony
  responsibilities <- gmm_result$weights[colony_idx] * 
    dnorm(plot_df$x, mu, sigma)
  
  # Normalize by total density
  total_density <- rowSums(outer(plot_df$x, gmm_result$means, function(x, mu_j) {
    idx <- match(mu_j, gmm_result$means)
    gmm_result$weights[idx] * dnorm(x, mu_j, gmm_result$sds[idx])
  }))
  
  responsibilities <- responsibilities / total_density
  responsibilities[is.na(responsibilities)] <- 0
  
  plot_df$responsibility <- responsibilities
  
  p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = adjusted), color = "blue") +
    ggplot2::geom_line(ggplot2::aes(y = baseline), color = "red") +
    ggplot2::geom_vline(xintercept = mu, color = "purple", linetype = "dashed") +
    ggplot2::labs(title = paste0("Colony ", colony_idx, " - Signal"),
                  x = "X coordinate", y = "Intensity") +
    ggplot2::theme_minimal()
  
  p2 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = weights), color = "green", size = 1) +
    ggplot2::geom_vline(xintercept = mu, color = "purple", linetype = "dashed") +
    ggplot2::labs(title = "Weights", x = "X coordinate", y = "Weight") +
    ggplot2::theme_minimal()
  
  p3 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = responsibility), color = "purple", size = 1) +
    ggplot2::geom_area(ggplot2::aes(y = responsibility), fill = "purple", alpha = 0.3) +
    ggplot2::labs(title = paste0("Responsibility (Area = ", round(area, 2), ")"),
                  x = "X coordinate", y = "Responsibility") +
    ggplot2::theme_minimal()
  
  gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
}