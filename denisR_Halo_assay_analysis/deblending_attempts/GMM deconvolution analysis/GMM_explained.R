#' GMM Scanline Analysis Visualization
#'
#' Creates quality assurance plots for GMM-based scanline colony quantification.
#' Generates an overview plot showing the full scanline with present colonies,
#' and detailed plots for each colony showing context, zoom, and responsibilities.
#'
#' @param washed_path `character(1)`. Path to the washed plate image.
#' @param background_path `character(1)`. Path to the background image.
#' @param scanline_csv `character(1)`. Path to per-colony CSV from gmm_scanline.
#' @param overview_csv `character(1)`. Path to overview CSV from gmm_scanline.
#' @param out_dir `character(1)`. Output directory for plot files.
#' @param show_detrended `logical(1)`. Whether to show detrended signal in plots.
#'
#' @return A list containing plot objects and metadata.
#' @export
gmm_explained <- function(
  washed_path,
  background_path,
  scanline_csv = "results/scanline/gmm_scanline_per_colony.csv",
  overview_csv = "results/scanline/gmm_scanline_overview.csv",
  out_dir = "results/figs/colony",
  show_detrended = TRUE
) {
  
  # Load required libraries
  if (!require(ggplot2, quietly = TRUE)) stop("ggplot2 package required")
  if (!require(readr, quietly = TRUE)) stop("readr package required")
  if (!require(dplyr, quietly = TRUE)) stop("dplyr package required")
  if (!require(gridExtra, quietly = TRUE)) stop("gridExtra package required")
  if (!require(imager, quietly = TRUE)) stop("imager package required")
  if (!require(mixtools, quietly = TRUE)) stop("mixtools package required")
  
  # Create output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Read CSV files
  colony_df <- readr::read_csv(scanline_csv, show_col_types = FALSE)
  overview_df <- readr::read_csv(overview_csv, show_col_types = FALSE)
  
  if (nrow(overview_df) == 0) {
    stop("No overview data found")
  }
  
  overview_row <- overview_df[1, ]  # Take first row
  
  # Rebuild scanline for plotting
  scanline_data <- rebuild_scanline_for_plotting(
    washed_path, background_path, overview_row,
    overview_row$detrend_k, overview_row
  )
  
  # Filter for present colonies only
  present_df <- colony_df[colony_df$present, ]
  
  # Create overview plot
  overview_plot <- plot_overview(scanline_data, present_df)
  
  # Save overview plot
  overview_file <- file.path(out_dir, "scanline_overview.png")
  ggplot2::ggsave(overview_file, overview_plot, width = 12, height = 6, dpi = 300)
  
  # Create individual colony plots
  colony_plots <- list()
  
  for (i in seq_len(nrow(present_df))) {
    colony_row <- present_df[i, ]
    colony_id <- colony_row$colony_id
    
    # Create colony plot
    colony_plot <- plot_colony(scanline_data, colony_row, overview_row, show_detrended)
    
    # Save colony plot
    colony_file <- file.path(out_dir, paste0(colony_id, ".png"))
    ggplot2::ggsave(colony_file, colony_plot, width = 12, height = 8, dpi = 300)
    
    colony_plots[[colony_id]] <- colony_plot
  }
  
  cat("Created overview plot:", overview_file, "\n")
  cat("Created", length(colony_plots), "colony plots in", out_dir, "\n")
  
  list(
    overview_plot = overview_plot,
    colony_plots = colony_plots,
    scanline_data = scanline_data,
    present_colonies = present_df
  )
}

#' Rebuild Scanline for Plotting
#'
#' Reconstructs the scanline data from images using the same parameters
#' as the original analysis for consistent plotting.
#'
#' @param washed_path `character(1)`. Path to washed image.
#' @param background_path `character(1)`. Path to background image.
#' @param overview_row `data.frame`. Single row from overview CSV.
#' @param detrend_k `integer(1)`. Detrending parameter.
#' @param crop_info `data.frame`. Crop information from overview.
#' @return A list containing scanline vectors and metadata.
rebuild_scanline_for_plotting <- function(washed_path, background_path, 
                                          overview_row, detrend_k, crop_info) {
  
  # Load images
  washed_img <- imager::load.image(washed_path)
  background_img <- imager::load.image(background_path)
  
  # Convert to grayscale if needed
  if (dim(washed_img)[4] > 1) {
    washed_gray <- 0.299 * washed_img[,,1,1] + 0.587 * washed_img[,,1,2] + 0.114 * washed_img[,,1,3]
  } else {
    washed_gray <- washed_img[,,1,1]
  }
  
  if (dim(background_img)[4] > 1) {
    background_gray <- 0.299 * background_img[,,1,1] + 0.587 * background_img[,,1,2] + 0.114 * background_img[,,1,3]
  } else {
    background_gray <- background_img[,,1,1]
  }
  
  # Apply crop
  washed_cropped <- washed_gray[crop_info$crop_left:crop_info$crop_right, 
                                crop_info$crop_top:crop_info$crop_bottom]
  background_cropped <- background_gray[crop_info$crop_left:crop_info$crop_right, 
                                        crop_info$crop_top:crop_info$crop_bottom]
  
  # Extract scanlines
  y_row_adj <- overview_row$y_row - crop_info$crop_top + 1
  washed_line <- as.numeric(washed_cropped[, y_row_adj])
  background_line <- as.numeric(background_cropped[, y_row_adj])
  
  # Background correction
  adjusted <- washed_line - background_line + mean(background_line)
  
  # Detrending
  baseline <- runmed(adjusted, detrend_k, endrule = "constant")
  detrended <- adjusted - baseline
  
  # X coordinates (adjusted back to original image coordinates)
  x_coords <- seq_along(adjusted) + crop_info$crop_left - 1
  
  list(
    adjusted = adjusted,
    baseline = baseline,
    detrended = detrended,
    washed_line = washed_line,
    background_line = background_line,
    x_coords = x_coords,
    y_row = overview_row$y_row
  )
}

#' Plot Overview
#'
#' Creates an overview plot showing the full scanline with present colonies
#' and their analysis windows highlighted.
#'
#' @param scanline_data `list`. Scanline data from rebuild_scanline_for_plotting.
#' @param present_df `data.frame`. Data for present colonies only.
#' @return A ggplot object.
plot_overview <- function(scanline_data, present_df) {
  
  # Create plotting data frame
  plot_df <- data.frame(
    x = scanline_data$x_coords,
    adjusted = scanline_data$adjusted,
    baseline = scanline_data$baseline,
    detrended = scanline_data$detrended
  )
  
  # Base plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = adjusted), color = "blue", size = 0.5, alpha = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = baseline), color = "red", size = 0.5, alpha = 0.7) +
    ggplot2::labs(
      title = paste("Scanline Overview - Y =", scanline_data$y_row),
      subtitle = paste("Present colonies:", nrow(present_df), "of", nrow(present_df)),
      x = "X coordinate (pixels)",
      y = "Intensity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12)
    )
  
  # Add colony windows
  if (nrow(present_df) > 0) {
    for (i in seq_len(nrow(present_df))) {
      colony <- present_df[i, ]
      
      # Add window rectangle
      p <- p + ggplot2::annotate("rect",
                                 xmin = colony$left, xmax = colony$right,
                                 ymin = -Inf, ymax = Inf,
                                 alpha = 0.2, fill = "green")
      
      # Add colony ID label
      label_x <- (colony$left + colony$right) / 2
      label_y <- max(plot_df$adjusted) * 0.9
      p <- p + ggplot2::annotate("text",
                                 x = label_x, y = label_y,
                                 label = colony$colony_id,
                                 size = 3, color = "darkgreen", fontface = "bold")
    }
  }
  
  return(p)
}

#' Plot Colony
#'
#' Creates a detailed 3-panel plot for a single colony showing context,
#' zoom view with GMM fit, and responsibility weights.
#'
#' @param scanline_data `list`. Scanline data from rebuild_scanline_for_plotting.
#' @param colony_row `data.frame`. Single row for the colony to plot.
#' @param overview_row `data.frame`. Overview metadata.
#' @param show_detrended `logical(1)`. Whether to show detrended signal.
#' @return A combined ggplot object with 3 panels.
plot_colony <- function(scanline_data, colony_row, overview_row, show_detrended = TRUE) {
  
  colony_id <- colony_row$colony_id
  
  # Panel A: Full scanline with window highlighted
  plot_df_full <- data.frame(
    x = scanline_data$x_coords,
    adjusted = scanline_data$adjusted,
    baseline = scanline_data$baseline
  )
  
  panel_a <- ggplot2::ggplot(plot_df_full, ggplot2::aes(x = x)) +
    ggplot2::geom_line(ggplot2::aes(y = adjusted), color = "blue", size = 0.5, alpha = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = baseline), color = "red", size = 0.5, alpha = 0.6) +
    ggplot2::annotate("rect",
                      xmin = colony_row$left, xmax = colony_row$right,
                      ymin = -Inf, ymax = Inf,
                      alpha = 0.3, fill = "orange") +
    ggplot2::labs(
      title = paste("A. Full Scanline Context -", colony_id),
      x = "X coordinate (pixels)",
      y = "Intensity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"))
  
  # Panel B: Window zoom with GMM fit
  window_indices <- which(scanline_data$x_coords >= colony_row$left & 
                          scanline_data$x_coords <= colony_row$right)
  
  if (length(window_indices) == 0) {
    # Fallback if no indices found
    panel_b <- ggplot2::ggplot() + ggplot2::labs(title = "B. No data in window")
    panel_c <- ggplot2::ggplot() + ggplot2::labs(title = "C. No data")
  } else {
    
    window_x <- scanline_data$x_coords[window_indices]
    window_adjusted <- scanline_data$adjusted[window_indices]
    window_baseline <- scanline_data$baseline[window_indices]
    
    # Reconstruct GMM for plotting
    gmm_data <- reconstruct_gmm_for_plotting(window_adjusted, colony_row)
    
    plot_df_window <- data.frame(
      x = window_x,
      adjusted = window_adjusted,
      baseline = window_baseline,
      responsibilities = gmm_data$gamma
    )
    
    # Create ribbon data for responsibilities
    plot_df_window$upper <- plot_df_window$adjusted
    plot_df_window$lower <- plot_df_window$baseline
    
    panel_b <- ggplot2::ggplot(plot_df_window, ggplot2::aes(x = x)) +
      ggplot2::geom_line(ggplot2::aes(y = adjusted), color = "blue", size = 0.8) +
      ggplot2::geom_line(ggplot2::aes(y = baseline), color = "red", size = 0.8) +
      ggplot2::geom_hline(yintercept = colony_row$gmm_mu_col, 
                          color = "violet", size = 1, alpha = 0.8) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, alpha = responsibilities),
                           fill = "purple", color = NA) +
      ggplot2::scale_alpha_identity() +
      ggplot2::labs(
        title = paste("B. Window Zoom -", colony_id),
        subtitle = paste("Colony mean =", round(colony_row$gmm_mu_col, 2)),
        x = "X coordinate (pixels)",
        y = "Intensity"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"))
    
    # Panel C: Responsibilities
    panel_c <- ggplot2::ggplot(plot_df_window, ggplot2::aes(x = x, y = responsibilities)) +
      ggplot2::geom_line(color = "purple", size = 1) +
      ggplot2::geom_area(fill = "purple", alpha = 0.3) +
      ggplot2::labs(
        title = paste("C. Colony Responsibilities -", colony_id),
        subtitle = paste("Area =", round(colony_row$area_weighted, 2), 
                         "| SNR =", round(colony_row$snr, 2)),
        x = "X coordinate (pixels)",
        y = "Responsibility"
      ) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"))
  }
  
  # Combine panels
  combined_plot <- gridExtra::grid.arrange(panel_a, panel_b, panel_c, 
                                           ncol = 1, heights = c(1, 1.2, 1))
  
  return(combined_plot)
}

#' Reconstruct GMM for Plotting
#'
#' Reconstructs the GMM fit for a colony window to extract responsibilities
#' for visualization purposes.
#'
#' @param window_y `numeric`. Intensity values in the window.
#' @param colony_row `data.frame`. Colony data with GMM parameters.
#' @return A list with gamma (responsibilities) and other GMM info.
reconstruct_gmm_for_plotting <- function(window_y, colony_row) {
  
  # Extract GMM parameters
  mu_bg <- colony_row$gmm_mu_bg
  mu_col <- colony_row$gmm_mu_col
  sd_bg <- colony_row$gmm_sd_bg
  sd_col <- colony_row$gmm_sd_col
  pi_bg <- colony_row$gmm_pi_bg
  pi_col <- colony_row$gmm_pi_col
  
  # Check for missing parameters
  if (any(is.na(c(mu_bg, mu_col, sd_bg, sd_col, pi_bg, pi_col)))) {
    # Return uniform responsibilities if GMM failed
    return(list(gamma = rep(0.5, length(window_y))))
  }
  
  # Compute posterior probabilities
  # P(colony | y) = pi_col * N(y; mu_col, sd_col) / [pi_bg * N(y; mu_bg, sd_bg) + pi_col * N(y; mu_col, sd_col)]
  
  likelihood_bg <- pi_bg * dnorm(window_y, mean = mu_bg, sd = sd_bg)
  likelihood_col <- pi_col * dnorm(window_y, mean = mu_col, sd = sd_col)
  
  total_likelihood <- likelihood_bg + likelihood_col
  
  # Avoid division by zero
  total_likelihood[total_likelihood < 1e-10] <- 1e-10
  
  gamma <- likelihood_col / total_likelihood
  
  # Clamp between 0 and 1
  gamma <- pmax(0, pmin(1, gamma))
  
  list(gamma = gamma)
}