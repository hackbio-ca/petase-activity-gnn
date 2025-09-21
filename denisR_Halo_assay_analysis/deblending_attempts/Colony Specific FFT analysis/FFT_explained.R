# FFT_explained.R - Robust plotting system for scanline analysis results
# Reads CSVs from robust FFT_scanline and produces overview plots showing only 
# present colonies and individual 3-panel analysis figures with context, 
# zoomed segments, and simplified FFT visualization.

library(ggplot2)
library(readr)
library(imager)
library(gridExtra)

#' FFT Explained - Robust plotting and QA for scanline analysis
#'
#' Creates comprehensive plots from robust FFT_scanline results including overview
#' showing only present colonies and detailed per-colony analysis with context,
#' segments, and single-line FFT visualization.
#'
#' @param washed_path Path to washed plate image
#' @param background_path Path to background image  
#' @param scanline_csv Path to per-colony CSV from FFT_scanline
#' @param overview_csv Path to overview CSV from FFT_scanline
#' @param out_dir Output directory for figure files
#' @param show_detrended Whether to show detrended signal in plots
#'
#' @return List with paths to generated figure files
fft_explained <- function(
  washed_path,
  background_path,
  scanline_csv = "results/scanline/fft_scanline_per_colony.csv",
  overview_csv = "results/scanline/fft_scanline_overview.csv",
  out_dir = "results/figs/colony",
  show_detrended = TRUE
) {
  
  # Ensure output directory exists
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("=== ROBUST FFT PLOTTING ===\n")
  
  # Read CSV results
  cat("Reading robust analysis results...\n")
  per_colony <- read_csv(scanline_csv, show_col_types = FALSE)
  overview <- read_csv(overview_csv, show_col_types = FALSE)
  
  # Filter to only present colonies
  present_colonies <- per_colony[per_colony$present, ]
  n_total <- nrow(per_colony)
  n_present <- nrow(present_colonies)
  
  cat(sprintf("Present colonies: %d / %d\n", n_present, n_total))
  
  if (n_present == 0) {
    cat("Warning: No present colonies found - creating overview only\n")
  }
  
  # Re-extract the same scanline data (accounting for robust cropping)
  cat("Re-extracting scanline data with robust cropping...\n")
  scanline_data <- extract_scanline_data_robust(
    washed_path, background_path, overview
  )
  
  # Create overview plot (shows only present colonies)
  cat("Creating robust overview plot...\n")
  overview_path <- file.path(out_dir, "scanline_overview.png")
  create_overview_plot_robust(scanline_data, present_colonies, overview_path, 
                             show_detrended, n_present, n_total)
  
  # Create individual colony plots (only for present colonies)  
  colony_paths <- c()
  if (n_present > 0) {
    cat("Creating individual colony plots...\n")
    colony_paths <- create_colony_plots_robust(scanline_data, present_colonies, out_dir, show_detrended)
  }
  
  # Save condensed summary CSV (only present colonies)
  summary_path <- file.path(out_dir, "fft_explained_summary.csv")
  if (n_present > 0) {
    summary_df <- present_colonies[, c("colony_id", "present", "center_x", "left", "right", 
                                      "width_px", "area_adjusted", "area_detrended", "snr", 
                                      "continuity", "touches_edge", "dom_freq_cyc_per_px", 
                                      "dom_period_px", "reason")]
  } else {
    summary_df <- data.frame()  # Empty dataframe if no present colonies
  }
  write_csv(summary_df, summary_path)
  
  cat(paste("✓ Saved overview plot to:", overview_path, "\n"))
  cat(paste("✓ Saved", length(colony_paths), "colony plots to:", out_dir, "\n"))
  cat(paste("✓ Saved summary CSV to:", summary_path, "\n"))
  cat("=== ROBUST PLOTTING COMPLETE ===\n")
  
  return(list(
    overview_plot = overview_path,
    colony_plots = colony_paths,
    summary_csv = summary_path,
    n_present = n_present,
    n_total = n_total
  ))
}

#' Extract scanline data for plotting with robust cropping support
extract_scanline_data_robust <- function(washed_path, background_path, overview) {
  # Load images
  washed_img <- load.image(washed_path)
  background_img <- load.image(background_path)
  
  # Convert to grayscale if needed
  if (dim(washed_img)[4] > 1) washed_img <- grayscale(washed_img)
  if (dim(background_img)[4] > 1) background_img <- grayscale(background_img)
  
  # Apply the same robust cropping as in analysis
  crop_info <- list(
    left = overview$crop_left[1],
    right = overview$crop_right[1], 
    top = overview$crop_top[1],
    bottom = overview$crop_bottom[1]
  )
  
  # Crop images
  washed_img <- washed_img[crop_info$left:crop_info$right, crop_info$top:crop_info$bottom, , ]
  background_img <- background_img[crop_info$left:crop_info$right, crop_info$top:crop_info$bottom, , ]
  
  # Extract scanlines
  y_row <- overview$y_row[1]
  detrend_k <- overview$detrend_k[1]
  img_dims <- dim(washed_img)
  
  # Handle potential dimension variations after cropping
  if (length(img_dims) == 4) {
    washed_line <- as.numeric(as.array(washed_img)[, y_row, 1, 1])
    background_line <- as.numeric(as.array(background_img)[, y_row, 1, 1])
  } else if (length(img_dims) == 3) {
    washed_line <- as.numeric(as.array(washed_img)[, y_row, 1])
    background_line <- as.numeric(as.array(background_img)[, y_row, 1])
  } else {
    washed_line <- as.numeric(as.array(washed_img)[, y_row])
    background_line <- as.numeric(as.array(background_img)[, y_row])
  }
  
  # Process signals exactly as in robust analysis
  adjusted <- washed_line - background_line + mean(background_line)
  detrended <- adjusted - runmed(adjusted, detrend_k)
  
  return(list(
    x = seq_along(adjusted),
    washed = washed_line,
    background = background_line,
    adjusted = adjusted,
    detrended = detrended,
    y_row = y_row
  ))
}

#' Create overview plot showing only present colonies on the scanline
create_overview_plot_robust <- function(scanline_data, present_colonies, output_path, 
                                       show_detrended, n_present, n_total) {
  # Prepare data frame for ggplot
  plot_data <- data.frame(
    x = scanline_data$x,
    adjusted = scanline_data$adjusted
  )
  
  if (show_detrended) {
    plot_data$detrended <- scanline_data$detrended
  }
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = .data$x)) +
    geom_line(aes(y = .data$adjusted, color = "Adjusted"), linewidth = 1) +
    labs(title = paste("Robust Scanline Analysis Overview (y =", scanline_data$y_row, 
                      ") - Present:", n_present, "/", n_total),
         x = "X position (pixels)", 
         y = "Intensity") +
    theme_minimal() +
    theme(legend.position = "top")
  
  if (show_detrended) {
    p <- p + geom_line(aes(y = .data$detrended, color = "Detrended"), linewidth = 0.8, alpha = 0.7)
  }
  
  # Add colony windows and annotations (only for present colonies)
  if (nrow(present_colonies) > 0) {
    for (i in seq_len(nrow(present_colonies))) {
      colony <- present_colonies[i, ]
      
      # Add shaded window
      y_min <- min(plot_data$adjusted, na.rm = TRUE)
      y_max <- max(plot_data$adjusted, na.rm = TRUE)
      
      p <- p + annotate("rect", 
                       xmin = colony$left, xmax = colony$right,
                       ymin = y_min, ymax = y_max,
                       alpha = 0.3, fill = "green")
      
      # Add center line
      p <- p + annotate("segment",
                       x = colony$center_x, xend = colony$center_x,
                       y = y_min, yend = y_max,
                       color = "darkgreen", linewidth = 1.5)
      
      # Add colony label with area
      p <- p + annotate("text", 
                       x = colony$center_x, y = y_max * 0.95,
                       label = paste("C", colony$colony_id, "\nA=", round(colony$area_adjusted, 0), sep = ""),
                       color = "darkgreen", size = 3, fontface = "bold")
    }
  }
  
  # Adjust colors
  color_values <- c("Adjusted" = "blue")
  if (show_detrended) {
    color_values <- c("Adjusted" = "blue", "Detrended" = "red")
  }
  
  p <- p + scale_color_manual(name = "Signal", values = color_values)
  
  # Save plot
  ggsave(output_path, p, width = 14, height = 6, dpi = 300)
  
  return(p)
}

#' Create individual colony analysis plots (only for present colonies)
create_colony_plots_robust <- function(scanline_data, present_colonies, out_dir, show_detrended) {
  colony_paths <- character(nrow(present_colonies))
  
  for (i in seq_len(nrow(present_colonies))) {
    colony <- present_colonies[i, ]
    
    # Create 3-panel plot for this colony
    output_path <- file.path(out_dir, paste0("C", colony$colony_id, ".png"))
    
    p <- create_single_colony_plot_robust(scanline_data, colony, show_detrended)
    
    ggsave(output_path, p, width = 12, height = 10, dpi = 300)
    colony_paths[i] <- output_path
  }
  
  return(colony_paths)
}

#' Create a single colony's 3-panel analysis plot (robust version)
create_single_colony_plot_robust <- function(scanline_data, colony, show_detrended) {
  
  # Panel A: Context - full scanline with this colony highlighted
  plot_data_full <- data.frame(
    x = scanline_data$x,
    adjusted = scanline_data$adjusted
  )
  
  p1 <- ggplot(plot_data_full, aes(x = .data$x, y = .data$adjusted)) +
    geom_line(color = "blue", linewidth = 1) +
    # Highlight this colony's window
    annotate("rect", 
             xmin = colony$left, xmax = colony$right,
             ymin = min(plot_data_full$adjusted, na.rm = TRUE), 
             ymax = max(plot_data_full$adjusted, na.rm = TRUE),
             alpha = 0.3, fill = "orange") +
    # Center line
    annotate("segment",
             x = colony$center_x, xend = colony$center_x,
             y = min(plot_data_full$adjusted, na.rm = TRUE), 
             yend = max(plot_data_full$adjusted, na.rm = TRUE),
             color = "red", linewidth = 2) +
    labs(title = paste("Panel A: Colony", colony$colony_id, 
                      "Context (Area =", round(colony$area_adjusted, 1), 
                      ", SNR =", round(colony$snr, 1), ")"),
         x = "X position (pixels)", y = "Intensity") +
    theme_minimal()
  
  # Panel B: Zoomed segment with narrowband component
  segment_indices <- colony$left:colony$right
  segment_adjusted <- scanline_data$adjusted[segment_indices]
  segment_x <- segment_indices
  
  plot_data_seg <- data.frame(
    x = segment_x,
    adjusted = segment_adjusted
  )
  
  if (show_detrended) {
    plot_data_seg$detrended <- scanline_data$detrended[segment_indices]
  }
  
  # Reconstruct narrowband component
  narrowband <- reconstruct_narrowband_component_robust(segment_adjusted, colony$dom_freq_cyc_per_px)
  plot_data_seg$narrowband <- narrowband
  
  p2 <- ggplot(plot_data_seg, aes(x = .data$x)) +
    geom_line(aes(y = .data$adjusted, color = "Adjusted"), linewidth = 1.5) +
    geom_line(aes(y = .data$narrowband, color = "Narrowband"), linewidth = 1.2, linetype = "dashed") +
    labs(title = paste("Panel B: Colony", colony$colony_id, "Segment (Period =", 
                      round(colony$dom_period_px, 1), "px)"),
         x = "X position (pixels)", y = "Intensity") +
    scale_color_manual(name = "Signal", 
                      values = c("Adjusted" = "blue", "Narrowband" = "purple")) +
    theme_minimal() +
    theme(legend.position = "top")
  
  if (show_detrended) {
    p2 <- p2 + geom_line(aes(y = .data$detrended, color = "Detrended"), linewidth = 0.8, alpha = 0.7)
    p2 <- p2 + scale_color_manual(name = "Signal", 
                                 values = c("Adjusted" = "blue", "Detrended" = "red", 
                                           "Narrowband" = "purple"))
  }
  
  # Panel C: FFT of the segment (simplified - single dashed line)
  p3 <- create_fft_plot_robust(segment_adjusted, colony)
  
  # Combine panels
  combined_plot <- grid.arrange(p1, p2, p3, ncol = 1, heights = c(1, 1, 1))
  
  return(combined_plot)
}

#' Reconstruct narrowband component from dominant frequency (robust version)
reconstruct_narrowband_component_robust <- function(segment, dom_freq) {
  n <- length(segment)
  if (n < 8 || dom_freq <= 0) {
    return(rep(mean(segment), n))  # Return flat baseline if no valid frequency
  }
  
  # Center and window the segment
  segment_centered <- segment - mean(segment)
  hann_window <- 0.5 * (1 - cos(2 * pi * (0:(n-1)) / (n-1)))
  segment_windowed <- segment_centered * hann_window
  
  # FFT
  fft_result <- fft(segment_windowed)
  
  # Create narrowband version by keeping only dominant frequency and neighbors
  fft_narrowband <- fft_result * 0  # Zero everything initially
  
  # Find the bin corresponding to dominant frequency
  freq_resolution <- 1 / n
  dom_bin <- round(dom_freq / freq_resolution) + 1  # +1 for 1-based indexing
  
  # Keep dominant frequency and ±1 neighbors for smoother reconstruction
  if (dom_bin > 1 && dom_bin <= n) {
    for (offset in -1:1) {
      bin_idx <- dom_bin + offset
      if (bin_idx > 1 && bin_idx <= n) {
        fft_narrowband[bin_idx] <- fft_result[bin_idx]
        
        # Add symmetric component for real signal
        symmetric_bin <- n - bin_idx + 2
        if (symmetric_bin > 1 && symmetric_bin <= n && symmetric_bin != bin_idx) {
          fft_narrowband[symmetric_bin] <- fft_result[symmetric_bin]
        }
      }
    }
  }
  
  # Inverse FFT to get narrowband signal
  narrowband_complex <- fft(fft_narrowband, inverse = TRUE) / n
  narrowband <- Re(narrowband_complex) + mean(segment)
  
  return(narrowband)
}

#' Create FFT plot for a segment with single dashed line at dominant frequency
create_fft_plot_robust <- function(segment, colony) {
  n <- length(segment)
  
  # Prepare segment for FFT
  segment_centered <- segment - mean(segment)
  hann_window <- 0.5 * (1 - cos(2 * pi * (0:(n-1)) / (n-1)))
  segment_windowed <- segment_centered * hann_window
  
  # FFT
  fft_result <- fft(segment_windowed)
  
  # One-sided amplitude spectrum
  one_sided_length <- ceiling(n / 2)
  amplitude <- abs(fft_result[1:one_sided_length])
  freq <- (0:(one_sided_length-1)) / n
  
  # Create plot data
  fft_data <- data.frame(
    frequency = freq,
    amplitude = amplitude
  )
  
  # Create simplified FFT plot with single dashed line
  p <- ggplot(fft_data, aes(x = .data$frequency, y = .data$amplitude)) +
    geom_line(color = "black", linewidth = 1) +
    # Mark dominant frequency with single dashed line only
    geom_vline(xintercept = colony$dom_freq_cyc_per_px, 
               color = "red", linewidth = 2, linetype = "dashed") +
    annotate("text", 
             x = colony$dom_freq_cyc_per_px * 1.2, 
             y = max(amplitude, na.rm = TRUE) * 0.8,
             label = paste("f =", round(colony$dom_freq_cyc_per_px, 4), "cyc/px\nT =", 
                          round(colony$dom_period_px, 1), "px"),
             color = "red", size = 3, hjust = 0) +
    labs(title = paste("Panel C: Colony", colony$colony_id, "FFT Spectrum"),
         x = "Frequency (cycles/pixel)", 
         y = "Amplitude") +
    theme_minimal() +
    xlim(0, min(0.15, max(freq, na.rm = TRUE)))  # Focus on relevant frequency range
  
  return(p)
}