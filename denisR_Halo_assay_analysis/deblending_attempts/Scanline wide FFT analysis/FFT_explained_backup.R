#' FFT Explanation for Colony Analysis
#'
#' This script visualizes how FFT improves colony analysis by examining a single
#' horizontal scanline through the first row of colonies on a washed plate.
#' It demonstrates the frequency domain analysis, fundamental frequency identification,
#' and colony center mapping using FFT-based signal processing.
#'
#' @param washed_path `character(1)`. Path to the washed plate image.
#' @param background_path `character(1)`. Path to the background image.
#' @param gitter_dat_path `character(1)` or `NULL`. Path to the gitter .dat file

# Simple Hilbert transform implementation
simple_hilbert <- function(x) {
  n <- length(x)
  X <- fft(x)
  # Create Hilbert multiplier
  h <- rep(0, n)
  if (n %% 2 == 0) {
    h[1] <- 1
    h[2:(n/2)] <- 2
    h[n/2 + 1] <- 1
  } else {
    h[1] <- 1
    h[2:((n+1)/2)] <- 2
  }
  # Apply multiplier and inverse FFT
  return(fft(X * h, inverse = TRUE) / n)
}
#'   for colony grid information (optional).
#' @param out_png `character(1)`. Output path for the explanation figure.
#' @param out_csv `character(1)`. Output path for the metrics CSV file.
#' @param colony_index_for_area `numeric(1)`. Which colony to use for area integration (default 2).
#' @param blur_sigma `numeric(1)` or `NA`. Sigma for Gaussian blur before FFT (default NA = no blur).
#' @param y_row `numeric(1)` or `NULL`. Force specific y-row for scanline (default NULL = auto-detect).
#'
#' @return A `list` containing analysis results and generated plots.
#' @export
fft_explained <- function(washed_path,
                         background_path,
                         gitter_dat_path = NULL,
                         out_png = "results/figs/fft_explained.png",
                         out_csv = "results/fft_explained_metrics.csv",
                         colony_index_for_area = 2,
                         blur_sigma = NA,
                         y_row = NULL) {
  
  # Load required libraries
  if (!require(imager, quietly = TRUE)) stop("Package 'imager' is required")
  if (!require(ggplot2, quietly = TRUE)) stop("Package 'ggplot2' is required")
  if (!require(gridExtra, quietly = TRUE)) stop("Package 'gridExtra' is required")
  if (!require(readr, quietly = TRUE)) stop("Package 'readr' is required")
  if (!require(dplyr, quietly = TRUE)) stop("Package 'dplyr' is required")
  if (!require(pracma, quietly = TRUE)) stop("Package 'pracma' is required for peak finding")
  
  cat("FFT Colony Analysis Explanation\n")
  cat("==============================\n")
  cat(sprintf("Washed image: %s\n", washed_path))
  cat(sprintf("Background image: %s\n", background_path))
  cat(sprintf("Gitter file: %s\n", ifelse(is.null(gitter_dat_path), "None (will auto-detect)", gitter_dat_path)))
  
  # Load images
  cat("\nLoading images...\n")
  washed_img <- imager::load.image(washed_path)
  background_img <- imager::load.image(background_path)
  
  # Convert to grayscale if needed
  if (dim(washed_img)[4] > 1) {
    washed_img <- imager::grayscale(washed_img)
  }
  if (dim(background_img)[4] > 1) {
    background_img <- imager::grayscale(background_img)
  }
  
  cat(sprintf("Image dimensions: %d x %d\n", width(washed_img), height(washed_img)))
  
  # Determine y-coordinate for scanline (center of first row)
  if (!is.null(y_row)) {
    cat(sprintf("Using forced y-row: %d\n", y_row))
    spacing_info <- NULL  # No gitter analysis when forcing y-row
  } else {
    y_row_result <- determine_scanline_y_enhanced(washed_img, gitter_dat_path)
    y_row <- y_row_result$y_row
    spacing_info <- y_row_result$spacing_info
  }
  
  cat(sprintf("Using scanline at y = %d\n", y_row))
  
  # Extract horizontal scanlines
  cat("\nExtracting scanlines...\n")
  washed_line <- extract_horizontal_line(washed_img, y_row)
  background_line <- extract_horizontal_line(background_img, y_row)
  
  # Apply optional blur
  if (!is.na(blur_sigma)) {
    cat(sprintf("Applying Gaussian blur with sigma = %.2f\n", blur_sigma))
    washed_line <- apply_blur_to_line(washed_line, blur_sigma)
    background_line <- apply_blur_to_line(background_line, blur_sigma)
  }
  
  # Background adjustment
  cat("Performing background adjustment...\n")
  adjusted_signal <- perform_background_adjustment(washed_line, background_line)
  
  # Enhanced FFT analysis with detrending and dual frequencies
  cat("Performing enhanced FFT analysis...\n")
  fft_results <- perform_fft_analysis(adjusted_signal, spacing_info, detrend_method = "runmed")
  
  # Validate FFT results and potentially fall back
  validation_result <- validate_fft_results(fft_results, y_row, washed_img, gitter_dat_path)
  if (!validation_result$valid) {
    cat("FFT validation failed, attempting fallback...\n")
    y_row <- validation_result$fallback_y
    cat(sprintf("Using fallback scanline at y = %d\n", y_row))
    
    # Re-extract scanlines
    washed_line <- extract_horizontal_line(washed_img, y_row)
    background_line <- extract_horizontal_line(background_img, y_row)
    if (!is.na(blur_sigma)) {
      washed_line <- apply_blur_to_line(washed_line, blur_sigma)
      background_line <- apply_blur_to_line(background_line, blur_sigma)
    }
    adjusted_signal <- perform_background_adjustment(washed_line, background_line)
    fft_results <- perform_fft_analysis(adjusted_signal, spacing_info = NULL, detrend_method = "runmed")
    spacing_info <- NULL  # Clear gitter info since we're using fallback
  }
  
  # Map frequencies to individual colony centers
  cat("Mapping frequency components to individual colony centers...\n")
  colony_mapping <- map_frequencies_to_individual_colonies(fft_results, length(adjusted_signal))
  
  # Calculate per-colony areas with local minima bounds
  cat("Calculating per-colony areas...\n")
  per_colony_results <- calculate_per_colony_areas(fft_results$signal_detrended, colony_mapping)
  
  # Select specific colony for detailed visualization
  area_results <- select_colony_for_visualization(per_colony_results, colony_index_for_area)
  
  # Create visualization
  cat("Creating visualization...\n")
  plots <- create_fft_visualization(
    washed_line = washed_line,
    background_line = background_line,
    adjusted_signal = adjusted_signal,
    fft_results = fft_results,
    colony_mapping = colony_mapping,
    area_results = area_results,
    y_row = y_row,
    blur_sigma = blur_sigma
  )
  
  # Save outputs
  cat("Saving outputs...\n")
  save_fft_outputs_enhanced(plots, fft_results, colony_mapping, area_results, per_colony_results, 
                           y_row, spacing_info, out_png, out_csv, blur_sigma)
  
  cat(sprintf("Results saved to:\n"))
  cat(sprintf("  Figure: %s\n", out_png))
  cat(sprintf("  Metrics: %s\n", out_csv))
  cat(sprintf("  Per-colony data: %s\n", gsub("\\.csv$", "_per_colony.csv", out_csv)))
  
  # Return results
  return(list(
    y_row = y_row,
    washed_line = washed_line,
    background_line = background_line,
    adjusted_signal = adjusted_signal,
    fft_results = fft_results,
    colony_mapping = colony_mapping,
    area_results = area_results,
    per_colony_results = per_colony_results,
    spacing_info = spacing_info,
    plots = plots
  ))
}

#' Determine Y-coordinate for Scanline
#'
#' Determines the y-coordinate for the horizontal scanline based on either
#' gitter data (first row center) or auto-detection via y-axis projection.
#'
#' @param img `cimg`. The image to analyze.
#' @param gitter_dat_path `character(1)` or `NULL`. Path to gitter .dat file.
#'
#' @return `numeric(1)`. Y-coordinate for the scanline.
determine_scanline_y <- function(img, gitter_dat_path) {
  if (!is.null(gitter_dat_path) && file.exists(gitter_dat_path)) {
    cat("Using gitter file to determine first row center...\n")
    
    # Load gitter data using the existing function pattern
    con <- file(gitter_dat_path, "r")
    lines_to_skip <- 0
    while(TRUE) {
      line <- readLines(con, n = 1)
      if(length(line) == 0) break
      if(startsWith(line, "#")) {
        lines_to_skip <- lines_to_skip + 1
      } else {
        break
      }
    }
    close(con)
    
    gitter_data <- readr::read_tsv(gitter_dat_path, skip = lines_to_skip - 1, show_col_types = FALSE)
    colnames(gitter_data)[1] <- "row"
    
    # Find first row (should be row 1) and calculate mean y-coordinate
    first_row_data <- gitter_data[gitter_data$row == 1, ]
    if (nrow(first_row_data) == 0) {
      warning("No row 1 found in gitter data, using fallback method")
      return(determine_scanline_y_fallback(img))
    }
    
    y_row <- round(mean(first_row_data$y, na.rm = TRUE))
    cat(sprintf("Found %d colonies in first row, mean y = %.1f\n", nrow(first_row_data), mean(first_row_data$y, na.rm = TRUE)))
    
    return(y_row)
  } else {
    cat("No gitter file available, using y-axis projection method...\n")
    return(determine_scanline_y_fallback(img))
  }
}

#' Fallback Method for Scanline Y-coordinate
#'
#' Uses y-axis projection and peak finding to estimate the first row center.
#'
#' @param img `cimg`. The image to analyze.
#'
#' @return `numeric(1)`. Y-coordinate for the scanline.
determine_scanline_y_fallback <- function(img) {
  # Project image onto y-axis (sum over x)
  y_projection <- apply(as.array(img)[,,1,1], 2, sum)
  
  # Simple peak finding: look for local maxima
  # This is a basic implementation - could be enhanced with proper peak detection
  smoothed <- stats::filter(y_projection, rep(1/5, 5), sides = 2)
  smoothed[is.na(smoothed)] <- y_projection[is.na(smoothed)]
  
  # Find the first significant peak (assuming it's the first row)
  # Look in the upper third of the image
  search_range <- 1:floor(length(smoothed) / 3)
  peak_idx <- which.max(smoothed[search_range])
  
  y_row <- search_range[peak_idx]
  
  cat(sprintf("Auto-detected first row at y = %d using y-axis projection\n", y_row))
  
  return(y_row)
}

#' Enhanced Scanline Y-coordinate Determination
#'
#' Determines y-coordinate with gitter spacing analysis and validation.
#'
#' @param img `cimg`. The image to analyze.
#' @param gitter_dat_path `character(1)` or `NULL`. Path to gitter .dat file.
#'
#' @return `list` with y_row and spacing_info.
determine_scanline_y_enhanced <- function(img, gitter_dat_path) {
  spacing_info <- NULL
  
  if (!is.null(gitter_dat_path) && file.exists(gitter_dat_path)) {
    cat("Using gitter file to determine first row center...\n")
    
    # Load gitter data using the existing function pattern
    con <- file(gitter_dat_path, "r")
    lines_to_skip <- 0
    while(TRUE) {
      line <- readLines(con, n = 1)
      if(length(line) == 0) break
      if(startsWith(line, "#")) {
        lines_to_skip <- lines_to_skip + 1
      } else {
        break
      }
    }
    close(con)
    
    gitter_data <- readr::read_tsv(gitter_dat_path, skip = lines_to_skip - 1, show_col_types = FALSE)
    colnames(gitter_data)[1] <- "row"
    
    # Check if row coordinates make sense (detect if row 1 is top or bottom)
    image_height <- height(img)
    first_row_data <- gitter_data[gitter_data$row == 1, ]
    
    if (nrow(first_row_data) == 0) {
      warning("No row 1 found in gitter data, using fallback method")
      return(list(y_row = determine_scanline_y_fallback(img), spacing_info = NULL))
    }
    
    y_row <- round(mean(first_row_data$y, na.rm = TRUE))
    
    # Sanity check: if y_row is in lower 10% of image, it might be correct
    # if y_row is in upper 10%, row numbering might be inverted
    if (y_row > 0.9 * image_height) {
      cat("Warning: First row appears at bottom of image - this may indicate inverted row numbering\n")
    }
    
    cat(sprintf("Found %d colonies in first row, mean y = %.1f\n", nrow(first_row_data), mean(first_row_data$y, na.rm = TRUE)))
    
    # Analyze spacing
    spacing_info <- analyze_gitter_spacing(gitter_data, target_row = 1)
    
    return(list(y_row = y_row, spacing_info = spacing_info))
  } else {
    cat("No gitter file available, using y-axis projection method...\n")
    return(list(y_row = determine_scanline_y_fallback(img), spacing_info = NULL))
  }
}

#' Extract Horizontal Line from Image
#'
#' Extracts a horizontal line of pixels from an image at the specified y-coordinate.
#'
#' @param img `cimg`. The image to extract from.
#' @param y_row `numeric(1)`. Y-coordinate of the line to extract.
#'
#' @return `numeric`. Vector of pixel intensities along the horizontal line.
extract_horizontal_line <- function(img, y_row) {
  # Ensure y_row is within image bounds
  y_row <- max(1, min(y_row, height(img)))
  
  # Extract the line
  line <- as.array(img)[, y_row, 1, 1]
  
  return(as.numeric(line))
}

#' Apply Blur to Line Signal
#'
#' Applies Gaussian blur to a 1D signal (simulating isoblur on the scanline).
#'
#' @param line `numeric`. The line signal to blur.
#' @param sigma `numeric(1)`. Gaussian blur sigma.
#'
#' @return `numeric`. Blurred line signal.
apply_blur_to_line <- function(line, sigma) {
  # Create Gaussian kernel for 1D blur
  kernel_size <- ceiling(6 * sigma)
  if (kernel_size %% 2 == 0) kernel_size <- kernel_size + 1
  
  x <- seq(-(kernel_size-1)/2, (kernel_size-1)/2, 1)
  kernel <- exp(-x^2 / (2 * sigma^2))
  kernel <- kernel / sum(kernel)
  
  # Apply convolution
  blurred <- stats::filter(line, kernel, sides = 2)
  
  # Handle edges by replicating the original values where filter returns NA
  blurred[is.na(blurred)] <- line[is.na(blurred)]
  
  return(as.numeric(blurred))
}

#' Detrend Line Signal
#'
#' Removes low-frequency baseline using running median or LOESS smoothing.
#'
#' @param line `numeric`. The line signal to detrend.
#' @param method `character(1)`. Detrending method: "runmed" or "loess" (default "runmed").
#' @param k `numeric(1)`. For runmed: window size (should be odd). For loess: span parameter.
#'
#' @return `list` with components `detrended` (baseline-removed signal) and `baseline` (estimated baseline).
detrend_line <- function(line, method = "runmed", k = NULL) {
  n <- length(line)
  
  if (method == "runmed") {
    # Default k to ~5x expected intra spacing (estimate ~n/20, make odd)
    if (is.null(k)) {
      k <- 2 * floor(n / 40) + 1  # Roughly n/20, forced to be odd
      k <- max(k, 3)  # Minimum window size
    }
    if (k %% 2 == 0) k <- k + 1  # Ensure odd
    
    # Running median for robust baseline estimation
    baseline <- stats::runmed(line, k)
    
  } else if (method == "loess") {
    # Default span to 0.1-0.2 of line length
    if (is.null(k)) {
      k <- 0.15  # Default span
    }
    
    # LOESS smoothing for baseline
    loess_fit <- stats::loess(line ~ seq_along(line), span = k)
    baseline <- predict(loess_fit)
    
  } else {
    stop("Detrending method must be 'runmed' or 'loess'")
  }
  
  detrended <- line - baseline
  
  return(list(
    detrended = detrended,
    baseline = baseline,
    method = method,
    parameter = k
  ))
}

#' Analyze Gitter Spacing
#'
#' Extracts colony spacing information from gitter data to seed frequency analysis.
#'
#' @param gitter_data `data.frame`. Parsed gitter data with columns row, col, x, y.
#' @param target_row `numeric(1)`. Which row to analyze for spacing.
#'
#' @return `list` with intra-pair and inter-pair spacing estimates.
analyze_gitter_spacing <- function(gitter_data, target_row = 1) {
  # Filter to target row and sort by column
  row_data <- gitter_data[gitter_data$row == target_row, ]
  if (nrow(row_data) < 2) {
    warning("Not enough colonies in target row for spacing analysis")
    return(list(
      dx_intra_seed = NA,
      dx_inter_seed = NA,
      f_intra_seed = NA,
      f_inter_seed = NA,
      n_colonies = 0
    ))
  }
  
  row_data <- row_data[order(row_data$col), ]
  
  # Calculate all adjacent distances
  x_positions <- row_data$x
  distances <- diff(x_positions)
  
  cat(sprintf("Found %d colonies in row %d\n", nrow(row_data), target_row))
  cat("Colony x-positions:", paste(round(x_positions[seq_len(min(8, length(x_positions)))]), collapse = ", "))
  if (length(x_positions) > 8) cat("...")
  cat("\n")
  cat("Adjacent distances:", paste(round(distances[seq_len(min(6, length(distances)))]), collapse = ", "))
  if (length(distances) > 6) cat("...")
  cat("\n")
  
  # Estimate intra and inter spacings
  # For tetramers: typically see alternating small-small-large pattern
  # or similar pattern depending on tetramer arrangement
  
  if (length(distances) >= 2) {
    # Simple approach: assume smallest distance is intra-pair, largest is inter-pair
    dx_intra_seed <- min(distances)
    dx_inter_seed <- max(distances)
    
    # If they're too similar, try a different approach
    if (dx_inter_seed / dx_intra_seed < 1.5) {
      # Look for pattern in distances - sort and take median of smaller vs larger half
      sorted_dist <- sort(distances)
      n_dist <- length(sorted_dist)
      dx_intra_seed <- median(sorted_dist[1:floor(n_dist/2)])
      dx_inter_seed <- median(sorted_dist[ceiling(n_dist/2):n_dist])
    }
  } else {
    dx_intra_seed <- distances[1]
    dx_inter_seed <- distances[1] * 2  # Rough estimate
  }
  
  # Convert to frequency seeds
  f_intra_seed <- 1 / dx_intra_seed
  f_inter_seed <- 1 / dx_inter_seed
  
  cat(sprintf("Estimated spacings - Intra: %.1f px (f=%.4f), Inter: %.1f px (f=%.4f)\n", 
              dx_intra_seed, f_intra_seed, dx_inter_seed, f_inter_seed))
  
  return(list(
    dx_intra_seed = dx_intra_seed,
    dx_inter_seed = dx_inter_seed,
    f_intra_seed = f_intra_seed,
    f_inter_seed = f_inter_seed,
    n_colonies = nrow(row_data),
    x_positions = x_positions,
    distances = distances
  ))
}

#' Find Two Principal Frequencies
#'
#' Identifies intra-tetramer and inter-tetramer frequencies using seeded search bands.
#'
#' @param frequencies `numeric`. Frequency vector (cycles/pixel).
#' @param power_spectrum `numeric`. Power spectrum values.
#' @param f_intra_seed `numeric(1)`. Seed frequency for intra-tetramer spacing.
#' @param f_inter_seed `numeric(1)`. Seed frequency for inter-tetramer spacing.
#' @param search_fraction `numeric(1)`. Search band as fraction of seed frequency (default 0.25).
#'
#' @return `list` with detected frequencies and their properties.
find_dual_frequencies <- function(frequencies, power_spectrum, f_intra_seed, f_inter_seed, 
                                search_fraction = 0.25) {
  
  # Define search bands
  f_intra_band <- c(f_intra_seed * (1 - search_fraction), f_intra_seed * (1 + search_fraction))
  f_inter_band <- c(f_inter_seed * (1 - search_fraction), f_inter_seed * (1 + search_fraction))
  
  # Clip to valid frequency range (0, 0.5]
  f_intra_band <- pmax(0.001, pmin(f_intra_band, 0.5))
  f_inter_band <- pmax(0.001, pmin(f_inter_band, 0.5))
  
  cat(sprintf("Search bands - Intra: [%.4f, %.4f], Inter: [%.4f, %.4f]\n", 
              f_intra_band[1], f_intra_band[2], f_inter_band[1], f_inter_band[2]))
  
  # Find peaks in each band
  f_intra <- find_peak_in_band(frequencies, power_spectrum, f_intra_band)
  f_inter <- find_peak_in_band(frequencies, power_spectrum, f_inter_band)
  
  # Check for harmonic redundancy
  if (!is.na(f_intra) && !is.na(f_inter)) {
    ratio <- max(f_inter, f_intra) / min(f_inter, f_intra)
    if (abs(ratio - round(ratio)) < 0.05) {  # Within 5% of integer ratio
      # Keep the one closer to its seed
      dist_intra <- abs(f_intra - f_intra_seed) / f_intra_seed
      dist_inter <- abs(f_inter - f_inter_seed) / f_inter_seed
      
      if (dist_intra <= dist_inter) {
        f_inter <- NA  # Remove inter frequency
        cat("Harmonic redundancy detected - keeping intra frequency\n")
      } else {
        f_intra <- NA  # Remove intra frequency  
        cat("Harmonic redundancy detected - keeping inter frequency\n")
      }
    }
  }
  
  # Calculate periods
  period_intra <- ifelse(is.na(f_intra), NA, 1 / f_intra)
  period_inter <- ifelse(is.na(f_inter), NA, 1 / f_inter)
  
  # Get amplitudes
  amp_intra <- if (!is.na(f_intra)) {
    idx <- which.min(abs(frequencies - f_intra))
    power_spectrum[idx]
  } else NA
  
  amp_inter <- if (!is.na(f_inter)) {
    idx <- which.min(abs(frequencies - f_inter))
    power_spectrum[idx]
  } else NA
  
  cat(sprintf("Detected frequencies - Intra: %.4f (period=%.1f), Inter: %.4f (period=%.1f)\n", 
              f_intra, period_intra, f_inter, period_inter))
  
  return(list(
    f_intra = f_intra,
    f_inter = f_inter,
    period_intra = period_intra,
    period_inter = period_inter,
    amp_intra = amp_intra,
    amp_inter = amp_inter,
    f_intra_band = f_intra_band,
    f_inter_band = f_inter_band
  ))
}

#' Find Peak in Frequency Band
#'
#' Helper function to find the strongest peak within a frequency band.
#'
#' @param frequencies `numeric`. Frequency vector.
#' @param power_spectrum `numeric`. Power spectrum values.
#' @param band `numeric(2)`. Frequency band [min, max].
#'
#' @return `numeric(1)`. Peak frequency or NA if no peak found.
find_peak_in_band <- function(frequencies, power_spectrum, band) {
  # Find indices within band
  in_band <- frequencies >= band[1] & frequencies <= band[2]
  
  if (sum(in_band) < 3) {
    return(NA)  # Not enough points for peak detection
  }
  
  band_freqs <- frequencies[in_band]
  band_power <- power_spectrum[in_band]
  
  # Simple peak finding - look for local maxima
  n <- length(band_power)
  if (n < 3) return(band_freqs[which.max(band_power)])
  
  # Find local maxima
  peaks <- c()
  for (i in 2:(n-1)) {
    if (band_power[i] > band_power[i-1] && band_power[i] > band_power[i+1]) {
      peaks <- c(peaks, i)
    }
  }
  
  if (length(peaks) == 0) {
    # No local maxima, return global maximum
    return(band_freqs[which.max(band_power)])
  }
  
  # Return frequency of strongest peak
  peak_powers <- band_power[peaks]
  strongest_peak <- peaks[which.max(peak_powers)]
  
  return(band_freqs[strongest_peak])
}

#' Perform Background Adjustment
#'
#' Performs background subtraction and adjustment on the scanline signal.
#'
#' @param washed_line `numeric`. Raw washed scanline.
#' @param background_line `numeric`. Background scanline.
#'
#' @return `numeric`. Background-adjusted signal.
perform_background_adjustment <- function(washed_line, background_line) {
  # Ensure both lines have the same length
  min_length <- min(length(washed_line), length(background_line))
  washed_line <- washed_line[1:min_length]
  background_line <- background_line[1:min_length]
  
  # Background adjustment: signal = washed - background + mean(background)
  bg_mean <- mean(background_line, na.rm = TRUE)
  adjusted_signal <- washed_line - background_line + bg_mean
  
  return(adjusted_signal)
}

#' Perform Enhanced FFT Analysis
#'
#' Performs FFT analysis with detrending and dual-frequency detection.
#'
#' @param signal `numeric`. The background-adjusted signal.
#' @param spacing_info `list`. Gitter spacing analysis results (optional).
#' @param detrend_method `character(1)`. Detrending method to use.
#'
#' @return `list`. Enhanced FFT analysis results.
perform_fft_analysis <- function(signal, spacing_info = NULL, detrend_method = "runmed") {
  N <- length(signal)
  
  # Step 1: Detrend the signal
  cat("Detrending signal...\n")
  detrend_result <- detrend_line(signal, method = detrend_method)
  signal_detrended <- detrend_result$detrended
  
  # Step 2: Remove DC offset and apply Hann window
  signal_dc_removed <- signal_detrended - mean(signal_detrended, na.rm = TRUE)
  hann_window <- 0.5 * (1 - cos(2 * pi * (0:(N-1)) / (N-1)))
  signal_windowed <- signal_dc_removed * hann_window
  
  # Step 3: Perform FFT and compute power spectrum
  fft_result <- fft(signal_windowed)
  power_spectrum <- Mod(fft_result)^2
  
  # For one-sided spectrum
  n_freq <- floor(N/2) + 1
  frequencies <- (0:(n_freq-1)) / N
  power_spectrum <- power_spectrum[1:n_freq]
  
  # Step 4: Find principal frequencies
  if (!is.null(spacing_info) && !is.na(spacing_info$f_intra_seed)) {
    cat("Using gitter-seeded frequency detection...\n")
    freq_results <- find_dual_frequencies(
      frequencies, power_spectrum, 
      spacing_info$f_intra_seed, spacing_info$f_inter_seed
    )
  } else {
    cat("No gitter info available, using fallback frequency detection...\n")
    # Fallback: find the strongest non-DC peak
    if (length(power_spectrum) > 1) {
      peak_idx <- which.max(power_spectrum[2:length(power_spectrum)]) + 1
      freq_results <- list(
        f_intra = frequencies[peak_idx],
        f_inter = NA,
        period_intra = 1 / frequencies[peak_idx],
        period_inter = NA,
        amp_intra = power_spectrum[peak_idx],
        amp_inter = NA
      )
    } else {
      freq_results <- list(
        f_intra = NA, f_inter = NA,
        period_intra = NA, period_inter = NA,
        amp_intra = NA, amp_inter = NA
      )
    }
  }
  
  # Step 5: Narrowband reconstruction
  components <- perform_narrowband_reconstruction(fft_result, frequencies, freq_results, N)
  
  return(list(
    # Original signals
    signal_original = signal,
    signal_detrended = signal_detrended,
    baseline = detrend_result$baseline,
    signal_dc_removed = signal_dc_removed,
    signal_windowed = signal_windowed,
    hann_window = hann_window,
    
    # Frequency domain
    frequencies = frequencies,
    power_spectrum = power_spectrum,
    fft_complex = fft_result,
    
    # Detected frequencies
    f_intra = freq_results$f_intra,
    f_inter = freq_results$f_inter,
    period_intra = freq_results$period_intra,
    period_inter = freq_results$period_inter,
    amp_intra = freq_results$amp_intra,
    amp_inter = freq_results$amp_inter,
    
    # Reconstructed components
    component_intra = components$component_intra,
    component_inter = components$component_inter,
    envelope_intra = components$envelope_intra,
    envelope_inter = components$envelope_inter,
    
    # Metadata
    detrend_method = detrend_method,
    detrend_parameter = detrend_result$parameter
  ))
}

#' Perform Narrowband Reconstruction
#'
#' Reconstructs narrowband components using inverse FFT and computes Hilbert envelopes.
#'
#' @param fft_result `complex`. Full FFT result.
#' @param frequencies `numeric`. Frequency vector.
#' @param freq_results `list`. Detected frequency information.
#' @param N `numeric(1)`. Signal length.
#' @param bandwidth `numeric(1)`. Number of bins to keep around each peak (default 2).
#'
#' @return `list`. Reconstructed components and envelopes.
perform_narrowband_reconstruction <- function(fft_result, frequencies, freq_results, N, bandwidth = 2) {
  
  # Initialize components
  component_intra <- rep(0, N)
  component_inter <- rep(0, N)
  envelope_intra <- rep(0, N)
  envelope_inter <- rep(0, N)
  
  # Reconstruct intra component
  if (!is.na(freq_results$f_intra)) {
    peak_idx <- which.min(abs(frequencies - freq_results$f_intra))
    clean_fft <- rep(0 + 0i, N)
    
    # Keep bins within bandwidth of the peak
    for (b in -bandwidth:bandwidth) {
      idx_pos <- peak_idx + b
      if (idx_pos >= 1 && idx_pos <= length(frequencies)) {
        clean_fft[idx_pos] <- fft_result[idx_pos]
        # Also set negative frequency counterpart
        idx_neg <- N - idx_pos + 2
        if (idx_neg >= 1 && idx_neg <= N && idx_pos > 1) {
          clean_fft[idx_neg] <- fft_result[idx_neg]
        }
      }
    }
    
    # Inverse FFT
    component_intra <- Re(fft(clean_fft, inverse = TRUE)) / N
    
    # Hilbert envelope
    if (requireNamespace("pracma", quietly = TRUE)) {
      analytic_signal <- simple_hilbert(component_intra)
      envelope_intra <- Mod(analytic_signal)
    } else {
      envelope_intra <- abs(component_intra)  # Fallback
    }
  }
  
  # Reconstruct inter component
  if (!is.na(freq_results$f_inter)) {
    peak_idx <- which.min(abs(frequencies - freq_results$f_inter))
    clean_fft <- rep(0 + 0i, N)
    
    # Keep bins within bandwidth of the peak
    for (b in -bandwidth:bandwidth) {
      idx_pos <- peak_idx + b
      if (idx_pos >= 1 && idx_pos <= length(frequencies)) {
        clean_fft[idx_pos] <- fft_result[idx_pos]
        # Also set negative frequency counterpart
        idx_neg <- N - idx_pos + 2
        if (idx_neg >= 1 && idx_neg <= N && idx_pos > 1) {
          clean_fft[idx_neg] <- fft_result[idx_neg]
        }
      }
    }
    
    # Inverse FFT
    component_inter <- Re(fft(clean_fft, inverse = TRUE)) / N
    
    # Hilbert envelope
    if (requireNamespace("pracma", quietly = TRUE)) {
      analytic_signal <- simple_hilbert(component_inter)
      envelope_inter <- Mod(analytic_signal)
    } else {
      envelope_inter <- abs(component_inter)  # Fallback
    }
  }
  
  return(list(
    component_intra = component_intra,
    component_inter = component_inter,
    envelope_intra = envelope_intra,
    envelope_inter = envelope_inter
  ))
}

#' Validate FFT Results
#'
#' Validates FFT results and provides fallback if needed.
#'
#' @param fft_results `list`. FFT analysis results.
#' @param y_row `numeric(1)`. Current y-row.
#' @param img `cimg`. Image for fallback calculation.
#' @param gitter_dat_path `character(1)` or `NULL`. Gitter file path.
#'
#' @return `list` with validation status and fallback information.
validate_fft_results <- function(fft_results, y_row, img, gitter_dat_path) {
  # Check if we have a valid intra frequency
  valid <- !is.na(fft_results$f_intra) && fft_results$f_intra > 0
  
  if (valid) {
    # Additional validation: check if power in the intra band is significant
    # compared to background noise level
    if (!is.na(fft_results$amp_intra)) {
      # Simple threshold: intra amplitude should be at least 2x the median power
      median_power <- median(fft_results$power_spectrum, na.rm = TRUE)
      if (fft_results$amp_intra < 2 * median_power) {
        cat("Warning: Intra frequency amplitude is low relative to background\n")
        valid <- FALSE
      }
    }
  }
  
  fallback_y <- if (!valid) {
    determine_scanline_y_fallback(img)
  } else {
    y_row
  }
  
  return(list(
    valid = valid,
    fallback_y = fallback_y
  ))
}

#' Map Frequencies to Individual Colony Centers
#'
#' Uses envelope detection to find individual colony centers.
#'
#' @param fft_results `list`. FFT analysis results.
#' @param signal_length `numeric(1)`. Length of the signal.
#'
#' @return `list` with individual colony centers.
map_frequencies_to_individual_colonies <- function(fft_results, signal_length) {
  
  if (is.na(fft_results$f_intra) || length(fft_results$envelope_intra) == 0) {
    warning("No valid intra frequency found for colony mapping")
    return(list(
      colony_centers = numeric(0),
      colony_labels = character(0),
      method = "none"
    ))
  }
  
  # Use envelope of intra component to find colony centers
  envelope <- fft_results$envelope_intra
  
  # Find local maxima in the envelope
  colony_centers <- find_local_maxima(envelope)
  
  # Filter centers that are too close to edges or each other
  min_separation <- fft_results$period_intra * 0.3  # Minimum 30% of period between centers
  if (!is.na(min_separation)) {
    colony_centers <- filter_close_centers(colony_centers, min_separation)
  }
  
  # Create labels
  colony_labels <- paste0("C", seq_along(colony_centers))
  
  cat(sprintf("Found %d individual colony centers using envelope detection\n", length(colony_centers)))
  
  return(list(
    colony_centers = colony_centers,
    colony_labels = colony_labels,
    method = "envelope_maxima",
    envelope = envelope
  ))
}

#' Find Local Maxima
#'
#' Finds local maxima in a signal.
#'
#' @param signal `numeric`. Input signal.
#' @param min_height `numeric(1)`. Minimum relative height (default 0.1).
#'
#' @return `numeric`. Indices of local maxima.
find_local_maxima <- function(signal, min_height = 0.1) {
  n <- length(signal)
  if (n < 3) return(numeric(0))
  
  # Normalize signal
  signal_range <- max(signal, na.rm = TRUE) - min(signal, na.rm = TRUE)
  threshold <- min(signal, na.rm = TRUE) + min_height * signal_range
  
  maxima <- c()
  for (i in 2:(n-1)) {
    if (signal[i] > signal[i-1] && signal[i] > signal[i+1] && signal[i] > threshold) {
      maxima <- c(maxima, i)
    }
  }
  
  return(maxima)
}

#' Filter Close Centers
#'
#' Removes centers that are too close to each other.
#'
#' @param centers `numeric`. Center positions.
#' @param min_separation `numeric(1)`. Minimum separation distance.
#'
#' @return `numeric`. Filtered centers.
filter_close_centers <- function(centers, min_separation) {
  if (length(centers) <= 1) return(centers)
  
  centers <- sort(centers)
  filtered <- centers[1]  # Keep first center
  
  for (i in 2:length(centers)) {
    if (centers[i] - filtered[length(filtered)] >= min_separation) {
      filtered <- c(filtered, centers[i])
    }
  }
  
  return(filtered)
}

#' Calculate Per-Colony Areas
#'
#' Calculates area for each detected colony using local minima bounds.
#'
#' @param signal `numeric`. Detrended signal for integration.
#' @param colony_mapping `list`. Colony center mapping results.
#'
#' @return `list` with per-colony area results.
calculate_per_colony_areas <- function(signal, colony_mapping) {
  colony_centers <- colony_mapping$colony_centers
  
  if (length(colony_centers) == 0) {
    return(list(
      colony_areas = data.frame(),
      total_colonies = 0
    ))
  }
  
  # Find integration bounds for each colony using local minima
  colony_areas <- data.frame(
    colony_id = seq_along(colony_centers),
    center_x = colony_centers,
    left_bound = NA,
    right_bound = NA,
    area = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(colony_centers)) {
    center <- colony_centers[i]
    
    # Find local minima to left and right of center
    left_bound <- find_nearest_minimum(signal, center, direction = "left")
    right_bound <- find_nearest_minimum(signal, center, direction = "right")
    
    # Ensure bounds are within signal range
    left_bound <- max(1, left_bound)
    right_bound <- min(length(signal), right_bound)
    
    # Calculate area using trapezoidal integration
    if (right_bound > left_bound) {
      x_vals <- left_bound:right_bound
      y_vals <- signal[x_vals]
      area <- sum(diff(x_vals) * (y_vals[-1] + y_vals[-length(y_vals)]) / 2)
    } else {
      area <- 0
    }
    
    colony_areas$left_bound[i] <- left_bound
    colony_areas$right_bound[i] <- right_bound
    colony_areas$area[i] <- area
  }
  
  cat(sprintf("Calculated areas for %d colonies\n", nrow(colony_areas)))
  cat(sprintf("Area range: %.1f to %.1f\n", min(colony_areas$area, na.rm = TRUE), max(colony_areas$area, na.rm = TRUE)))
  
  return(list(
    colony_areas = colony_areas,
    total_colonies = nrow(colony_areas)
  ))
}

#' Find Nearest Minimum
#'
#' Finds the nearest local minimum in a specified direction from a center point.
#'
#' @param signal `numeric`. Input signal.
#' @param center `numeric(1)`. Center position.
#' @param direction `character(1)`. Search direction: "left" or "right".
#'
#' @return `numeric(1)`. Position of nearest minimum.
find_nearest_minimum <- function(signal, center, direction) {
  center <- round(center)
  n <- length(signal)
  
  if (direction == "left") {
    # Search from center towards beginning
    search_range <- max(1, center - round(n/4)):max(1, center - 1)
    search_range <- search_range[search_range >= 1]
  } else {
    # Search from center towards end
    search_range <- min(n, center + 1):min(n, center + round(n/4))
    search_range <- search_range[search_range <= n]
  }
  
  if (length(search_range) < 2) {
    return(ifelse(direction == "left", 1, n))
  }
  
  # Find local minima in search range
  local_mins <- c()
  for (i in 2:(length(search_range)-1)) {
    idx <- search_range[i]
    if (signal[idx] < signal[idx-1] && signal[idx] < signal[idx+1]) {
      local_mins <- c(local_mins, idx)
    }
  }
  
  if (length(local_mins) == 0) {
    # No local minimum found, return boundary
    return(ifelse(direction == "left", search_range[1], search_range[length(search_range)]))
  }
  
  # Return the nearest minimum to center
  distances <- abs(local_mins - center)
  nearest_idx <- which.min(distances)
  return(local_mins[nearest_idx])
}

#' Select Colony for Visualization
#'
#' Selects a specific colony for detailed visualization.
#'
#' @param per_colony_results `list`. Per-colony area results.
#' @param colony_index `numeric(1)`. Index of colony to visualize.
#'
#' @return `list` with selected colony information.
select_colony_for_visualization <- function(per_colony_results, colony_index) {
  colony_areas <- per_colony_results$colony_areas
  
  if (nrow(colony_areas) == 0) {
    return(list(
      colony_index = colony_index,
      start_px = NA,
      end_px = NA,
      area_under_one_period = NA,
      integration_x = numeric(0),
      integration_y = numeric(0)
    ))
  }
  
  # Validate colony index
  if (colony_index > nrow(colony_areas)) {
    cat(sprintf("Colony index %d not available, using colony 1\n", colony_index))
    colony_index <- 1
  }
  
  selected_colony <- colony_areas[colony_index, ]
  
  return(list(
    colony_index = colony_index,
    start_px = selected_colony$left_bound,
    end_px = selected_colony$right_bound,
    area_under_one_period = selected_colony$area,
    integration_x = selected_colony$left_bound:selected_colony$right_bound,
    integration_y = rep(0, length(selected_colony$left_bound:selected_colony$right_bound))  # Will be filled by visualization
  ))
}

#' Create Enhanced FFT Visualization
#'
#' Creates the enhanced 3-panel visualization with detrended signals and dual frequencies.
#'
#' @param washed_line `numeric`. Raw washed scanline.
#' @param background_line `numeric`. Background scanline.
#' @param adjusted_signal `numeric`. Background-adjusted signal.
#' @param fft_results `list`. Enhanced FFT analysis results.
#' @param colony_mapping `list`. Individual colony center mapping.
#' @param area_results `list`. Selected colony area results.
#' @param y_row `numeric(1)`. Y-coordinate of the scanline.
#' @param blur_sigma `numeric(1)` or `NA`. Blur sigma used.
#'
#' @return `list`. List of ggplot objects for each panel.
create_fft_visualization <- function(washed_line, background_line, adjusted_signal, 
                                   fft_results, colony_mapping, area_results, 
                                   y_row, blur_sigma) {
  
  # Create data frames for plotting
  x_vals <- seq_along(adjusted_signal)
  
  # Panel 1: Raw, adjusted, and detrended signals with individual colony centers
  signal_df <- data.frame(
    x = rep(x_vals, 4),
    intensity = c(
      washed_line[seq_along(adjusted_signal)], 
      background_line[seq_along(adjusted_signal)], 
      adjusted_signal,
      fft_results$signal_detrended
    ),
    signal_type = rep(c("Raw Washed", "Background", "Background-Adjusted", "Detrended"), 
                      each = length(adjusted_signal)),
    stringsAsFactors = FALSE
  )
  
  panel1 <- ggplot(signal_df, aes(x = .data$x, y = .data$intensity, color = .data$signal_type)) +
    geom_line(size = 0.7) +
    scale_color_manual(values = c("Raw Washed" = "blue", 
                                  "Background" = "gray", 
                                  "Background-Adjusted" = "red",
                                  "Detrended" = "darkred")) +
    labs(title = sprintf("Enhanced Scanline Analysis (y = %d)", y_row),
         subtitle = sprintf("Detrend: %s | Blur: %s", 
                           fft_results$detrend_method,
                           ifelse(is.na(blur_sigma), "none", sprintf("σ=%.2f", blur_sigma))),
         x = "Pixel Position", 
         y = "Intensity",
         color = "Signal Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Add individual colony center vertical lines
  if (length(colony_mapping$colony_centers) > 0) {
    for (i in seq_along(colony_mapping$colony_centers)) {
      panel1 <- panel1 + 
        geom_vline(xintercept = colony_mapping$colony_centers[i], 
                   linetype = "dashed", alpha = 0.7, color = "darkgreen") +
        annotate("text", x = colony_mapping$colony_centers[i], 
                 y = max(signal_df$intensity) * 0.95,
                 label = colony_mapping$colony_labels[i], 
                 color = "darkgreen", size = 3, hjust = 0.5)
    }
  }
  
  # Add shaded integration area for selected colony
  if (!is.na(area_results$area_under_one_period) && length(area_results$integration_x) > 0) {
    # Get the actual signal values for the integration region
    integration_signal <- fft_results$signal_detrended[area_results$integration_x]
    integration_df <- data.frame(
      x = area_results$integration_x,
      y = integration_signal
    )
    panel1 <- panel1 + 
      geom_ribbon(data = integration_df, aes(x = .data$x, ymin = min(signal_df$intensity), ymax = .data$y),
                  alpha = 0.3, fill = "orange", inherit.aes = FALSE) +
      annotate("text", x = mean(area_results$integration_x), 
               y = min(signal_df$intensity) + 0.1 * diff(range(signal_df$intensity)),
               label = sprintf("C%d Area = %.1f", area_results$colony_index, area_results$area_under_one_period),
               color = "darkorange", size = 3, fontface = "bold")
  }
  
  # Panel 2: Pre-FFT signal with dual frequency components and envelopes
  prefft_df <- data.frame(
    x = x_vals,
    windowed_signal = fft_results$signal_windowed,
    component_intra = fft_results$component_intra[seq_along(x_vals)],
    component_inter = fft_results$component_inter[seq_along(x_vals)],
    envelope_intra = fft_results$envelope_intra[seq_along(x_vals)],
    envelope_inter = fft_results$envelope_inter[seq_along(x_vals)],
    stringsAsFactors = FALSE
  )
  
  panel2 <- ggplot(prefft_df, aes(x = .data$x)) +
    geom_line(aes(y = .data$windowed_signal, color = "Windowed Signal"), size = 0.7) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Add intra component if available
  if (!all(is.na(fft_results$f_intra))) {
    panel2 <- panel2 +
      geom_line(aes(y = .data$component_intra, color = sprintf("Intra f=%.4f", fft_results$f_intra)), size = 0.5) +
      geom_line(aes(y = .data$envelope_intra, color = sprintf("Intra Envelope")), size = 0.3, linetype = "dotted")
  }
  
  # Add inter component if available
  if (!all(is.na(fft_results$f_inter))) {
    panel2 <- panel2 +
      geom_line(aes(y = .data$component_inter, color = sprintf("Inter f=%.4f", fft_results$f_inter)), size = 0.5) +
      geom_line(aes(y = .data$envelope_inter, color = sprintf("Inter Envelope")), size = 0.3, linetype = "dotted")
  }
  
  # Create dynamic color mapping
  color_values <- c("Windowed Signal" = "black")
  if (!all(is.na(fft_results$f_intra))) {
    intra_label <- sprintf("Intra f=%.4f", fft_results$f_intra)
    color_values[intra_label] <- "red"
    color_values["Intra Envelope"] <- "red"
  }
  if (!all(is.na(fft_results$f_inter))) {
    inter_label <- sprintf("Inter f=%.4f", fft_results$f_inter)
    color_values[inter_label] <- "orange"
    color_values["Inter Envelope"] <- "orange"
  }
  
  panel2 <- panel2 +
    scale_color_manual(values = color_values) +
    labs(title = "Frequency Component Analysis",
         subtitle = sprintf("Intra: %.4f cyc/px (%.1fpx) | Inter: %.4f cyc/px (%.1fpx)", 
                           fft_results$f_intra, fft_results$period_intra,
                           fft_results$f_inter, fft_results$period_inter),
         x = "Pixel Position", 
         y = "Normalized Intensity",
         color = "Component")
  
  # Panel 3: Power spectrum with dual frequency peaks
  spectrum_df <- data.frame(
    frequency = fft_results$frequencies,
    power = fft_results$power_spectrum,
    stringsAsFactors = FALSE
  )
  
  panel3 <- ggplot(spectrum_df, aes(x = .data$frequency, y = .data$power)) +
    geom_line(color = "blue", size = 0.7) +
    labs(title = "Power Spectrum Analysis",
         subtitle = sprintf("Dual-frequency detection: Intra & Inter tetramer spacings"),
         x = "Spatial Frequency (cycles/pixel)", 
         y = "Power") +
    theme_minimal()
  
  # Mark DC component
  panel3 <- panel3 + 
    geom_point(data = data.frame(x = 0, y = spectrum_df$power[1]), 
               aes(x = .data$x, y = .data$y), color = "gray", size = 3) +
    annotate("text", x = 0, y = spectrum_df$power[1] * 1.1, 
             label = "DC", color = "gray", size = 3)
  
  # Mark intra frequency
  if (!is.na(fft_results$f_intra) && fft_results$f_intra > 0) {
    panel3 <- panel3 + 
      geom_point(data = data.frame(x = fft_results$f_intra, y = fft_results$amp_intra), 
                 aes(x = .data$x, y = .data$y), color = "red", size = 4) +
      annotate("text", x = fft_results$f_intra, y = fft_results$amp_intra * 1.15, 
               label = sprintf("f_intra\n%.1fpx", fft_results$period_intra), 
               color = "red", size = 3, hjust = 0.5)
  }
  
  # Mark inter frequency
  if (!is.na(fft_results$f_inter) && fft_results$f_inter > 0) {
    panel3 <- panel3 + 
      geom_point(data = data.frame(x = fft_results$f_inter, y = fft_results$amp_inter), 
                 aes(x = .data$x, y = .data$y), color = "orange", size = 4) +
      annotate("text", x = fft_results$f_inter, y = fft_results$amp_inter * 1.15, 
               label = sprintf("f_inter\n%.1fpx", fft_results$period_inter), 
               color = "orange", size = 3, hjust = 0.5)
  }
  
  return(list(
    panel1 = panel1,
    panel2 = panel2,
    panel3 = panel3
  ))
}

#' Save Enhanced FFT Outputs
#'
#' Saves the enhanced visualization and metrics to files.
#'
#' @param plots `list`. List of ggplot objects.
#' @param fft_results `list`. Enhanced FFT analysis results.
#' @param colony_mapping `list`. Individual colony center mapping.
#' @param area_results `list`. Selected colony area results.
#' @param per_colony_results `list`. Per-colony area results.
#' @param y_row `numeric(1)`. Y-coordinate of scanline.
#' @param spacing_info `list` or `NULL`. Gitter spacing information.
#' @param out_png `character(1)`. Output PNG path.
#' @param out_csv `character(1)`. Output CSV path.
#' @param blur_sigma `numeric(1)` or `NA`. Blur sigma used.
save_fft_outputs_enhanced <- function(plots, fft_results, colony_mapping, area_results, 
                                     per_colony_results, y_row, spacing_info, 
                                     out_png, out_csv, blur_sigma) {
  
  # Create output directories if they don't exist
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  
  # Save the 3-panel figure
  combined_plot <- gridExtra::grid.arrange(
    plots$panel1, 
    plots$panel2, 
    plots$panel3,
    ncol = 1,
    heights = c(1, 1, 1)
  )
  
  ggplot2::ggsave(out_png, combined_plot, width = 12, height = 10, dpi = 300)
  
  # Create main metrics data frame
  metrics_df <- data.frame(
    y_row = y_row,
    f_intra = ifelse(is.na(fft_results$f_intra), NA, fft_results$f_intra),
    period_intra = ifelse(is.na(fft_results$period_intra), NA, fft_results$period_intra),
    f_inter = ifelse(is.na(fft_results$f_inter), NA, fft_results$f_inter),
    period_inter = ifelse(is.na(fft_results$period_inter), NA, fft_results$period_inter),
    n_colonies_detected = ifelse(is.null(per_colony_results$total_colonies), 
                                 length(colony_mapping$colony_centers), 
                                 per_colony_results$total_colonies),
    colony_index_used = area_results$colony_index,
    x_start = area_results$start_px,
    x_end = area_results$end_px,
    area = area_results$area_under_one_period,
    detrend_method = fft_results$detrend_method,
    blur_sigma = ifelse(is.na(blur_sigma), "none", as.character(blur_sigma)),
    notes = sprintf("gitter=%s; spacing_method=%s; envelope_method=%s", 
                   ifelse(is.null(spacing_info), "none", "used"),
                   ifelse(is.null(spacing_info), "fallback", "gitter_seeded"),
                   colony_mapping$method),
    stringsAsFactors = FALSE
  )
  
  # Save main metrics CSV
  readr::write_csv(metrics_df, out_csv)
  
  # Save per-colony CSV
  per_colony_csv <- gsub("\\.csv$", "_per_colony.csv", out_csv)
  if (nrow(per_colony_results$colony_areas) > 0) {
    readr::write_csv(per_colony_results$colony_areas, per_colony_csv)
  } else {
    # Create empty per-colony file with headers
    empty_df <- data.frame(
      colony_id = integer(0),
      center_x = numeric(0),
      left_bound = numeric(0),
      right_bound = numeric(0),
      area = numeric(0)
    )
    readr::write_csv(empty_df, per_colony_csv)
  }
}