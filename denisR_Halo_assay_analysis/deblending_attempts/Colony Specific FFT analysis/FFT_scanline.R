# FFT_scanline.R - Scanline analyzer for colony intensity measurement (no gitter)
# Implements the analysis pipeline without gitter dependency:
# A) Full image processing (no gitter-based cropping)
# B) Manual colony center specification
# C) Adaptive detrending and baseline computation 
# D) Minima-based windowing with background-aware tail expansion
# E) SNR/area/continuity filtering for present colonies
# F) Per-colony FFT frequency detection
# G) Trapezoid area integration
# H) CSV outputs with comprehensive metrics

library(imager)
library(readr)
library(dplyr)

#' Estimate robust scanline background from regions outside colony neighborhoods
estimate_scanline_background <- function(adjusted, centers, pad_px = 80) {
  keep <- rep(TRUE, length(adjusted))
  for (cx in centers) {
    L <- max(1, round(cx - pad_px))
    R <- min(length(adjusted), round(cx + pad_px))
    keep[L:R] <- FALSE
  }
  # Fallback if everything got masked (very dense line)
  if (!any(keep)) keep[] <- TRUE
  list(mu = median(adjusted[keep], na.rm = TRUE),
       mad = mad(adjusted[keep], na.rm = TRUE))
}

#' Extend bounds until adjusted signal looks like background for a sustained run
extend_until_background <- function(adjusted, L, R, bg_mu, bg_tol, 
                                    limit_left = 1, limit_right = length(adjusted),
                                    quiet_run = 30) {
  # Left extension
  i <- L; run <- 0
  while (i > limit_left && run < quiet_run) {
    i <- i - 1
    if (abs(adjusted[i] - bg_mu) < bg_tol) run <- run + 1 else run <- 0
  }
  newL <- i + 1

  # Right extension
  j <- R; run <- 0
  while (j < limit_right && run < quiet_run) {
    j <- j + 1
    if (abs(adjusted[j] - bg_mu) < bg_tol) run <- run + 1 else run <- 0
  }
  newR <- j - 1

  c(newL, newR)
}

#' FFT Scanline Analysis - Robust implementation with comprehensive filtering
#'
#' Performs robust scanline analysis with global cropping, signal-minima based
#' windowing, comprehensive present-colony filtering, and accurate area measurement.
#'
#' @param washed_path Path to washed plate image
#' @param background_path Path to background image
#' @param colony_centers Vector of x-coordinates for colony centers (required)
#' @param y_row Y-coordinate for the scanline to analyze (required)
#' @param out_dir Output directory for CSV files
#' @param snr_min Present colony SNR threshold (default 3)
#' @param area_min Present colony area threshold (default 50)
#' @param eps_baseline Baseline return threshold multiplier (default 0.5)
#' @param r_intra_min Minimum intra-colony search radius (default 40)
#' @param r_intra_max Maximum intra-colony search radius (default 220)
#' @param continuity_q Continuity fraction requirement (default 0.6)
#' @param continuity_k Continuity threshold multiplier (default 0.5)
#' @param background_k Background tolerance multiplier (default 2.5)
#' @param quiet_run Consecutive quiet pixels required for tail termination (default 30)
#' @param pad_px_bg Padding around centers when estimating background (default 80)
#' @param save_debug Save debug plots (default TRUE)
#'
#' @return List containing processed data and file paths
fft_scanline <- function(
  washed_path,
  background_path,
  colony_centers,
  y_row,
  out_dir = "results/scanline",
  snr_min = 3,
  area_min = 50,
  eps_baseline = 0.5,
  r_intra_min = 40,
  r_intra_max = 220,
  continuity_q = 0.6,
  continuity_k = 0.5,
  background_k = 2.5,
  quiet_run = 30,
  pad_px_bg = 80,
  save_debug = TRUE
) {
  
  # Ensure output directory exists
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Input validation
  if (is.null(colony_centers) || length(colony_centers) == 0) {
    stop("Colony centers are required for scanline analysis")
  }
  if (is.null(y_row) || !is.numeric(y_row)) {
    stop("Y-row coordinate is required for scanline analysis")
  }
  
  cat("=== FFT SCANLINE ANALYSIS (No Gitter) ===\n")
  
  # Load images
  cat("Loading images...\n")
  washed_img <- load.image(washed_path)
  background_img <- load.image(background_path)
  
  # Convert to grayscale if needed
  if (dim(washed_img)[4] > 1) washed_img <- grayscale(washed_img)
  if (dim(background_img)[4] > 1) background_img <- grayscale(background_img)
  
  # Use full image dimensions (no cropping)
  cat("Using full image dimensions (no gitter-based cropping)...\n")
  
  # No cropping - record full image dimensions for reference
  img_dims <- dim(washed_img)
  crop_info <- list(
    left = 1,
    right = img_dims[1], 
    top = 1,
    bottom = img_dims[2]
  )
  
  cat(sprintf("Full image size: %dx%d\n", img_dims[1], img_dims[2]))
  
  # B) SCANLINE & BASELINE - Extract and process signals
  cat("Step B: Extracting scanline and computing baseline...\n")
  
  cat(paste("Using y_row:", y_row, "\n"))
  
  # Extract scanlines
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
  
  # Background adjustment 
  adjusted <- washed_line - background_line + mean(background_line)
  
  # Use provided colony centers
  centers <- sort(colony_centers)
  cat(paste("Using", length(centers), "provided colony centers\n"))
  
  # Adaptive detrending
  if (length(centers) > 1) {
    median_spacing <- median(diff(centers))
    detrend_k <- round(3.5 * median_spacing)
    if (detrend_k %% 2 == 0) detrend_k <- detrend_k + 1
    detrend_k <- max(21, min(detrend_k, 201))
  } else {
    detrend_k <- 101
  }
  
  detrended <- adjusted - runmed(adjusted, detrend_k)
  cat(paste("Applied detrending with window size:", detrend_k, "\n"))
  
  
  # C) WINDOWING - Derive bounds from signal minima + tail expansion
  cat("Step C: Computing colony bounds from signal minima...\n")
  windows_and_flags <- derive_bounds_from_minima(
    centers, adjusted, detrended, 
    r_intra_min, r_intra_max, eps_baseline,
    background_k, quiet_run, pad_px_bg
  )
  
  # D) PRESENT COLONY FILTERING - SNR, area, continuity tests
  cat("Step D: Filtering for present colonies...\n")
  colony_results <- list()
  n_present <- 0
  
  for (i in seq_along(centers)) {
    L <- windows_and_flags$left[i]
    R <- windows_and_flags$right[i]
    cx <- centers[i]
    
    # Check edge touching
    touches_edge <- (L <= 5) || (R >= (length(adjusted) - 5))
    
    # Extract segments
    seg_adjusted <- adjusted[L:R]
    seg_detrended <- detrended[L:R]
    
    # Compute baseline reference from flanks
    bg_ref <- compute_baseline_ref(adjusted, detrended, L, R, length(centers))
    
    # D1) Excess area
    excess_values <- pmax(seg_adjusted - bg_ref, 0)
    area_adjusted <- trapezoid_area(excess_values)
    area_detrended <- trapezoid_area(pmax(seg_detrended, 0))
    
    # D2) SNR calculation
    top_10_pct <- quantile(seg_detrended, 0.9, na.rm = TRUE)
    median_detrended <- median(seg_detrended, na.rm = TRUE)
    mad_detrended <- mad(seg_detrended, na.rm = TRUE)
    snr <- (top_10_pct - median_detrended) / (mad_detrended + 1e-8)
    
    # D3) Continuity test
    continuity_thresh <- bg_ref + continuity_k * mad_detrended
    above_thresh <- sum(seg_adjusted > continuity_thresh)
    continuity <- (above_thresh / length(seg_adjusted)) >= continuity_q
    
    # Determine presence and rejection reason
    present <- FALSE
    reason <- ""
    
    if (touches_edge) {
      reason <- "edge"
    } else if (area_adjusted < area_min) {
      reason <- "low_area"
    } else if (snr < snr_min) {
      reason <- "low_snr" 
    } else if (!continuity) {
      reason <- "discontinuous"
    } else {
      present <- TRUE
      n_present <- n_present + 1
    }
    
    # E) FFT ANALYSIS - Only for present colonies
    fft_result <- if (present) {
      analyze_segment_fft_robust(seg_adjusted, L, R)
    } else {
      list(dom_freq = 0, dom_period = Inf)
    }
    
    # Store results
    colony_results[[i]] <- data.frame(
      colony_id = i,
      present = present,
      center_x = cx,
      left = L,
      right = R,
      width_px = R - L + 1,
      area_adjusted = area_adjusted,
      area_detrended = area_detrended,
      snr = snr,
      continuity = continuity,
      touches_edge = touches_edge,
      dom_freq_cyc_per_px = fft_result$dom_freq,
      dom_period_px = fft_result$dom_period,
      reason = reason,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results
  per_colony_df <- do.call(rbind, colony_results)
  
  cat(sprintf("Present colonies: %d / %d\n", n_present, length(centers)))
  
  
  # Apply guardrails - retry with lower SNR if no present colonies found
  if (n_present == 0) {
    cat("Warning: No present colonies found, retrying with 20% lower SNR threshold...\n")
    lower_snr <- snr_min * 0.8
    
    # Re-evaluate presence with lower threshold
    for (i in seq_along(centers)) {
      if (per_colony_df$reason[i] == "low_snr" && 
          per_colony_df$snr[i] >= lower_snr &&
          !per_colony_df$touches_edge[i] &&
          per_colony_df$area_adjusted[i] >= area_min) {
        
        per_colony_df$present[i] <- TRUE
        per_colony_df$reason[i] <- ""
        n_present <- n_present + 1
        
        # Recompute FFT for this colony
        L <- per_colony_df$left[i]
        R <- per_colony_df$right[i]
        seg_adjusted <- adjusted[L:R]
        fft_result <- analyze_segment_fft_robust(seg_adjusted, L, R)
        per_colony_df$dom_freq_cyc_per_px[i] <- fft_result$dom_freq
        per_colony_df$dom_period_px[i] <- fft_result$dom_period
      }
    }
    
    if (n_present > 0) {
      cat(sprintf("Recovered %d colonies with lower SNR threshold (%.2f)\n", n_present, lower_snr))
    }
  }
  
  # H) OUTPUT CSV FILES
  cat("Step H: Writing output files...\n")
  
  # Overview metadata
  overview_df <- data.frame(
    y_row = y_row,
    n_colonies_total = length(centers),
    n_colonies_present = n_present,
    crop_left = crop_info$left,
    crop_right = crop_info$right,
    crop_top = crop_info$top,
    crop_bottom = crop_info$bottom,
    detrend_k = detrend_k,
    eps_baseline = eps_baseline,
    area_min = area_min,
    snr_min = snr_min,
    continuity_q = continuity_q,
    continuity_k = continuity_k,
    washed_path = washed_path,
    background_path = background_path,
    stringsAsFactors = FALSE
  )
  
  # Write CSV files
  per_colony_path <- file.path(out_dir, "fft_scanline_per_colony.csv")
  overview_path <- file.path(out_dir, "fft_scanline_overview.csv")
  
  write_csv(per_colony_df, per_colony_path)
  write_csv(overview_df, overview_path)
  
  cat(paste("Saved per-colony results to:", per_colony_path, "\n"))
  cat(paste("Saved overview metadata to:", overview_path, "\n"))
  
  # Save debug plot if requested
  if (save_debug) {
    debug_path <- file.path(out_dir, "scanline_overview.png")
    save_debug_plot_robust(adjusted, detrended, centers, windows_and_flags, 
                          per_colony_df, debug_path, y_row, n_present, length(centers))
    cat(paste("Saved debug plot to:", debug_path, "\n"))
  }
  
  cat("=== ROBUST ANALYSIS COMPLETE ===\n")
  cat(sprintf("Present colonies: %d / %d\n", n_present, length(centers)))
  
  # Return processed data
  return(list(
    per_colony = per_colony_df,
    overview = overview_df,
    adjusted = adjusted,
    detrended = detrended,
    centers = centers,
    windows = windows_and_flags,
    y_row = y_row,
    crop_info = crop_info
  ))
}

#' Derive bounds from signal minima with robust background-aware tail expansion
derive_bounds_from_minima <- function(centers, adjusted, detrended, 
                                     r_intra_min, r_intra_max, eps_baseline,
                                     background_k = 2.5, quiet_run = 30, pad_px_bg = 80) {
  n_colonies <- length(centers)
  left_bounds <- numeric(n_colonies)
  right_bounds <- numeric(n_colonies)
  
  # Compute intra-colony search radius
  if (length(centers) > 1) {
    r_intra <- 0.6 * median(diff(sort(centers)))
    r_intra <- max(r_intra_min, min(r_intra_max, r_intra))
  } else {
    r_intra <- (r_intra_min + r_intra_max) / 2
  }
  
  # --- NEW: robust background-aware expansion ---
  # 1) Get background stats on adjusted signal (outside neighborhoods)
  bg <- estimate_scanline_background(adjusted, centers, pad_px = pad_px_bg)
  bg_mu  <- bg$mu
  bg_tol <- background_k * bg$mad
  
  for (i in seq_along(centers)) {
    cx <- centers[i]
    
    # initial minima within intra range (computed above)
    left_bound  <- left_bounds[i]
    right_bound <- right_bounds[i]
    
    # hard limits so we never invade a neighbor's core region:
    left_limit  <- if (i == 1) 1 else floor((centers[i-1] + cx) / 2)
    right_limit <- if (i == length(centers)) length(adjusted) else ceiling((cx + centers[i+1]) / 2)
    
    # extend until we've seen a sustained run of background-like pixels
    ext <- extend_until_background(
      adjusted, left_bound, right_bound,
      bg_mu = bg_mu, bg_tol = bg_tol,
      limit_left = left_limit, limit_right = right_limit,
      quiet_run = quiet_run
    )
    left_bounds[i]  <- max(left_limit,  ext[1])
    right_bounds[i] <- min(right_limit, ext[2])
  }
  # --- END NEW ---
  
  # Enforce non-overlap between adjacent colonies (keep existing logic)
  for (i in 2:n_colonies) {
    if (left_bounds[i] <= right_bounds[i-1]) {
      midpoint <- round((left_bounds[i] + right_bounds[i-1]) / 2)
      right_bounds[i-1] <- midpoint - 1
      left_bounds[i] <- midpoint + 1
    }
  }
  
  return(list(left = left_bounds, right = right_bounds))
}

#' Compute baseline reference from flanks
compute_baseline_ref <- function(adjusted, detrended, L, R, n_colonies) {
  # Use flanks to estimate baseline: L-2W:L-W and R+W:R+2W
  W <- round((R - L + 1) / 3)  # Window width based on colony size
  
  # Left flank
  left_flank_start <- max(1, L - 2 * W)
  left_flank_end <- max(1, L - W)
  if (left_flank_start < left_flank_end) {
    left_baseline <- median(detrended[left_flank_start:left_flank_end], na.rm = TRUE)
  } else {
    left_baseline <- 0
  }
  
  # Right flank  
  right_flank_start <- min(length(adjusted), R + W)
  right_flank_end <- min(length(adjusted), R + 2 * W)
  if (right_flank_start < right_flank_end) {
    right_baseline <- median(detrended[right_flank_start:right_flank_end], na.rm = TRUE)
  } else {
    right_baseline <- 0
  }
  
  # Average flanks and map back to adjusted space
  baseline_detrended <- mean(c(left_baseline, right_baseline), na.rm = TRUE)
  bg_ref <- mean(adjusted[L:R], na.rm = TRUE) - baseline_detrended
  
  return(bg_ref)
}

#' Calculate trapezoid area under curve  
trapezoid_area <- function(y) {
  if (length(y) < 2) return(0)
  x <- seq_along(y)
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
}

#' Analyze segment with robust FFT within search band
analyze_segment_fft_robust <- function(segment, L, R) {
  n <- length(segment)
  if (n < 8) {
    return(list(dom_freq = 0, dom_period = Inf))
  }
  
  # Remove DC and apply Hann window
  segment_centered <- segment - mean(segment)
  hann_window <- 0.5 * (1 - cos(2 * pi * (0:(n-1)) / (n-1)))
  segment_windowed <- segment_centered * hann_window
  
  # FFT
  fft_result <- fft(segment_windowed)
  
  # One-sided amplitude spectrum
  one_sided_length <- ceiling(n / 2)
  amplitude <- abs(fft_result[1:one_sided_length])
  freq <- (0:(one_sided_length-1)) / n
  
  # Define search band based on bump width
  W <- R - L + 1
  f_center <- 1 / W  # Expected fundamental
  freq_min <- 0.7 / W
  freq_max <- 1.3 / W
  
  # Find dominant frequency within search band
  freq_mask <- freq >= freq_min & freq <= freq_max & freq > 0
  if (any(freq_mask)) {
    valid_amp <- amplitude[freq_mask]
    valid_freq <- freq[freq_mask]
    
    if (length(valid_amp) > 0 && max(valid_amp) > 0) {
      max_idx <- which.max(valid_amp)
      dom_freq <- valid_freq[max_idx]
    } else {
      dom_freq <- f_center  # Fallback to expected
    }
  } else {
    # No valid frequencies in band
    dom_freq <- 0
  }
  
  dom_period <- if (dom_freq > 0) 1 / dom_freq else Inf
  
  return(list(dom_freq = dom_freq, dom_period = dom_period))
}

#' Save robust debug plot showing scanline with present colonies only
save_debug_plot_robust <- function(adjusted, detrended, centers, windows_and_flags, 
                                  per_colony_df, file_path, y_row, n_present, n_total) {
  png(file_path, width = 1200, height = 400)
  
  x <- seq_along(adjusted)
  
  # Plot adjusted and detrended signals
  plot(x, adjusted, type = "l", col = "blue", lwd = 2,
       xlab = "X position (pixels)", ylab = "Intensity",
       main = paste("Robust Scanline Analysis Overview (y =", y_row, ")"))
  lines(x, detrended, col = "red", lwd = 1)
  
  # Add colony windows and centers (only for present colonies)
  present_colonies <- per_colony_df[per_colony_df$present, ]
  
  if (nrow(present_colonies) > 0) {
    for (i in seq_len(nrow(present_colonies))) {
      colony <- present_colonies[i, ]
      
      # Shade window
      polygon(c(colony$left, colony$right, colony$right, colony$left),
               c(par("usr")[3], par("usr")[3], par("usr")[4], par("usr")[4]),
               col = rgb(0, 1, 0, 0.2), border = NA)
      
      # Mark center
      abline(v = colony$center_x, col = "green", lwd = 2)
      
      # Label colony
      text(colony$center_x, par("usr")[4] * 0.9, paste("C", colony$colony_id, sep = ""),
           pos = 3, col = "darkgreen", font = 2)
    }
  }
  
  # Add legend with summary info
  legend("topright", 
         legend = c("Adjusted", "Detrended", "Present windows", 
                   paste("Present:", n_present, "/", n_total)),
         col = c("blue", "red", rgb(0, 1, 0, 0.2), "black"),
         lty = c(1, 1, NA, NA),
         pch = c(NA, NA, 15, NA),
         pt.cex = c(1, 1, 2, 1),
         text.col = c("black", "black", "darkgreen", "darkred"))
  
  dev.off()
}