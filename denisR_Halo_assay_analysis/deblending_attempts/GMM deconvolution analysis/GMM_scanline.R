#' GMM-based Scanline Colony Quantification
#'
#' Performs scanline analysis using Gaussian Mixture Models to quantify
#' colony intensity across a single row of colonies. Background-corrects,
#' detrends, defines colony windows, and fits 2-component GMMs to separate
#' background from colony signal.
#'
#' @param washed_path `character(1)`. Path to the washed plate image.
#' @param background_path `character(1)`. Path to the background image.
#' @param gitter_dat_path `character(1)` or `NULL`. Path to gitter .dat file.
#'   If NULL, will try to infer from washed_path.
#' @param out_dir `character(1)`. Output directory for CSV files.
#' @param row_index `integer(1)`. Which gitter row to analyze (1-based).
#' @param y_row `integer(1)` or `NULL`. Explicit y-coordinate for scanline.
#'   If NULL, uses mean y of gitter row after cropping.
#' @param use_gitter_centers `logical(1)`. Whether to use gitter centers for
#'   window definition.
#' @param margin `integer(1)`. Extra margin around gitter bounding box for crop.
#' @param snap_radius `integer(1)`. Radius for snapping window bounds to local minima.
#' @param r_intra_min `integer(1)`. Minimum distance for intra-colony minima search.
#' @param r_intra_max `integer(1)`. Maximum distance for intra-colony minima search.
#' @param detrend_k `integer(1)` or `NULL`. Window size for running median detrend.
#'   If NULL, auto-computed from median center spacing.
#' @param area_min `numeric(1)`. Minimum weighted area for colony presence.
#' @param snr_min `numeric(1)`. Minimum signal-to-noise ratio for presence.
#' @param resp_min `numeric(1)`. Minimum mean responsibility for presence.
#' @param save_debug `logical(1)`. Whether to save debug information.
#'
#' @return A list containing scanline data, colony metrics, and metadata.
#' @export
gmm_scanline <- function(
  washed_path,
  background_path,
  gitter_dat_path = NULL,
  out_dir = "results/scanline",
  row_index = 1,
  y_row = NULL,
  use_gitter_centers = TRUE,
  margin = 20,
  snap_radius = 30,
  r_intra_min = 80,
  r_intra_max = 300,
  detrend_k = NULL,
  area_min = 2.0,
  snr_min = 0.8,
  resp_min = 0.25,
  save_debug = TRUE
) {
  
  # Load required libraries
  if (!require(imager, quietly = TRUE)) stop("imager package required")
  if (!require(mixtools, quietly = TRUE)) stop("mixtools package required")
  if (!require(readr, quietly = TRUE)) stop("readr package required")
  if (!require(dplyr, quietly = TRUE)) stop("dplyr package required")
  
  # Create output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Infer gitter path if not provided
  if (is.null(gitter_dat_path)) {
    gitter_dat_path <- paste0(washed_path, ".dat")
  }
  
  # Read gitter data
  gitter_df <- read_gitter_dat(gitter_dat_path)
  cat("Read", nrow(gitter_df), "gitter entries\n")
  cat("Gitter bounds: xl =", range(gitter_df$xl, na.rm = TRUE), 
      "xr =", range(gitter_df$xr, na.rm = TRUE), "\n")
  cat("Gitter bounds: yt =", range(gitter_df$yt, na.rm = TRUE), 
      "yb =", range(gitter_df$yb, na.rm = TRUE), "\n")
  
  # Load images
  washed_img <- imager::load.image(washed_path)
  background_img <- imager::load.image(background_path)
  
  cat("Image dimensions: washed =", dim(washed_img)[1:2], 
      "background =", dim(background_img)[1:2], "\n")
  
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
  
  # Compute global crop
  img_dim <- c(width = dim(washed_img)[1], height = dim(washed_img)[2])
  cat("Image dimensions: width =", img_dim["width"], "height =", img_dim["height"], "\n")
  crop_info <- compute_global_crop(gitter_df, img_dim, margin)
  
  # Apply crop with safe bounds checking
  cat("Applying crop to images...\n")
  
  # Ensure bounds are within image dimensions
  left <- max(1, min(crop_info$left, dim(washed_gray)[1]))
  right <- min(dim(washed_gray)[1], max(crop_info$right, 1))
  top <- max(1, min(crop_info$top, dim(washed_gray)[2]))
  bottom <- min(dim(washed_gray)[2], max(crop_info$bottom, 1))
  
  # Double-check bounds are reasonable
  if (left >= right || top >= bottom) {
    stop(paste("Invalid crop boundaries: left=", left, "right=", right, "top=", top, "bottom=", bottom))
  }
  
  crop_width <- right - left + 1
  crop_height <- bottom - top + 1
  
  if (crop_width > 50000 || crop_height > 50000) {
    stop(paste("Crop dimensions too large: width=", crop_width, "height=", crop_height))
  }
  
  cat("Final crop bounds: [", left, ":", right, "] x [", top, ":", bottom, "]\n")
  cat("Crop size:", crop_width, "x", crop_height, "pixels\n")
  
  # Apply the crop
  washed_cropped <- washed_gray[left:right, top:bottom]
  background_cropped <- background_gray[left:right, top:bottom]
  
  # Update crop_info with final bounds
  crop_info$left <- left
  crop_info$right <- right  
  crop_info$top <- top
  crop_info$bottom <- bottom
  
  cat("Cropping completed successfully\n")
  
  # Adjust gitter coordinates for crop
  gitter_df$x_adj <- gitter_df$x - crop_info$left + 1
  gitter_df$y_adj <- gitter_df$y - crop_info$top + 1
  
  # Get target row data
  row_data <- gitter_df[gitter_df$row == row_index, ]
  if (nrow(row_data) == 0) {
    stop(paste("No colonies found for row_index", row_index))
  }
  
  # Determine y_row for scanline
  if (is.null(y_row)) {
    y_row <- round(mean(row_data$y_adj))
  } else {
    y_row <- y_row - crop_info$top + 1  # Adjust for crop
  }
  
  # Extract scanlines
  washed_line <- extract_scanline(washed_cropped, y_row)
  background_line <- extract_scanline(background_cropped, y_row)
  
  # Background correction
  adjusted <- washed_line - background_line + mean(background_line)
  
  # Detrending
  if (is.null(detrend_k)) {
    centers_sorted <- sort(row_data$x_adj)
    if (length(centers_sorted) > 1) {
      median_spacing <- median(diff(centers_sorted))
      detrend_k <- as.integer(3.5 * median_spacing)
      detrend_k <- max(41, min(201, detrend_k))
      if (detrend_k %% 2 == 0) detrend_k <- detrend_k + 1  # Ensure odd
    } else {
      detrend_k <- 101
    }
  }
  
  baseline <- detrend_runmed(adjusted, detrend_k)
  detrended <- adjusted - baseline
  
  # Define windows from centers
  windows <- define_windows_from_centers(
    centers = row_data$x_adj,
    n = length(adjusted),
    detrended = detrended,
    snap_radius = snap_radius,
    r_intra_min = r_intra_min,
    r_intra_max = r_intra_max
  )
  
  # Process each window
  colony_results <- list()
  
  for (i in seq_along(windows)) {
    window <- windows[[i]]
    colony_id <- paste0("C", i)
    
    # Extract window data
    window_indices <- window$left:window$right
    window_adjusted <- adjusted[window_indices]
    window_baseline <- baseline[window_indices]
    window_x <- window_indices
    
    # Check if window touches edges
    touches_edge <- (window$left <= 1) || (window$right >= length(adjusted))
    
    # Fit GMM
    gmm_result <- fit_gmm_window(window_adjusted, k = 2)
    
    if (is.null(gmm_result)) {
      # Failed GMM fit
      colony_results[[i]] <- data.frame(
        colony_id = colony_id,
        present = FALSE,
        center_x = window$center + crop_info$left - 1,
        left = window$left + crop_info$left - 1,
        right = window$right + crop_info$left - 1,
        width_px = window$right - window$left + 1,
        area_weighted = NA,
        mu_lift = NA,
        snr = NA,
        resp_mean = NA,
        gmm_mu_bg = NA,
        gmm_sd_bg = NA,
        gmm_pi_bg = NA,
        gmm_mu_col = NA,
        gmm_sd_col = NA,
        gmm_pi_col = NA,
        touches_edge = touches_edge,
        notes = "GMM_fit_failed",
        y_row = y_row + crop_info$top - 1,
        row_index = row_index,
        source_centers = paste(row_data$x, collapse = ";")
      )
      next
    }
    
    # Identify colony component (higher mean)
    if (gmm_result$mu[1] > gmm_result$mu[2]) {
      col_idx <- 1
      bg_idx <- 2
    } else {
      col_idx <- 2
      bg_idx <- 1
    }
    
    # Get posteriors for colony component
    gamma <- gmm_result$posteriors[, col_idx]
    
    # Compute metrics
    area_weighted <- compute_weighted_area(window_adjusted, window_baseline, gamma)
    mu_lift <- sum(gamma * (window_adjusted - window_baseline)) / sum(gamma)
    
    # SNR calculation
    residuals <- window_adjusted - window_baseline
    snr <- mu_lift / mad(residuals, na.rm = TRUE)
    
    resp_mean <- mean(gamma)
    
    # Presence check
    present <- presence_check(area_weighted, snr, resp_mean, touches_edge,
                              area_min, snr_min, resp_min)
    
    # Store results
    colony_results[[i]] <- data.frame(
      colony_id = colony_id,
      present = present,
      center_x = window$center + crop_info$left - 1,
      left = window$left + crop_info$left - 1,
      right = window$right + crop_info$left - 1,
      width_px = window$right - window$left + 1,
      area_weighted = area_weighted,
      mu_lift = mu_lift,
      snr = snr,
      resp_mean = resp_mean,
      gmm_mu_bg = gmm_result$mu[bg_idx],
      gmm_sd_bg = gmm_result$sd[bg_idx],
      gmm_pi_bg = gmm_result$pi[bg_idx],
      gmm_mu_col = gmm_result$mu[col_idx],
      gmm_sd_col = gmm_result$sd[col_idx],
      gmm_pi_col = gmm_result$pi[col_idx],
      touches_edge = touches_edge,
      notes = ifelse(present, "OK", "not_present"),
      y_row = y_row + crop_info$top - 1,
      row_index = row_index,
      source_centers = paste(row_data$x, collapse = ";")
    )
  }
  
  # Combine results
  colony_df <- do.call(rbind, colony_results)
  
  # Overview data
  overview_df <- data.frame(
    y_row = y_row + crop_info$top - 1,
    row_index = row_index,
    n_colonies_total = nrow(colony_df),
    n_colonies_present = sum(colony_df$present),
    crop_left = crop_info$left,
    crop_right = crop_info$right,
    crop_top = crop_info$top,
    crop_bottom = crop_info$bottom,
    detrend_k = detrend_k,
    area_min = area_min,
    snr_min = snr_min,
    resp_min = resp_min,
    snap_radius = snap_radius,
    r_intra_min = r_intra_min,
    r_intra_max = r_intra_max,
    margin = margin,
    washed_path = washed_path,
    background_path = background_path,
    gitter_dat_path = gitter_dat_path
  )
  
  # Save CSV files
  readr::write_csv(colony_df, file.path(out_dir, "gmm_scanline_per_colony.csv"))
  readr::write_csv(overview_df, file.path(out_dir, "gmm_scanline_overview.csv"))
  
  # Return results
  list(
    colony_data = colony_df,
    overview_data = overview_df,
    scanline = list(
      adjusted = adjusted,
      baseline = baseline,
      detrended = detrended,
      x = seq_along(adjusted) + crop_info$left - 1,
      y_row = y_row + crop_info$top - 1
    ),
    windows = windows,
    crop_info = crop_info
  )
}

#' Read Gitter DAT File
#'
#' Reads a gitter .dat file, skipping comment lines and parsing colony data.
#'
#' @param path `character(1)`. Path to the .dat file.
#' @return A data.frame with colony information.
read_gitter_dat <- function(path) {
  if (!file.exists(path)) {
    stop(paste("Gitter file not found:", path))
  }
  
  # Read all lines to inspect the file structure
  all_lines <- readLines(path)
  cat("Total lines in file:", length(all_lines), "\n")
  
  # Show first few lines for debugging
  cat("First 6 lines:\n")
  for (i in 1:min(6, length(all_lines))) {
    cat(sprintf("Line %d: %s\n", i, all_lines[i]))
  }
  
  # Count comment lines
  comment_lines <- startsWith(all_lines, "#")
  num_comments <- sum(comment_lines)
  cat("Number of comment lines:", num_comments, "\n")
  
  # The gitter format appears to be:
  # Comment lines starting with #
  # Then header line with column names  
  # Then data lines
  
  # Try reading with different approaches
  
  # Approach 1: Assume header is immediately after comments
  lines_to_skip <- num_comments
  cat("Trying to read with skip =", lines_to_skip, "\n")
  
  tryCatch({
    # First try with readr
    data <- readr::read_tsv(path, skip = lines_to_skip, show_col_types = FALSE)
    cat("Successfully read", nrow(data), "rows with readr\n")
    cat("Column names:", paste(names(data), collapse = ", "), "\n")
    
    # Check if we have the required columns
    required_cols <- c("xl", "xr", "yt", "yb")
    if (all(required_cols %in% names(data))) {
      cat("All required columns found!\n")
      return(data)
    } else {
      missing <- setdiff(required_cols, names(data))
      cat("Missing columns with readr:", paste(missing, collapse = ", "), "\n")
    }
  }, error = function(e) {
    cat("readr failed:", conditionMessage(e), "\n")
  })
  
  # Approach 2: Try base R
  cat("Trying base R read.delim...\n")
  tryCatch({
    data <- read.delim(path, skip = lines_to_skip, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
    cat("Successfully read", nrow(data), "rows with base R\n")
    cat("Column names:", paste(names(data), collapse = ", "), "\n")
    
    # Check if we have the required columns
    required_cols <- c("xl", "xr", "yt", "yb")
    if (all(required_cols %in% names(data))) {
      cat("All required columns found!\n")
      return(data)
    } else {
      missing <- setdiff(required_cols, names(data))
      cat("Missing columns with base R:", paste(missing, collapse = ", "), "\n")
    }
  }, error = function(e) {
    cat("Base R failed:", conditionMessage(e), "\n")
  })
  
  # Approach 3: Manual parsing - assume fixed column structure
  cat("Trying manual parsing with known column structure...\n")
  
  # Based on the gitter file format, the columns should be:
  # row col size circularity flags x y xl xr yt yb
  col_names <- c("row", "col", "size", "circularity", "flags", "x", "y", "xl", "xr", "yt", "yb")
  
  tryCatch({
    # Read data lines (skip comments)
    data_lines <- all_lines[!comment_lines]
    
    # Skip the header line if it exists (first non-comment line might be header)
    if (length(data_lines) > 0 && !grepl("^\\d", data_lines[1])) {
      cat("Skipping apparent header line:", data_lines[1], "\n")
      data_lines <- data_lines[-1]
    }
    
    # Parse data manually
    data_list <- list()
    for (col_name in col_names) {
      data_list[[col_name]] <- character(length(data_lines))
    }
    
    for (i in seq_along(data_lines)) {
      line <- data_lines[i]
      if (nchar(trimws(line)) == 0) next  # Skip empty lines
      
      parts <- strsplit(line, "\t")[[1]]
      
      # Handle missing values or extra columns
      for (j in seq_along(col_names)) {
        if (j <= length(parts)) {
          data_list[[col_names[j]]][i] <- parts[j]
        } else {
          data_list[[col_names[j]]][i] <- NA
        }
      }
    }
    
    # Convert to data frame
    data <- data.frame(data_list, stringsAsFactors = FALSE)
    
    # Convert numeric columns
    numeric_cols <- c("row", "col", "size", "circularity", "x", "y", "xl", "xr", "yt", "yb")
    for (col in numeric_cols) {
      if (col %in% names(data)) {
        data[[col]] <- as.numeric(data[[col]])
      }
    }
    
    # Remove empty rows
    data <- data[!is.na(data$row), ]
    
    cat("Manual parsing successful:", nrow(data), "rows\n")
    cat("Column names:", paste(names(data), collapse = ", "), "\n")
    
    # Check required columns
    required_cols <- c("xl", "xr", "yt", "yb")
    if (all(required_cols %in% names(data))) {
      cat("All required columns found!\n")
      return(data)
    } else {
      missing <- setdiff(required_cols, names(data))
      cat("Still missing columns:", paste(missing, collapse = ", "), "\n")
    }
    
  }, error = function(e) {
    cat("Manual parsing failed:", conditionMessage(e), "\n")
  })
  
  stop("Could not read gitter file with any method")
}

#' Compute Global Crop from Gitter Data
#'
#' Computes crop boundaries from the union of all gitter bounding boxes
#' plus a margin.
#'
#' @param gitter_df `data.frame`. Gitter data with xl, xr, yt, yb columns.
#' @param img_dim `named numeric`. Image dimensions with width and height.
#' @param margin `integer(1)`. Margin to add around bounding box.
#' @return A list with left, right, top, bottom crop coordinates.
compute_global_crop <- function(gitter_df, img_dim, margin) {
  # Check for valid gitter data
  if (nrow(gitter_df) == 0) {
    stop("No gitter data provided")
  }
  
  # Check for required columns
  required_cols <- c("xl", "xr", "yt", "yb")
  missing_cols <- setdiff(required_cols, names(gitter_df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Compute boundaries with safety checks
  left <- max(1, min(gitter_df$xl, na.rm = TRUE) - margin)
  right <- min(img_dim["width"], max(gitter_df$xr, na.rm = TRUE) + margin)
  top <- max(1, min(gitter_df$yt, na.rm = TRUE) - margin)
  bottom <- min(img_dim["height"], max(gitter_df$yb, na.rm = TRUE) + margin)
  
  # Ensure left < right and top < bottom
  if (left >= right) {
    stop(paste("Invalid crop boundaries: left (", left, ") >= right (", right, ")"))
  }
  if (top >= bottom) {
    stop(paste("Invalid crop boundaries: top (", top, ") >= bottom (", bottom, ")"))
  }
  
  # Check for reasonable crop size (prevent extremely large crops)
  crop_width <- right - left + 1
  crop_height <- bottom - top + 1
  
  if (crop_width > 10000 || crop_height > 10000) {
    stop(paste("Crop size too large: width =", crop_width, "height =", crop_height, 
               "\nImage dimensions:", img_dim["width"], "x", img_dim["height"],
               "\nGitter bounds: xl=[", min(gitter_df$xl, na.rm = TRUE), ",", max(gitter_df$xr, na.rm = TRUE), 
               "] yt=[", min(gitter_df$yt, na.rm = TRUE), ",", max(gitter_df$yb, na.rm = TRUE), "]"))
  }
  
  crop_info <- list(left = as.integer(left), right = as.integer(right), 
                    top = as.integer(top), bottom = as.integer(bottom))
  
  # Debug output
  cat("Crop boundaries: left =", crop_info$left, "right =", crop_info$right, 
      "top =", crop_info$top, "bottom =", crop_info$bottom, "\n")
  cat("Crop size:", crop_width, "x", crop_height, "pixels\n")
  
  return(crop_info)
}

#' Extract Scanline from Image
#'
#' Extracts a horizontal scanline at the specified y-coordinate.
#'
#' @param img `matrix`. Grayscale image matrix.
#' @param y `integer(1)`. Y-coordinate for scanline extraction.
#' @return A numeric vector representing the scanline intensities.
extract_scanline <- function(img, y) {
  if (y < 1 || y > dim(img)[2]) {
    stop(paste("y coordinate", y, "out of bounds"))
  }
  return(as.numeric(img[, y]))
}

#' Detrend with Running Median
#'
#' Applies running median smoothing to remove slow trends.
#'
#' @param y `numeric`. Signal to detrend.
#' @param k `integer(1)`. Window size for running median (must be odd).
#' @return A numeric vector of the same length as y representing the baseline.
detrend_runmed <- function(y, k) {
  if (k %% 2 == 0) k <- k + 1  # Ensure odd
  return(runmed(y, k, endrule = "constant"))
}

#' Define Windows from Centers
#'
#' Defines analysis windows around colony centers using midpoints and
#' optional snapping to local minima.
#'
#' @param centers `numeric`. X-coordinates of colony centers.
#' @param n `integer(1)`. Length of the scanline.
#' @param detrended `numeric`. Detrended signal for minima finding.
#' @param snap_radius `integer(1)`. Radius for snapping to local minima.
#' @param r_intra_min `integer(1)`. Minimum distance for minima search.
#' @param r_intra_max `integer(1)`. Maximum distance for minima search.
#' @return A list of windows, each with left, right, center elements.
define_windows_from_centers <- function(centers, n, detrended, snap_radius,
                                        r_intra_min, r_intra_max) {
  centers_sorted <- sort(centers)
  windows <- list()
  
  for (i in seq_along(centers_sorted)) {
    center <- centers_sorted[i]
    
    # Determine initial bounds using midpoints
    if (i == 1) {
      left_bound <- 1
    } else {
      left_bound <- round((centers_sorted[i-1] + center) / 2)
    }
    
    if (i == length(centers_sorted)) {
      right_bound <- n
    } else {
      right_bound <- round((center + centers_sorted[i+1]) / 2)
    }
    
    # Snap to local minima within snap_radius
    left_snap <- snap_to_minimum(detrended, left_bound, snap_radius)
    right_snap <- snap_to_minimum(detrended, right_bound, snap_radius)
    
    # Apply intra-colony distance constraints
    left_final <- max(left_snap, center - r_intra_max)
    left_final <- min(left_final, center - r_intra_min)
    left_final <- max(1, left_final)
    
    right_final <- min(right_snap, center + r_intra_max)
    right_final <- max(right_final, center + r_intra_min)
    right_final <- min(n, right_final)
    
    windows[[i]] <- list(
      left = as.integer(left_final),
      right = as.integer(right_final),
      center = as.integer(center)
    )
  }
  
  return(windows)
}

#' Snap to Local Minimum
#'
#' Finds the nearest local minimum within a specified radius.
#'
#' @param signal `numeric`. Signal to search for minima.
#' @param pos `integer(1)`. Initial position.
#' @param radius `integer(1)`. Search radius.
#' @return Integer position of the nearest minimum.
snap_to_minimum <- function(signal, pos, radius) {
  n <- length(signal)
  left <- max(1, pos - radius)
  right <- min(n, pos + radius)
  
  search_region <- signal[left:right]
  min_idx <- which.min(search_region)
  
  return(left + min_idx - 1)
}

#' Fit GMM to Window
#'
#' Fits a 2-component Gaussian mixture model to intensity values in a window.
#'
#' @param y `numeric`. Intensity values to fit.
#' @param k `integer(1)`. Number of components (fixed at 2).
#' @return A list with mu, sd, pi, and posteriors, or NULL if fit fails.
fit_gmm_window <- function(y, k = 2) {
  if (length(y) < 10 || var(y) < 1e-6) {
    return(NULL)  # Insufficient data or no variation
  }
  
  # Initialize with quantiles
  init_mu <- c(quantile(y, 0.2), quantile(y, 0.8))
  init_lambda <- c(0.5, 0.5)
  
  tryCatch({
    fit <- mixtools::normalmixEM(y, k = k, lambda = init_lambda, mu = init_mu,
                                 maxit = 1000, epsilon = 1e-8)
    
    # Extract parameters
    mu <- fit$mu
    sigma <- fit$sigma
    lambda <- fit$lambda
    
    # Compute posteriors
    posteriors <- fit$posterior
    
    list(
      mu = mu,
      sd = sigma,
      pi = lambda,
      posteriors = posteriors
    )
  }, error = function(e) {
    # Fallback: try with mclust if available
    if (requireNamespace("mclust", quietly = TRUE)) {
      tryCatch({
        fit <- mclust::Mclust(y, G = 2, modelNames = "V")
        if (!is.null(fit)) {
          posteriors <- fit$z
          list(
            mu = fit$parameters$mean,
            sd = sqrt(fit$parameters$variance$sigmasq),
            pi = fit$parameters$pro,
            posteriors = posteriors
          )
        } else {
          NULL
        }
      }, error = function(e2) NULL)
    } else {
      NULL
    }
  })
}

#' Compute Weighted Area
#'
#' Computes the responsibility-weighted area above baseline.
#'
#' @param adjusted `numeric`. Background-corrected signal.
#' @param baseline `numeric`. Baseline signal.
#' @param gamma `numeric`. Posterior responsibilities for colony component.
#' @return Numeric value representing the weighted area.
compute_weighted_area <- function(adjusted, baseline, gamma) {
  lift <- adjusted - baseline
  area <- sum(gamma * lift)
  return(area)
}

#' Presence Check
#'
#' Determines if a colony is present based on multiple criteria.
#'
#' @param area `numeric(1)`. Weighted area.
#' @param snr `numeric(1)`. Signal-to-noise ratio.
#' @param resp_mean `numeric(1)`. Mean responsibility.
#' @param touches_edge `logical(1)`. Whether window touches image edge.
#' @param area_min `numeric(1)`. Minimum area threshold.
#' @param snr_min `numeric(1)`. Minimum SNR threshold.
#' @param resp_min `numeric(1)`. Minimum responsibility threshold.
#' @return Logical indicating presence.
presence_check <- function(area, snr, resp_mean, touches_edge,
                           area_min, snr_min, resp_min) {
  if (is.na(area) || is.na(snr) || is.na(resp_mean)) {
    return(FALSE)
  }
  
  area_ok <- area >= area_min
  snr_ok <- snr >= snr_min
  resp_ok <- resp_mean >= resp_min
  edge_ok <- !touches_edge
  
  return(area_ok && snr_ok && resp_ok && edge_ok)
}