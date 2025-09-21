#' FFT-Based Bleed Suppression for Halo Assay Images
#'
#' Removes periodic bleed from tetramer colony lattices and low-frequency
#' illumination artifacts using 2D Fourier Transform filtering. This function
#' is specifically designed for halo assay plates where colonies are arranged
#' in regular tetramers and halos bleed into each other.
#'
#' @param washed_img `cimg`. The washed image to clean.
#' @param background_img `cimg`. Background image for subtraction (optional).
#' @param hp_radius `numeric(1)`. High-pass filter radius to remove low-frequency
#'   haze (default 6). Smaller values remove more low frequencies.
#' @param hp_slope `numeric(1)`. High-pass filter slope/softness (default 3).
#'   Higher values create sharper cutoffs.
#' @param notch_radius `numeric(1)`. Radius of notch filters for colony lattice
#'   peaks (default 2). Larger values remove more around each peak.
#' @param harmonics `numeric`. Which harmonics of the lattice to notch (default 1:3).
#'   Higher harmonics correspond to finer periodic structures.
#' @param spacing_hint `numeric(2)`. Manual hint for colony spacing [x_spacing, y_spacing]
#'   in pixels. If NULL, auto-detects from FFT peaks.
#' @param auto_detect_lattice `logical(1)`. Whether to automatically detect
#'   lattice peaks in the FFT (default TRUE).
#' @param background_subtract `logical(1)`. Whether to subtract background before
#'   FFT filtering (default FALSE).
#'
#' @return `list` containing:
#'   \itemize{
#'     \item `cleaned_img`: FFT-filtered image ready for HAIP analysis
#'     \item `detected_peaks`: Detected lattice peaks in frequency domain
#'     \item `filter_mask`: The frequency domain filter that was applied
#'     \item `fft_magnitude`: Magnitude of the FFT for visualization
#'   }
#' @export
fft_bleed_suppression <- function(washed_img,
                                 background_img = NULL,
                                 hp_radius = 6,
                                 hp_slope = 3,
                                 notch_radius = 2,
                                 harmonics = 1:3,
                                 spacing_hint = NULL,
                                 auto_detect_lattice = TRUE,
                                 background_subtract = FALSE) {

  # Convert to grayscale if needed
  if (dim(washed_img)[4] > 1) {
    gray_img <- grayscale(washed_img)
  } else {
    gray_img <- washed_img
  }

  # Optional background subtraction
  if (background_subtract && !is.null(background_img)) {
    if (dim(background_img)[4] > 1) {
      bg_gray <- grayscale(background_img)
    } else {
      bg_gray <- background_img
    }
    # Subtract background with offset to avoid negative values
    gray_img <- gray_img - bg_gray + mean(bg_gray)
  }

  # Convert to matrix for FFT processing
  img_matrix <- as.matrix(gray_img[,,1,1])
  img_dims <- dim(img_matrix)

  # Apply window function to reduce edge artifacts
  window_x <- outer(rep(1, img_dims[1]),
                    0.5 * (1 - cos(2 * pi * (0:(img_dims[2]-1)) / (img_dims[2]-1))))
  window_y <- outer(0.5 * (1 - cos(2 * pi * (0:(img_dims[1]-1)) / (img_dims[1]-1))),
                    rep(1, img_dims[2]))
  window <- window_x * window_y

  windowed_img <- img_matrix * window

  # Calculate 2D FFT
  img_fft <- fft(windowed_img)
  fft_magnitude <- abs(img_fft)

  # Create frequency coordinate grids
  freq_x <- seq(0, img_dims[1]-1) - floor(img_dims[1]/2)
  freq_y <- seq(0, img_dims[2]-1) - floor(img_dims[2]/2)
  freq_grid_x <- outer(freq_x, rep(1, img_dims[2]))
  freq_grid_y <- outer(rep(1, img_dims[1]), freq_y)
  freq_radius <- sqrt(freq_grid_x^2 + freq_grid_y^2)

  # Shift FFT to center zero frequency
  img_fft_shifted <- fftshift(img_fft)
  fft_mag_shifted <- fftshift(fft_magnitude)

  # Initialize filter mask (starts as all-pass)
  filter_mask <- array(1, dim = img_dims)

  # 1. High-pass filter to remove low-frequency haze
  hp_filter <- 1 / (1 + (hp_radius / (freq_radius + 1e-8))^(2 * hp_slope))
  hp_filter[freq_radius == 0] <- 0  # Remove DC component
  filter_mask <- filter_mask * hp_filter

  # 2. Auto-detect lattice peaks or use manual spacing
  detected_peaks <- list()

  if (auto_detect_lattice) {
    # Find peaks in the FFT magnitude (excluding center)
    # Smooth the FFT magnitude to reduce noise
    smoothed_fft <- isoblur(as.cimg(fft_mag_shifted), sigma = 1)
    smoothed_matrix <- as.matrix(smoothed_fft[,,1,1])

    # Mask out center region
    center_mask <- freq_radius > (max(img_dims) * 0.05)
    masked_fft <- smoothed_matrix * center_mask

    # Find local maxima
    threshold <- quantile(masked_fft[masked_fft > 0], 0.95)
    peak_candidates <- which(masked_fft > threshold, arr.ind = TRUE)

    # Select strongest peaks that could represent lattice structure
    if (nrow(peak_candidates) > 0) {
      peak_values <- masked_fft[peak_candidates]
      sorted_indices <- order(peak_values, decreasing = TRUE)

      # Take top peaks that are not too close to each other
      selected_peaks <- peak_candidates[sorted_indices[1:min(10, length(sorted_indices))], , drop = FALSE]

      detected_peaks <- list(
        positions = selected_peaks,
        frequencies = data.frame(
          fx = freq_x[selected_peaks[,1]],
          fy = freq_y[selected_peaks[,2]],
          magnitude = masked_fft[selected_peaks]
        )
      )
    }
  }

  # 3. Apply notch filters for lattice suppression
  if (!is.null(spacing_hint)) {
    # Use manual spacing to place notches
    fundamental_fx <- img_dims[1] / spacing_hint[1]
    fundamental_fy <- img_dims[2] / spacing_hint[2]

    for (h in harmonics) {
      for (sx in c(-1, 1)) {
        for (sy in c(-1, 1)) {
          notch_fx <- sx * h * fundamental_fx
          notch_fy <- sy * h * fundamental_fy

          # Find corresponding indices
          fx_idx <- which.min(abs(freq_x - notch_fx))
          fy_idx <- which.min(abs(freq_y - notch_fy))

          # Create notch filter
          notch_dist <- sqrt((freq_grid_x - freq_x[fx_idx])^2 +
                           (freq_grid_y - freq_y[fy_idx])^2)
          notch_filter <- 1 / (1 + (notch_radius / (notch_dist + 1e-8))^8)
          filter_mask <- filter_mask * notch_filter
        }
      }
    }
  } else if (auto_detect_lattice && length(detected_peaks) > 0) {
    # Use detected peaks to place notches
    for (i in 1:min(nrow(detected_peaks$positions), 6)) {  # Limit to strongest peaks
      peak_pos <- detected_peaks$positions[i,]
      center_x <- freq_x[peak_pos[1]]
      center_y <- freq_y[peak_pos[2]]

      # Create notch at this frequency and its harmonics
      for (h in harmonics) {
        notch_dist <- sqrt((freq_grid_x - h * center_x)^2 +
                         (freq_grid_y - h * center_y)^2)
        notch_filter <- 1 / (1 + (notch_radius / (notch_dist + 1e-8))^8)
        filter_mask <- filter_mask * notch_filter

        # Also notch the negative frequency
        notch_dist_neg <- sqrt((freq_grid_x + h * center_x)^2 +
                             (freq_grid_y + h * center_y)^2)
        notch_filter_neg <- 1 / (1 + (notch_radius / (notch_dist_neg + 1e-8))^8)
        filter_mask <- filter_mask * notch_filter_neg
      }
    }
  }

  # Apply filter in frequency domain
  filtered_fft <- img_fft_shifted * filter_mask

  # Shift back and inverse FFT
  filtered_fft_unshifted <- ifftshift(filtered_fft)
  cleaned_matrix <- Re(fft(filtered_fft_unshifted, inverse = TRUE))

  # Convert back to cimg format
  cleaned_img <- as.cimg(cleaned_matrix)

  # Restore original image dimensions and properties
  if (dim(washed_img)[4] > 1) {
    # For color images, apply the same filtering to each channel
    cleaned_color <- washed_img
    for (ch in 1:dim(washed_img)[4]) {
      channel_matrix <- as.matrix(washed_img[,,1,ch])
      channel_windowed <- channel_matrix * window
      channel_fft <- fft(channel_windowed)
      channel_fft_shifted <- fftshift(channel_fft)
      channel_filtered <- channel_fft_shifted * filter_mask
      channel_filtered_unshifted <- ifftshift(channel_filtered)
      channel_cleaned <- Re(fft(channel_filtered_unshifted, inverse = TRUE))
      cleaned_color[,,1,ch] <- channel_cleaned
    }
    cleaned_img <- cleaned_color
  }

  return(list(
    cleaned_img = cleaned_img,
    detected_peaks = detected_peaks,
    filter_mask = filter_mask,
    fft_magnitude = fft_mag_shifted
  ))
}

#' Helper function to shift FFT for centered display
#' @param x Complex matrix to shift
#' @return Shifted matrix
fftshift <- function(x) {
  dims <- dim(x)
  shift_x <- floor(dims[1]/2)
  shift_y <- floor(dims[2]/2)

  result <- x
  if (shift_x > 0) {
    result <- result[c((shift_x+1):dims[1], 1:shift_x), ]
  }
  if (shift_y > 0) {
    result <- result[, c((shift_y+1):dims[2], 1:shift_y)]
  }
  return(result)
}

#' Helper function to reverse FFT shift
#' @param x Complex matrix to unshift
#' @return Unshifted matrix
ifftshift <- function(x) {
  dims <- dim(x)
  shift_x <- ceiling(dims[1]/2)
  shift_y <- ceiling(dims[2]/2)

  result <- x
  if (shift_x > 0) {
    result <- result[c((shift_x+1):dims[1], 1:shift_x), ]
  }
  if (shift_y > 0) {
    result <- result[, c((shift_y+1):dims[2], 1:shift_y)]
  }
  return(result)
}
