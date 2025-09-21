#' Align Images Using Phase Correlation
#'
#' Aligns two images using phase correlation to find the optimal translation.
#' This function is essential for ensuring that background subtraction and
#' FFT filtering work correctly across different images of the same plate.
#'
#' @param reference_img `cimg`. Reference image (typically background).
#' @param target_img `cimg`. Image to be aligned to the reference.
#'
#' @return `list` containing:
#'   \itemize{
#'     \item `aligned_img`: The target image aligned to the reference
#'     \item `shift_x`: Horizontal shift applied (pixels)
#'     \item `shift_y`: Vertical shift applied (pixels)
#'   }
#' @export
align_images_phase_correlation <- function(reference_img, target_img) {
  # Convert to grayscale if needed
  if (dim(reference_img)[4] > 1) {
    ref_gray <- grayscale(reference_img)
  } else {
    ref_gray <- reference_img
  }

  if (dim(target_img)[4] > 1) {
    target_gray <- grayscale(target_img)
  } else {
    target_gray <- target_img
  }

  # Ensure images are the same size
  ref_dims <- dim(ref_gray)
  target_dims <- dim(target_gray)

  if (!all(ref_dims[1:2] == target_dims[1:2])) {
    # Crop to minimum common size
    min_width <- min(ref_dims[1], target_dims[1])
    min_height <- min(ref_dims[2], target_dims[2])

    ref_gray <- imsub(ref_gray, x %inr% c(1, min_width), y %inr% c(1, min_height))
    target_gray <- imsub(target_gray, x %inr% c(1, min_width), y %inr% c(1, min_height))
  }

  # Convert to matrices for FFT
  ref_matrix <- as.matrix(ref_gray[,,1,1])
  target_matrix <- as.matrix(target_gray[,,1,1])

  # Apply window function to reduce edge effects
  window_x <- outer(rep(1, nrow(ref_matrix)),
                    0.5 * (1 - cos(2 * pi * (0:(ncol(ref_matrix)-1)) / (ncol(ref_matrix)-1))))
  window_y <- outer(0.5 * (1 - cos(2 * pi * (0:(nrow(ref_matrix)-1)) / (nrow(ref_matrix)-1))),
                    rep(1, ncol(ref_matrix)))
  window <- window_x * window_y

  ref_windowed <- ref_matrix * window
  target_windowed <- target_matrix * window

  # Calculate FFTs
  ref_fft <- fft(ref_windowed)
  target_fft <- fft(target_windowed)

  # Phase correlation
  cross_power_spectrum <- (ref_fft * Conj(target_fft)) / (abs(ref_fft * Conj(target_fft)) + 1e-8)
  correlation <- Re(fft(cross_power_spectrum, inverse = TRUE))

  # Find peak correlation
  peak_pos <- which(correlation == max(correlation), arr.ind = TRUE)[1,]

  # Calculate shifts (accounting for wraparound)
  shift_x <- peak_pos[1] - 1
  shift_y <- peak_pos[2] - 1

  # Adjust for wraparound (shifts > image_size/2 are negative shifts)
  if (shift_x > nrow(ref_matrix) / 2) {
    shift_x <- shift_x - nrow(ref_matrix)
  }
  if (shift_y > ncol(ref_matrix) / 2) {
    shift_y <- shift_y - ncol(ref_matrix)
  }

  # Apply the shift to the original target image (ensure integer shifts)
  aligned_img <- imshift(target_img, delta_x = as.integer(shift_x), delta_y = as.integer(shift_y),
                        boundary = 1)

  return(list(
    aligned_img = aligned_img,
    shift_x = shift_x,
    shift_y = shift_y
  ))
}

#' Simple Image Alignment Using Cross-Correlation
#'
#' A simpler alternative alignment method using imager's built-in functions.
#' May be more robust for images with large illumination differences.
#'
#' @param reference_img `cimg`. Reference image (typically background).
#' @param target_img `cimg`. Image to be aligned to the reference.
#' @param search_radius `numeric(1)`. Maximum search radius for alignment (default 50).
#'
#' @return `list` containing aligned image and shift information.
#' @export
align_images_simple <- function(reference_img, target_img, search_radius = 50) {
  # Convert to grayscale
  if (dim(reference_img)[4] > 1) {
    ref_gray <- grayscale(reference_img)
  } else {
    ref_gray <- reference_img
  }

  if (dim(target_img)[4] > 1) {
    target_gray <- grayscale(target_img)
  } else {
    target_gray <- target_img
  }

  # Simple grid search for best alignment
  best_correlation <- -Inf
  best_shift_x <- 0
  best_shift_y <- 0

  for (dx in -search_radius:search_radius) {
    for (dy in -search_radius:search_radius) {
      # Apply shift
      shifted <- imshift(target_gray, delta_x = dx, delta_y = dy, boundary = "neumann")

      # Calculate correlation (using a central crop to avoid edge effects)
      crop_size <- min(width(ref_gray), height(ref_gray)) * 0.8
      center_x <- width(ref_gray) / 2
      center_y <- height(ref_gray) / 2

      ref_crop <- imsub(ref_gray,
                       x %inr% c(center_x - crop_size/2, center_x + crop_size/2),
                       y %inr% c(center_y - crop_size/2, center_y + crop_size/2))
      shifted_crop <- imsub(shifted,
                           x %inr% c(center_x - crop_size/2, center_x + crop_size/2),
                           y %inr% c(center_y - crop_size/2, center_y + crop_size/2))

      # Calculate correlation coefficient
      correlation <- cor(as.vector(ref_crop), as.vector(shifted_crop))

      if (correlation > best_correlation) {
        best_correlation <- correlation
        best_shift_x <- dx
        best_shift_y <- dy
      }
    }
  }

  # Apply best shift to original image (ensure integer shifts)
  aligned_img <- imshift(target_img, delta_x = as.integer(best_shift_x), delta_y = as.integer(best_shift_y),
                        boundary = "neumann")

  return(list(
    aligned_img = aligned_img,
    shift_x = best_shift_x,
    shift_y = best_shift_y,
    correlation = best_correlation
  ))
}
