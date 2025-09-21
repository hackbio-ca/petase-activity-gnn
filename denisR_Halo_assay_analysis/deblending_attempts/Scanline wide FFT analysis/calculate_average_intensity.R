#' Calculate Average Pixel Intensity in Bounding Boxes
#'
#' For each bounding box defined by coordinates (new_xl, new_xr, new_yt, new_yb),
#' calculates the mean, median, and standard deviation of pixel intensities within
#' that area using the imager package for efficient computation.
#'
#' Requires `create_grid_boxes()` to be called first to establish
#' bounding box coordinates.
#'
#' Outputs a warning message if there are no pixels within the bounding box
#' created
#'
#' @param img `cimg`. Image object from the imager package, typically loaded
#'   with `load.image()` or created with other imager functions.
#' @param colony_data `data.frame`. Colony data with bounding box coordinates
#'   (new_xl, new_xr, new_yt, new_yb) created by `create_grid_boxes()`.
#' @param prefix `character(1)`. Prefix to add to the new intensity columns.
#'   Default is "bg" (background).
#'
#' @return A `data.frame` containing the original colony data plus three new
#'   columns with the specified prefix: mean_intensity, median_intensity, and
#'   sd_intensity, containing the statistical measures of pixel intensity
#'   within each bounding box.
#' @export
calculate_average_intensity <- function(img, colony_data, prefix = "bg") {

  # Convert to grayscale if image is color
  if (dim(img)[4] > 1) {
    gray_img <- imager::grayscale(img)
  } else {
    gray_img <- img
  }

  # Initialize columns
  colony_data$mean_intensity <- NA
  colony_data$median_intensity <- NA
  colony_data$sd_intensity <- NA

  # Calculate statistics for each bounding box
  for(i in 1:nrow(colony_data)) {

    # Get bounding box coordinates
    xl <- colony_data$new_xl[i]
    xr <- colony_data$new_xr[i]
    yt <- colony_data$new_yt[i]
    yb <- colony_data$new_yb[i]

    # Check if any coordinates are NA
    if (any(is.na(c(xl, xr, yt, yb)))) {
      next
    }

    # Round to ensure integer indices
    xl <- round(xl)
    xr <- round(xr)
    yt <- round(yt)
    yb <- round(yb)

    # Ensure coordinates are within image bounds
    xl <- max(1, xl)
    xr <- min(width(gray_img), xr)
    yt <- max(1, yt)
    yb <- min(height(gray_img), yb)

    # Check if valid bounding box exists
    if (xl <= xr && yt <= yb) {

      # Extract the region of interest using imager's imsub function
      roi <- imsub(gray_img, x %inr% c(xl, xr), y %inr% c(yt, yb))

      # Convert to vector for statistical calculations
      roi_values <- as.vector(roi)

      if(length(roi_values) > 0) {
        colony_data$mean_intensity[i] <- mean(roi_values, na.rm = TRUE)
        colony_data$median_intensity[i] <- median(roi_values, na.rm = TRUE)
        colony_data$sd_intensity[i] <- sd(roi_values, na.rm = TRUE)
      } else {
        cat(sprintf("Warning: No pixels found in box %d\n", i))
      }
    } else {
      cat(sprintf("Warning: Invalid bounding box coordinates for box %d\n", i))
    }
  }

  # Rename columns with prefix
  colony_data <- colony_data %>%
    rename_with(~ paste0(prefix, "_", .x),
                c(mean_intensity, median_intensity, sd_intensity))

  return(colony_data)
}
