# Quick Start Example for FFT Bleed Suppression
# =============================================
# This script shows the simplest way to use the FFT enhancement

library(imager)
library(dplyr)
library(ggplot2)

# Source the FFT-enhanced HAIP functions
source("R/align_images.R")
source("R/fft_bleed_suppression.R")
source("R/fft_config.R")

# Source original HAIP functions
source("R/initialize_colony_data.R")
source("R/create_grid_boxes.R")
source("R/calculate_average_intensity.R")
source("R/calculate_edge_intensity.R")
source("R/classify_pixels.R")
source("R/rgb_to_gray.R")

# EXAMPLE 1: Basic usage with standard settings
# ==============================================

# Load your images
bg_img     <- imager::load.image("background/background/background_BHET25_5_1.JPG")
unwashed   <- imager::load.image("all_plate_images/all_plate_images/colony_images/unwashed/BHET25_8h_5_1.JPG")
washed_img <- imager::load.image("all_plate_images/all_plate_images/wash_images/washed_BHET25_8h_5_1.JPG")


# Create standard FFT configuration
config <- fft_preset_config("standard")

# Align images
washed_aligned <- align_images_phase_correlation(bg_img, washed_img)
unwashed_aligned <- align_images_phase_correlation(bg_img, unwashed)

# Apply FFT cleaning
fft_result <- fft_bleed_suppression(
  washed_img = washed_aligned$aligned_img,
  hp_radius = config$hp_radius,
  hp_slope = config$hp_slope,
  notch_radius = config$notch_radius,
  harmonics = config$harmonics
)

# Use cleaned image in HAIP pipeline
washed_clean <- fft_result$cleaned_img

# Continue with standard HAIP analysis...
# colony_data <- initialize_colony_data("path/to/your/data.dat", unwashed_aligned$aligned_img)
# colony_data <- create_grid_boxes(colony_data)
# etc.

# EXAMPLE 2: Custom parameters for heavily blended halos
# =======================================================

# Use aggressive settings
aggressive_config <- fft_preset_config("aggressive")

# Or create completely custom settings
custom_config <- create_fft_config(
  hp_radius = 4,        # Remove more low frequencies
  hp_slope = 4,         # Sharper cutoff
  notch_radius = 3,     # Larger notches
  harmonics = 1:4,      # More harmonics
  spacing_hint = c(150, 150)  # Manual spacing if known
)

# EXAMPLE 3: Comparing before and after
# =====================================

# Function to quickly compare original vs cleaned
compare_fft_cleaning <- function(original_img, cleaned_img, title = "FFT Comparison") {
  # Convert to grayscale for comparison
  orig_gray <- if(dim(original_img)[4] > 1) grayscale(original_img) else original_img
  clean_gray <- if(dim(cleaned_img)[4] > 1) grayscale(cleaned_img) else cleaned_img

  # Convert to data frames (sample for faster plotting)
  orig_df <- as.data.frame(orig_gray)
  clean_df <- as.data.frame(clean_gray)

  sample_size <- min(20000, nrow(orig_df))
  sample_indices <- sample(nrow(orig_df), sample_size)

  orig_df <- orig_df[sample_indices,]
  clean_df <- clean_df[sample_indices,]

  # Create comparison plot
  p1 <- ggplot(orig_df, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.2) + scale_y_reverse() +
    scale_color_gradient(low = "black", high = "white") +
    theme_minimal() + labs(title = "Original", x = "", y = "") +
    theme(legend.position = "none", axis.text = element_blank())

  p2 <- ggplot(clean_df, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.2) + scale_y_reverse() +
    scale_color_gradient(low = "black", high = "white") +
    theme_minimal() + labs(title = "FFT Cleaned", x = "", y = "") +
    theme(legend.position = "none", axis.text = element_blank())

  gridExtra::grid.arrange(p1, p2, ncol = 2, top = title)
}

# Use the comparison function - UNCOMMENT TO RUN
compare_fft_cleaning(washed_img, washed_clean, "Before vs After FFT")

# EXAMPLE 4: Parameter optimization workflow
# ==========================================

# Function to test different parameter combinations
test_parameters <- function(washed_img, param_list) {
  results <- list()

  for (i in seq_along(param_list)) {
    params <- param_list[[i]]

    cat("Testing configuration", i, "\n")

    fft_result <- fft_bleed_suppression(
      washed_img = washed_img,
      hp_radius = params$hp_radius,
      hp_slope = params$hp_slope,
      notch_radius = params$notch_radius,
      harmonics = params$harmonics
    )

    results[[i]] <- list(
      config = params,
      cleaned_img = fft_result$cleaned_img,
      detected_peaks = length(fft_result$detected_peaks)
    )
  }

  return(results)
}

# Example parameter sweep
# param_combinations <- list(
#   list(hp_radius = 4, hp_slope = 3, notch_radius = 2, harmonics = 1:3),
#   list(hp_radius = 6, hp_slope = 3, notch_radius = 2, harmonics = 1:3),
#   list(hp_radius = 8, hp_slope = 3, notch_radius = 2, harmonics = 1:3)
# )
#
# test_results <- test_parameters(washed_img, param_combinations)

cat("Quick start examples ready!\n")
cat("Uncomment and modify the example sections above to use with your data.\n")
