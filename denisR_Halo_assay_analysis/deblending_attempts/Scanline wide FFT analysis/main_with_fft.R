# HAIP with FFT Bleed Suppression Pipeline
# ===========================================
# This script integrates FFT-based bleed suppression into the HAIP workflow
# to improve halo detection accuracy in tetramer colony plates.

# Install & load required packages
install.packages(c("imager","dplyr","ggplot2","readr","gridExtra"))
library(imager); library(dplyr); library(ggplot2); library(readr); library(gridExtra)

# Install fftwtools if not available (for FFT operations)
if (!require(fftwtools, quietly = TRUE)) {
  install.packages("fftwtools")
  library(fftwtools)
}

# Source all HAIP functions
source("R/initialize_colony_data.R")
source("R/create_grid_boxes.R")
source("R/calculate_average_intensity.R")
source("R/calculate_edge_intensity.R")
source("R/classify_pixels.R")
source("R/plot_pixel_classification.R")  # optional
source("R/rgb_to_gray.R")

# Source new FFT and alignment functions
source("R/align_images.R")
source("R/fft_bleed_suppression.R")

# Configuration Parameters for FFT Bleed Suppression
# ===================================================
# Tuning parameters - adjust these based on your plate characteristics

# High-pass filter parameters (removes low-frequency illumination)
hp_radius <- 6        # Lower values remove more low frequencies (try 4-10)
hp_slope <- 3         # Higher values create sharper cutoffs (try 2-5)

# Notch filter parameters (removes tetramer lattice bleed)
notch_radius <- 2     # Size of notches around lattice peaks (try 1-4)
harmonics <- 1:3      # Which harmonics to suppress (try 1:2 or 1:4)

# Colony spacing hint (pixels) - set if you know your tetramer spacing
# Leave as NULL for automatic detection
spacing_hint <- NULL  # e.g., c(150, 150) for 150-pixel spacing

# Alignment parameters
use_phase_correlation <- TRUE  # Use phase correlation (TRUE) or simple alignment (FALSE)
alignment_search_radius <- 50  # For simple alignment only

# ===========================================
# MAIN PIPELINE
# ===========================================

cat("Loading images...\n")
# Load images (use forward slashes on Windows)
bg_img     <- imager::load.image("background/background/background_BHET25_5_1.JPG")
unwashed   <- imager::load.image("all_plate_images/all_plate_images/colony_images/unwashed/BHET25_8h_5_1.JPG")
washed_img <- imager::load.image("all_plate_images/all_plate_images/wash_images/washed_BHET25_8h_5_1.JPG")

cat("Image dimensions:\n")
cat("  Background:", dim(bg_img)[1:2], "\n")
cat("  Unwashed:  ", dim(unwashed)[1:2], "\n")
cat("  Washed:    ", dim(washed_img)[1:2], "\n")

# ===========================================
# STEP 1: Image Alignment
# ===========================================
cat("\nAligning images to background...\n")

if (use_phase_correlation) {
  # Align using phase correlation (more accurate)
  unwashed_aligned <- align_images_phase_correlation(bg_img, unwashed)
  washed_aligned <- align_images_phase_correlation(bg_img, washed_img)
  
  cat("  Unwashed shift: (", unwashed_aligned$shift_x, ",", unwashed_aligned$shift_y, ")\n")
  cat("  Washed shift:   (", washed_aligned$shift_x, ",", washed_aligned$shift_y, ")\n")
  
  unwashed <- unwashed_aligned$aligned_img
  washed_img <- washed_aligned$aligned_img
  
} else {
  # Align using simple cross-correlation
  unwashed_aligned <- align_images_simple(bg_img, unwashed, alignment_search_radius)
  washed_aligned <- align_images_simple(bg_img, washed_img, alignment_search_radius)
  
  cat("  Unwashed shift: (", unwashed_aligned$shift_x, ",", unwashed_aligned$shift_y, 
      ") correlation:", round(unwashed_aligned$correlation, 3), "\n")
  cat("  Washed shift:   (", washed_aligned$shift_x, ",", washed_aligned$shift_y, 
      ") correlation:", round(washed_aligned$correlation, 3), "\n")
  
  unwashed <- unwashed_aligned$aligned_img
  washed_img <- washed_aligned$aligned_img
}

# ===========================================
# STEP 2: FFT Bleed Suppression
# ===========================================
cat("\nApplying FFT bleed suppression...\n")

# Apply FFT cleaning to the washed image
fft_result <- fft_bleed_suppression(
  washed_img = washed_img,
  background_img = bg_img,
  hp_radius = hp_radius,
  hp_slope = hp_slope,
  notch_radius = notch_radius,
  harmonics = harmonics,
  spacing_hint = spacing_hint,
  auto_detect_lattice = is.null(spacing_hint),
  background_subtract = FALSE  # We handle this in HAIP pipeline
)

# Use the cleaned image for further processing
washed_clean <- fft_result$cleaned_img

cat("  FFT cleaning complete.\n")
if (length(fft_result$detected_peaks) > 0 && !is.null(fft_result$detected_peaks$frequencies)) {
  cat("  Detected", nrow(fft_result$detected_peaks$frequencies), "lattice peaks in FFT\n")
}

# ===========================================
# STEP 3: HAIP Classification Pipeline
# ===========================================
cat("\nRunning HAIP classification pipeline...\n")

# Prepare colony data (expects a gitter .dat file next to the unwashed image)
colony_dat_path <- "all_plate_images/all_plate_images/colony_images/unwashed/BHET12.5_8h_5_1.JPG.dat"

if (file.exists(colony_dat_path)) {
  # Initialize colony data
  colony_data <- initialize_colony_data(colony_dat_path, unwashed)
  colony_data <- create_grid_boxes(colony_data)
  colony_data <- calculate_average_intensity(bg_img, colony_data)
  
  # Compute edge normalization using the FFT-cleaned washed image
  edge_norm <- calculate_edge_intensity(bg_img, washed_clean, edge_width = 150)$difference
  cat("  Edge normalization factor:", round(edge_norm, 3), "\n")
  
  # Run classification on the FFT-cleaned image
  res <- classify_pixels(washed_clean, colony_data, edge_norm)
  
  cat("  Classification complete.\n")
  cat("  Total colonies processed:", nrow(res$colony_summary), "\n")
  
  # ===========================================
  # STEP 4: Results and Visualization
  # ===========================================
  cat("\nResults summary:\n")
  print(res$total_stats)
  
  cat("\nFirst few colony summaries:\n")
  print(head(res$colony_summary))
  
  # Quick visualization of classified pixels
  cat("\nGenerating classification plot...\n")
  classification_plot <- ggplot(res$classified_data, aes(x = x, y = y, color = label)) +
    geom_point(size = 0.3, alpha = 0.7) + 
    scale_y_reverse() + 
    theme_minimal() +
    labs(title = "Pixel Classification (FFT-cleaned)",
         subtitle = paste("Processed", nrow(res$colony_summary), "colonies"),
         x = "X position", y = "Y position", color = "Classification") +
    theme(legend.position = "bottom")
  
  print(classification_plot)
  
  # Optional: Create comparison plots
  cat("\nCreating before/after comparison...\n")
  
  # Convert images to data frames for plotting
  washed_df <- as.data.frame(grayscale(washed_img))
  washed_clean_df <- as.data.frame(grayscale(washed_clean))
  
  # Sample pixels for faster plotting
  sample_size <- min(50000, nrow(washed_df))
  sample_indices <- sample(nrow(washed_df), sample_size)
  
  comparison_plot <- grid.arrange(
    ggplot(washed_df[sample_indices,], aes(x = x, y = y, color = value)) +
      geom_point(size = 0.2) + scale_y_reverse() + scale_color_gradient(low = "black", high = "white") +
      theme_minimal() + labs(title = "Original Washed Image") + theme(legend.position = "none"),
    
    ggplot(washed_clean_df[sample_indices,], aes(x = x, y = y, color = value)) +
      geom_point(size = 0.2) + scale_y_reverse() + scale_color_gradient(low = "black", high = "white") +
      theme_minimal() + labs(title = "FFT-Cleaned Image") + theme(legend.position = "none"),
    
    ncol = 2
  )
  
} else {
  cat("Warning: Colony data file not found at:", colony_dat_path, "\n")
  cat("Please ensure the .dat file exists next to your unwashed image.\n")
  cat("Continuing with FFT cleaning demonstration only...\n")
  
  # Still show the FFT cleaning results
  washed_df <- as.data.frame(grayscale(washed_img))
  washed_clean_df <- as.data.frame(grayscale(washed_clean))
  
  sample_size <- min(50000, nrow(washed_df))
  sample_indices <- sample(nrow(washed_df), sample_size)
  
  comparison_plot <- grid.arrange(
    ggplot(washed_df[sample_indices,], aes(x = x, y = y, color = value)) +
      geom_point(size = 0.2) + scale_y_reverse() + scale_color_gradient(low = "black", high = "white") +
      theme_minimal() + labs(title = "Original Washed Image") + theme(legend.position = "none"),
    
    ggplot(washed_clean_df[sample_indices,], aes(x = x, y = y, color = value)) +
      geom_point(size = 0.2) + scale_y_reverse() + scale_color_gradient(low = "black", high = "white") +
      theme_minimal() + labs(title = "FFT-Cleaned Image") + theme(legend.position = "none"),
    
    ncol = 2
  )
}

cat("\nPipeline complete!\n")
cat("\n===========================================\n")
cat("TUNING TIPS:\n")
cat("- If halos look over-smoothed: decrease hp_radius (try 4-5)\n")
cat("- If lattice bleed remains: increase notch_radius (try 3-4)\n")
cat("- If artifacts appear: decrease hp_slope or notch harmonics\n")
cat("- For manual tuning: set spacing_hint to your tetramer spacing\n")
cat("===========================================\n")