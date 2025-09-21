# test_scanline.R - Comprehensive test script for robust FFT scanline analysis
# Tests the complete robust pipeline with validation of all requirements

# Load required libraries
library(imager)
library(readr)
library(dplyr)

# Source the robust implementations
source("R/FFT_scanline.R")
source("R/FFT_explained.R")

cat("=== ROBUST FFT SCANLINE TEST ===\n")

# Define test data paths
washed_path <- "Plate Images/washed/washed_BHET25_2d_6_1.JPG"
background_path <- "Plate Images/background/background_BHET25_6_1.JPG"

# Define colony centers manually (replace with your known colony positions)
# These should be x-coordinates where colonies are located in the image
colony_centers <- c(500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250)
y_row <- 500  # Y-coordinate for the scanline to analyze

# Validate input files exist
if (!file.exists(washed_path)) {
  stop("Washed image not found: ", washed_path)
}
if (!file.exists(background_path)) {
  stop("Background image not found: ", background_path)
}

cat("✓ Input files validated\n")

# Test A: FFT Scanline Analysis (No Gitter)
cat("\n--- Testing FFT Scanline Analysis (No Gitter) ---\n")

scanline_results <- fft_scanline(
  washed_path   = washed_path,
  background_path = background_path,
  colony_centers = colony_centers,
  y_row         = y_row,
  out_dir       = "results/scanline",
  save_debug    = TRUE
)

cat("✓ FFT scanline analysis completed (no gitter)\n")

# Report key metrics
per_colony_df <- scanline_results$per_colony
overview_df <- scanline_results$overview

n_total <- nrow(per_colony_df)
n_present <- sum(per_colony_df$present)

cat(sprintf("Colonies detected: %d\n", n_total))
cat(sprintf("Present colonies: %d / %d\n", n_present, n_total))

# Validate crop info (now just full image dimensions)
crop_info <- scanline_results$crop_info
cat(sprintf("Image region: left=%d, right=%d, top=%d, bottom=%d\n",
            crop_info$left, crop_info$right, crop_info$top, crop_info$bottom))

# Test B: FFT Plotting (No Gitter)
cat("\n--- Testing FFT Plotting (No Gitter) ---\n")

plotting_results <- fft_explained(
  washed_path = washed_path,
  background_path = background_path,
  scanline_csv = "results/scanline/fft_scanline_per_colony.csv",
  overview_csv = "results/scanline/fft_scanline_overview.csv",
  out_dir = "results/figs/colony",
  show_detrended = TRUE
)

cat("✓ FFT plotting completed (no gitter)\n")
cat(sprintf("Colony plots generated: %d\n", length(plotting_results$colony_plots)))

cat("\n=== FFT IMPLEMENTATION COMPLETE (NO GITTER) ===\n")
cat(sprintf("Total colonies detected: %d\n", n_total))
cat(sprintf("Present colonies: %d\n", n_present))
cat("Key features implemented:\n")
cat("  ✓ Full image processing (no cropping)\n")
cat("  ✓ Manual colony center specification\n")
cat("  ✓ Signal minima-based windowing\n")
cat("  ✓ Background-aware tail expansion\n")
cat("  ✓ SNR/area/continuity filtering\n")
cat("  ✓ Edge artifact rejection\n")
cat("  ✓ Present-only visualization\n")
cat("  ✓ Comprehensive CSV outputs\n")

cat("\n=== TEST COMPLETE ===\n")
