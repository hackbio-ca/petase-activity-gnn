# Enhanced test script for FFT_explained.R
# This script tests the enhanced FFT explanation functionality with dual-frequency detection

# Load required libraries
library(imager)
library(ggplot2)
library(gridExtra)
library(readr)
library(dplyr)

# Check for pracma package
if (!requireNamespace("pracma", quietly = TRUE)) {
  cat("Installing pracma package for enhanced functionality...\n")
  install.packages("pracma")
}
library(pracma)

# Source the enhanced FFT explanation script
source("R/FFT_explained.R")

# Set up paths for testing
washed_path <- "spwsh/washed_BHET25_2d_6_1.JPG"
background_path <- "spbg/background_BHET25_6_1.JPG"
gitter_dat_path <- "spuwsh/BHET25_2d_6_1.JPG.dat"

# Create results directory if it doesn't exist
dir.create("results", showWarnings = FALSE)
dir.create("results/figs", showWarnings = FALSE)

# Run the enhanced FFT explanation
cat("Testing Enhanced FFT_explained function...\n")
cat("==========================================\n")

# Test 1: Enhanced version with gitter file
cat("\nTest 1: Enhanced analysis with gitter file\n")
try({
  result1 <- fft_explained(
    washed_path = washed_path,
    background_path = background_path,
    gitter_dat_path = gitter_dat_path,
    out_png = "results/figs/fft_explained_enhanced_test1.png",
    out_csv = "results/fft_explained_enhanced_metrics_test1.csv",
    colony_index_for_area = 2,
    blur_sigma = NA
  )
  
  cat("Test 1 completed successfully!\n")
  cat(sprintf("  - Detected %d colonies\n", length(result1$colony_mapping$colony_centers)))
  cat(sprintf("  - Intra frequency: %.4f cyc/px (period: %.1f px)\n", 
              result1$fft_results$f_intra, result1$fft_results$period_intra))
  cat(sprintf("  - Inter frequency: %.4f cyc/px (period: %.1f px)\n", 
              result1$fft_results$f_inter, result1$fft_results$period_inter))
  cat(sprintf("  - Per-colony areas calculated: %d\n", result1$per_colony_results$total_colonies))
})

# Test 2: Without gitter file (auto-detection fallback)
cat("\nTest 2: Without gitter file (enhanced auto-detection)\n")
try({
  result2 <- fft_explained(
    washed_path = washed_path,
    background_path = background_path,
    gitter_dat_path = NULL,
    out_png = "results/figs/fft_explained_enhanced_test2.png",
    out_csv = "results/fft_explained_enhanced_metrics_test2.csv",
    colony_index_for_area = 1,
    blur_sigma = NA
  )
  cat("Test 2 completed successfully!\n")
  cat(sprintf("  - Fallback y-row detection: %d\n", result2$y_row))
  cat(sprintf("  - Detected colonies: %d\n", length(result2$colony_mapping$colony_centers)))
})

# Test 3: With detrending and blur
cat("\nTest 3: Enhanced analysis with blur and detrending\n")
try({
  result3 <- fft_explained(
    washed_path = washed_path,
    background_path = background_path,
    gitter_dat_path = gitter_dat_path,
    out_png = "results/figs/fft_explained_enhanced_test3.png",
    out_csv = "results/fft_explained_enhanced_metrics_test3.csv",
    colony_index_for_area = 3,
    blur_sigma = 1.0
  )
  cat("Test 3 completed successfully!\n")
  cat(sprintf("  - Detrending method: %s\n", result3$fft_results$detrend_method))
  cat(sprintf("  - Blur applied: σ = %.1f\n", 1.0))
})

# Test 4: Different colony for area calculation
cat("\nTest 4: Testing different colony selection for area calculation\n")
try({
  result4 <- fft_explained(
    washed_path = washed_path,
    background_path = background_path,
    gitter_dat_path = gitter_dat_path,
    out_png = "results/figs/fft_explained_enhanced_test4.png",
    out_csv = "results/fft_explained_enhanced_metrics_test4.csv",
    colony_index_for_area = 5,  # Try a different colony
    blur_sigma = NA
  )
  cat("Test 4 completed successfully!\n")
  if (result4$area_results$colony_index <= length(result4$colony_mapping$colony_centers)) {
    cat(sprintf("  - Selected colony %d for area calculation\n", result4$area_results$colony_index))
    cat(sprintf("  - Colony area: %.1f\n", result4$area_results$area_under_one_period))
  }
})

# Validation tests
cat("\nValidation Summary:\n")
cat("==================\n")

# Check if per-colony CSV files were created
test_files <- list.files("results", pattern = "per_colony\\.csv$", full.names = TRUE)
cat(sprintf("Per-colony CSV files created: %d\n", length(test_files)))

# Check if figures were created
fig_files <- list.files("results/figs", pattern = "enhanced.*\\.png$", full.names = TRUE)
cat(sprintf("Enhanced figure files created: %d\n", length(fig_files)))

# Basic validation of one result
if (exists("result1") && !is.null(result1)) {
  cat("\nResult 1 validation:\n")
  cat(sprintf("  - Y-row used: %d\n", result1$y_row))
  cat(sprintf("  - Has dual frequencies: %s\n", 
              ifelse(!is.na(result1$fft_results$f_intra) && !is.na(result1$fft_results$f_inter), "YES", "NO")))
  cat(sprintf("  - Individual colonies detected: %d\n", length(result1$colony_mapping$colony_centers)))
  cat(sprintf("  - Detrending applied: %s\n", result1$fft_results$detrend_method))
  
  if (length(result1$per_colony_results$colony_areas) > 0) {
    areas <- result1$per_colony_results$colony_areas$area
    cat(sprintf("  - Colony area range: %.1f to %.1f\n", min(areas, na.rm = TRUE), max(areas, na.rm = TRUE)))
  }
}

cat("\nAll enhanced tests completed! Check the results/ directory for outputs.\n")
cat("Key improvements:\n")
cat("  - Dual-frequency detection (intra & inter tetramer)\n")
cat("  - Individual colony center identification\n")
cat("  - Robust detrending before FFT\n")
cat("  - Per-colony area calculation with local minima bounds\n")
cat("  - Enhanced validation and fallback mechanisms\n")