#!/usr/bin/env Rscript

# Test script for the simplified per-colony FFT analysis
source("R/FFT_explained.R")

cat("Testing Simple Per-Colony FFT Analysis\n")
cat("=====================================\n")

# Test with y_row = 341
result <- fft_explained(
  washed_path = "spwsh/washed_BHET25_2d_6_1.JPG",
  background_path = "spbg/background_BHET25_6_1.JPG", 
  gitter_dat_path = "spuwsh/BHET25_2d_6_1.JPG.dat",
  out_csv = "results/fft_explained_per_colony.csv",
  out_overview_png = "results/figs/fft_overview.png",
  colony_index_for_area = 2,
  blur_sigma = NA,
  y_row = 341  # Force specific y-row
)

cat("\nTest completed!\n")
cat(sprintf("Y-row used: %d\n", result$y_row))
cat(sprintf("Colonies detected: %d\n", length(result$colony_centers)))
if (nrow(result$colony_data) > 0) {
  cat(sprintf("Colony areas: %.1f to %.1f\n", 
              min(result$colony_data$area, na.rm = TRUE),
              max(result$colony_data$area, na.rm = TRUE)))
  
  cat("\nFirst few colonies:\n")
  print(head(result$colony_data, 3))
}

cat("\nFiles generated:\n")
cat("  Per-colony CSV: results/fft_explained_per_colony.csv\n")
cat("  Overview plot: results/figs/fft_overview.png\n")
cat("  Individual plots: results/figs/colony/C*.png\n")