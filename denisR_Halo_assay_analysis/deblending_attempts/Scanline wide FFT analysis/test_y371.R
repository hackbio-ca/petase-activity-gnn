#!/usr/bin/env Rscript

# Test script for FFT_explained with forced y_row = 371
source("R/FFT_explained.R")

cat("Testing FFT_explained with forced y_row = 371\n")
cat("============================================\n")

# Run with forced y_row = 371
result <- fft_explained(
  washed_path = "spwsh/washed_BHET25_2d_6_1.JPG",
  background_path = "spbg/background_BHET25_6_1.JPG",
  gitter_dat_path = "spuwsh/BHET25_2d_6_1.JPG.dat",
  out_png = "results/figs/fft_explained_y371.png",
  out_csv = "results/fft_explained_y371_metrics.csv",
  colony_index_for_area = 2,
  blur_sigma = NA,
  y_row = 371
)

cat("\nTest completed successfully!\n")
cat(sprintf("Y-row used: %d\n", result$y_row))
cat(sprintf("Individual colonies detected: %d\n", length(result$colony_mapping$colony_centers)))
cat(sprintf("Intra frequency: %.4f cyc/px (period: %.1f px)\n",
            result$fft_results$f_intra, result$fft_results$period_intra))
cat(sprintf("Inter frequency: %.4f cyc/px (period: %.1f px)\n",
            result$fft_results$f_inter, result$fft_results$period_inter))

if (result$per_colony_results$total_colonies > 0) {
  cat(sprintf("Per-colony areas: %.1f to %.1f\n",
              min(result$per_colony_results$colony_areas$area, na.rm = TRUE),
              max(result$per_colony_results$colony_areas$area, na.rm = TRUE)))
}

cat("\nOutput files generated:\n")
cat(sprintf("  Figure: results/figs/fft_explained_y371.png\n"))
cat(sprintf("  Metrics: results/fft_explained_y371_metrics.csv\n"))
cat(sprintf("  Per-colony: results/fft_explained_y371_metrics_per_colony.csv\n"))
