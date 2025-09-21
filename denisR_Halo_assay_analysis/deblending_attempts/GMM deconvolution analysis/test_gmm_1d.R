# Test script for 1D Weighted GMM approach
# This implements the new algorithm: GMM over x with intensity weights

# Load the functions  
source("R/GMM_scanline.R")  # For gitter reading function
source("R/GMM_1D_weighted.R")

# Set up paths
washed_path <- "Plate Images/washed/washed_BHET25_2d_6_1.JPG"
background_path <- "Plate Images/background/background_BHET25_6_1.JPG"
gitter_dat_path <- "Plate Images/unwashed/BHET25_2d_6_1.JPG.dat"

cat("=== Testing 1D Weighted GMM Approach ===\n")

# Check files exist
if (!all(file.exists(c(washed_path, background_path, gitter_dat_path)))) {
  stop("Required files not found")
}

tryCatch({
  # Run the 1D weighted GMM analysis
  result <- gmm_1d_weighted(
    washed_path = washed_path,
    background_path = background_path,
    gitter_dat_path = gitter_dat_path,
    out_dir = "results/gmm_1d",
    row_index = 1,
    baseline_window = 500,      # Wide baseline window
    baseline_quantile = 0.15,   # 15th percentile baseline
    margin = 25,
    min_mass = 0.01,           # Minimum component mass
    min_snr = 2.0,             # Minimum SNR
    save_plots = TRUE
  )
  
  cat("\n=== Analysis Complete! ===\n")
  cat("Found", nrow(result$results_df), "colonies\n")
  cat("Total weighted area:", sum(result$colony_areas), "\n")
  
  # Show results
  print(result$results_df[, c("colony_id", "center_x", "area_weighted")])
  
  cat("\nFiles saved to:", file.path(getwd(), "results/gmm_1d"), "\n")
  cat("- gmm_1d_colony_areas.csv\n")
  cat("- plots/overview.png\n")
  cat("- plots/colony_XX.png (one per colony)\n")
  
}, error = function(e) {
  cat("Error:", conditionMessage(e), "\n")
  traceback()
})