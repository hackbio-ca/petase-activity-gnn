# Test script for GMM scanline analysis
# This script tests the implementation with BHET25 row 1 data

# Check if files exist
washed_path <- "Plate Images/washed/washed_BHET25_2d_6_1.JPG"
background_path <- "Plate Images/background/background_BHET25_6_1.JPG"
gitter_dat_path <- "Plate Images/unwashed/BHET25_2d_6_1.JPG.dat"

cat("=== Checking file existence ===\n")
cat("Washed image:", ifelse(file.exists(washed_path), "EXISTS", "MISSING"), "\n")
cat("Background image:", ifelse(file.exists(background_path), "EXISTS", "MISSING"), "\n")
cat("Gitter data:", ifelse(file.exists(gitter_dat_path), "EXISTS", "MISSING"), "\n")

if (!all(file.exists(c(washed_path, background_path, gitter_dat_path)))) {
  cat("ERROR: Some required files are missing. Please check file paths.\n")
  stop("Missing files")
}

# Load the functions
cat("\n=== Loading functions ===\n")
tryCatch({
  source("R/GMM_scanline.R")
  source("R/GMM_explained.R")
  cat("Functions loaded successfully\n")
}, error = function(e) {
  cat("Error loading functions:", conditionMessage(e), "\n")
  stop("Function loading failed")
})

# Test with the specified parameters for BHET25 row 1
cat("\n=== Running GMM scanline analysis ===\n")

tryCatch({
  # Run the analysis
  result <- gmm_scanline(
    washed_path = washed_path,
    background_path = background_path,
    gitter_dat_path = gitter_dat_path,
    out_dir = "results/scanline",
    row_index = 1,
    margin = 25,
    snap_radius = 30,
    r_intra_min = 110,
    r_intra_max = 360,
    detrend_k = NULL,     # auto from spacing
    area_min = 2.0,
    snr_min = 0.6,
    resp_min = 0.25,
    save_debug = TRUE
  )
  
  cat("\n=== Analysis completed successfully! ===\n")
  cat("Found", sum(result$colony_data$present), "present colonies out of", 
      nrow(result$colony_data), "total\n")
  
  # Display some results
  present_colonies <- result$colony_data[result$colony_data$present, ]
  if (nrow(present_colonies) > 0) {
    cat("\nPresent colonies:\n")
    print(present_colonies[, c("colony_id", "area_weighted", "snr", "resp_mean", "center_x")])
  }
  
  # Run visualization if analysis succeeded
  cat("\n=== Creating visualization plots ===\n")
  
  viz_result <- gmm_explained(
    washed_path = washed_path,
    background_path = background_path,
    scanline_csv = "results/scanline/gmm_scanline_per_colony.csv",
    overview_csv = "results/scanline/gmm_scanline_overview.csv",
    out_dir = "results/figs/colony",
    show_detrended = TRUE
  )
  
  cat("Visualization completed successfully!\n")
  cat("Plots saved to results/figs/colony/\n")
  
}, error = function(e) {
  cat("\n=== ERROR OCCURRED ===\n")
  cat("Error message:", conditionMessage(e), "\n")
  cat("\nFull error details:\n")
  print(e)
  cat("\nThis error suggests there may be an issue with:\n")
  cat("1. Image file formats or corruption\n")
  cat("2. Gitter data format\n") 
  cat("3. Required R packages not installed\n")
  cat("4. Memory limitations\n")
})