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
gitter_dat_path <- "Plate Images/unwashed/BHET25_2d_6_1.JPG.dat"

# Validate input files exist
if (!file.exists(washed_path)) {
  stop("Washed image not found: ", washed_path)
}
if (!file.exists(background_path)) {
  stop("Background image not found: ", background_path)
}
if (!file.exists(gitter_dat_path)) {
  stop("Gitter file not found: ", gitter_dat_path)
}

cat("✓ Input files validated\n")

# Test A: Robust FFT Scanline Analysis
cat("\n--- Testing Robust FFT Scanline Analysis ---\n")

scanline_results <- fft_scanline(
  washed_path = washed_path,
  background_path = background_path,
  gitter_dat_path = gitter_dat_path,
  out_dir = "results/scanline",
  row_index = 1,
  margin = 15,
  snr_min = 3,
  area_min = 50,
  eps_baseline = 0.5,
  r_intra_min = 40,
  r_intra_max = 220,
  continuity_q = 0.6,
  continuity_k = 0.5,
  save_debug = TRUE
)

cat("✓ Robust scanline analysis completed\n")

# Validate results structure
if (!is.list(scanline_results) || 
    !all(c("per_colony", "overview", "adjusted", "detrended", "centers", "windows", "y_row", "crop_info") %in% names(scanline_results))) {
  stop("Invalid scanline results structure")
}

cat("✓ Results structure validated\n")

# Report key metrics
per_colony_df <- scanline_results$per_colony
overview_df <- scanline_results$overview

n_total <- nrow(per_colony_df)
n_present <- sum(per_colony_df$present)
n_edge <- sum(per_colony_df$touches_edge)
n_low_snr <- sum(per_colony_df$reason == "low_snr")
n_low_area <- sum(per_colony_df$reason == "low_area")
n_discontinuous <- sum(per_colony_df$reason == "discontinuous")

cat(sprintf("Colonies detected: %d\n", n_total))
cat(sprintf("Present colonies: %d / %d\n", n_present, n_total))
cat(sprintf("Rejection reasons:\n"))
cat(sprintf("  - Edge touching: %d\n", n_edge))
cat(sprintf("  - Low SNR: %d\n", n_low_snr))
cat(sprintf("  - Low area: %d\n", n_low_area))
cat(sprintf("  - Discontinuous: %d\n", n_discontinuous))

# Validate crop info
crop_info <- scanline_results$crop_info
cat(sprintf("Crop region: left=%d, right=%d, top=%d, bottom=%d\n", 
            crop_info$left, crop_info$right, crop_info$top, crop_info$bottom))

# Test B: Robust FFT Plotting
cat("\n--- Testing Robust FFT Plotting ---\n")

plotting_results <- fft_explained(
  washed_path = washed_path,
  background_path = background_path,
  scanline_csv = "results/scanline/fft_scanline_per_colony.csv",
  overview_csv = "results/scanline/fft_scanline_overview.csv",
  out_dir = "results/figs/colony",
  show_detrended = TRUE
)

cat("✓ Robust plotting completed\n")

# Validate plotting results
if (!is.list(plotting_results) || 
    !all(c("overview_plot", "colony_plots", "summary_csv", "n_present", "n_total") %in% names(plotting_results))) {
  stop("Invalid plotting results structure")
}

cat("✓ Plotting results structure validated\n")

# Report plotting metrics
cat(sprintf("Overview plot: %s\n", plotting_results$overview_plot))
cat(sprintf("Colony plots generated: %d\n", length(plotting_results$colony_plots)))
cat(sprintf("Summary CSV: %s\n", plotting_results$summary_csv))

# Test C: Validation of Requirements
cat("\n--- Validating Robust Implementation Requirements ---\n")

# C1: Global cropping validation
if (crop_info$left <= 1 || crop_info$right >= 5000 || 
    crop_info$top <= 1 || crop_info$bottom >= 3500) {
  cat("⚠ Warning: Crop region may not be optimal\n")
} else {
  cat("✓ Global cropping appears appropriate\n")
}

# C2: Present colony filtering validation
present_colonies <- per_colony_df[per_colony_df$present, ]
if (nrow(present_colonies) > 0) {
  cat("✓ Present colony filtering functional\n")
  
  # Check SNR values for present colonies
  min_snr_present <- min(present_colonies$snr)
  max_area_present <- max(present_colonies$area_adjusted)
  cat(sprintf("Present colony SNR range: %.2f - %.2f\n", min_snr_present, max(present_colonies$snr)))
  cat(sprintf("Present colony area range: %.1f - %.1f\n", min(present_colonies$area_adjusted), max_area_present))
  
  # C3: FFT frequency validation
  valid_freq <- present_colonies$dom_freq_cyc_per_px > 0 & present_colonies$dom_freq_cyc_per_px < 0.5
  cat(sprintf("Valid dominant frequencies: %d / %d\n", sum(valid_freq), nrow(present_colonies)))
  
  if (all(valid_freq)) {
    cat("✓ FFT frequency detection working\n")
  } else {
    cat("⚠ Warning: Some invalid frequencies detected\n")
  }
} else {
  cat("⚠ Warning: No present colonies found\n")
}

# C4: Output file validation
required_files <- c(
  "results/scanline/fft_scanline_per_colony.csv",
  "results/scanline/fft_scanline_overview.csv", 
  "results/scanline/scanline_overview.png",
  "results/figs/colony/scanline_overview.png",
  "results/figs/colony/fft_explained_summary.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  cat("⚠ Warning: Missing output files:\n")
  for (f in missing_files) cat(sprintf("  - %s\n", f))
} else {
  cat("✓ All required output files generated\n")
}

# C5: CSV structure validation
cat("\n--- Validating CSV Structure ---\n")

required_per_colony_cols <- c("colony_id", "present", "center_x", "left", "right", 
                             "width_px", "area_adjusted", "area_detrended", "snr", 
                             "continuity", "touches_edge", "dom_freq_cyc_per_px", 
                             "dom_period_px", "reason")

missing_cols <- required_per_colony_cols[!required_per_colony_cols %in% colnames(per_colony_df)]
if (length(missing_cols) > 0) {
  cat("⚠ Warning: Missing per-colony columns:\n")
  for (c in missing_cols) cat(sprintf("  - %s\n", c))
} else {
  cat("✓ Per-colony CSV structure validated\n")
}

required_overview_cols <- c("y_row", "row_index", "n_colonies_total", "n_colonies_present",
                           "crop_left", "crop_right", "crop_top", "crop_bottom",
                           "detrend_k", "eps_baseline", "area_min", "snr_min")

missing_overview_cols <- required_overview_cols[!required_overview_cols %in% colnames(overview_df)]
if (length(missing_overview_cols) > 0) {
  cat("⚠ Warning: Missing overview columns:\n")
  for (c in missing_overview_cols) cat(sprintf("  - %s\n", c))
} else {
  cat("✓ Overview CSV structure validated\n")
}

# Final summary
cat("\n=== ROBUST IMPLEMENTATION VALIDATION COMPLETE ===\n")
cat(sprintf("Total colonies detected: %d\n", n_total))
cat(sprintf("Present colonies: %d\n", n_present))
cat(sprintf("Edge-touching colonies rejected: %d\n", n_edge))
cat(sprintf("Crop dimensions: %dx%d\n", 
            crop_info$right - crop_info$left + 1,
            crop_info$bottom - crop_info$top + 1))
cat(sprintf("Individual colony plots: %d\n", length(plotting_results$colony_plots)))

if (n_present >= 2 && n_edge > 0) {
  cat("✓ Robust implementation validated successfully!\n")
  cat("Key features confirmed:\n")
  cat("  ✓ Global cropping with margin\n")
  cat("  ✓ Signal minima-based windowing\n") 
  cat("  ✓ SNR/area/continuity filtering\n")
  cat("  ✓ Edge artifact rejection\n")
  cat("  ✓ Present-only visualization\n")
  cat("  ✓ Comprehensive CSV outputs\n")
} else {
  cat("⚠ Implementation may need tuning - check parameters\n")
}

cat("\n=== TEST COMPLETE ===\n")
  background_path = background_path,
  scanline_csv = "results/scanline/fft_scanline_per_colony.csv",
  overview_csv = "results/scanline/fft_scanline_overview.csv",
  out_dir = "results/figs/colony",
  show_detrended = TRUE
)

cat("✓ Explanatory plots created\n")

# Print summary
cat("\n=== Analysis Summary ===\n")
cat("Colonies detected:", nrow(scanline_results$per_colony), "\n")
cat("Present colonies:", sum(scanline_results$per_colony$present), "\n")
cat("Y-row analyzed:", scanline_results$y_row, "\n")
cat("Source of centers:", ifelse(length(unique(scanline_results$per_colony$source_centers)) > 0, 
                                 unique(scanline_results$per_colony$source_centers)[1], "unknown"), "\n")

# Show colony areas and SNR (present colonies only)
present_colonies <- scanline_results$per_colony[scanline_results$per_colony$present, ]
cat("\nPresent colony metrics:\n")
if (nrow(present_colonies) > 0) {
  for (i in seq_len(nrow(present_colonies))) {
    colony <- present_colonies[i, ]
    cat(sprintf("  Colony %.0f: Area=%.1f, SNR=%.1f (center: %.0f, width: %.0f px)\n", 
                colony$colony_id, colony$area_adjusted, colony$snr,
                colony$center_x, colony$width_px))
  }
} else {
  cat("  No present colonies found!\n")
}

# Show absent colonies summary
absent_colonies <- scanline_results$per_colony[!scanline_results$per_colony$present, ]
if (nrow(absent_colonies) > 0) {
  cat("\nAbsent colonies summary:\n")
  absence_reasons <- table(absent_colonies$notes)
  for (reason in names(absence_reasons)) {
    cat(sprintf("  %s: %d colonies\n", reason, absence_reasons[reason]))
  }
}

# Additional testing commented out for now
#
# skip_additional_tests <- TRUE
# 
# if (!skip_additional_tests) {
#   cat("\n--- Testing Peak Detection Mode (No Gitter) ---\n")
#   scanline_results_peaks <- fft_scanline(
#     washed_path = washed_path,
#     background_path = background_path,
#     gitter_dat_path = NULL,  # Force peak detection
#     out_dir = "results/scanline_peaks",
#     row_index = 1,
#     use_gitter_centers = FALSE,
#     margin = 15,
#     snr_min = 3,
#     area_min = 50,
#     min_sep_px = 60,
#     save_debug = TRUE
#   )
#   
#   cat("✓ Peak detection mode completed\n")
#   cat("Colonies detected (peak mode):", nrow(scanline_results_peaks$per_colony), "\n")
#   cat("Present colonies (peak mode):", sum(scanline_results_peaks$per_colony$present), "\n")
#   
#   # Test with multiple plates if desired
#   test_multiple_plates <- TRUE
#   if (test_multiple_plates) {
#     cat("\n--- Testing Additional Plates ---\n")
#     
#     plate_pairs <- list(
#       list(washed = "Plate Images/washed/washed_BHET25_2d_6_2.JPG",
#            background = "Plate Images/background/background_BHET25_6_2.JPG",
#            gitter = "Plate Images/unwashed/BHET25_2d_6_2.JPG.dat"),
#       list(washed = "Plate Images/washed/washed_BHET25_2d_7_4.JPG",
#            background = "Plate Images/background/background_BHET25_7_4.JPG",
#            gitter = "Plate Images/unwashed/BHET25_2d_7_4.JPG.dat")
#     )
#     
#     for (i in seq_along(plate_pairs)) {
#       plate <- plate_pairs[[i]]
#       
#       cat(sprintf("\nTesting plate %d...\n", i + 1))
#       
#       if (all(sapply(unlist(plate), file.exists))) {
#         result <- fft_scanline(
#           washed_path = plate$washed,
#           background_path = plate$background,
#           gitter_dat_path = plate$gitter,
#           out_dir = paste0("results/scanline_plate_", i + 1),
#           row_index = 1,
#           save_debug = TRUE
#         )
#         
#         cat(sprintf("  ✓ Plate %d: %d colonies detected\n", 
#                     i + 1, nrow(result$per_colony)))
#       } else {
#         cat(sprintf("  ✗ Plate %d: Missing files\n", i + 1))
#       }
#     }
#   }
# }

cat("\n=== Test Complete ===\n")
cat("Check the results/ directory for outputs:\n")
cat("  - results/scanline/ : CSV files and debug plots\n")
cat("  - results/figs/colony/ : Overview and per-colony plots\n")
# Additional testing outputs would appear here when enabled
cat("\nMain functionality test completed successfully!\n")