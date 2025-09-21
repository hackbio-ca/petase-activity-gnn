# Quick Test of Fixed Alignment Function
# =====================================
# This script tests that the alignment error is fixed

cat("Testing fixed alignment function...\n")

# Load required libraries
suppressMessages({
  library(imager)
})

# Source the fixed alignment function
source("R/align_images.R")

cat("Testing with provided image files...\n")

# Test with the specific files mentioned
tryCatch({
  # Load images
  bg_img <- imager::load.image("background/background/background_BHET25_5_1.JPG")
  washed_img <- imager::load.image("all_plate_images/all_plate_images/wash_images/washed_BHET25_8h_5_1.JPG")
  
  cat("✓ Images loaded successfully\n")
  cat("  Background dimensions:", dim(bg_img)[1:2], "\n")
  cat("  Washed dimensions:    ", dim(washed_img)[1:2], "\n")
  
  # Test phase correlation alignment (the one that was failing)
  cat("Testing phase correlation alignment...\n")
  result <- align_images_phase_correlation(bg_img, washed_img)
  
  cat("✓ Phase correlation alignment successful!\n")
  cat("  Detected shift: (", result$shift_x, ",", result$shift_y, ")\n")
  cat("  Aligned image dimensions:", dim(result$aligned_img)[1:2], "\n")
  
  # Test simple alignment as backup
  cat("Testing simple alignment...\n")
  result2 <- align_images_simple(bg_img, washed_img, search_radius = 20)
  
  cat("✓ Simple alignment successful!\n")
  cat("  Detected shift: (", result2$shift_x, ",", result2$shift_y, ")\n")
  cat("  Correlation:", round(result2$correlation, 3), "\n")
  
  cat("\n🎉 Alignment functions are working correctly!\n")
  cat("The type error has been fixed. You can now run the full pipeline.\n")
  
}, error = function(e) {
  cat("✗ Error:", e$message, "\n")
  cat("Make sure the image files exist at the specified paths.\n")
})

cat("\nNext steps:\n")
cat("1. Run source('R/colony_summary.R') for complete analysis\n")
cat("2. Or run source('examples_fft.R') for step-by-step examples\n")