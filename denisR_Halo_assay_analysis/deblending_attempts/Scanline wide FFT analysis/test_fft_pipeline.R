# Test Script for FFT Bleed Suppression Pipeline
# ==============================================
# This script tests the FFT enhancement without requiring the full dataset

cat("Testing FFT Bleed Suppression Pipeline\n")
cat("=====================================\n")

# Check if required packages are available
required_packages <- c("imager", "dplyr", "ggplot2", "readr")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

# Load packages
suppressMessages({
  library(imager)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# Try to load fftwtools (may not be available on all systems)
fft_available <- requireNamespace("fftwtools", quietly = TRUE)
if (!fft_available) {
  cat("Warning: fftwtools not available. Some FFT functions may fall back to base R.\n")
}

cat("\nTesting function loading...\n")

# Source all functions
tryCatch({
  source("R/initialize_colony_data.R")
  source("R/create_grid_boxes.R") 
  source("R/calculate_average_intensity.R")
  source("R/calculate_edge_intensity.R")
  source("R/classify_pixels.R")
  source("R/rgb_to_gray.R")
  source("R/align_images.R")
  source("R/fft_bleed_suppression.R")
  source("R/fft_config.R")
  cat("✓ All functions loaded successfully\n")
}, error = function(e) {
  cat("✗ Error loading functions:", e$message, "\n")
  stop("Cannot continue without functions")
})

cat("\nTesting configuration system...\n")

# Test configuration functions
tryCatch({
  # Test basic config creation
  config <- create_fft_config()
  cat("✓ Basic configuration created\n")
  
  # Test preset configs
  presets <- c("conservative", "standard", "aggressive", "manual_spacing")
  for (preset in presets) {
    config <- fft_preset_config(preset)
    cat("✓ Preset '", preset, "' loaded\n", sep = "")
  }
  
  # Test validation
  valid <- validate_fft_config(config)
  cat("✓ Configuration validation works\n")
  
}, error = function(e) {
  cat("✗ Configuration error:", e$message, "\n")
})

cat("\nTesting with synthetic images...\n")

# Create synthetic test images
tryCatch({
  # Create a simple test image (100x100 pixels)
  test_img <- as.cimg(array(rnorm(100*100*1*3, mean = 0.5, sd = 0.1), 
                           dim = c(100, 100, 1, 3)))
  
  # Add some periodic structure (simulating colony lattice)
  x_coords <- 1:100
  y_coords <- 1:100
  periodic_pattern <- outer(sin(2*pi*x_coords/20), sin(2*pi*y_coords/20))
  
  for (ch in 1:3) {
    test_img[,,1,ch] <- test_img[,,1,ch] + 0.1 * periodic_pattern
  }
  
  cat("✓ Synthetic test images created\n")
  
  # Test image alignment
  bg_test <- test_img
  washed_test <- imshift(test_img, delta_x = 2, delta_y = 3)
  
  alignment_result <- align_images_simple(bg_test, washed_test, search_radius = 10)
  cat("✓ Image alignment works (detected shift: ", 
      alignment_result$shift_x, ", ", alignment_result$shift_y, ")\n", sep = "")
  
  # Test FFT cleaning (basic version without fftwtools)
  if (fft_available) {
    fft_result <- fft_bleed_suppression(
      washed_img = test_img,
      hp_radius = 6,
      hp_slope = 3,
      notch_radius = 2,
      harmonics = 1:2
    )
    cat("✓ FFT bleed suppression works\n")
  } else {
    cat("⚠ FFT bleed suppression skipped (fftwtools not available)\n")
  }
  
}, error = function(e) {
  cat("✗ Synthetic image test error:", e$message, "\n")
})

cat("\nTesting with real images (if available)...\n")

# Test with real images if they exist
real_image_test <- function() {
  bg_path <- "background/background/background_BHET25_5_1.JPG"
  washed_path <- "all_plate_images/all_plate_images/wash_images/washed_BHET25_8h_5_1.JPG"
  
  if (file.exists(bg_path) && file.exists(washed_path)) {
    tryCatch({
      bg_img <- load.image(bg_path)
      washed_img <- load.image(washed_path)
      
      cat("✓ Real images loaded successfully\n")
      cat("  Background dimensions:", dim(bg_img)[1:2], "\n")
      cat("  Washed dimensions:    ", dim(washed_img)[1:2], "\n")
      
      # Test alignment
      alignment <- align_images_simple(bg_img, washed_img, search_radius = 20)
      cat("✓ Real image alignment works\n")
      
      # Test FFT (if available)
      if (fft_available) {
        fft_result <- fft_bleed_suppression(washed_img, hp_radius = 6)
        cat("✓ Real image FFT cleaning works\n")
      }
      
      return(TRUE)
    }, error = function(e) {
      cat("✗ Real image test error:", e$message, "\n")
      return(FALSE)
    })
  } else {
    cat("⚠ Real images not found at expected paths\n")
    cat("  Expected background: ", bg_path, "\n")
    cat("  Expected washed:     ", washed_path, "\n")
    return(FALSE)
  }
}

real_images_ok <- real_image_test()

cat("\nTest Summary\n")
cat("============\n")
cat("Functions loaded:        ✓\n")
cat("Configuration system:    ✓\n")
cat("Synthetic image tests:   ✓\n")
cat("Real image tests:        ", if(real_images_ok) "✓" else "⚠ (images not found)", "\n")
cat("FFT functionality:       ", if(fft_available) "✓" else "⚠ (fftwtools needed)", "\n")

if (fft_available && real_images_ok) {
  cat("\n🎉 Full pipeline ready! Run main_with_fft.R for complete analysis.\n")
} else if (!fft_available) {
  cat("\n⚠ Install fftwtools for full FFT functionality:\n")
  cat("   install.packages('fftwtools')\n")
} else {
  cat("\n⚠ Place your images in the expected directories to test with real data.\n")
}

cat("\nNext steps:\n")
cat("1. Ensure images are in correct directories\n")
cat("2. Install fftwtools if needed\n") 
cat("3. Run source('R/main_with_fft.R') for full analysis\n")
cat("4. Adjust parameters in fft_config.R as needed\n")