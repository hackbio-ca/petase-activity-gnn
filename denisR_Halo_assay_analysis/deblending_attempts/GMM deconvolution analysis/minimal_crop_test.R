# Minimal test to isolate the cropping issue
# This will help identify where exactly the "vector too long" error occurs

cat("=== Minimal Crop Test ===\n")

# Check if required packages are available
required_packages <- c("imager", "readr")
missing_packages <- character(0)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  cat("Missing packages:", paste(missing_packages, collapse = ", "), "\n")
  cat("Please install with: install.packages(c(", 
      paste(paste0('"', missing_packages, '"'), collapse = ", "), "))\n")
  stop("Required packages missing")
}

# Load just the helper functions we need
source("R/GMM_scanline.R")

# Set up paths
washed_path <- "Plate Images/washed/washed_BHET25_2d_6_1.JPG"
gitter_dat_path <- "Plate Images/unwashed/BHET25_2d_6_1.JPG.dat"

cat("\n=== Step 1: Check file existence ===\n")
cat("Washed image exists:", file.exists(washed_path), "\n")
cat("Gitter data exists:", file.exists(gitter_dat_path), "\n")

if (!file.exists(washed_path) || !file.exists(gitter_dat_path)) {
  stop("Required files not found")
}

cat("\n=== Step 2: Read gitter data ===\n")
tryCatch({
  gitter_df <- read_gitter_dat(gitter_dat_path)
  cat("Gitter data loaded:", nrow(gitter_df), "rows\n")
  cat("Coordinate ranges:\n")
  cat("  xl:", min(gitter_df$xl), "to", max(gitter_df$xl), "\n")
  cat("  xr:", min(gitter_df$xr), "to", max(gitter_df$xr), "\n")
  cat("  yt:", min(gitter_df$yt), "to", max(gitter_df$yt), "\n")
  cat("  yb:", min(gitter_df$yb), "to", max(gitter_df$yb), "\n")
}, error = function(e) {
  cat("Error reading gitter data:", conditionMessage(e), "\n")
  stop("Gitter reading failed")
})

cat("\n=== Step 3: Load image ===\n")
tryCatch({
  washed_img <- imager::load.image(washed_path)
  cat("Image loaded successfully\n")
  cat("Image dimensions:", dim(washed_img), "\n")
  cat("Width x Height:", dim(washed_img)[1], "x", dim(washed_img)[2], "\n")
}, error = function(e) {
  cat("Error loading image:", conditionMessage(e), "\n")
  stop("Image loading failed")
})

cat("\n=== Step 4: Convert to grayscale ===\n")
tryCatch({
  if (dim(washed_img)[4] > 1) {
    washed_gray <- 0.299 * washed_img[,,1,1] + 0.587 * washed_img[,,1,2] + 0.114 * washed_img[,,1,3]
    cat("Converted RGB to grayscale\n")
  } else {
    washed_gray <- washed_img[,,1,1]
    cat("Image already grayscale\n")
  }
  cat("Grayscale dimensions:", dim(washed_gray), "\n")
}, error = function(e) {
  cat("Error converting to grayscale:", conditionMessage(e), "\n")
  stop("Grayscale conversion failed")
})

cat("\n=== Step 5: Compute crop boundaries ===\n")
tryCatch({
  img_dim <- c(width = dim(washed_img)[1], height = dim(washed_img)[2])
  cat("Image dimension array:", img_dim, "\n")
  
  # Test with small margin first
  margin <- 25
  crop_info <- compute_global_crop(gitter_df, img_dim, margin)
  cat("Crop computation successful\n")
}, error = function(e) {
  cat("Error computing crop:", conditionMessage(e), "\n")
  stop("Crop computation failed")
})

cat("\n=== Step 6: Test range creation ===\n")
tryCatch({
  left_range <- crop_info$left:crop_info$right
  top_range <- crop_info$top:crop_info$bottom
  
  cat("Left range length:", length(left_range), "\n")
  cat("Top range length:", length(top_range), "\n")
  cat("Total crop pixels:", length(left_range) * length(top_range), "\n")
  
  if (length(left_range) > 1000000 || length(top_range) > 1000000) {
    stop("Range too large - this would cause the 'vector too long' error")
  }
  
  cat("Range creation successful\n")
}, error = function(e) {
  cat("Error creating ranges:", conditionMessage(e), "\n")
  cat("This is likely the source of the 'vector too long' error\n")
  stop("Range creation failed")
})

cat("\n=== Step 7: Test actual cropping ===\n")
tryCatch({
  # Just test a small subset first
  test_left <- max(1, crop_info$left)
  test_right <- min(test_left + 100, crop_info$right, dim(washed_gray)[1])
  test_top <- max(1, crop_info$top)
  test_bottom <- min(test_top + 100, crop_info$bottom, dim(washed_gray)[2])
  
  cat("Testing crop with bounds: [", test_left, ":", test_right, "] x [", test_top, ":", test_bottom, "]\n")
  test_crop <- washed_gray[test_left:test_right, test_top:test_bottom]
  cat("Test crop successful, size:", dim(test_crop), "\n")
  
  cat("\n=== All tests passed! ===\n")
  cat("The cropping should work with the corrected bounds checking.\n")
  
}, error = function(e) {
  cat("Error in test cropping:", conditionMessage(e), "\n")
  stop("Test cropping failed")
})