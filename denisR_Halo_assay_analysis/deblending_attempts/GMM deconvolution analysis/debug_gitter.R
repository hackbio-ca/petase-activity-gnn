# Debug script to check gitter data and image dimensions
source("R/GMM_scanline.R")

# Set up paths
gitter_dat_path <- "Plate Images/unwashed/BHET25_2d_6_1.JPG.dat"
washed_path <- "Plate Images/washed/washed_BHET25_2d_6_1.JPG"

cat("=== Debug: Reading gitter data ===\n")
gitter_df <- read_gitter_dat(gitter_dat_path)
cat("Gitter data shape:", nrow(gitter_df), "rows,", ncol(gitter_df), "columns\n")
cat("Column names:", paste(names(gitter_df), collapse = ", "), "\n")

if (nrow(gitter_df) > 0) {
  cat("First few rows:\n")
  print(head(gitter_df))
  
  cat("\nGitter coordinate ranges:\n")
  cat("xl: min =", min(gitter_df$xl, na.rm = TRUE), "max =", max(gitter_df$xl, na.rm = TRUE), "\n")
  cat("xr: min =", min(gitter_df$xr, na.rm = TRUE), "max =", max(gitter_df$xr, na.rm = TRUE), "\n")
  cat("yt: min =", min(gitter_df$yt, na.rm = TRUE), "max =", max(gitter_df$yt, na.rm = TRUE), "\n")
  cat("yb: min =", min(gitter_df$yb, na.rm = TRUE), "max =", max(gitter_df$yb, na.rm = TRUE), "\n")
  
  # Check for any extreme or invalid values
  cat("\nChecking for problematic values:\n")
  cat("Any NA in xl:", any(is.na(gitter_df$xl)), "\n")
  cat("Any NA in xr:", any(is.na(gitter_df$xr)), "\n") 
  cat("Any NA in yt:", any(is.na(gitter_df$yt)), "\n")
  cat("Any NA in yb:", any(is.na(gitter_df$yb)), "\n")
  
  cat("Any negative xl:", any(gitter_df$xl < 0, na.rm = TRUE), "\n")
  cat("Any negative yt:", any(gitter_df$yt < 0, na.rm = TRUE), "\n")
}

cat("\n=== Debug: Loading image ===\n")
if (file.exists(washed_path)) {
  tryCatch({
    washed_img <- imager::load.image(washed_path)
    cat("Image loaded successfully\n")
    cat("Image dimensions:", dim(washed_img), "\n")
    cat("Width x Height:", dim(washed_img)[1], "x", dim(washed_img)[2], "\n")
  }, error = function(e) {
    cat("Error loading image:", conditionMessage(e), "\n")
  })
} else {
  cat("Image file not found:", washed_path, "\n")
}

cat("\n=== Computing crop with margin = 25 ===\n")
if (exists("washed_img") && nrow(gitter_df) > 0) {
  img_dim <- c(width = dim(washed_img)[1], height = dim(washed_img)[2])
  cat("Image dimensions array:", img_dim, "\n")
  
  tryCatch({
    crop_info <- compute_global_crop(gitter_df, img_dim, margin = 25)
    cat("Crop computed successfully\n")
  }, error = function(e) {
    cat("Error computing crop:", conditionMessage(e), "\n")
  })
}