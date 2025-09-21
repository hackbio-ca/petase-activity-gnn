# Test script to verify background image matching
unwashed_files <- c("BHET25_2d_6_1.JPG", "BHET25_2d_6_2.JPG", "BHET25_2d_6_3.JPG", 
                    "BHET25_2d_7_4.JPG", "BHET25_2d_7_5.JPG", "BHET25_2d_7_6.JPG")

background_files <- c("background_BHET25_6_1.JPG", "background_BHET25_6_2.JPG", 
                     "background_BHET25_6_3.JPG", "background_BHET25_7_4.JPG", 
                     "background_BHET25_7_5.JPG", "background_BHET25_7_6.JPG")

cat("Testing background image matching:\n")
for (img_file in unwashed_files) {
  # Apply the same logic as in the main script
  bg_basename <- paste0("background_", gsub("_2d_", "_", img_file))
  
  if (bg_basename %in% background_files) {
    cat("✓", img_file, "->", bg_basename, "(FOUND)\n")
  } else {
    cat("✗", img_file, "->", bg_basename, "(NOT FOUND)\n")
  }
}

cat("\nBackground files available:\n")
for (bg_file in background_files) {
  cat(" -", bg_file, "\n")
}