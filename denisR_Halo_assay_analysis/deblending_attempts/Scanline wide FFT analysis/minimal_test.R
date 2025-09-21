# Minimal test for FFT_explained.R functionality
# This checks if the basic functions can load and parse files

# Source the script
cat("Testing FFT_explained.R basic functionality...\n")

# Check if required packages are available
required_packages <- c("imager", "ggplot2", "gridExtra", "readr", "dplyr")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Missing required packages:", paste(missing_packages, collapse = ", "), "\n")
  cat("Please install them using: install.packages(c(", paste(paste0('"', missing_packages, '"'), collapse = ", "), "))\n")
} else {
  cat("All required packages are available.\n")
}

# Test file existence
washed_path <- "spwsh/washed_BHET25_2d_6_1.JPG"
background_path <- "spbg/background_BHET25_6_1.JPG"
gitter_dat_path <- "spuwsh/BHET25_2d_6_1.JPG.dat"

cat("\nChecking file availability:\n")
cat("Washed image:", ifelse(file.exists(washed_path), "EXISTS", "NOT FOUND"), "\n")
cat("Background image:", ifelse(file.exists(background_path), "EXISTS", "NOT FOUND"), "\n")
cat("Gitter file:", ifelse(file.exists(gitter_dat_path), "EXISTS", "NOT FOUND"), "\n")

# Test gitter file parsing
if (file.exists(gitter_dat_path)) {
  cat("\nTesting gitter file parsing...\n")
  
  # Simple parsing test
  con <- file(gitter_dat_path, "r")
  lines_to_skip <- 0
  while(TRUE) {
    line <- readLines(con, n = 1)
    if(length(line) == 0) break
    if(startsWith(line, "#")) {
      lines_to_skip <- lines_to_skip + 1
    } else {
      break
    }
  }
  close(con)
  
  cat("Header lines to skip:", lines_to_skip, "\n")
  
  if (requireNamespace("readr", quietly = TRUE)) {
    gitter_data <- readr::read_tsv(gitter_dat_path, skip = lines_to_skip - 1, show_col_types = FALSE)
    colnames(gitter_data)[1] <- "row"
    
    cat("Gitter data dimensions:", dim(gitter_data), "\n")
    cat("First few rows:\n")
    print(head(gitter_data, 3))
    
    # Check first row
    first_row_data <- gitter_data[gitter_data$row == 1, ]
    cat("First row colonies:", nrow(first_row_data), "\n")
    cat("Mean y-coordinate of first row:", mean(first_row_data$y, na.rm = TRUE), "\n")
  }
}

cat("\nBasic functionality test completed.\n")