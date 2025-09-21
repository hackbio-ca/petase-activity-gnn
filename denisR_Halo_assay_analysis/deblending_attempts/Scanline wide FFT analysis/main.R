# Install & load
install.packages(c("imager","dplyr","ggplot2"))
library(imager); library(dplyr); library(ggplot2)

install.packages("readr")
library(readr)
install.packages("gridExtra")
library(gridExtra)

# Source repo functions (or load package if this is a package)
source("R/initialize_colony_data.R")
source("R/create_grid_boxes.R")
source("R/calculate_average_intensity.R")
source("R/calculate_edge_intensity.R")
source("R/classify_pixels.R")
source("R/plot_pixel_classification.R")  # optional
source("R/rgb_to_gray.R")

# Source FFT bleed suppression functions (new)
source("R/align_images.R")
source("R/fft_bleed_suppression.R")
source("R/fft_config.R")

# Load images (use forward slashes on Windows)
bg_img     <- imager::load.image("background/background/background_BHET25_5_1.JPG")
unwashed   <- imager::load.image("all_plate_images/all_plate_images/colony_images/unwashed/BHET25_8h_5_1.JPG")
washed_img <- imager::load.image("all_plate_images/all_plate_images/wash_images/washed_BHET25_8h_5_1.JPG")

# Prepare colony data (expects a gitter .dat file next to the unwashed image)
colony_dat_path <- "all_plate_images/all_plate_images/colony_images/unwashed/BHET12.5_8h_5_1.JPG.dat"
colony_data <- initialize_colony_data(colony_dat_path, unwashed)
colony_data <- create_grid_boxes(colony_data)
colony_data <- calculate_average_intensity(bg_img, colony_data)

# Compute edge normalization (returns a list; use $difference)
edge_norm <- calculate_edge_intensity(bg_img, washed_img, edge_width = 150)$difference

# Run classification
res <- classify_pixels(washed_img, colony_data, edge_norm)

# Inspect outputs
str(res$classified_data)   # data.frame with x,y,intensity,label,colony_id
head(res$colony_summary)   # per-colony counts / medians
res$total_stats

# Quick ggplot of classified pixels (example)
ggplot(res$classified_data, aes(x = x, y = y, color = label)) +
  geom_point(size = 0.4) + scale_y_reverse() + theme_minimal()
