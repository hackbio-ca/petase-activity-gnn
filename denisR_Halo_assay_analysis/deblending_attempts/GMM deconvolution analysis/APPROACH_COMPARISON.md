# Comparison of GMM Approaches
# 
# This document explains the difference between the two GMM approaches

cat("=== GMM Approach Comparison ===\n\n")

cat("ORIGINAL APPROACH (gmm_scanline.R):\n")
cat("- Scanline analysis along y-row\n")
cat("- Defines windows around colonies\n") 
cat("- Fits 2-component GMM on INTENSITIES within each window\n")
cat("- Components: background vs colony signal\n")
cat("- Area = sum of responsibilities × (adjusted - baseline)\n")
cat("- One GMM per colony window\n\n")

cat("NEW APPROACH (gmm_1d_weighted.R):\n")
cat("- 1D GMM over X-COORDINATES across entire scanline\n")
cat("- Uses intensity weights = max(adjusted - baseline, 0)\n")
cat("- Fits multi-component GMM on X positions, weighted by intensity\n")
cat("- Components: spatial locations of colonies\n")
cat("- Initialize K components from gitter centers\n")
cat("- Area = responsibility-weighted sum of positive signal\n")
cat("- Wide rolling quantile baseline (≈500 px, q=0.15)\n")
cat("- Prunes components by mass and SNR\n")
cat("- Captures full colony tails naturally\n\n")

cat("KEY DIFFERENCES:\n")
cat("1. GMM DOMAIN:\n")
cat("   - Original: GMM on intensities (y-values)\n") 
cat("   - New: GMM on positions (x-values)\n\n")

cat("2. INITIALIZATION:\n")
cat("   - Original: 2 components (bg + colony) per window\n")
cat("   - New: K components from gitter centers\n\n")

cat("3. BASELINE:\n")
cat("   - Original: Running median detrend\n")
cat("   - New: Wide rolling quantile (q=0.15)\n\n")

cat("4. AREA CALCULATION:\n")
cat("   - Original: Window-based with presence checks\n")
cat("   - New: Global responsibility-weighted integration\n\n")

cat("5. TAIL HANDLING:\n")
cat("   - Original: Limited by window bounds + snapping\n")
cat("   - New: Full tails captured by Gaussian components\n\n")

cat("The new approach should be more robust for:\n")
cat("- Overlapping colony tails\n")
cat("- Variable colony sizes\n")
cat("- Weak colonies near strong ones\n")
cat("- Automatic component pruning\n\n")

# Example usage:
cat("USAGE EXAMPLE:\n")
cat('result <- gmm_1d_weighted(\n')
cat('  washed_path = "path/to/washed.jpg",\n')
cat('  background_path = "path/to/background.jpg", \n')
cat('  gitter_dat_path = "path/to/gitter.dat",\n')
cat('  row_index = 1,\n')
cat('  baseline_window = 500,     # Wide baseline\n')
cat('  baseline_quantile = 0.15,  # 15th percentile\n')
cat('  min_mass = 0.01,          # Prune tiny components\n')
cat('  min_snr = 2.0             # Prune low SNR\n')
cat(')\n')