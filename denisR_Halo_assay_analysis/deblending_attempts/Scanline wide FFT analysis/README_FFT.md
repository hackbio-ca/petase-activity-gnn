# FFT Bleed Suppression for Halo Assay Analysis

## Overview

This enhanced HAIP pipeline uses **2D Fourier Transform (FFT) filtering** to remove problematic artifacts from halo assay images, specifically:

1. **Low-frequency illumination haze** that makes it hard to distinguish true halos from lighting variations
2. **Periodic bleed from tetramer colony lattices** where halos from neighboring colonies interfere with each other

The FFT cleaning step is strategically placed **after image alignment** but **before HAIP's classification functions**, ensuring that per-colony measurements are accurate while preserving all of HAIP's analytical capabilities.

## Why FFT Filtering?

### The Problem
- **Tetramer colony plates** have colonies arranged in regular 4-colony groups across the plate
- **Halo bleed**: On washed images, inhibition halos from neighboring colonies blur together due to the regular lattice structure
- **Illumination artifacts**: Uneven lighting creates gradual intensity changes that interfere with halo detection
- **Standard approaches fail**: Simple background subtraction can't handle the periodic interference patterns

### The Solution
FFT filtering works by:
1. **Converting images to frequency domain** where periodic patterns appear as distinct peaks
2. **Removing specific frequencies** corresponding to the tetramer lattice spacing
3. **Suppressing low frequencies** that correspond to gradual illumination changes
4. **Preserving halo frequencies** that correspond to the actual biological signal

This is like using noise-canceling headphones that specifically target the frequencies of unwanted sounds while preserving the music you want to hear.

## Pipeline Overview

```
1. Load background, unwashed, and washed images
         ↓
2. Align images using phase correlation
         ↓
3. FFT bleed suppression on washed image ← NEW STEP
         ↓
4. HAIP classification (calculate_edge_intensity, classify_pixels)
         ↓
5. Results with improved per-colony accuracy
```

## Quick Start

### Basic Usage
```r
# Source the enhanced pipeline
source("R/main_with_fft.R")

# The script will automatically:
# - Load your images
# - Align them
# - Apply FFT cleaning
# - Run HAIP analysis
# - Show before/after comparison
```

### Using Preset Configurations
```r
# Load configuration system
source("R/fft_config.R")

# Try different presets
config_conservative <- fft_preset_config("conservative")  # Light filtering
config_standard <- fft_preset_config("standard")         # Balanced (default)
config_aggressive <- fft_preset_config("aggressive")     # Strong filtering

# Print configuration details
print(config_standard)
```

### Manual Parameter Tuning
```r
# Create custom configuration
my_config <- create_fft_config(
  hp_radius = 5,        # Remove more low frequencies
  notch_radius = 3,     # Larger notches around lattice peaks
  spacing_hint = c(120, 120)  # Manual colony spacing
)

# Validate configuration
validate_fft_config(my_config)
```

## Parameter Tuning Guide

### High-Pass Filter Parameters

**`hp_radius`** (default: 6)
- **What it does**: Removes low-frequency illumination gradients
- **Lower values (3-5)**: Remove more low frequencies, stronger haze removal
- **Higher values (8-12)**: More conservative, preserves larger-scale features
- **Too low**: May remove genuine large halos
- **Too high**: Won't remove illumination artifacts

**`hp_slope`** (default: 3)
- **What it does**: Controls how sharply the filter cuts off frequencies
- **Lower values (1-2)**: Gradual transition, softer filtering
- **Higher values (4-6)**: Sharp cutoff, more aggressive filtering
- **Too low**: May not cleanly separate wanted/unwanted frequencies
- **Too high**: Can create ringing artifacts

### Notch Filter Parameters

**`notch_radius`** (default: 2)
- **What it does**: Size of "holes" punched around lattice frequency peaks
- **Smaller values (1-1.5)**: Precise removal of just the peak frequencies
- **Larger values (3-4)**: Remove broader range around each peak
- **Too small**: May not fully suppress lattice bleed
- **Too large**: May remove legitimate halo frequencies

**`harmonics`** (default: 1:3)
- **What it does**: Which multiples of the fundamental lattice frequency to suppress
- **Fewer harmonics (1:2)**: Conservative, only removes main lattice frequencies
- **More harmonics (1:4)**: Aggressive, removes finer periodic structures
- **Too few**: Subtle lattice patterns may remain
- **Too many**: May over-filter and remove legitimate features

### Lattice Detection

**`spacing_hint`** (default: NULL for auto-detection)
- **When to use**: If auto-detection fails or for reproducible results
- **How to find**: Measure distance between colony centers in pixels
- **Format**: `c(x_spacing, y_spacing)` in pixels
- **Example**: `c(150, 150)` for 150-pixel spacing in both directions

**`auto_detect_lattice`** (default: TRUE)
- **TRUE**: Automatically find lattice peaks in the FFT
- **FALSE**: Use manual `spacing_hint` instead
- **Auto-detection may fail if**: Lattice is irregular, image is very noisy, or colonies are sparse

## Troubleshooting

### Common Issues and Solutions

**Problem**: Halos look over-smoothed or washed out
- **Solution**: Decrease `hp_radius` (try 4-5) or increase `hp_slope` (try 4-5)
- **Cause**: High-pass filter is removing legitimate halo frequencies

**Problem**: Lattice bleed still visible between colonies
- **Solution**: Increase `notch_radius` (try 3-4) or add more `harmonics` (try 1:4)
- **Cause**: Notch filters aren't fully suppressing the periodic patterns

**Problem**: Strange artifacts or ringing around colonies
- **Solution**: Decrease `hp_slope` (try 2) or reduce `harmonics` (try 1:2)
- **Cause**: Filters are too aggressive and creating edge artifacts

**Problem**: Auto-detection isn't finding the right lattice spacing
- **Solution**: Set `spacing_hint` manually by measuring your colony spacing
- **How to measure**: Open image in ImageJ, measure distance between colony centers

**Problem**: Images don't align properly
- **Solution**: Try `use_phase_correlation = FALSE` for simpler alignment
- **Alternative**: Check that your images are actually the same plate

### Optimization Workflow

1. **Start with "standard" preset** - works for most cases
2. **Check if lattice bleed is gone** - if not, try "aggressive" preset
3. **Check if halos look natural** - if over-filtered, try "conservative" preset
4. **Fine-tune individual parameters** based on the guide above
5. **Test on multiple plates** to ensure settings work consistently

## Advanced Usage

### Processing Multiple Plates
```r
# Define your optimal configuration once
optimal_config <- create_fft_config(
  hp_radius = 5,
  notch_radius = 2.5,
  harmonics = 1:3
)

# Apply to multiple plates
plates <- c("plate1", "plate2", "plate3")
for (plate in plates) {
  # Load images for this plate
  # Apply FFT with optimal_config
  # Save results
}
```

### Comparing Different Settings
```r
# Test multiple configurations
configs <- list(
  conservative = fft_preset_config("conservative"),
  standard = fft_preset_config("standard"),
  aggressive = fft_preset_config("aggressive")
)

# Apply each and compare results
results <- lapply(configs, function(config) {
  # Apply FFT with this config
  # Return classification results
})
```

## Technical Details

### Image Alignment
- **Phase correlation** finds optimal translation between images using FFT
- **Cross-correlation** uses brute-force search (slower but more robust)
- **Why alignment matters**: Ensures background subtraction and FFT work correctly

### FFT Processing
- **Window function** applied to reduce edge artifacts
- **Zero-padding** may be used to improve frequency resolution
- **Filter design** uses smooth transitions to avoid ringing

### Integration with HAIP
- **Preserves all HAIP outputs**: colony_summary, pixel classifications, etc.
- **Compatible with existing workflows**: just replace `washed_img` with `washed_clean`
- **Maintains coordinate systems**: no transformation of colony positions needed

## Files Added

- `R/align_images.R` - Image alignment functions
- `R/fft_bleed_suppression.R` - Core FFT filtering
- `R/fft_config.R` - Parameter management and presets
- `R/main_with_fft.R` - Complete pipeline with FFT
- `README_FFT.md` - This documentation

## Requirements

```r
# Additional packages needed
install.packages(c("fftwtools"))  # For FFT operations
```

All other dependencies are the same as standard HAIP (imager, dplyr, ggplot2, readr, gridExtra).

## Citation

If you use this FFT enhancement in your research, please cite both:
1. The original HAIP package/paper
2. This FFT enhancement (include GitHub repository or relevant publication)

## Support

For questions about:
- **Parameter tuning**: Start with presets, then follow the tuning guide
- **Technical issues**: Check the troubleshooting section
- **Integration problems**: Ensure all dependencies are installed and images are properly formatted