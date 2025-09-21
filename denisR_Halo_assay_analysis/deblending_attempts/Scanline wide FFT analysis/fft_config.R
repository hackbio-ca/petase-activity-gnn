#' FFT Bleed Suppression Configuration
#'
#' This file contains configurable parameters for the FFT bleed suppression
#' pipeline. Modify these values to tune the cleaning for your specific
#' plate characteristics.
#'
#' @export

#' Create FFT Configuration Object
#'
#' Creates a configuration list with all tunable parameters for FFT bleed
#' suppression. This makes it easy to experiment with different settings
#' and maintain consistent parameters across analyses.
#'
#' @param hp_radius `numeric(1)`. High-pass filter radius for removing 
#'   low-frequency illumination haze (default 6). Smaller values remove
#'   more low frequencies. Try 4-10.
#' @param hp_slope `numeric(1)`. High-pass filter slope/softness (default 3).
#'   Higher values create sharper frequency cutoffs. Try 2-5.
#' @param notch_radius `numeric(1)`. Radius of notch filters around detected
#'   lattice peaks (default 2). Larger values remove more frequency content
#'   around each peak. Try 1-4.
#' @param harmonics `numeric`. Which harmonics of the fundamental lattice
#'   frequency to suppress (default 1:3). Higher harmonics correspond to
#'   finer periodic structures. Try 1:2 for lighter filtering or 1:4 for
#'   stronger suppression.
#' @param spacing_hint `numeric(2)` or `NULL`. Manual hint for colony spacing
#'   as [x_spacing, y_spacing] in pixels. If NULL, auto-detects from FFT.
#'   Use this if auto-detection fails or for consistent results.
#' @param auto_detect_lattice `logical(1)`. Whether to automatically detect
#'   lattice peaks in the FFT (default TRUE). Set FALSE if using spacing_hint.
#' @param background_subtract `logical(1)`. Whether to subtract background
#'   before FFT filtering (default FALSE). Usually handled by HAIP pipeline.
#' @param use_phase_correlation `logical(1)`. Use phase correlation (TRUE)
#'   or simple cross-correlation (FALSE) for image alignment (default TRUE).
#' @param alignment_search_radius `numeric(1)`. Search radius for simple
#'   alignment method (default 50). Only used if use_phase_correlation = FALSE.
#'
#' @return `list` containing all configuration parameters.
#' @export
create_fft_config <- function(hp_radius = 6,
                             hp_slope = 3,
                             notch_radius = 2,
                             harmonics = 1:3,
                             spacing_hint = NULL,
                             auto_detect_lattice = TRUE,
                             background_subtract = FALSE,
                             use_phase_correlation = TRUE,
                             alignment_search_radius = 50) {
  
  config <- list(
    # High-pass filter parameters
    hp_radius = hp_radius,
    hp_slope = hp_slope,
    
    # Notch filter parameters  
    notch_radius = notch_radius,
    harmonics = harmonics,
    
    # Lattice detection
    spacing_hint = spacing_hint,
    auto_detect_lattice = auto_detect_lattice,
    
    # Processing options
    background_subtract = background_subtract,
    
    # Alignment parameters
    use_phase_correlation = use_phase_correlation,
    alignment_search_radius = alignment_search_radius
  )
  
  class(config) <- "fft_config"
  return(config)
}

#' Print FFT Configuration
#' @param x FFT configuration object
#' @param ... Additional arguments (ignored)
#' @export
print.fft_config <- function(x, ...) {
  cat("FFT Bleed Suppression Configuration\n")
  cat("=====================================\n")
  cat("High-pass filter:\n")
  cat("  hp_radius:", x$hp_radius, "(lower = remove more low freq)\n")
  cat("  hp_slope: ", x$hp_slope, "(higher = sharper cutoff)\n")
  cat("\nNotch filter:\n")
  cat("  notch_radius:", x$notch_radius, "(larger = remove more around peaks)\n")
  cat("  harmonics:   ", paste(x$harmonics, collapse = ", "), "\n")
  cat("\nLattice detection:\n")
  if (is.null(x$spacing_hint)) {
    cat("  spacing_hint: AUTO-DETECT\n")
  } else {
    cat("  spacing_hint:", paste(x$spacing_hint, collapse = ", "), "pixels\n")
  }
  cat("  auto_detect: ", x$auto_detect_lattice, "\n")
  cat("\nAlignment:\n")
  cat("  method:", ifelse(x$use_phase_correlation, "phase correlation", "cross-correlation"), "\n")
  if (!x$use_phase_correlation) {
    cat("  search_radius:", x$alignment_search_radius, "\n")
  }
}

#' Preset Configurations for Common Scenarios
#'
#' Pre-defined parameter sets for typical halo assay scenarios.
#' These can serve as starting points for optimization.
#'
#' @param preset `character(1)`. One of:
#'   - "conservative": Light filtering, preserves fine details
#'   - "standard": Balanced filtering (default)
#'   - "aggressive": Strong filtering for heavily blended halos
#'   - "manual_spacing": Template for manual spacing specification
#'
#' @return FFT configuration object
#' @export
fft_preset_config <- function(preset = "standard") {
  
  configs <- list(
    conservative = create_fft_config(
      hp_radius = 8,
      hp_slope = 2,
      notch_radius = 1.5,
      harmonics = 1:2
    ),
    
    standard = create_fft_config(
      hp_radius = 6,
      hp_slope = 3,
      notch_radius = 2,
      harmonics = 1:3
    ),
    
    aggressive = create_fft_config(
      hp_radius = 4,
      hp_slope = 4,
      notch_radius = 3,
      harmonics = 1:4
    ),
    
    manual_spacing = create_fft_config(
      hp_radius = 6,
      hp_slope = 3,
      notch_radius = 2,
      harmonics = 1:3,
      spacing_hint = c(150, 150),  # Adjust these values!
      auto_detect_lattice = FALSE
    )
  )
  
  if (!preset %in% names(configs)) {
    stop("Unknown preset: ", preset, 
         ". Available presets: ", paste(names(configs), collapse = ", "))
  }
  
  config <- configs[[preset]]
  attr(config, "preset") <- preset
  return(config)
}

#' Validate FFT Configuration
#'
#' Checks that configuration parameters are within reasonable ranges
#' and warns about potentially problematic combinations.
#'
#' @param config FFT configuration object
#' @return `logical(1)` indicating if configuration is valid
#' @export
validate_fft_config <- function(config) {
  valid <- TRUE
  warnings <- character(0)
  
  # Check parameter ranges
  if (config$hp_radius < 1 || config$hp_radius > 20) {
    warnings <- c(warnings, "hp_radius should typically be between 1-20")
  }
  
  if (config$hp_slope < 1 || config$hp_slope > 10) {
    warnings <- c(warnings, "hp_slope should typically be between 1-10")  
  }
  
  if (config$notch_radius < 0.5 || config$notch_radius > 10) {
    warnings <- c(warnings, "notch_radius should typically be between 0.5-10")
  }
  
  if (max(config$harmonics) > 5) {
    warnings <- c(warnings, "harmonics above 5 may cause over-filtering")
  }
  
  # Check logical combinations
  if (!is.null(config$spacing_hint) && config$auto_detect_lattice) {
    warnings <- c(warnings, "auto_detect_lattice should be FALSE when using spacing_hint")
  }
  
  if (length(warnings) > 0) {
    cat("Configuration warnings:\n")
    for (w in warnings) {
      cat("  -", w, "\n")
    }
    cat("\n")
  }
  
  return(length(warnings) == 0)
}