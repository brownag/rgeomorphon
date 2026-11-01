#' Calculate Focal (Moving Window) PTRM
#'
#' Calculates Perceived Terrain Realism Metric (PTRM) for each pixel using
#' a circular moving window. Produces a spatially-varying measure of local
#' terrain complexity and realism.
#'
#' @section Algorithm:
#' Uses a two-pass approach for efficiency:
#'
#' 1. **Calculate geomorphons globally** (if `mode = "dem"`)
#'    - Applies the full geomorphon classification algorithm once
#'    - Produces classification raster with classes 1-10
#'
#' 2. **Apply focal PTRM calculation** to each pixel's neighborhood
#'    - For each interior pixel, extracts geomorphon classes within circular window
#'    - Calculates proportions of each geomorphon form
#'    - Applies PTRM formula using empirical coefficients
#'    - Returns spatially-varying PTRM scores
#'
#' This two-pass approach is much more efficient than recalculating geomorphons
#' for each window (~38x speedup for typical parameters).
#'
#' @section Window Type:
#' Uses a circular window based on Euclidean distance. For an integer
#' `radius = r`, includes all cells where `sqrt(dx^2 + dy^2) <= r`.
#'
#' @section Border Handling:
#' Cells within `radius` distance of raster edges are set to NA, as their
#' neighborhoods are incomplete. This prevents edge artifacts in the output.
#'
#' @section Memory Requirements:
#' For small/medium rasters, uses the fast C++ implementation requiring
#' approximately 4 bytes per input cell plus 8 bytes per output cell (about
#' 12 bytes per cell) in RAM. For large rasters that exceed half of the
#' available memory, automatically falls back to `terra::focal()`, which is
#' slower but memory-efficient.
#'
#' @param x `SpatRaster`. Digital Elevation Model (when `mode = "dem"`) or
#'   pre-computed geomorphon classification (when `mode = "class"`).
#' @param radius Numeric. Integer window radius in cells (not map units).
#'   Default: 3. Larger values compute PTRM from larger neighborhoods.
#' @param mode Character. Input type: `"dem"` to calculate geomorphons internally,
#'   or `"class"` if input is pre-computed geomorphon classification.
#'   Default: `"dem"`.
#' @param coefficients Named numeric vector. Custom PTRM coefficients, either
#'   named for the requested form system (length k) or a named ten-form vector,
#'   which is reduced to the requested system. If `NULL`, uses the paper weights
#'   from Rajasekaran et al. (2022), derived for reduced form systems.
#'   Default: `NULL`.
#' @param forms Character. Form system to use: `"forms10"`, `"forms6"`,
#'   `"forms5"`, or `"forms4"`. Default: `"forms10"`.
#' @param intercept Numeric. PTRM intercept. Default: `-38.02`.
#' @param normalizer Numeric. PTRM normalizer; must be finite and non-zero.
#'   Default: `69.96`.
#' @param search Numeric. Search radius for geomorphon calculation (only used
#'   when `mode = "dem"`). Default: 3.
#' @param skip Numeric. Skip radius for geomorphon calculation (only used when
#'   `mode = "dem"`). Default: 0.
#' @param flat_angle_deg Numeric. Flatness angle threshold in degrees (only used
#'   when `mode = "dem"`). Default: 1.0.
#' @param dist Numeric. Flatness distance parameter (only used when `mode = "dem"`).
#'   Default: 0.
#' @param filename Character. Output filename. If provided, the result is written
#'   to disk instead of returned in memory. Default: `NULL`.
#' @param overwrite Logical. Overwrite existing file? Default: `FALSE`.
#' @param ... Additional arguments passed to [geomorphons()] when `mode = "dem"`.
#'
#' @return
#' A `SpatRaster` with raw focal PTRM scores. Values are not normalized to
#' `[0, 1]` and may fall outside that range. Border cells are NA.
#'
#' @details
#' The PTRM score measures local terrain realism on a continuous scale. Higher
#' values indicate more realistic terrain patterns according to the perception
#' model of Rajasekaran et al. (2022). The focal version applies this metric
#' locally to map spatial variation in terrain complexity.
#'
#' **Computational Complexity:**
#' - Geomorphons: O(N × 8R) where R = search radius
#' - Focal PTRM: O(N × W) where W ≈ π × radius²
#' - Total: O(N × (8R + W))
#'
#' For a 1000×1000 DEM with radius=5, search=7:
#' - Expected runtime: 3-6 seconds on modern multi-core system
#' - Memory: ~20 MB
#'
#' @references
#' Rajasekaran, S.D., Kang, H., Čadík, M., Galin, E., Guérin, E., Peytavie, A.,
#' Slavík, P., & Mech, R. (2022). PTRM: Perceived Terrain Realism Metric.
#' ACM Transactions on Applied Perception, 19(2), Article 5.
#' \doi{10.1145/3514244}
#'
#' Jasiewicz, J., Stepinski, T. (2013). Geomorphons - a pattern recognition
#' approach to classification and mapping of landforms. Geomorphology, 182,
#' 147-156. \doi{10.1016/j.geomorph.2012.11.005}
#'
#' @seealso [geomorphon_ptrm()], [geomorphons()], [terra::focal()]
#'
#' @export
#' @examplesIf requireNamespace("terra")
#' library(terra)
#' library(rgeomorphon)
#'
#' # Load example DEM
#' data("volcano", package = "datasets")
#' dem <- rast(volcano)
#' crs(dem) <- crs("EPSG:2193")
#' ext(dem) <- c(1756968, 1757578, 5917000, 5917870)
#'
#' # Calculate focal PTRM with 5-cell radius
#' ptrm_focal <- geomorphon_ptrm_focal(dem, radius = 5, mode = "dem")
#'
#' # Compare different window sizes
#' ptrm_r3 <- geomorphon_ptrm_focal(dem, radius = 3)
#' ptrm_r7 <- geomorphon_ptrm_focal(dem, radius = 7)
#'
geomorphon_ptrm_focal <- function(
  x,
  radius = 3,
  mode = c("dem", "class"),
  coefficients = NULL,
  forms = "forms10",
  intercept = -38.02,
  normalizer = 69.96,
  search = 3,
  skip = 0,
  flat_angle_deg = 1.0,
  dist = 0,
  filename = NULL,
  overwrite = FALSE,
  ...) {

  # Validate terra package availability
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for this function", call. = FALSE)
  }

  # Input validation
  if (!inherits(x, "SpatRaster")) {
    stop("Input 'x' must be a SpatRaster object", call. = FALSE)
  }

  if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) ||
      !is.finite(radius) || radius < 0) {
    stop("radius must be a non-negative number", call. = FALSE)
  }

  if (radius != as.integer(radius)) {
    stop("radius must be a whole number of cells", call. = FALSE)
  }

  max_radius <- floor((min(terra::nrow(x), terra::ncol(x)) - 1) / 2)
  if (radius > max_radius) {
    stop(
      "radius is too large for this raster: maximum is ", max_radius,
      " (the C++ kernel requires radius <= floor((min(nrow, ncol) - 1) / 2))",
      call. = FALSE
    )
  }

  # Match and validate mode parameter
  mode <- match.arg(mode)

  # Shared PTRM model: coefficients, intercept, normalizer, labels and codes
  model <- .get_ptrm_model(
    forms = forms,
    coefficients = coefficients,
    intercept = intercept,
    normalizer = normalizer
  )

  # Step 1: Get geomorphon classification
  if (mode == "dem") {
    # Calculate geomorphons from DEM
    geo_class <- geomorphons(
      x,
      search = search,
      skip = skip,
      flat_angle_deg = flat_angle_deg,
      dist = dist,
      forms = forms,
      ...
    )
  } else {
    # Use provided classification directly
    geo_class <- x
  }

  # Step 2: Apply focal PTRM calculation
  result <- .geomorphon_ptrm_focal_impl(
    geo_class,
    radius,
    model,
    force_cpp = NULL
  )

  # Write to file if requested
  if (!is.null(filename)) {
    result <- terra::writeRaster(
      result,
      filename = filename,
      overwrite = overwrite,
      datatype = "FLT4S"
    )
  }

  return(result)
}

#' @noRd
.geomorphon_ptrm_focal_impl <- function(geo_class, radius, model,
                                        force_cpp = NULL) {
  # Memory-aware implementation selector.
  # force_cpp = NULL picks automatically; TRUE/FALSE is used by tests.
  remapped <- .remap_ptrm_classes(geo_class, model$classes)

  use_cpp <- if (is.null(force_cpp)) {
    .focal_memory_ok(geo_class, bytes_per_cell = 12)
  } else {
    isTRUE(force_cpp)
  }

  if (use_cpp) {
    # Fast path: use C++ implementation
    mat <- terra::as.matrix(remapped, wide = TRUE)

    ptrm_result <- geomorphon_ptrm_focal_cpp(
      geomorphon_class = mat,
      radius = as.integer(radius),
      coefficients = model$coefficients,
      intercept = model$intercept,
      normalizer = model$normalizer,
      na_value = NA_real_
    )

    result <- .new_ptrm_raster(geo_class)
    terra::values(result) <- as.vector(t(ptrm_result))
  } else {
    # Fallback path: use terra::focal() for memory efficiency
    if (is.null(force_cpp)) {
      warning(
        "Raster too large for fast C++ implementation (",
        format(terra::ncell(geo_class), big.mark = ","), " pixels). ",
        "Using terra::focal() instead (slower but memory-efficient).",
        call. = FALSE
      )
    }

    k <- length(model$coefficients)

    # Define PTRM calculation function for moving window
    ptrm_window_fun <- function(v, ...) {
      v <- as.integer(v[!is.na(v)])

      # Skip all-NA windows
      if (length(v) == 0) {
        return(NA_real_)
      }

      # Count only contiguous class codes 1..k, matching the C++ kernel
      counts <- tabulate(v, nbins = k)
      n_total <- sum(counts)
      if (n_total == 0) {
        return(NA_real_)
      }

      weighted_sum <- sum((counts / n_total) * model$coefficients)
      (model$intercept + weighted_sum) / model$normalizer
    }

    # Explicit circular window using sqrt(dx^2 + dy^2) <= radius
    window_size <- 2L * as.integer(radius) + 1L
    offsets <- expand.grid(
      dr = -as.integer(radius):as.integer(radius),
      dc = -as.integer(radius):as.integer(radius)
    )
    window_mat <- matrix(0, nrow = window_size, ncol = window_size)
    in_window <- sqrt(offsets$dr^2 + offsets$dc^2) <= radius
    window_mat[cbind(offsets$dr + radius + 1L, offsets$dc + radius + 1L)] <-
      as.integer(in_window)

    focal_result <- terra::focal(
      remapped,
      w = window_mat,
      fun = ptrm_window_fun,
      na.rm = TRUE,
      na.policy = "all"
    )

    # Match C++ border handling: cells within radius of an edge are NA
    result <- .new_ptrm_raster(geo_class)
    out_mat <- terra::as.matrix(focal_result, wide = TRUE)
    if (radius > 0) {
      r <- as.integer(radius)
      out_mat[seq_len(r), ] <- NA_real_
      out_mat[nrow(out_mat) - seq_len(r) + 1L, ] <- NA_real_
      out_mat[, seq_len(r)] <- NA_real_
      out_mat[, ncol(out_mat) - seq_len(r) + 1L] <- NA_real_
    }
    terra::values(result) <- as.vector(t(out_mat))
  }

  names(result) <- "ptrm"
  return(result)
}

#' @noRd
.focal_memory_ok <- function(x, bytes_per_cell) {
  # Shared memory gate for focal C++ paths.
  # Returns TRUE when the estimated working set fits in half of available RAM.
  required_gb <- terra::ncell(x) * bytes_per_cell / 1e9
  available_gb <- terra::mem_info(x, print = FALSE)["available"]
  required_gb < available_gb / 2
}

#' @noRd
.new_ptrm_raster <- function(template) {
  # Fresh raster with geometry only: never inherit values, levels or colours.
  terra::rast(
    terra::ext(template),
    nrows = terra::nrow(template),
    ncols = terra::ncol(template),
    crs = terra::crs(template)
  )
}

#' @noRd
.remap_ptrm_classes <- function(geo_class, classes) {
  # Map system-specific codes (e.g. forms4: 1, 3, 6, 9) to contiguous 1..k.
  # Values outside the requested system become NA and are ignored downstream.
  remapped <- match(terra::values(geo_class, mat = FALSE), classes)
  result <- .new_ptrm_raster(geo_class)
  terra::values(result) <- remapped
  result
}
