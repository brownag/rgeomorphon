# Terrain Analysis and Metrics Functions

#' @noRd
.calculate_entropy <- function(x) {
  # Remove NAs
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(0)
  }

  # Calculate frequency table and proportions
  freq <- table(x)
  prob <- freq / sum(freq)

  # Remove zero probabilities to avoid log(0)
  prob <- prob[prob > 0]

  # Shannon entropy: H(X) = -sum(p(x) * log2(p(x)))
  -sum(prob * log2(prob))
}

#' @noRd
.calculate_mutual_information <- function(x, y) {
  # Remove NAs (paired)
  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]
  y <- y[valid]

  if (length(x) == 0) {
    return(0)
  }

  # Calculate individual entropies
  hx <- .calculate_entropy(x)
  hy <- .calculate_entropy(y)

  # Calculate joint entropy
  # Create joint categorical variable
  xy <- paste(x, y, sep = "_")
  hxy <- .calculate_entropy(xy)

  # Mutual information: I(X;Y) = H(X) + H(Y) - H(X,Y)
  mi <- hx + hy - hxy

  # Ensure non-negative (due to numerical precision)
  max(0, mi)
}

#' @noRd
.calculate_nmi <- function(x, y) {
  # Remove NAs (paired)
  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]
  y <- y[valid]

  if (length(x) == 0) {
    return(0)
  }

  # Calculate entropies
  hx <- .calculate_entropy(x)
  hy <- .calculate_entropy(y)

  # If both are constant (entropy = 0), return 0
  if (hx + hy == 0) {
    return(0)
  }

  # Calculate mutual information
  mi <- .calculate_mutual_information(x, y)

  # Normalized Mutual Information: NMI = 2 * I(X;Y) / (H(X) + H(Y))
  nmi <- 2 * mi / (hx + hy)

  # Clamp to [0, 1] to handle numerical precision issues
  min(1, max(0, nmi))
}

#' @noRd
.detect_edges <- function(x) {
  # Memory-aware edge detection: uses fast C++ for small/medium rasters,
  # falls back to terra::boundaries() for large rasters that don't fit in RAM

  # Estimate memory requirements for C++ implementation
  # Formula: 2x raster size (input + output matrices as doubles)
  raster_bytes <- terra::ncell(x) * 8
  total_bytes_needed <- raster_bytes * 2
  total_gb_needed <- total_bytes_needed / 1e9

  # Get available system memory (in GB)
  mem_info <- terra::mem_info(x, print = FALSE)
  available_gb <- mem_info["available"]

  # Use C++ if raster fits comfortably in available memory (with 50% safety margin)
  memory_threshold_gb <- available_gb / 2

  if (total_gb_needed < memory_threshold_gb) {
    # Fast path: use C++ implementation
    mat <- terra::as.matrix(x, wide = TRUE)
    edge_mat <- focal_edge_cpp(mat, na_value = NA_real_)
    result <- x
    terra::values(result) <- edge_mat
  } else {
    # Fallback path: use terra::boundaries() for large rasters
    # This is slower but memory-efficient and prevents out-of-memory crashes
    warning(
      "Raster is too large for fast C++ edge detection (",
      format(terra::ncell(x), big.mark = ","), " pixels). ",
      "Using terra::boundaries() instead, which will be slower. ",
      "For better performance, consider processing the raster in tiles.",
      call. = FALSE
    )

    # terra::boundaries() includes border cells, so we need to set them to NA
    # to match the behavior of focal_edge_cpp()
    result <- terra::boundaries(x, classes = TRUE, directions = 8)

    # Set border cells to NA to match focal_edge_cpp behavior
    # (focal_edge_cpp excludes borders to prevent artifacts in concordance analysis)
    # Create a simple mask: set borders to NA, keep interior as is
    border_vals <- terra::values(result)
    nrows <- terra::nrow(x)
    ncols <- terra::ncol(x)
    ncells <- terra::ncell(x)

    # Identify border cell indices
    # Row 1 (cells 1 to ncols) and last row (cells ncells - ncols + 1 to ncells)
    # Column 1 (every ncols-th cell starting at 1) and last column (every ncols-th cell ending at ncols)
    border_indices <- c(
      1:ncols,  # First row
      (ncells - ncols + 1):ncells,  # Last row
      seq(1, ncells, ncols),  # First column
      seq(ncols, ncells, ncols)  # Last column
    )
    border_indices <- unique(border_indices)

    # Set border cells to NA
    border_vals[border_indices] <- NA
    terra::values(result) <- border_vals
  }

  return(result)
}

#' @noRd
.get_ptrm_coefficients <- function() {
  # Coefficients from Rajasekaran et al. (2022), Table 1
  c(
    flat = 25.12,
    peak = 1.75,
    ridge = 7.59,
    shoulder = 7.31,
    spur = 9.02,
    slope = 28.95,
    hollow = 6.71,
    footslope = 7.63,
    valley = 9.61,
    depression = 3.55
  )
}

#' Calculate Perceived Terrain Realism Metric (PTRM)
#'
#' Calculates the Perceived Terrain Realism Metric based on the distribution
#' of geomorphon landform features. The PTRM provides a quantitative estimate
#' of human-perceived terrain realism by relating the proportions of terrain
#' features to their perceived realism.
#'
#' # Distance Calculation and Coordinate Reference Systems
#'
#' The geomorphon classification algorithm assumes planar distances, so it is
#' strongly recommended that elevation data be in a projected coordinate system.
#'
#' # Tiled Processing
#'
#' When `mode = "dem"`, geomorphons are calculated internally. For large DEMs,
#' this may trigger the automatic tiled processing workflow. See [geomorphons()]
#' for details.
#'
#' @param x `SpatRaster`. Either a Digital Elevation Model (when `mode = "dem"`)
#'   or a geomorphon classification raster (when `mode = "class"`).
#' @param mode Character. One of `"dem"` or `"class"`. When `"dem"`, geomorphons
#'   are calculated from the input DEM. When `"class"`, the input is assumed to
#'   be a pre-computed geomorphon classification. Default: `"dem"`.
#' @param search Numeric. Search radius for geomorphon calculation (only used
#'   when `mode = "dem"`). Units depend on the CRS of the input. Default: `3`.
#' @param skip Numeric. Skip radius for geomorphon calculation (only used when
#'   `mode = "dem"`). Default: `0`.
#' @param flat_angle_deg Numeric. Flatness angle threshold in degrees (only
#'   used when `mode = "dem"`). Default: `1.0`.
#' @param dist Numeric. Flatness distance (only used when `mode = "dem"`).
#'   Default: `0`.
#' @param forms Character. Form system to use: `"forms10"`, `"forms6"`,
#'   `"forms5"`, or `"forms4"`. Default: `"forms10"`.
#' @param coefficients Named numeric vector. Custom PTRM coefficients for each
#'   geomorphon class. If `NULL`, uses default values from Rajasekaran et al.
#'   (2022). Names should match geomorphon form names. Default: `NULL`.
#' @param ... Additional arguments passed to [geomorphons()] when
#'   `mode = "dem"`.
#'
#' @return A list containing:
#'   \describe{
#'     \item{`ptrm_score`}{Numeric. PTRM value on a scale from 0 to 1. Higher
#'       values indicate more realistic terrain.}
#'     \item{`raw_score`}{Numeric. Unnormalized score before division.}
#'     \item{`proportions`}{Named numeric vector. Proportion of each geomorphon
#'       form in the input.}
#'     \item{`n_cells`}{Integer. Total number of valid (non-NA) cells analyzed.}
#'   }
#'
#' @details
#' The PTRM is calculated using the formula:
#'
#' \deqn{PTRM = \frac{c_0 + \sum_{i=1}^{n} c_i \times p_i}{D}}
#'
#' where \eqn{c_0} is the intercept (-38.02), \eqn{c_i} are the coefficients
#' for each geomorphon form, \eqn{p_i} are the proportions of each form in
#' the terrain, and \eqn{D} is the divisor (69.96).
#'
#' Higher PTRM scores indicate more realistic terrain from a human perception
#' perspective. Rajasekaran et al. (2022) found that real terrain from Arizona
#' achieved a PTRM of 0.76, ranking in the 78th percentile for perceived
#' realism.
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
#' @seealso [geomorphons()], [geomorphon_categories()]
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
#' # Calculate PTRM from DEM
#' ptrm1 <- geomorphon_ptrm(dem, mode = "dem", search = 7)
#' print(ptrm1$ptrm_score)
#'
#' # Or calculate geomorphons first, then PTRM
#' geo <- geomorphons(dem, search = 7, simplify = TRUE)
#' ptrm2 <- geomorphon_ptrm(geo, mode = "class")
#'
#' # Results should be identical
#' all.equal(ptrm1$ptrm_score, ptrm2$ptrm_score)
#'
geomorphon_ptrm <- function(x,
                           mode = c("dem", "class"),
                           search = 3,
                           skip = 0,
                           flat_angle_deg = 1.0,
                           dist = 0,
                           forms = "forms10",
                           coefficients = NULL,
                           ...) {

  # Validate terra package availability
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for this function", call. = FALSE)
  }

  # Input validation
  if (!inherits(x, "SpatRaster")) {
    stop("Input 'x' must be a SpatRaster object", call. = FALSE)
  }

  # Validate mode parameter
  mode <- match.arg(mode)

  # Get geomorphon classification
  if (mode == "dem") {
    # Calculate geomorphons from DEM
    geo <- geomorphons(
      x,
      search = search,
      skip = skip,
      flat_angle_deg = flat_angle_deg,
      dist = dist,
      forms = forms,
      ...
    )
  } else if (mode == "class") {
    # Use input as geomorphon classification
    geo <- x
  }

  # Extract values and remove NAs
  vals <- terra::values(geo, mat = FALSE, na.rm = TRUE)
  n_cells <- length(vals)

  if (n_cells == 0) {
    stop("No valid (non-NA) cells found in input raster", call. = FALSE)
  }

  # Calculate proportions from frequency table
  freq_table <- table(vals)
  proportions <- as.numeric(freq_table) / n_cells

  # Get form categories to map numeric values to form names
  form_cats <- geomorphon_categories()

  # Map proportions to form names
  form_indices <- as.integer(names(freq_table))
  form_names <- form_cats$form[form_indices]
  names(proportions) <- form_names

  # Get or validate coefficients
  if (is.null(coefficients)) {
    coefficients <- .get_ptrm_coefficients()
  }

  # Validate coefficient structure
  if (any(is.na(coefficients))) {
    stop("Coefficients contain NA values. Please provide valid coefficients",
         call. = FALSE)
  }

  # Ensure all geomorphon forms are represented in proportions
  full_proportions <- setNames(rep(0, length(coefficients)), names(coefficients))
  for (i in seq_along(proportions)) {
    form_name <- names(proportions)[i]
    if (form_name %in% names(full_proportions)) {
      full_proportions[form_name] <- proportions[i]
    }
  }

  # Calculate weighted sum
  weighted_sum <- sum(full_proportions * coefficients, na.rm = TRUE)

  # Apply PTRM formula: (intercept + weighted_sum) / divisor
  intercept <- -38.02
  normalizer <- 69.96
  raw_score <- intercept + weighted_sum
  ptrm_score <- raw_score / normalizer

  # Return structured result
  list(
    ptrm_score = ptrm_score,
    raw_score = raw_score,
    proportions = proportions,
    n_cells = n_cells
  )
}

#' Calculate Structural Concordance Between Geomorphon and Categorical Maps
#'
#' Measures the degree of structural agreement between a geomorphon
#' classification and another categorical map (e.g., soil types, vegetation
#' classes) using Normalized Mutual Information and boundary alignment metrics.
#'
#' # Information-Theoretic Concordance
#'
#' Normalized Mutual Information (NMI) measures the mutual dependence between
#' two categorical variables, independent of the specific class mappings. An
#' NMI of 1 indicates perfect association, while 0 indicates independence.
#'
#' # Boundary Alignment
#'
#' Boundary alignment measures physical proximity between geomorphon and
#' categorical boundaries. The `edge_offset` is the median distance from
#' categorical boundaries to the nearest geomorphon boundary. The
#' `edge_concordance` is the proportion of categorical boundaries within a
#' threshold distance.
#'
#' @param geo_rast `SpatRaster`. Geomorphon classification raster (typically
#'   output from [geomorphons()]).
#' @param cat_rast `SpatRaster`. Categorical comparison raster (e.g., soil map,
#'   vegetation map). Must have matching extent and resolution to `geo_rast`.
#' @param method Character. Which concordance metrics to calculate. One of
#'   `"nmi"` (Normalized Mutual Information only), `"boundary"` (boundary
#'   alignment only), or `"both"` (all metrics). Default: `"both"`.
#' @param threshold_cells Integer. Buffer distance in cells for boundary
#'   concordance calculation. Categorical boundaries within this distance of
#'   geomorphon boundaries are counted as concordant. Default: `2`.
#' @param ... Reserved for future extensions.
#'
#' @return A list with the following components (depending on `method`):
#'   \describe{
#'     \item{`nmi`}{Numeric. Normalized Mutual Information (0-1). Higher values
#'       indicate stronger statistical association between classifications.
#'       Only included when method includes "nmi".}
#'     \item{`edge_offset`}{Numeric. Mean offset distance (in map units) from
#'       categorical boundaries to the nearest geomorphon boundary. Lower values
#'       indicate better boundary alignment. Only included when method includes
#'       "boundary".}
#'     \item{`edge_concordance`}{Numeric. Proportion (0-1) of categorical
#'       boundaries within `threshold_cells` of geomorphon boundaries. Higher
#'       values indicate better alignment. Only included when method includes
#'       "boundary".}
#'     \item{`n_cells`}{Integer. Total number of valid cells analyzed.}
#'     \item{`n_edges_geo`}{Integer. Number of geomorphon boundary cells.
#'       Only included when method includes "boundary".}
#'     \item{`n_edges_cat`}{Integer. Number of categorical boundary cells.
#'       Only included when method includes "boundary".}
#'     \item{`contingency_table`}{Matrix. Cross-tabulation of geomorphon vs.
#'       categorical classes. Only included when method includes "nmi".}
#'   }
#'
#' @details
#' For rasters to be comparable, they must have matching extents and
#' resolutions. Use [terra::crop()], [terra::resample()], or
#' [terra::disagg()] to align rasters before calling this function.
#'
#' The NMI is symmetric and normalized to the range 0-1, with higher values
#' indicating stronger association regardless of how classes map to each other. The
#' boundary metrics are useful for validating whether geomorphological
#' transitions correspond to ecological or pedological boundaries.
#'
#' @note
#' **Memory Requirements:** This function performs edge detection on both input
#' rasters. For small to medium rasters (< 1GB each), a fast C++ implementation
#' is used. For larger rasters, the function automatically falls back to
#' `terra::boundaries()` which is slower but more memory-efficient.
#'
#' Since both rasters are processed simultaneously, the total memory requirement
#' can approach 4x the size of a single raster in the worst case. For rasters
#' larger than available RAM:
#' - Process in tiles using [terra::makeTiles()]
#' - Increase available system memory
#' - Consider splitting the analysis into smaller regions
#'
#' See `misc/MEMORY_CONSTRAINTS_ANALYSIS.md` for technical details.
#'
#' @seealso [geomorphons()], [terra::boundaries()]
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
#' # Calculate geomorphons
#' geo <- geomorphons(dem, search = 7, simplify = TRUE)
#'
#' # Create synthetic categorical map for demonstration
#' cat_map <- dem
#' values(cat_map) <- cut(values(dem), breaks = 5, labels = FALSE)
#'
#' # Calculate concordance
#' conc <- geomorphon_concordance(geo, cat_map)
#' print(conc$nmi)
#' print(conc$edge_concordance)
#'
geomorphon_concordance <- function(geo_rast,
                                  cat_rast,
                                  method = c("both", "nmi", "boundary"),
                                  threshold_cells = 2,
                                  ...) {

  # Validate terra package availability
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for this function", call. = FALSE)
  }

  # Input validation
  if (!inherits(geo_rast, "SpatRaster")) {
    stop("Input 'geo_rast' must be a SpatRaster object", call. = FALSE)
  }

  if (!inherits(cat_rast, "SpatRaster")) {
    stop("Input 'cat_rast' must be a SpatRaster object", call. = FALSE)
  }

  # Validate method parameter
  method <- match.arg(method)

  # Check extents match
  if (!isTRUE(all.equal(terra::ext(geo_rast), terra::ext(cat_rast)))) {
    stop("Rasters must have matching extents. ",
         "Use terra::crop() or terra::resample() to align them.",
         call. = FALSE)
  }

  # Check resolutions match (within tolerance)
  res_geo <- terra::res(geo_rast)
  res_cat <- terra::res(cat_rast)
  if (!isTRUE(all.equal(res_geo, res_cat, tolerance = 1e-6))) {
    stop("Rasters must have matching resolutions. ",
         "Use terra::resample() to align them.",
         call. = FALSE)
  }

  # Initialize result list
  result <- list()

  # Calculate NMI if requested
  if (method %in% c("nmi", "both")) {
    geo_vals <- terra::values(geo_rast, mat = FALSE)
    cat_vals <- terra::values(cat_rast, mat = FALSE)

    # Remove NA pairs
    valid <- !is.na(geo_vals) & !is.na(cat_vals)
    geo_vals_clean <- geo_vals[valid]
    cat_vals_clean <- cat_vals[valid]

    n_cells <- length(geo_vals_clean)

    if (n_cells == 0) {
      stop("No overlapping valid cells found between rasters",
           call. = FALSE)
    }

    # Calculate NMI
    result$nmi <- .calculate_nmi(geo_vals_clean, cat_vals_clean)
    result$n_cells <- n_cells

    # Add contingency table
    result$contingency_table <- table(
      Geomorphon = geo_vals_clean,
      Category = cat_vals_clean
    )
  }

  # Calculate boundary alignment if requested
  if (method %in% c("boundary", "both")) {

    # Detect edges
    geo_edges <- .detect_edges(geo_rast)
    cat_edges <- .detect_edges(cat_rast)

    # Get edge cell values
    geo_edge_vals <- terra::values(geo_edges, mat = FALSE)
    cat_edge_vals <- terra::values(cat_edges, mat = FALSE)

    # Count edges
    n_edges_geo <- sum(geo_edge_vals == 1, na.rm = TRUE)
    n_edges_cat <- sum(cat_edge_vals == 1, na.rm = TRUE)

    result$n_edges_geo <- n_edges_geo
    result$n_edges_cat <- n_edges_cat

    # Calculate distance metrics if there are edges
    if (n_edges_geo > 0 && n_edges_cat > 0) {
      # Create a raster with NA where there are no cat edges, and non-NA elsewhere
      # Then compute distance from non-NA to NA cells
      cat_edge_dist_rast <- cat_edges

      # Set NA where there are NO edges, and keep values where there ARE edges
      terra::values(cat_edge_dist_rast)[cat_edge_vals != 1] <- NA

      # Compute distance from NA cells to nearest non-NA cell
      # This gives us the distance from each cell to the nearest cat edge
      cat_edge_dist <- terra::distance(cat_edge_dist_rast)

      # Extract distances at geomorphon edge cells
      cat_edge_dist_vals <- terra::values(cat_edge_dist, mat = FALSE)
      geo_edge_dists <- cat_edge_dist_vals[geo_edge_vals == 1]
      geo_edge_dists <- geo_edge_dists[!is.na(geo_edge_dists)]

      if (length(geo_edge_dists) > 0) {
        # Convert to map units
        cell_size <- terra::res(geo_rast)[1]
        edge_dists_mapunits <- geo_edge_dists * cell_size

        # Calculate metrics
        result$edge_offset <- mean(edge_dists_mapunits, na.rm = TRUE)
        threshold_mapunits <- threshold_cells * cell_size
        result$edge_concordance <- sum(edge_dists_mapunits <= threshold_mapunits,
                                       na.rm = TRUE) / length(geo_edge_dists)
      } else {
        result$edge_offset <- NA_real_
        result$edge_concordance <- NA_real_
      }
    } else {
      result$edge_offset <- NA_real_
      result$edge_concordance <- NA_real_
    }
  }

  return(result)
}
