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
.detect_edges <- function(x, force_fallback = FALSE) {
  # Memory-aware edge detection: uses fast C++ for small/medium rasters,
  # falls back to terra::boundaries() for large rasters that don't fit in RAM

  # C++ working set: input + output matrices as doubles.
  # force_fallback exists so tests can exercise the terra path on small rasters.
  if (!force_fallback && .focal_memory_ok(x, bytes_per_cell = 16)) {
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
    # to match the behavior of focal_edge_cpp(). ignoreNA = TRUE skips NA
    # neighbors; differing valid neighbors still make a cell an edge.
    result <- terra::boundaries(x, classes = TRUE, directions = 8,
                                ignoreNA = TRUE)

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

#' @noRd
.get_ptrm_model <- function(forms = "forms10",
                            coefficients = NULL,
                            intercept = -38.02,
                            normalizer = 69.96) {
  form_systems <- c("forms10", "forms6", "forms5", "forms4")

  if (!is.character(forms) || length(forms) != 1L || !forms %in% form_systems) {
    stop(
      "'forms' must be one of: ",
      paste(form_systems, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(intercept) || length(intercept) != 1L ||
      !is.finite(intercept)) {
    stop("'intercept' must be a length-1 finite numeric value", call. = FALSE)
  }

  if (!is.numeric(normalizer) || length(normalizer) != 1L ||
      !is.finite(normalizer) || normalizer == 0) {
    stop(
      "'normalizer' must be a length-1 finite, non-zero numeric value",
      call. = FALSE
    )
  }

  # Paper weights for the ten-form system are the single source of truth
  paper <- .get_ptrm_coefficients()
  paper_labels <- names(paper)

  num_forms <- as.integer(sub("^forms", "", forms))
  forms_table <- get_forms_matrix_cpp(num_forms)
  ten_form_table <- get_forms_matrix_cpp(10L)

  # Reachable codes are the non-zero (non-G_NONE) entries of the forms matrix
  classes <- sort(unique(as.vector(forms_table)[as.vector(forms_table) != 0]))
  k <- length(classes)
  if (k < 1L) {
    stop("The requested form system has no reachable classes", call. = FALSE)
  }
  labels <- paper_labels[classes]

  # Average the ten-form paper weights over every reachable (negative,
  # positive) pair that maps to each reduced class, repeats included
  reduce_weights <- function(weights) {
    vapply(
      classes,
      function(cl) {
        ten_classes <- ten_form_table[forms_table == cl]
        mean(weights[ten_classes])
      },
      numeric(1)
    )
  }

  if (is.null(coefficients)) {
    coefficients <- if (num_forms == 10L) {
      paper
    } else {
      reduce_weights(paper)
    }
  } else {
    if (!is.numeric(coefficients) || is.null(names(coefficients))) {
      stop("'coefficients' must be a named numeric vector", call. = FALSE)
    }

    coeff_names <- tolower(names(coefficients))
    if (anyDuplicated(coeff_names)) {
      stop("'coefficients' must not contain duplicate names", call. = FALSE)
    }

    if (length(coefficients) == k && all(labels %in% coeff_names)) {
      coefficients <- coefficients[match(labels, coeff_names)]
    } else if (length(coefficients) == length(paper_labels) &&
               all(paper_labels %in% coeff_names)) {
      ten_weights <- coefficients[match(paper_labels, coeff_names)]
      names(ten_weights) <- paper_labels
      coefficients <- if (num_forms == 10L) {
        ten_weights
      } else {
        reduce_weights(ten_weights)
      }
    } else {
      stop(
        "'coefficients' must be named for the ",
        forms,
        " classes (",
        paste(labels, collapse = ", "),
        ") or for all ten geomorphon forms",
        call. = FALSE
      )
    }

    if (any(!is.finite(coefficients))) {
      stop("'coefficients' must contain only finite values", call. = FALSE)
    }
  }

  names(coefficients) <- labels

  list(
    classes = classes,
    labels = labels,
    coefficients = coefficients,
    intercept = intercept,
    normalizer = normalizer
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
#'   This is a convenience default, not a value calibrated against the PTRM
#'   paper, so it should not be read as a validated choice.
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
#'   (2022). Names should match geomorphon form names. A ten-form vector is
#'   reduced to the requested `forms` system. Default: `NULL`.
#' @param intercept Numeric. PTRM intercept. Default: `-38.02`.
#' @param normalizer Numeric. PTRM divisor. Must be finite and non-zero.
#'   Default: `69.96`.
#' @param ... Additional arguments passed to [geomorphons()] when
#'   `mode = "dem"`. Ignored when `mode = "class"`.
#'
#' @return A list containing:
#'   \describe{
#'     \item{`ptrm_score`}{Numeric. PTRM value obtained by dividing the raw
#'       score by the normalizer. Scores are raw and may fall outside
#'       `[0, 1]`; higher values indicate more realistic terrain.}
#'     \item{`raw_score`}{Numeric. Unnormalized score before division.}
#'     \item{`proportions`}{Named numeric vector. Proportion of each observed
#'       geomorphon form in the input, named with the canonical system labels.}
#'     \item{`n_cells`}{Integer. Total number of valid (non-NA) cells analyzed.}
#'   }
#'
#' @details
#' The PTRM is calculated using the formula:
#'
#' \deqn{PTRM = \frac{c_0 + \sum_{i=1}^{n} c_i \times p_i}{D}}
#'
#' where \eqn{c_0} is the intercept (-38.02 by default), \eqn{c_i} are the
#' coefficients for each geomorphon form, \eqn{p_i} are the proportions of each
#' form in the terrain, and \eqn{D} is the divisor (69.96 by default). The
#' intercept, normalizer, and coefficients can all be overridden.
#'
#' The returned score is not clamped or rescaled, so it may fall outside
#' `[0, 1]` depending on the input distribution and the model parameters.
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
                           intercept = -38.02,
                           normalizer = 69.96,
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

  # Resolve the PTRM model (classes, labels, coefficients, constants)
  model <- .get_ptrm_model(
    forms = forms,
    coefficients = coefficients,
    intercept = intercept,
    normalizer = normalizer
  )

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
  } else {
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
  observed <- as.integer(names(freq_table))

  # Codes must belong to the requested form system
  invalid <- setdiff(observed, model$classes)
  if (length(invalid) > 0) {
    stop(
      "Invalid geomorphon codes for ",
      forms,
      ": ",
      paste(sort(invalid), collapse = ", "),
      " (expected one of: ",
      paste(model$classes, collapse = ", "),
      ")",
      call. = FALSE
    )
  }

  proportions <- as.numeric(freq_table) / n_cells
  names(proportions) <- model$labels[match(observed, model$classes)]

  # Ensure every class in the model is represented, in class-code order
  full_proportions <- stats::setNames(rep(0, length(model$coefficients)),
                               model$labels)
  full_proportions[names(proportions)] <- proportions

  # Calculate weighted sum
  weighted_sum <- sum(full_proportions * model$coefficients)

  # Apply PTRM formula: (intercept + weighted_sum) / normalizer
  raw_score <- model$intercept + weighted_sum
  ptrm_score <- raw_score / model$normalizer

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
#' categorical boundaries. The `edge_offset` is the mean distance from each
#' categorical boundary cell to the nearest geomorphon boundary cell, in map
#' units of the input CRS. The `edge_concordance` is the proportion of
#' categorical boundary cells within a threshold distance of the nearest
#' geomorphon boundary cell.
#'
#' @param geo_rast `SpatRaster`. Geomorphon classification raster (typically
#'   output from [geomorphons()]).
#' @param cat_rast `SpatRaster`. Categorical comparison raster (e.g., soil map,
#'   vegetation map). Must have matching extent and resolution to `geo_rast`.
#' @param method Character. Which concordance metrics to calculate. One of
#'   `"nmi"` (Normalized Mutual Information only), `"boundary"` (boundary
#'   alignment only), or `"both"` (all metrics). Default: `"both"`.
#' @param threshold_cells Integer. Buffer distance in cells for boundary
#'   concordance calculation, converted to map units with the cell size.
#'   Categorical boundary cells within this distance of the nearest
#'   geomorphon boundary cell are counted as concordant. Must be a finite,
#'   non-negative whole number when `method` includes boundary metrics.
#'   Geographic (longitude/latitude) CRS data trigger a warning because
#'   boundary distances are calculated in CRS map units. Default: `2`.
#'
#' @return A list with the following components (depending on `method`):
#'   \describe{
#'     \item{`nmi`}{Numeric. Normalized Mutual Information (0-1). Higher values
#'       indicate stronger statistical association between classifications.
#'       Only included when method includes "nmi".}
#'     \item{`edge_offset`}{Numeric. Mean distance (in map units) from each
#'       categorical boundary cell to the nearest geomorphon boundary cell.
#'       Lower values indicate better boundary alignment. Only included when
#'       method includes "boundary".}
#'     \item{`edge_concordance`}{Numeric. Proportion (0-1) of categorical
#'       boundary cells within `threshold_cells` of the nearest geomorphon
#'       boundary cell. Higher values indicate better alignment. Only included
#'       when method includes "boundary".}
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
                                  threshold_cells = 2) {

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

  # Validate threshold only when boundary distances are requested.
  if (method %in% c("boundary", "both")) {
    if (!is.numeric(threshold_cells) || length(threshold_cells) != 1L ||
        !is.finite(threshold_cells) || threshold_cells < 0 ||
        threshold_cells != floor(threshold_cells)) {
      stop("'threshold_cells' must be a finite, non-negative whole number",
           call. = FALSE)
    }

    if (terra::is.lonlat(geo_rast)) {
      warning(
        "Geographic CRS detected; boundary distances are in CRS map units. ",
        "Use a projected CRS for meaningful geomorphon distance calculations.",
        call. = FALSE
      )
    }
  }

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
      # Keep only geomorphon edge cells as non-NA targets, then measure the
      # distance from every cell to the nearest geomorphon boundary.
      geo_edge_dist_rast <- geo_edges
      terra::values(geo_edge_dist_rast)[geo_edge_vals != 1] <- NA

      geo_edge_dist <- terra::distance(geo_edge_dist_rast)

      # Sample those map-unit distances at categorical edge cells
      geo_edge_dist_vals <- terra::values(geo_edge_dist, mat = FALSE)
      cat_edge_dists <- geo_edge_dist_vals[cat_edge_vals == 1]
      cat_edge_dists <- cat_edge_dists[!is.na(cat_edge_dists)]

      if (length(cat_edge_dists) > 0) {
        # terra::distance() already returns map units; do not rescale again
        cell_size <- terra::res(geo_rast)[1]
        threshold_mapunits <- threshold_cells * cell_size

        result$edge_offset <- mean(cat_edge_dists)
        result$edge_concordance <- sum(cat_edge_dists <= threshold_mapunits) /
          length(cat_edge_dists)
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
