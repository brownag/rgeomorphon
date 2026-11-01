library(tinytest)
library(rgeomorphon)

if (requireNamespace("terra", quietly = TRUE)) {
  # Setup test data
  data("volcano", package = "datasets")
  dem <- terra::rast(volcano)
  terra::crs(dem) <- terra::crs("EPSG:2193")
  terra::ext(dem) <- c(1756968, 1757578, 5917000, 5917870)
  names(dem) <- "elevation"

  # ===========================================================================
  # PTRM Function Tests
  # ===========================================================================

  # Test 1: PTRM from DEM returns correct structure
  ptrm_dem <- geomorphon_ptrm(dem, mode = "dem", search = 3)

  expect_true(inherits(ptrm_dem, "list"), "PTRM should return a list")
  expect_true("ptrm_score" %in% names(ptrm_dem), "Result should have ptrm_score")
  expect_true("raw_score" %in% names(ptrm_dem), "Result should have raw_score")
  expect_true("proportions" %in% names(ptrm_dem), "Result should have proportions")
  expect_true("n_cells" %in% names(ptrm_dem), "Result should have n_cells")

  # Test 2: PTRM score is a numeric value
  expect_true(is.numeric(ptrm_dem$ptrm_score), "PTRM score should be numeric")
  expect_true(is.finite(ptrm_dem$ptrm_score), "PTRM score should be finite")

  # Test 3: Proportions sum to 1.0
  expect_equal(sum(ptrm_dem$proportions), 1.0, tolerance = 1e-6,
               info = "Proportions should sum to 1.0")

  # Test 4: n_cells is positive integer
  expect_true(ptrm_dem$n_cells > 0, "n_cells should be positive")
  expect_equal(ptrm_dem$n_cells, as.integer(ptrm_dem$n_cells),
               info = "n_cells should be an integer")

  # Test 5: mode="dem" and mode="class" give identical results
  geo <- geomorphons(dem, search = 3)
  ptrm_class <- geomorphon_ptrm(geo, mode = "class")

  expect_equal(ptrm_dem$ptrm_score, ptrm_class$ptrm_score,
               tolerance = 1e-10,
               info = "mode='dem' and mode='class' should give same score")
  expect_equal(ptrm_dem$proportions, ptrm_class$proportions,
               tolerance = 1e-10,
               info = "Proportions should be identical for both modes")

  # Test 6: PTRM with different search radius
  ptrm_search5 <- geomorphon_ptrm(dem, mode = "dem", search = 5)

  expect_true(inherits(ptrm_search5, "list"), "Should return list with search=5")
  # Different search radius may give different proportions
  expect_true(is.finite(ptrm_search5$ptrm_score),
              info = "Score with search=5 should be finite")

  # Test 7: Custom coefficients
  custom_coeffs <- c(
    flat = 1, peak = 1, ridge = 1, shoulder = 1, spur = 1,
    slope = 1, hollow = 1, footslope = 1, valley = 1, depression = 1
  )
  ptrm_custom <- geomorphon_ptrm(geo, mode = "class",
                                coefficients = custom_coeffs)

  expect_true(inherits(ptrm_custom, "list"), "Should work with custom coefficients")
  expect_true(is.finite(ptrm_custom$ptrm_score),
              info = "Custom coefficients should produce finite score")

  # Test 8: Error handling - non-SpatRaster input
  expect_error(geomorphon_ptrm(volcano, mode = "dem"),
               pattern = "SpatRaster",
               info = "Should error on non-SpatRaster input")

  # Test 9: Error handling - invalid mode
  expect_error(geomorphon_ptrm(dem, mode = "invalid"),
               info = "Should error on invalid mode")

  # Test 10: Error handling - empty raster (all NA)
  empty_rast <- dem
  terra::values(empty_rast) <- NA
  expect_error(geomorphon_ptrm(empty_rast, mode = "class"),
               pattern = "valid",
               info = "Should error on all-NA raster")

  # Test 10a: Known two-class proportions match the hand-computed formula
  two_class <- terra::rast(nrows = 2, ncols = 2,
                           xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::values(two_class) <- c(1, 1, 1, 6)  # flat: 3/4, slope: 1/4

  ptrm_two <- geomorphon_ptrm(two_class, mode = "class")
  expected_two <- (-38.02 + 0.75 * 25.12 + 0.25 * 28.95) / 69.96

  expect_equal(ptrm_two$ptrm_score, expected_two, tolerance = 1e-10,
               info = "Two-class score should match the hand-computed PTRM")
  expect_equal(ptrm_two$raw_score, expected_two * 69.96, tolerance = 1e-10,
               info = "Raw score should be the unnormalized PTRM numerator")
  expect_equal(unname(ptrm_two$proportions), c(0.75, 0.25),
               tolerance = 1e-10,
               info = "Two-class proportions should match the cell counts")
  expect_identical(names(ptrm_two$proportions), c("flat", "slope"),
                   info = "Proportions should use canonical system labels")

  # Test 10b: Single-class raster uses that class's coefficient, not -0.543
  single_ridge <- terra::rast(nrows = 2, ncols = 2,
                              xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::values(single_ridge) <- 3  # ridge

  ptrm_single <- geomorphon_ptrm(single_ridge, mode = "class")
  expect_equal(ptrm_single$ptrm_score, (-38.02 + 7.59) / 69.96,
               tolerance = 1e-10,
               info = "Single-class score should include that class's weight")
  expect_equal(ptrm_single$ptrm_score, -0.434962835906232,
               tolerance = 1e-12,
               info = "Single ridge score should match the literal expectation")
  expect_false(isTRUE(all.equal(ptrm_single$ptrm_score, -0.543453401943968)),
               info = "Single-class score must not collapse to the intercept-only value")

  # Test 10c: Reduced form systems return finite scores and system labels
  reduced_systems <- list(
    forms4 = c(1, 3, 6, 9),
    forms5 = c(1, 2, 3, 6, 9),
    forms6 = c(1, 3, 4, 6, 8, 9)
  )
  expected_labels <- list(
    forms4 = c("flat", "ridge", "slope", "valley"),
    forms5 = c("flat", "peak", "ridge", "slope", "valley"),
    forms6 = c("flat", "ridge", "shoulder", "slope", "footslope",
               "valley")
  )

  for (system in names(reduced_systems)) {
    codes <- reduced_systems[[system]]
    r <- terra::rast(nrows = 1, ncols = length(codes),
                     xmin = 0, xmax = length(codes), ymin = 0, ymax = 1)
    terra::values(r) <- codes

    ptrm_reduced <- geomorphon_ptrm(r, mode = "class", forms = system)

    expect_true(is.finite(ptrm_reduced$ptrm_score),
                info = paste(system, "score should be finite"))
    expect_identical(names(ptrm_reduced$proportions), expected_labels[[system]],
                     info = paste(system, "proportions should be named by system"))
    expect_equal(sum(ptrm_reduced$proportions), 1.0, tolerance = 1e-10,
                 info = paste(system, "proportions should sum to 1"))
  }

  # Test 10d: A ten-form custom vector reduces to the derived defaults
  paper_coefficients <- c(
    flat = 25.12, peak = 1.75, ridge = 7.59, shoulder = 7.31, spur = 9.02,
    slope = 28.95, hollow = 6.71, footslope = 7.63, valley = 9.61,
    depression = 3.55
  )

  for (system in c("forms4", "forms5", "forms6")) {
    codes <- reduced_systems[[system]]
    r <- terra::rast(nrows = 1, ncols = length(codes),
                     xmin = 0, xmax = length(codes), ymin = 0, ymax = 1)
    terra::values(r) <- codes

    ptrm_default <- geomorphon_ptrm(r, mode = "class", forms = system)
    ptrm_ten <- geomorphon_ptrm(r, mode = "class", forms = system,
                               coefficients = paper_coefficients)

    expect_equal(ptrm_ten$ptrm_score, ptrm_default$ptrm_score,
                 tolerance = 1e-10,
                 info = paste(system, "ten-form custom vector should match defaults"))
  }

  # Test 10e: Intercept and normalizer overrides change the score as expected
  ptrm_override <- geomorphon_ptrm(two_class, mode = "class",
                                   intercept = -10, normalizer = 10)
  expected_override <- (-10 + 0.75 * 25.12 + 0.25 * 28.95) / 10

  expect_equal(ptrm_override$ptrm_score, expected_override,
               tolerance = 1e-10,
               info = "Override score should use the supplied constants")
  expect_false(isTRUE(all.equal(ptrm_override$ptrm_score,
                                ptrm_two$ptrm_score)),
               info = "Overrides should change the score")

  # Test 10f: Invalid model parameters and coefficients error
  expect_error(geomorphon_ptrm(two_class, mode = "class", normalizer = 0),
               pattern = "non-zero",
               info = "normalizer = 0 should error")
  expect_error(
    geomorphon_ptrm(two_class, mode = "class", forms = "forms4",
                    coefficients = c(flat = 1, ridge = 2, slope = 3,
                                     peak = 4)),
    pattern = "must be named for",
    info = "Wrong coefficient names should error"
  )
  expect_error(
    geomorphon_ptrm(two_class, mode = "class",
                    coefficients = c(flat = 1, peak = 2, ridge = 3)),
    pattern = "must be named for",
    info = "Wrong coefficient length should error"
  )

  # Test 10g: Codes outside the requested system error informatively
  out_of_system <- terra::rast(nrows = 1, ncols = 2,
                               xmin = 0, xmax = 2, ymin = 0, ymax = 1)
  terra::values(out_of_system) <- c(1, 2)

  expect_error(
    geomorphon_ptrm(out_of_system, mode = "class", forms = "forms4"),
    pattern = "Invalid geomorphon codes for forms4: 2",
    info = "Out-of-system codes should name the invalid codes and system"
  )

  # ===========================================================================
  # Concordance Function Tests
  # ===========================================================================

  # Create categorical test raster
  cat_map <- dem
  terra::values(cat_map) <- as.numeric(cut(terra::values(dem), breaks = 5, labels = FALSE))

  # Test 11: Concordance returns correct structure (both method)
  conc_both <- geomorphon_concordance(geo, cat_map, method = "both")

  expect_true(inherits(conc_both, "list"), "Should return a list")
  expect_true("nmi" %in% names(conc_both), "Should have nmi with method='both'")
  expect_true("edge_offset" %in% names(conc_both),
              info = "Should have edge_offset with method='both'")
  expect_true("edge_concordance" %in% names(conc_both),
              info = "Should have edge_concordance with method='both'")

  # Test 12: NMI is in valid range [0, 1]
  expect_true(conc_both$nmi >= 0, "NMI should be >= 0")
  expect_true(conc_both$nmi <= 1, "NMI should be <= 1")

  # Test 13: Edge concordance is in valid range [0, 1]
  expect_true(conc_both$edge_concordance >= 0, "Edge concordance should be >= 0")
  expect_true(conc_both$edge_concordance <= 1, "Edge concordance should be <= 1")

  # Test 14: Perfect concordance (identical rasters)
  conc_perfect <- geomorphon_concordance(geo, geo, method = "nmi")

  expect_true(conc_perfect$nmi >= 0.99, "Perfect concordance NMI should be ~1.0")
  expect_true("contingency_table" %in% names(conc_perfect),
              info = "Should have contingency_table with NMI method")

  # Test 15: NMI only method
  conc_nmi <- geomorphon_concordance(geo, cat_map, method = "nmi")

  expect_false("edge_offset" %in% names(conc_nmi),
               info = "Should not have edge metrics with method='nmi'")
  expect_true("nmi" %in% names(conc_nmi), "Should have nmi with method='nmi'")

  # Test 16: Boundary only method
  conc_boundary <- geomorphon_concordance(geo, cat_map, method = "boundary")

  expect_false("nmi" %in% names(conc_boundary),
               info = "Should not have nmi with method='boundary'")
  expect_true("edge_offset" %in% names(conc_boundary),
              info = "Should have edge_offset with method='boundary'")

  # Test 17: Different threshold_cells values
  conc_thresh1 <- geomorphon_concordance(geo, cat_map, method = "boundary",
                                        threshold_cells = 1)
  conc_thresh3 <- geomorphon_concordance(geo, cat_map, method = "boundary",
                                        threshold_cells = 3)

  expect_true(conc_thresh3$edge_concordance >= conc_thresh1$edge_concordance,
              info = "Larger threshold should have >= concordance")

  # Threshold must be a finite, non-negative whole number for boundary methods.
  for (boundary_method in c("boundary", "both")) {
    for (invalid_threshold in list(NULL, numeric(0), c(1, 2), NA_real_, NaN,
                                   Inf, -Inf, -1, 1.5, "2", TRUE)) {
      expect_error(
        geomorphon_concordance(geo, cat_map, method = boundary_method,
                               threshold_cells = invalid_threshold),
        pattern = "threshold_cells",
        info = paste("Invalid threshold should error for", boundary_method)
      )
    }
    for (valid_threshold in list(0, 0L, 2, 2L)) {
      expect_true(is.list(geomorphon_concordance(
        geo, cat_map, method = boundary_method,
        threshold_cells = valid_threshold
      )), info = "Zero and positive whole-number doubles and integers are valid")
    }
  }
  expect_true(is.list(geomorphon_concordance(
    geo, cat_map, method = "nmi", threshold_cells = -1
  )), info = "NMI-only method does not validate boundary threshold")

  # Test 18: Error handling - non-SpatRaster inputs
  expect_error(geomorphon_concordance(volcano, cat_map),
               pattern = "SpatRaster",
               info = "Should error on non-SpatRaster geo_rast")

  expect_error(geomorphon_concordance(geo, volcano),
               pattern = "SpatRaster",
               info = "Should error on non-SpatRaster cat_rast")

  # Test 19: Error handling - mismatched extents
  dem_small <- terra::crop(dem, terra::ext(dem) / 2)
  geo_small <- geomorphons(dem_small, search = 3, simplify = TRUE)

  expect_error(geomorphon_concordance(geo, geo_small),
               pattern = "matching extents",
               info = "Should error on mismatched extents")

  # Test 20: Error handling - mismatched resolutions (also changes extent)
  geo_coarse <- terra::aggregate(geo, fact = 2)

  expect_error(geomorphon_concordance(geo, geo_coarse),
               pattern = "matching",
               info = "Should error on mismatched resolutions/extents")

  # Test 21: Error handling - all NA rasters
  empty_cat <- cat_map
  terra::values(empty_cat) <- NA

  expect_error(geomorphon_concordance(geo, empty_cat, method = "nmi"),
               pattern = "valid",
               info = "Should error on all-NA categorical raster")

  # Test 22: Contingency table structure
  expect_true(inherits(conc_both$contingency_table, "table"),
              info = "Contingency table should be a table object")
  expect_true(all(dim(conc_both$contingency_table) > 0),
              info = "Contingency table should have non-zero dimensions")

  # Test 23: n_cells reported correctly
  expect_true(conc_both$n_cells > 0, "n_cells should be positive")
  expect_equal(conc_both$n_cells, sum(conc_both$contingency_table),
               info = "n_cells should equal sum of contingency table")

  # Test 24: Edge count consistency
  expect_true(conc_both$n_edges_geo >= 0, "n_edges_geo should be non-negative")
  expect_true(conc_both$n_edges_cat >= 0, "n_edges_cat should be non-negative")

  # Test 25: Edge offset uses map units and is not double-scaled
  make_concordance_raster <- function(mat, res = 1) {
    r <- terra::rast(
      nrows = nrow(mat), ncols = ncol(mat),
      xmin = 0, xmax = ncol(mat) * res,
      ymin = 0, ymax = nrow(mat) * res,
      crs = "EPSG:3857"
    )
    terra::values(r) <- as.vector(t(mat))
    r
  }

  geo_a <- matrix(1, nrow = 7, ncol = 7)
  geo_a[, 4:7] <- 2
  cat_a <- matrix(1, nrow = 7, ncol = 7)
  cat_a[, 6:7] <- 2

  conc_a <- geomorphon_concordance(
    make_concordance_raster(geo_a, res = 30),
    make_concordance_raster(cat_a, res = 30),
    method = "boundary"
  )

  # Categorical edges are 1 and 2 cells from the nearest geomorphon edge
  expect_equal(conc_a$edge_offset, 45, tolerance = 1e-8,
               info = "edge_offset should be in map units, without double scaling")

  # Test 26: Edge offset is a mean, not a median, on asymmetric offsets
  geo_b <- matrix(1, nrow = 9, ncol = 9)
  geo_b[, 5:9] <- 2
  cat_b <- matrix(1, nrow = 9, ncol = 9)
  cat_b[2, 6:9] <- 2
  cat_b[3:8, 8:9] <- 2

  conc_b <- geomorphon_concordance(
    make_concordance_raster(geo_b),
    make_concordance_raster(cat_b),
    method = "boundary",
    threshold_cells = 2
  )

  # Distance counts 0:2, 1:2, 2:7, 3:7; mean = 37/18, median = 2
  expect_equal(conc_b$edge_offset, 37 / 18, tolerance = 1e-8,
               info = "edge_offset should be the mean, not the median")
  expect_false(isTRUE(all.equal(conc_b$edge_offset, 2)),
               info = "Asymmetric offsets should distinguish mean from median")

  # Test 27: Edge concordance matches a manually computed proportion
  expect_equal(conc_b$edge_concordance, 11 / 18, tolerance = 1e-8,
               info = "edge_concordance should be the proportion within threshold")
  expect_true(conc_b$edge_concordance >= 0 && conc_b$edge_concordance <= 1,
              info = "edge_concordance should lie in [0, 1]")

  # Test 28: Fallback edge detection matches the C++ path
  # Rectangular dimensions expose row/column transposition mistakes.
  edge_mat <- matrix(1, nrow = 7, ncol = 9)
  edge_mat[3, 2] <- NA
  edge_mat[4, 6] <- 2
  edge_mat[5, 7] <- NA
  edge_rast <- make_concordance_raster(edge_mat)

  edges_cpp <- .detect_edges(edge_rast)
  edges_fallback <- suppressWarnings(
    .detect_edges(edge_rast, force_fallback = TRUE)
  )

  expect_equal(
    terra::values(edges_fallback, mat = FALSE),
    terra::values(edges_cpp, mat = FALSE),
    info = "Fallback and C++ edge masks should be identical"
  )

  # NA is not a class: (4,2) has an NA neighbor and only class-1 valid neighbors.
  same_class_neighborhood <- edge_mat[3:5, 1:3]
  expect_true(anyNA(same_class_neighborhood))
  expect_true(all(same_class_neighborhood[!is.na(same_class_neighborhood)] ==
                    edge_mat[4, 2]))
  # Use [row, column], not flattened R matrix indices (terra cells are row-major).
  edge_values <- terra::as.matrix(edges_cpp, wide = TRUE)
  expect_equal(edge_values[4, 2], 0,
               info = "NA alone does not make a valid cell an edge")
  expect_equal(edge_values[4, 6], 1,
               info = "Differing valid neighbors remain edges even next to NA")
  expect_equal(edge_values[3, 5], 1,
               info = "A diagonal differing neighbor counts with 8-connectivity")

  # Only the class-2 cell and its eight neighbors are edges, excluding NA centers.
  expected_edges <- matrix(0, nrow = 7, ncol = 9)
  expected_edges[3:5, 5:7] <- 1
  expected_edges[is.na(edge_mat)] <- NA
  expected_edges[c(1, 7), ] <- NA
  expected_edges[, c(1, 9)] <- NA
  expect_equal(edge_values, expected_edges,
               info = "Complete edge mask includes NA centers and all four borders")
  expect_equal(terra::as.matrix(edges_fallback, wide = TRUE), expected_edges,
               info = "Fallback also matches the independently specified mask")

  # Geographic CRS warning is limited to methods that request boundary metrics.
  geo_ll <- make_concordance_raster(geo_a)
  terra::crs(geo_ll) <- "EPSG:4326"
  cat_ll <- geo_ll
  terra::values(cat_ll) <- as.vector(t(cat_a))
  expect_warning(
    geomorphon_concordance(geo_ll, cat_ll, method = "boundary"),
    pattern = "Geographic CRS"
  )
  expect_warning(
    geomorphon_concordance(geo_ll, cat_ll, method = "both"),
    pattern = "Geographic CRS"
  )
  expect_silent(
    geomorphon_concordance(geo_ll, cat_ll, method = "nmi")
  )
  expect_silent(
    geomorphon_concordance(make_concordance_raster(geo_a),
                           make_concordance_raster(cat_a), method = "both")
  )

  # ===========================================================================
  # Integration Tests
  # ===========================================================================

  # Test 29: Full workflow with PTRM and concordance
  e <- terra::ext(dem)
  dem_subset <- terra::crop(dem, terra::ext(e[1] + 100, e[1] + 400, e[3] + 100, e[3] + 400))
  geo_subset <- geomorphons(dem_subset, search = 5)
  ptrm_subset <- geomorphon_ptrm(geo_subset, mode = "class")
  cat_subset <- dem_subset
  terra::values(cat_subset) <- as.numeric(cut(terra::values(dem_subset), breaks = 3,
                                              labels = FALSE))
  conc_subset <- geomorphon_concordance(geo_subset, cat_subset)

  expect_true(inherits(ptrm_subset, "list"), "PTRM should work on subset")
  expect_true(is.finite(ptrm_subset$ptrm_score),
              info = "PTRM score should be finite")
  expect_true(inherits(conc_subset, "list"), "Concordance should work on subset")
  expect_true(conc_subset$nmi >= 0 && conc_subset$nmi <= 1,
              info = "Concordance NMI should be valid")
}
