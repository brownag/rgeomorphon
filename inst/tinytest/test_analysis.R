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

  frac_geo <- geo
  terra::values(frac_geo)[1:2] <- c(1.5, 6.2)
  expect_error(geomorphon_ptrm(frac_geo, mode = "class"),
               pattern = "Invalid geomorphon codes")

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
  expect_true(all(c("homogeneity", "completeness") %in% names(conc_both)),
              info = "NMI components should be returned with method='both'")
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
  expect_true(all(c("homogeneity", "completeness") %in% names(conc_perfect)),
              info = "NMI components should be returned with method='nmi'")

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
  expect_false(any(c("homogeneity", "completeness") %in%
                   names(conc_boundary)),
               info = "Boundary-only output should omit NMI components")

  # Spatial-null opt-out preserves observed scores and zeroes null summaries.
  for (null_method in c("both", "nmi", "boundary")) {
    enabled <- geomorphon_concordance(geo, cat_map, method = null_method)
    disabled <- geomorphon_concordance(geo, cat_map, method = null_method,
                                       spatial_null = FALSE)
    observed_fields <- intersect(c("nmi", "homogeneity", "completeness",
                                   "edge_offset", "edge_concordance", "edge_f1",
                                   "edge_mhd", "edge_chamfer"), names(enabled))
    expect_equal(disabled[observed_fields], enabled[observed_fields])
    null_fields <- intersect(c("mean_null_nmi", "delta_nmi",
                               "mean_null_edge_f1", "delta_edge_f1"),
                             names(disabled))
    expect_true(all(vapply(disabled[null_fields], is.na, logical(1))))
    count_fields <- intersect(c("nmi_null_shifts", "edge_f1_null_shifts"),
                              names(disabled))
    expect_true(all(vapply(disabled[count_fields], identical, logical(1), 0L)))
    expect_identical(disabled$null_shifts_evaluated, 0L)
  }
  for (invalid_spatial_null in list(NA, "yes", c(TRUE, FALSE))) {
    expect_error(geomorphon_concordance(geo, cat_map,
                                        spatial_null = invalid_spatial_null),
                 pattern = "'spatial_null' must be TRUE or FALSE")
  }

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

  # Contingency table [[4, 2], [0, 2]] (geomorphon rows, category columns):
  # H(cat)=1, H(geo)=H(1/4, 3/4), H(geo,cat)=1.5, so
  # I=H(geo)-0.5. The unequal marginal entropies give distinct components.
  asymmetric_geo <- matrix(c(1, 1, 1, 1, 1, 1, 2, 2), nrow = 1)
  asymmetric_cat <- matrix(c(1, 1, 1, 1, 2, 2, 2, 2), nrow = 1)
  asymmetric_geo_entropy <- -(0.75 * log2(0.75) + 0.25 * log2(0.25))
  asymmetric_mi <- asymmetric_geo_entropy - 0.5
  expected_homogeneity <- asymmetric_mi
  expected_completeness <- asymmetric_mi / asymmetric_geo_entropy
  expected_nmi <- 2 * asymmetric_mi / (1 + asymmetric_geo_entropy)
  asymmetric_result <- geomorphon_concordance(
    make_concordance_raster(asymmetric_geo),
    make_concordance_raster(asymmetric_cat), method = "nmi"
  )
  expect_equal(asymmetric_result$homogeneity, expected_homogeneity,
               tolerance = 1e-10)
  expect_equal(asymmetric_result$completeness, expected_completeness,
               tolerance = 1e-10)
  expect_equal(asymmetric_result$nmi, expected_nmi, tolerance = 1e-10)
  expect_true(abs(asymmetric_result$homogeneity -
                    asymmetric_result$completeness) > 1e-3,
              info = "Unequal marginal entropies should yield distinct components")
  expect_equal(2 * asymmetric_result$homogeneity *
                 asymmetric_result$completeness /
                 (asymmetric_result$homogeneity +
                    asymmetric_result$completeness),
               asymmetric_result$nmi, tolerance = 1e-10)

  asymmetric_swapped <- geomorphon_concordance(
    make_concordance_raster(asymmetric_cat),
    make_concordance_raster(asymmetric_geo), method = "nmi"
  )
  expect_equal(asymmetric_result$homogeneity,
               asymmetric_swapped$completeness, tolerance = 1e-10)
  expect_equal(asymmetric_result$completeness,
               asymmetric_swapped$homogeneity, tolerance = 1e-10)
  expect_equal(asymmetric_result$nmi, asymmetric_swapped$nmi,
               tolerance = 1e-10)

  # Constant reference/target maps exercise zero-entropy conventions and
  # paired-NA filtering without changing existing NMI behavior.
  constant_map <- matrix(1, nrow = 1, ncol = 8)
  varying_map <- matrix(c(1, 1, 1, 1, 2, 2, 2, 2), nrow = 1)
  constant_cat <- geomorphon_concordance(
    make_concordance_raster(varying_map), make_concordance_raster(constant_map),
    method = "nmi"
  )
  expect_identical(constant_cat$homogeneity, NA_real_)
  expect_identical(constant_cat$completeness, 0)
  expect_identical(constant_cat$nmi, 0)
  constant_geo <- geomorphon_concordance(
    make_concordance_raster(constant_map), make_concordance_raster(varying_map),
    method = "nmi"
  )
  expect_identical(constant_geo$homogeneity, 0)
  expect_identical(constant_geo$completeness, NA_real_)
  expect_identical(constant_geo$nmi, 0)
  both_constant <- geomorphon_concordance(
    make_concordance_raster(constant_map), make_concordance_raster(constant_map),
    method = "nmi"
  )
  expect_identical(both_constant$homogeneity, NA_real_)
  expect_identical(both_constant$completeness, NA_real_)
  expect_identical(both_constant$nmi, 0)

  na_geo <- asymmetric_geo
  na_cat <- asymmetric_cat
  na_geo[1, 1] <- NA_real_
  na_cat[1, 2] <- NA_real_
  na_result <- geomorphon_concordance(
    make_concordance_raster(na_geo), make_concordance_raster(na_cat),
    method = "nmi"
  )
  na_valid <- !is.na(na_geo) & !is.na(na_cat)
  na_geo_valid <- na_geo[na_valid]
  na_cat_valid <- na_cat[na_valid]
  na_entropy <- function(x) {
    p <- as.numeric(table(x)) / length(x)
    -sum(p * log2(p))
  }
  na_mi <- na_entropy(na_geo_valid) + na_entropy(na_cat_valid) -
    na_entropy(paste(na_geo_valid, na_cat_valid, sep = "_"))
  expect_equal(na_result$homogeneity, na_mi / na_entropy(na_cat_valid),
               tolerance = 1e-10)
  expect_equal(na_result$completeness, na_mi / na_entropy(na_geo_valid),
               tolerance = 1e-10)

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

  edges_cpp <- rgeomorphon:::.detect_edges(edge_rast)
  edges_fallback <- suppressWarnings(
    rgeomorphon:::.detect_edges(edge_rast, force_fallback = TRUE)
  )

  expect_equal(
    terra::values(edges_fallback, mat = FALSE),
    terra::values(edges_cpp, mat = FALSE),
    info = "Fallback and C++ edge masks should be identical"
  )

  nan_edge_mat <- matrix(1, 5, 5)
  nan_edge_mat[3, 3] <- NaN
  nan_edge_rast <- make_concordance_raster(nan_edge_mat)
  nan_cpp <- rgeomorphon:::focal_edge_cpp(nan_edge_mat)
  na_edge_mat <- nan_edge_mat
  na_edge_mat[3, 3] <- NA_real_
  na_cpp <- rgeomorphon:::focal_edge_cpp(na_edge_mat)
  expect_identical(nan_cpp[3, 3], NA_real_)
  expect_equal(nan_cpp[2:4, 2:4], na_cpp[2:4, 2:4])
  nan_detected <- rgeomorphon:::.detect_edges(nan_edge_rast)
  detected_values <- terra::as.matrix(nan_detected, wide = TRUE)
  expect_equal(unname(detected_values), unname(na_cpp))

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

  # Task 3.1.1: tiny coordinate fixtures and independent nearest-center oracle.
  # Matrices use [row, column]; coordinates are cell centers in EPSG:3857.
  boundary_coordinates <- function(mask, xres = 1, yres = 1) {
    ij <- which(mask, arr.ind = TRUE)
    if (nrow(ij) == 0L) return(matrix(numeric(), ncol = 2L))
    cbind(x = (ij[, 2] - 0.5) * xres,
          y = (nrow(mask) - ij[, 1] + 0.5) * yres)
  }
  brute_nearest <- function(from, to) {
    if (!nrow(from) || !nrow(to)) return(numeric())
    vapply(seq_len(nrow(from)), function(i) {
      min(sqrt(rowSums((to - matrix(from[i, ], nrow(to), 2L,
                                    byrow = TRUE))^2)))
    }, numeric(1))
  }
  reference_mask <- function(mat, valid = !is.na(mat)) {
    mask <- matrix(FALSE, nrow(mat), ncol(mat))
    for (i in 2:(nrow(mat) - 1L)) {
      for (j in 2:(ncol(mat) - 1L)) {
        if (!valid[i, j]) next
        neighbors <- mat[(i - 1L):(i + 1L), (j - 1L):(j + 1L)]
        neighbor_valid <- valid[(i - 1L):(i + 1L), (j - 1L):(j + 1L)]
        mask[i, j] <- any(neighbor_valid & neighbors != mat[i, j])
      }
    }
    mask
  }
  oracle_scores <- function(geo_mat, cat_mat, xres = 1, yres = 1,
                            threshold_cells = 1) {
    shared <- !is.na(geo_mat) & !is.na(cat_mat)
    geo_xy <- boundary_coordinates(reference_mask(geo_mat, shared), xres, yres)
    cat_xy <- boundary_coordinates(reference_mask(cat_mat, shared), xres, yres)
    if (!nrow(geo_xy) || !nrow(cat_xy)) return(c(NA_real_, NA_real_, NA_real_))
    g2c <- brute_nearest(geo_xy, cat_xy)
    c2g <- brute_nearest(cat_xy, geo_xy)
    p <- mean(c2g <= threshold_cells * xres)
    r <- mean(g2c <= threshold_cells * xres)
    c(if (p + r == 0) 0 else 2 * p * r / (p + r),
      max(mean(g2c), mean(c2g)), (mean(g2c) + mean(c2g)) / 2)
  }

  # Manual checks: coincident centers are 0; adjacent centers are 2 by 3.
  point_a <- matrix(c(0.5, 1.5), nrow = 1)
  point_b <- matrix(c(2.5, 4.5), nrow = 1)
  expect_equal(brute_nearest(point_a, point_a), 0)
  expect_equal(brute_nearest(point_a, point_b), 3.605551275,
               tolerance = 1e-8)
  expect_true(brute_nearest(point_a, point_b) > 3)
  expect_true(brute_nearest(point_a, point_b) <= sqrt(13))
  expect_false(brute_nearest(point_a, point_b) <= 3)
  expect_equal(brute_nearest(point_a, point_a) <= 0, TRUE)
  expect_equal(brute_nearest(point_a, point_b) <= 0, FALSE)

  # Shared-domain fixture expectation is specified from class contrasts, not
  # production edge detection; differing NA masks exclude those cells first.
  fixture_geo <- matrix(1, 7, 9)
  fixture_geo[, 5:9] <- 2
  fixture_cat <- fixture_geo
  fixture_cat[, 6:9] <- 2
  expected_geo_fixture <- reference_mask(fixture_geo)
  expected_cat_fixture <- reference_mask(fixture_cat)
  expect_equal(nrow(which(expected_geo_fixture, arr.ind = TRUE)), 10)
  expect_equal(nrow(which(expected_cat_fixture, arr.ind = TRUE)), 10)
  shared_valid <- !is.na(fixture_geo) & !is.na(fixture_cat)
  expected_shared_geo <- reference_mask(fixture_geo, shared_valid)
  expected_shared_cat <- reference_mask(fixture_cat, shared_valid)
  expect_equal(expected_shared_geo, expected_geo_fixture)
  expect_equal(expected_shared_cat, expected_cat_fixture)
  geo_fixture_rast <- make_concordance_raster(fixture_geo)
  cat_fixture_rast <- make_concordance_raster(fixture_cat)
  fixture_legacy <- geomorphon_concordance(
    geo_fixture_rast, cat_fixture_rast, method = "boundary", threshold_cells = 1
  )
  expect_equal(fixture_legacy$n_edges_geo, 10)
  expect_equal(fixture_legacy$n_edges_cat, 10)
  expect_equal(fixture_legacy$edge_offset, 0)
  expect_equal(fixture_legacy$edge_concordance, 1)

  # With asymmetric NA masks, legacy edge counts remain per-input; pin outputs.
  geo_asym <- fixture_geo
  cat_asym <- fixture_cat
  geo_asym[4, 3] <- NA
  cat_asym[4, 7] <- NA
  asym_result <- geomorphon_concordance(
    make_concordance_raster(geo_asym), make_concordance_raster(cat_asym),
    method = "boundary", threshold_cells = 1
  )
  expect_equal(asym_result$n_edges_geo, 10)
  expect_equal(asym_result$n_edges_cat, 10)
  expect_equal(asym_result$edge_offset, 0)
  expect_equal(asym_result$edge_concordance, 1)

  asym_oracle <- oracle_scores(geo_asym, cat_asym, threshold_cells = 1)
  expect_equal(unname(unlist(asym_result[c("edge_f1", "edge_mhd", "edge_chamfer")])),
               asym_oracle, tolerance = 1e-8)

  # Symmetric shared-domain summaries match the independent coordinate oracle.
  fixture_geo_xy <- boundary_coordinates(expected_shared_geo)
  fixture_cat_xy <- boundary_coordinates(expected_shared_cat)
  fixture_geo_to_cat <- brute_nearest(fixture_geo_xy, fixture_cat_xy)
  fixture_cat_to_geo <- brute_nearest(fixture_cat_xy, fixture_geo_xy)
  fixture_p <- mean(fixture_cat_to_geo <= 1)
  fixture_r <- mean(fixture_geo_to_cat <= 1)
  expect_equal(fixture_legacy$edge_f1, 2 * fixture_p * fixture_r /
                 (fixture_p + fixture_r), tolerance = 1e-8)
  expect_equal(fixture_legacy$edge_mhd,
               max(mean(fixture_geo_to_cat), mean(fixture_cat_to_geo)),
               tolerance = 1e-8)
  expect_equal(fixture_legacy$edge_chamfer,
               mean(c(mean(fixture_geo_to_cat), mean(fixture_cat_to_geo))),
               tolerance = 1e-8)

  # Rectangular cells use x/y center spacing in distance, but threshold uses x.
  rectangular_geo <- fixture_geo
  rectangular_cat <- fixture_cat
  rect_geo_r <- terra::rast(nrows = 7, ncols = 9, xmin = 0, xmax = 18,
                            ymin = 0, ymax = 21, crs = "EPSG:3857")
  rect_cat_r <- terra::rast(nrows = 7, ncols = 9, xmin = 0, xmax = 18,
                            ymin = 0, ymax = 21, crs = "EPSG:3857")
  terra::values(rect_geo_r) <- as.vector(t(rectangular_geo))
  terra::values(rect_cat_r) <- as.vector(t(rectangular_cat))
  rect_result <- geomorphon_concordance(rect_geo_r, rect_cat_r,
                                        method = "boundary", threshold_cells = 2)
  expect_equal(unname(unlist(rect_result[c("edge_f1", "edge_mhd", "edge_chamfer")])),
               oracle_scores(rectangular_geo, rectangular_cat, 2, 3, 2),
               tolerance = 1e-8)

  asym_swapped <- geomorphon_concordance(
    make_concordance_raster(cat_asym), make_concordance_raster(geo_asym),
    method = "boundary", threshold_cells = 1
  )
  expect_equal(asym_result$edge_f1, asym_swapped$edge_f1, tolerance = 1e-8)
  expect_equal(asym_result$edge_mhd, asym_swapped$edge_mhd, tolerance = 1e-8)
  expect_equal(asym_result$edge_chamfer, asym_swapped$edge_chamfer,
               tolerance = 1e-8)
  expect_equal(asym_result$n_edges_geo, 10)
  expect_equal(asym_result$n_edges_cat, 10)
  expect_equal(asym_result$edge_offset, 0)
  expect_equal(asym_result$edge_concordance, 1)

  # New fields are boundary-only outputs; NMI retains its established shape.
  expect_false(any(c("edge_f1", "edge_mhd", "edge_chamfer") %in%
                   names(geomorphon_concordance(
                     geo_fixture_rast, cat_fixture_rast, method = "nmi"))))
  nmi_null_result <- geomorphon_concordance(
    geo_fixture_rast, cat_fixture_rast, method = "nmi"
  )
  expect_true(all(c("mean_null_nmi", "delta_nmi", "nmi_null_shifts",
                    "null_shifts_evaluated") %in% names(nmi_null_result)))
  expect_false("mean_null_edge_f1" %in% names(nmi_null_result))
  expect_equal(nmi_null_result$null_shifts_evaluated, 62L)
  # The capped shift schedule is independently specified by flattened ranks.
  large_nr <- 13L
  large_nc <- 13L
  large_n <- large_nr * large_nc
  # Sparse perpendicular boundaries keep threshold-one scores below saturation.
  large_geo <- large_cat <- matrix(1, large_nr, large_nc)
  large_geo[, 7:13] <- 2
  large_cat[7:13, ] <- 2
  large_null <- rgeomorphon:::spatial_null_cpp(
    as.vector(t(large_geo)), as.vector(t(large_cat)), large_nr, large_nc,
    1, FALSE, TRUE
  )
  expect_equal(large_null$evaluated, 128L)
  direct_shift_f1 <- function(dr, dc, threshold) {
    shifted <- large_cat
    for (r in seq_len(large_nr)) for (cc in seq_len(large_nc)) {
      shifted[r, cc] <- large_cat[(r - dr - 1L) %% large_nr + 1L,
                                  (cc - dc - 1L) %% large_nc + 1L]
    }
    valid <- matrix(TRUE, large_nr, large_nc)
    ge <- reference_mask(large_geo, valid)
    ce <- reference_mask(shifted, valid)
    hit_rate <- function(from, to) {
      ij_from <- which(from, arr.ind = TRUE)
      ij_to <- which(to, arr.ind = TRUE)
      if (!nrow(ij_from) || !nrow(ij_to)) return(NA_real_)
      mean(vapply(seq_len(nrow(ij_from)), function(k) {
        any(vapply(seq_len(nrow(ij_to)), function(j) {
          rr <- abs(ij_from[k, 1] - ij_to[j, 1])
          cc <- abs(ij_from[k, 2] - ij_to[j, 2])
          min(rr, large_nr - rr)^2 + min(cc, large_nc - cc)^2 <= threshold^2
        }, logical(1)))
      }, logical(1)))
    }
    p <- hit_rate(ce, ge); r <- hit_rate(ge, ce)
    if (is.na(p) || is.na(r)) return(NA_real_)
    if (p + r == 0) 0 else 2 * p * r / (p + r)
  }
  direct_scores <- vapply(0:127, function(i) {
    rank <- floor((i + 0.5) * (large_n - 1) / 128)
    q <- rank + 1
    direct_shift_f1(q %/% large_nc, q %% large_nc, 1)
  }, numeric(1))
  expect_equal(large_null$edge_mean, mean(direct_scores), tolerance = 1e-12)
  expect_true(large_null$edge_mean < 1)

  # Independently verify the capped NMI null schedule on a rectangular raster
  # with distinct geo/category missing-value masks.
  null_nr <- 5L
  null_nc <- 7L
  null_n <- null_nr * null_nc
  null_geo <- matrix(rep(c(1, 1, 2, 3), length.out = null_n), null_nr, null_nc)
  null_cat <- matrix(rep(c(1, 2, 2, 3, 1), length.out = null_n), null_nr, null_nc)
  null_geo[c(2, 6, 19)] <- NA_real_
  null_cat[c(4, 11, 28)] <- NA_real_
  null_cpp <- function() rgeomorphon:::spatial_null_cpp(
    as.vector(t(null_geo)), as.vector(t(null_cat)), null_nr, null_nc,
    0, TRUE, FALSE
  )
  null_entropy <- function(x) {
    p <- as.numeric(table(x)) / length(x)
    -sum(p * log2(p))
  }
  null_shifts <- min(128L, null_n - 1L)
  direct_nmi <- vapply(seq_len(null_shifts) - 1L, function(i) {
    rank <- floor((i + 0.5) * (null_n - 1) / null_shifts)
    q <- rank + 1L
    dr <- q %/% null_nc
    dc <- q %% null_nc
    shifted_cat <- matrix(NA_real_, null_nr, null_nc)
    for (r in seq_len(null_nr)) for (cc in seq_len(null_nc)) {
      shifted_cat[r, cc] <- null_cat[(r - dr - 1L) %% null_nr + 1L,
                                      (cc - dc - 1L) %% null_nc + 1L]
    }
    valid <- !is.na(null_geo) & !is.na(shifted_cat)
    if (!any(valid)) return(NA_real_)
    gv <- null_geo[valid]
    cv <- shifted_cat[valid]
    hg <- null_entropy(gv)
    hc <- null_entropy(cv)
    if (hg + hc == 0) return(0)
    hj <- null_entropy(paste(gv, cv, sep = "_"))
    max(0, min(1, 2 * max(0, hg + hc - hj) / (hg + hc)))
  }, numeric(1))
  direct_eligible <- direct_nmi[!is.na(direct_nmi)]
  null_first <- null_cpp()
  null_repeat <- null_cpp()
  expect_equal(null_first$nmi_mean, mean(direct_eligible), tolerance = 1e-12)
  expect_equal(null_first$nmi_eligible, length(direct_eligible))
  expect_identical(null_repeat$nmi_mean, null_first$nmi_mean)
  expect_identical(null_repeat$nmi_eligible, null_first$nmi_eligible)

  expect_equal(nmi_null_result$delta_nmi,
               nmi_null_result$nmi - nmi_null_result$mean_null_nmi,
               tolerance = 1e-12)
  expect_true(nmi_null_result$nmi_null_shifts <=
              nmi_null_result$null_shifts_evaluated)
  boundary_null_result <- geomorphon_concordance(
    geo_fixture_rast, cat_fixture_rast, method = "boundary",
    threshold_cells = 1
  )
  expect_true(all(c("mean_null_edge_f1", "delta_edge_f1",
                    "edge_f1_null_shifts") %in% names(boundary_null_result)))
  expect_false("mean_null_nmi" %in% names(boundary_null_result))
  expect_equal(boundary_null_result$delta_edge_f1,
               boundary_null_result$edge_f1 - boundary_null_result$mean_null_edge_f1,
               tolerance = 1e-12)
  expect_true(all(c("edge_f1", "edge_mhd", "edge_chamfer") %in%
                   names(geomorphon_concordance(
                     geo_fixture_rast, cat_fixture_rast, method = "both"))))

  # Empty shared boundaries produce typed missing scores, but retain legacy
  # no-valid-pairs behavior for methods that calculate NMI.
  empty_shared <- fixture_geo
  empty_shared[] <- NA_real_
  empty_boundary_result <- geomorphon_concordance(
    make_concordance_raster(empty_shared),
    make_concordance_raster(empty_shared), method = "boundary"
  )
  expect_identical(empty_boundary_result$edge_f1, NA_real_)
  expect_identical(empty_boundary_result$edge_mhd, NA_real_)
  expect_identical(empty_boundary_result$edge_chamfer, NA_real_)
  constant_fixture <- matrix(1, 7, 9)
  one_empty_result <- geomorphon_concordance(
    geo_fixture_rast, make_concordance_raster(constant_fixture),
    method = "boundary"
  )
  expect_identical(one_empty_result$edge_f1, NA_real_)
  expect_identical(one_empty_result$edge_mhd, NA_real_)
  expect_identical(one_empty_result$edge_chamfer, NA_real_)
  expect_error(geomorphon_concordance(
    make_concordance_raster(empty_shared), make_concordance_raster(empty_shared),
    method = "both"), pattern = "No overlapping valid cells")

  # Zero-match F1 and positive-distance summaries on a displaced boundary.
  shifted_fixture <- matrix(1, 7, 9)
  shifted_fixture[, 7:9] <- 2
  shifted_result <- geomorphon_concordance(
    geo_fixture_rast, make_concordance_raster(shifted_fixture),
    method = "boundary", threshold_cells = 0
  )
  expect_equal(shifted_result$edge_f1, 0)
  expect_true(shifted_result$edge_mhd > 0)
  expect_true(shifted_result$edge_chamfer > 0)
  expect_equal(geomorphon_concordance(
    geo_fixture_rast, cat_fixture_rast, method = "boundary",
    threshold_cells = 10)$edge_f1, 1)
  expect_equal(geomorphon_concordance(
    geo_fixture_rast, cat_fixture_rast, method = "boundary",
    threshold_cells = 1e9)$edge_f1,
    oracle_scores(fixture_geo, fixture_cat, threshold_cells = 1e9)[1],
    tolerance = 1e-8)

  # The C++ helper checks dimensions before allocating; products larger than
  # signed-int range must not wrap to a misleading cell count.
  expect_error(rgeomorphon:::spatial_null_cpp(
    1, 1, .Machine$integer.max, .Machine$integer.max, 0, TRUE, FALSE
  ), pattern = "value vectors must match")

  # Match R's is.na semantics: exclude NA and NaN, retain infinities as classes.
  nonfinite_geo <- matrix(c(1, 1, NaN, Inf, 2, 2, -Inf, NA_real_), 2, 4,
                          byrow = TRUE)
  nonfinite_cat <- matrix(c(1, 2, Inf, 3, 2, 2, -Inf, NA_real_), 2, 4,
                          byrow = TRUE)
  nonfinite_null <- rgeomorphon:::spatial_null_cpp(
    as.vector(t(nonfinite_geo)), as.vector(t(nonfinite_cat)), 2L, 4L,
    0, TRUE, FALSE
  )
  entropy_r <- function(x) {
    p <- as.numeric(table(x)) / length(x)
    -sum(p * log2(p))
  }
  direct_nonfinite <- vapply(0:6, function(i) {
    rank <- floor((i + 0.5) * 7 / 7)
    q <- rank + 1L
    dr <- q %/% 4L
    dc <- q %% 4L
    shifted <- matrix(NA_real_, 2, 4)
    for (r in 1:2) for (cc in 1:4) {
      shifted[r, cc] <- nonfinite_cat[(r - dr - 1L) %% 2L + 1L,
                                       (cc - dc - 1L) %% 4L + 1L]
    }
    valid <- !is.na(nonfinite_geo) & !is.na(shifted)
    if (!any(valid)) return(NA_real_)
    gv <- nonfinite_geo[valid]
    cv <- shifted[valid]
    hg <- entropy_r(gv); hc <- entropy_r(cv)
    if (hg + hc == 0) return(0)
    hj <- entropy_r(paste(gv, cv, sep = "_"))
    2 * max(0, hg + hc - hj) / (hg + hc)
  }, numeric(1))
  expect_equal(nonfinite_null$nmi_mean, mean(direct_nonfinite, na.rm = TRUE),
               tolerance = 1e-12)
  expect_equal(nonfinite_null$nmi_eligible, sum(!is.na(direct_nonfinite)))

  # A shared-domain hole is not an edge; compare expected masks independently.
  hole_geo <- fixture_geo
  hole_cat <- fixture_cat
  hole_geo[4, 4] <- NA_real_
  hole_cat[4, 4] <- NA_real_
  hole_valid <- !is.na(hole_geo) & !is.na(hole_cat)
  expect_equal(reference_mask(hole_geo, hole_valid),
               reference_mask(fixture_geo, hole_valid))
  hole_result <- geomorphon_concordance(
    make_concordance_raster(hole_geo), make_concordance_raster(hole_cat),
    method = "boundary", threshold_cells = 1
  )
  expect_true(is.finite(hole_result$edge_f1))

  # No-overlap and empty-boundary coordinate sets are explicit oracle fixtures.
  empty_mask <- matrix(FALSE, 5, 7)
  expect_equal(nrow(boundary_coordinates(empty_mask)), 0L)
  no_overlap_a <- matrix(c(0.5, 0.5), nrow = 1)
  no_overlap_b <- matrix(c(20.5, 20.5), nrow = 1)
  expect_true(brute_nearest(no_overlap_a, no_overlap_b) > 1)
  expect_equal(brute_nearest(no_overlap_a, matrix(numeric(), ncol = 2)),
               numeric())

  # Explicit diagonal and border fixtures: 8-neighbor mask expected by
  # coordinates, with raster borders excluded from edge candidates.
  diagonal_fixture <- matrix(1, 7, 7)
  diagonal_fixture[4, 4] <- 2
  diagonal_expected <- matrix(FALSE, 7, 7)
  diagonal_expected[3:5, 3:5] <- TRUE
  diagonal_expected[c(1, 7), ] <- FALSE
  diagonal_expected[, c(1, 7)] <- FALSE
  expect_equal(reference_mask(diagonal_fixture), diagonal_expected)
  diagonal_result <- geomorphon_concordance(
    make_concordance_raster(diagonal_fixture),
    make_concordance_raster(diagonal_fixture), method = "boundary",
    threshold_cells = 0
  )
  expect_equal(unname(unlist(diagonal_result[
    c("edge_f1", "edge_mhd", "edge_chamfer")])),
    oracle_scores(diagonal_fixture, diagonal_fixture, threshold_cells = 0),
    tolerance = 1e-8)
  border_fixture <- matrix(1, 7, 7)
  border_fixture[1, 4:7] <- 2
  border_expected <- reference_mask(border_fixture)
  expect_false(any(border_expected[c(1, 7), ]))
  expect_false(any(border_expected[, c(1, 7)]))
  border_result <- geomorphon_concordance(
    make_concordance_raster(border_fixture),
    make_concordance_raster(fixture_geo[, 1:7]), method = "boundary",
    threshold_cells = 1
  )
  expect_equal(unname(unlist(border_result[
    c("edge_f1", "edge_mhd", "edge_chamfer")])),
    oracle_scores(border_fixture, fixture_geo[, 1:7]), tolerance = 1e-8)

  # Controlled extra boundary keeps recall high but lowers symmetric F1.
  excess_fixture <- fixture_geo
  excess_fixture[3:5, 3:4] <- 2
  excess_masks_valid <- !is.na(fixture_geo) & !is.na(excess_fixture)
  excess_geo_xy <- boundary_coordinates(
    reference_mask(fixture_geo, excess_masks_valid))
  excess_cat_xy <- boundary_coordinates(
    reference_mask(excess_fixture, excess_masks_valid))
  excess_recall <- mean(brute_nearest(excess_geo_xy, excess_cat_xy) <= 1)
  excess_result <- geomorphon_concordance(
    geo_fixture_rast, make_concordance_raster(excess_fixture),
    method = "boundary", threshold_cells = 1
  )
  expect_equal(excess_result$edge_f1,
               oracle_scores(fixture_geo, excess_fixture)[1], tolerance = 1e-8)
  expect_true(excess_recall >= 0.8)
  expect_true(excess_result$edge_f1 < fixture_legacy$edge_f1)

  # Rectangular cells use both resolutions in the independent oracle; a
  # threshold at 2 map units includes the vertical displacement, not 3 units.
  anisotropic_a <- matrix(c(0.5, 1.5), nrow = 1)
  anisotropic_b <- matrix(c(0.5, 4.5), nrow = 1)
  expect_equal(brute_nearest(anisotropic_a, anisotropic_b), 3)
  expect_true(brute_nearest(anisotropic_a, anisotropic_b) > 2)

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
