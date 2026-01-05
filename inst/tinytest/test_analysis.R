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

  # ===========================================================================
  # Integration Tests
  # ===========================================================================

  # Test 25: Full workflow with PTRM and concordance
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
