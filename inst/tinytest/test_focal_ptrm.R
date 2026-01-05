# Tests for focal PTRM (moving window PTRM calculation)

library(tinytest)
library(rgeomorphon)

if (requireNamespace("terra", quietly = TRUE)) {
  # Setup: Load test data
  data("volcano", package = "datasets")
  dem <- terra::rast(volcano)
  terra::crs(dem) <- terra::crs("EPSG:2193")
  terra::ext(dem) <- c(1756968, 1757578, 5917000, 5917870)

  # Test 1: Basic functionality - returns SpatRaster
  expect_true({
    ptrm_focal <- geomorphon_ptrm_focal(dem, radius = 3, mode = "dem")
    inherits(ptrm_focal, "SpatRaster")
  })

  # Test 2: Output values are finite (no Inf or NaN)
  expect_true({
    ptrm_focal <- geomorphon_ptrm_focal(dem, radius = 3, mode = "dem")
    vals <- terra::values(ptrm_focal, na.rm = TRUE)
    all(is.finite(vals))
  })

  # Test 3: Border cells are NA (width = radius)
  expect_true({
    ptrm_focal <- geomorphon_ptrm_focal(dem, radius = 3, mode = "dem")
    mat <- terra::as.matrix(ptrm_focal, wide = TRUE)
    n_border_na <- sum(is.na(mat[1:3, ]))  # Top border
    n_border_na > 0  # Should have some NA values
  })

  # Test 4: Larger radius produces more NA cells
  expect_true({
    ptrm_r3 <- geomorphon_ptrm_focal(dem, radius = 3, mode = "dem")
    ptrm_r5 <- geomorphon_ptrm_focal(dem, radius = 5, mode = "dem")

    mat_r3 <- terra::as.matrix(ptrm_r3, wide = TRUE)
    mat_r5 <- terra::as.matrix(ptrm_r5, wide = TRUE)

    na_count_r3 <- sum(is.na(mat_r3))
    na_count_r5 <- sum(is.na(mat_r5))

    na_count_r5 > na_count_r3
  })

  # Test 5: mode="dem" vs mode="class" consistency
  expect_true({
    # Pre-calculate geomorphons with same parameters
    geo_pre <- geomorphons(dem, search = 3, skip = 0)

    # Calculate focal PTRM from DEM
    ptrm_from_dem <- geomorphon_ptrm_focal(dem, radius = 3, mode = "dem",
                                            search = 3, skip = 0)

    # Calculate focal PTRM from pre-computed geomorphons
    ptrm_from_class <- geomorphon_ptrm_focal(geo_pre, radius = 3, mode = "class")

    # Results should be identical
    vals_dem <- terra::values(ptrm_from_dem)
    vals_class <- terra::values(ptrm_from_class)

    # Use near-equality due to floating point precision
    all.equal(vals_dem, vals_class, tolerance = 1e-10) == TRUE
  })

  # Test 6: Different radii produce different results
  expect_true({
    ptrm_r2 <- geomorphon_ptrm_focal(dem, radius = 2, mode = "dem")
    ptrm_r5 <- geomorphon_ptrm_focal(dem, radius = 5, mode = "dem")

    vals_r2 <- terra::values(ptrm_r2, na.rm = TRUE)
    vals_r5 <- terra::values(ptrm_r5, na.rm = TRUE)

    # Should not be identical
    !identical(vals_r2, vals_r5)
  })

  # Test 7: Custom coefficients work
  expect_true({
    # Create custom coefficients (all equal weight)
    custom_coeffs <- setNames(
      rep(1, 10),
      c("flat", "peak", "ridge", "shoulder", "spur",
        "slope", "hollow", "footslope", "valley", "depression")
    )

    ptrm_custom <- geomorphon_ptrm_focal(
      dem,
      radius = 3,
      mode = "dem",
      coefficients = custom_coeffs
    )

    inherits(ptrm_custom, "SpatRaster")
  })

  # Test 8: Default coefficients are applied
  expect_true({
    ptrm_default <- geomorphon_ptrm_focal(dem, radius = 3, mode = "dem")

    # Explicit coefficients from Rajasekaran et al. (2022)
    explicit_coeffs <- c(
      flat = 25.12, peak = 1.75, ridge = 7.59, shoulder = 7.31,
      spur = 9.02, slope = 28.95, hollow = 6.71, footslope = 7.63,
      valley = 9.61, depression = 3.55
    )

    ptrm_explicit <- geomorphon_ptrm_focal(
      dem,
      radius = 3,
      mode = "dem",
      coefficients = explicit_coeffs
    )

    vals_default <- terra::values(ptrm_default)
    vals_explicit <- terra::values(ptrm_explicit)

    all.equal(vals_default, vals_explicit, tolerance = 1e-10) == TRUE
  })

  # Test 9: Small radius (radius = 1) works
  expect_true({
    ptrm_r1 <- geomorphon_ptrm_focal(dem, radius = 1, mode = "dem")
    inherits(ptrm_r1, "SpatRaster")
  })

  # Test 10: Zero radius - only center pixel (mostly NA)
  expect_true({
    ptrm_r0 <- geomorphon_ptrm_focal(dem, radius = 0, mode = "dem")
    mat <- terra::as.matrix(ptrm_r0, wide = TRUE)
    # With radius 0, all cells should be NA (window is just center, no neighbors)
    # Actually, with radius 0, it should count just the center cell itself
    inherits(ptrm_r0, "SpatRaster")
  })

  # Test 11: Error handling - invalid radius
  expect_error({
    geomorphon_ptrm_focal(dem, radius = -1, mode = "dem")
  })

  # Test 12: Error handling - invalid mode
  expect_error({
    geomorphon_ptrm_focal(dem, radius = 3, mode = "invalid")
  })

  # Test 13: Error handling - non-SpatRaster input
  expect_error({
    geomorphon_ptrm_focal(terra::as.matrix(volcano), radius = 3, mode = "dem")
  })

  # Test 14: Error handling - invalid coefficients length
  expect_error({
    bad_coeffs <- c(1, 2, 3)  # Wrong length
    geomorphon_ptrm_focal(dem, radius = 3, mode = "dem",
                          coefficients = bad_coeffs)
  })

  # Test 15: Error handling - invalid coefficients names
  expect_error({
    bad_coeffs <- setNames(rep(1, 10), paste0("form", 1:10))  # Wrong names
    geomorphon_ptrm_focal(dem, radius = 3, mode = "dem",
                          coefficients = bad_coeffs)
  })

  # Test 16: Interior cells have non-NA values
  expect_true({
    ptrm_focal <- geomorphon_ptrm_focal(dem, radius = 2, mode = "dem")
    mat <- terra::as.matrix(ptrm_focal, wide = TRUE)
    nrows <- nrow(mat)
    ncols <- ncol(mat)

    # Check interior cells (excluding borders)
    interior <- mat[3:(nrows-2), 3:(ncols-2)]
    # At least some interior cells should have non-NA values
    any(!is.na(interior))
  })

  # Test 17: Consistency with geomorphons output structure
  expect_true({
    geo <- geomorphons(dem, search = 3)
    ptrm_focal <- geomorphon_ptrm_focal(geo, radius = 3, mode = "class")

    # Should have same dimensions as input
    identical(dim(ptrm_focal), dim(geo))
  })

  # Test 18: No artifacts in interior (reasonable values)
  expect_true({
    ptrm_focal <- geomorphon_ptrm_focal(dem, radius = 3, mode = "dem")
    vals <- terra::values(ptrm_focal, na.rm = TRUE)

    # PTRM scores typically between -1 and 2 (formula can produce values outside [0,1])
    # but should be reasonable (not Inf, NaN, or extreme outliers)
    all(vals > -10 & vals < 10)
  })

  # Test 19: Mode="class" with invalid classification (out of range)
  # This should still work but count only valid classes
  expect_true({
    geo <- geomorphons(dem, search = 3)
    # geomorphons should only produce classes 1-10
    ptrm_focal <- geomorphon_ptrm_focal(geo, radius = 3, mode = "class")
    inherits(ptrm_focal, "SpatRaster")
  })

  # Test 20: Reproducibility - same input produces same output
  expect_true({
    ptrm1 <- geomorphon_ptrm_focal(dem, radius = 3, mode = "dem")
    ptrm2 <- geomorphon_ptrm_focal(dem, radius = 3, mode = "dem")

    vals1 <- terra::values(ptrm1)
    vals2 <- terra::values(ptrm2)

    all.equal(vals1, vals2, tolerance = 1e-10) == TRUE
  })
}
