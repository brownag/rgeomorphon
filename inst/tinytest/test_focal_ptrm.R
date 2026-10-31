# Tests for focal PTRM (moving window PTRM calculation)

library(tinytest)
library(rgeomorphon)

if (requireNamespace("terra", quietly = TRUE)) {
  make_raster <- function(values, nrow = NULL, ncol = NULL) {
    if (is.null(nrow)) {
      nrow <- base::nrow(values)
      ncol <- base::ncol(values)
    }
    r <- terra::rast(
      nrows = nrow,
      ncols = ncol,
      xmin = 0,
      xmax = ncol,
      ymin = 0,
      ymax = nrow
    )
    vals <- if (is.matrix(values)) as.vector(t(values)) else values
    terra::values(r) <- vals
    r
  }

  flat_score <- (25.12 - 38.02) / 69.96

  # Hand-computed score: constant 3x3, radius 1 (cross of 5 flat cells)
  r3 <- make_raster(matrix(1, nrow = 3, ncol = 3))
  ptrm3 <- geomorphon_ptrm_focal(r3, radius = 1, mode = "class")
  vals3 <- as.vector(terra::values(ptrm3))
  expect_equal(vals3[5], flat_score, tolerance = 1e-12)
  expect_true(all(is.na(vals3[-5])))

  # Border width equals radius
  r7 <- make_raster(matrix(1, nrow = 7, ncol = 7))
  ptrm_r2 <- geomorphon_ptrm_focal(r7, radius = 2, mode = "class")
  mat_r2 <- terra::as.matrix(ptrm_r2, wide = TRUE)
  expect_true(all(is.na(mat_r2[1:2, ])))
  expect_true(all(is.na(mat_r2[6:7, ])))
  expect_true(all(is.na(mat_r2[, 1:2])))
  expect_true(all(is.na(mat_r2[, 6:7])))
  expect_true(all(!is.na(mat_r2[3:5, 3:5])))

  # Radius 0 is valid and has no border NAs (centre cell only)
  ptrm_r0 <- geomorphon_ptrm_focal(r3, radius = 0, mode = "class")
  vals_r0 <- as.vector(terra::values(ptrm_r0))
  expect_true(all(!is.na(vals_r0)))
  expect_equal(vals_r0, rep(flat_score, 9), tolerance = 1e-12)

  # C++ vs fallback parity on identical remapped input
  pattern <- matrix(((seq_len(49) - 1) %% 10) + 1, nrow = 7, ncol = 7,
                    byrow = TRUE)
  pattern[2, 3] <- NA
  pattern[4, 4] <- NA
  pattern[1, 1] <- 99  # outside the forms10 system, mapped to NA
  r_pattern <- make_raster(pattern)
  model <- rgeomorphon:::.get_ptrm_model("forms10", NULL, -38.02, 69.96)
  ptrm_cpp <- rgeomorphon:::.geomorphon_ptrm_focal_impl(
    r_pattern, 2, model, force_cpp = TRUE
  )
  ptrm_fallback <- rgeomorphon:::.geomorphon_ptrm_focal_impl(
    r_pattern, 2, model, force_cpp = FALSE
  )
  expect_equal(
    terra::values(ptrm_cpp),
    terra::values(ptrm_fallback),
    tolerance = 1e-10
  )

  # Output metadata: fresh raster named ptrm, no levels, numeric values
  expect_equal(names(ptrm_cpp), "ptrm")
  no_levels <- function(r) {
    lv <- terra::levels(r)
    is.null(lv) || length(lv) == 0 ||
      all(vapply(
        lv,
        function(z) is.null(z) || length(z) == 0 ||
          (is.character(z) && all(!nzchar(z))),
        logical(1)
      ))
  }
  expect_true(no_levels(ptrm_cpp))
  expect_true(is.numeric(terra::values(ptrm_cpp)))
  expect_equal(terra::ncell(ptrm_cpp), 49)

  # Reduced form systems produce finite scores
  for (form_system in c("forms4", "forms5", "forms6")) {
    ptrm_form <- geomorphon_ptrm_focal(
      r3, radius = 1, mode = "class", forms = form_system
    )
    expect_true(is.finite(as.vector(terra::values(ptrm_form))[5]))
  }

  # mode = "dem" forwards forms to geomorphons and keeps ... pass-through
  dem <- make_raster(matrix(seq_len(49), nrow = 7, ncol = 7))
  ptrm_dem <- geomorphon_ptrm_focal(
    dem, radius = 1, mode = "dem", search = 1, forms = "forms4"
  )
  expect_true(any(is.finite(terra::values(ptrm_dem))))

  # Class codes outside the requested system are not counted
  r_outside <- make_raster(matrix(2, nrow = 3, ncol = 3))
  ptrm_outside <- geomorphon_ptrm_focal(
    r_outside, radius = 1, mode = "class", forms = "forms4"
  )
  expect_true(all(is.na(terra::values(ptrm_outside))))

  # Radius validation
  expect_error(
    geomorphon_ptrm_focal(r3, radius = -1, mode = "class"),
    "non-negative"
  )
  expect_error(
    geomorphon_ptrm_focal(r3, radius = 1.5, mode = "class"),
    "whole number"
  )
  expect_error(
    geomorphon_ptrm_focal(r3, radius = 2, mode = "class"),
    "too large"
  )

  # Invalid coefficient names are rejected by the shared model helper
  bad_coeffs <- setNames(rep(1, 10), paste0("form", 1:10))
  expect_error(
    geomorphon_ptrm_focal(
      r3, radius = 1, mode = "class", coefficients = bad_coeffs
    )
  )

  # Shared memory gate helper is callable with a bytes-per-cell argument
  expect_true(is.logical(rgeomorphon:::.focal_memory_ok(r3, 12)))
}
