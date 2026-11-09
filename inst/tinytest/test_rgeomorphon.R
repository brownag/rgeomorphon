library(tinytest)
library(rgeomorphon)

SEARCH <- 7        # outer search radius (cells)
SKIP <- 1          # inner skip radius (cells)
DIST <- 2          # flatness distance (cells)
FLAT <- 1          # flat angle threshold
MODE <- "anglev1"  # comparison mode

## classic volcano elevation dataset
data("volcano", package = "datasets")

# matrix interface
res0 <- geomorphons(volcano, xres = 10, simplify = TRUE)
expect_equivalent(dim(res0), dim(volcano))

# matrix interface: multi-form count reuse and single-form backward compatibility
res_mat_multi <- geomorphons(volcano, xres = 10, forms = c("forms10", "forms6"), search = SEARCH, skip = SKIP)
res_mat_10 <- geomorphons(volcano, xres = 10, forms = "forms10", search = SEARCH, skip = SKIP)
res_mat_6 <- geomorphons(volcano, xres = 10, forms = "forms6", search = SEARCH, skip = SKIP)
expect_true(is.list(res_mat_multi))
expect_true(all(c("forms10", "forms6", "forms") %in% names(res_mat_multi)))
expect_identical(res_mat_multi[["forms10"]], res_mat_10[["forms"]])
expect_identical(res_mat_multi[["forms6"]], res_mat_6[["forms"]])

res_mat_single <- geomorphons(volcano, xres = 10, forms = TRUE, search = SEARCH, skip = SKIP, simplify = TRUE)
expect_true(is.matrix(res_mat_single))
expect_equivalent(dim(res_mat_single), dim(volcano))

# SpatRaster interface requires terra
if (requireNamespace("terra", quietly = TRUE)) {
    # construct and georeference a SpatRaster object
    dem <- terra::flip(terra::rast(volcano))
    terra::crs(dem) <- terra::crs("EPSG:27200")
    terra::ext(dem) <- c(2667400, 2668010, 6478700, 6479570)
    names(dem) <- "Elevation (meters)"

    expect_silent({
        rg <- geomorphons(
            dem,
            search = SEARCH,
            skip = SKIP,
            dist = DIST,
            flat = FLAT,
            comparison_mode = MODE
        )
    })

    expect_true(inherits(rg, 'SpatRaster'))
    expect_equivalent(nrow(terra::unique(rg)), 10)
    expect_equivalent(nrow(terra::cats(rg)[[1]]), 10)

    expect_silent({
        res <- geomorphons(
            dem,
            search = SEARCH,
            skip = SKIP,
            dist = DIST,
            flat = FLAT,
            comparison_mode = MODE,
            forms = TRUE,
            ternary = TRUE,
            positive = TRUE,
            negative = TRUE
        )
    })

    res2 <- terra::rast(lapply(c(4, 5, 6), function(n) {
      geomorphon_theme(
        forms_matrix_apply(
            x = res[[c("positive", "negative")]],
            rcl = forms_matrix_get(n)
        )
      )
    }))
    names(res2) <- c("forms4", "forms5", "forms6")

    expect_equivalent(terra::nlyr(c(res, res2)), 7)

    # simplified SpatRaster output must keep its values and levels
    geo_simplified <- geomorphons(dem, search = 3, simplify = TRUE)
    geo_full <- geomorphons(dem, search = 3)

    expect_true(inherits(geo_simplified, 'SpatRaster'))
    expect_equivalent(terra::nlyr(geo_simplified), 1L)
    expect_true(all(!is.na(terra::values(geo_simplified))))
    expect_equivalent(
        as.vector(terra::values(geo_simplified)),
        as.vector(terra::values(geo_full[[1]]))
    )
    expect_equivalent(names(geo_simplified), names(geo_full[[1]]))
    expect_equivalent(nrow(terra::cats(geo_simplified)[[1]]), 10)
    expect_equivalent(
        terra::cats(geo_simplified)[[1]]$form,
        geomorphon_categories()$form
    )

    expect_silent({
        rg <- geomorphons(
            dem,
            search = SEARCH,
            skip = SKIP,
            dist = DIST,
            flat = FLAT,
            comparison_mode = MODE,
            nchunk = 2
        )
    })

    # Pixel-for-pixel identity between vector forms and individual calls
    res_multi <- geomorphons(dem, search = SEARCH, skip = SKIP, forms = c("forms10", "forms6"))
    res_10 <- geomorphons(dem, search = SEARCH, skip = SKIP, forms = "forms10")
    res_6 <- geomorphons(dem, search = SEARCH, skip = SKIP, forms = "forms6")

    expect_true(inherits(res_multi, 'SpatRaster'))
    expect_equal(terra::nlyr(res_multi), 2L)
    expect_equal(names(res_multi), c("forms10", "forms6"))
    expect_equal(as.vector(terra::values(res_multi[["forms10"]])), as.vector(terra::values(res_10)))
    expect_equal(as.vector(terra::values(res_multi[["forms6"]])), as.vector(terra::values(res_6)))
    expect_equal(terra::cats(res_multi[["forms10"]])[[1]]$form, geomorphon_categories()$form)
    expect_equal(terra::cats(res_multi[["forms6"]])[[1]]$form, geomorphon_categories()$form)

    # Pixel-for-pixel identity for forms5 and forms4 against individual calls
    res_multi_all <- geomorphons(dem, search = SEARCH, skip = SKIP,
                                 forms = c("forms10", "forms6", "forms5", "forms4"))
    res_5 <- geomorphons(dem, search = SEARCH, skip = SKIP, forms = "forms5")
    res_4 <- geomorphons(dem, search = SEARCH, skip = SKIP, forms = "forms4")
    expect_equal(names(res_multi_all), c("forms10", "forms6", "forms5", "forms4"))
    expect_equal(as.vector(terra::values(res_multi_all[["forms5"]])), as.vector(terra::values(res_5)))
    expect_equal(as.vector(terra::values(res_multi_all[["forms4"]])), as.vector(terra::values(res_4)))

    # Tiled execution (nchunk = 2) with multiple forms matches the untiled multi-form run
    res_tiled <- geomorphons(dem, search = SEARCH, skip = SKIP,
                             forms = c("forms10", "forms6"), nchunk = 2)
    expect_equal(names(res_tiled), c("forms10", "forms6"))
    expect_equal(as.vector(terra::values(res_tiled[["forms10"]])), as.vector(terra::values(res_multi[["forms10"]])))
    expect_equal(as.vector(terra::values(res_tiled[["forms6"]])), as.vector(terra::values(res_multi[["forms6"]])))

    # Single-form backward compatibility: returns 1-layer SpatRaster named "form"
    res_s_true <- geomorphons(dem, search = SEARCH, skip = SKIP, forms = TRUE)
    res_s_10 <- geomorphons(dem, search = SEARCH, skip = SKIP, forms = "forms10")
    res_s_6 <- geomorphons(dem, search = SEARCH, skip = SKIP, forms = "forms6")

    expect_equal(terra::nlyr(res_s_true), 1L)
    expect_equal(names(res_s_true), "form")
    expect_equal(terra::nlyr(res_s_10), 1L)
    expect_equal(names(res_s_10), "form")
    expect_equal(terra::nlyr(res_s_6), 1L)
    expect_equal(names(res_s_6), "form")

    # Equivalence between forms_matrix_apply() C++ backend and legacy results
    rcl6 <- forms_matrix_get(6)
    out_cpp <- forms_matrix_apply(res[[c("positive", "negative")]], rcl = rcl6)

    d <- nrow(rcl6)
    p_lyr <- res[["positive"]]
    n_lyr <- res[["negative"]]
    r_id <- n_lyr * d + p_lyr
    lut <- expand.grid(negative = seq(0, d - 1), positive = seq(0, d - 1))
    lut$form <- unlist(apply(rcl6, MARGIN = 2, function(a) a, simplify = FALSE))
    lut <- cbind(ID = lut$negative * d + lut$positive, lut)
    lut <- lut[order(lut$ID), ]
    out_legacy <- terra::subst(r_id, lut$ID, lut$form)
    names(out_legacy) <- "form"

    expect_equal(as.vector(terra::values(out_cpp)), as.vector(terra::values(out_legacy)))

    # forms_matrix_apply on matrix / list input
    pos_mat <- terra::as.matrix(p_lyr, wide = TRUE)
    neg_mat <- terra::as.matrix(n_lyr, wide = TRUE)
    mat_cpp <- forms_matrix_apply(list(positive = pos_mat, negative = neg_mat), rcl = rcl6)
    expect_equal(as.vector(t(mat_cpp)), as.vector(terra::values(out_legacy)))
}

# Out-of-range and infinite counts classify as Flat without an unsafe int cast
expect_equal(
    rgeomorphon:::forms_matrix_apply_cpp(matrix(c(Inf, 1e20, 9, -1e20), 2), matrix(0, 2, 2), 10L),
    matrix(1, 2, 2)
)

# Fast reclassification returns themed output: named "form", factor levels, color table
if (requireNamespace("terra", quietly = TRUE)) {
    p_fast <- terra::rast(nrows = 3, ncols = 3, vals = c(0, 1, 2, 0, 3, 1, 0, 0, 2))
    n_fast <- terra::rast(nrows = 3, ncols = 3, vals = c(0, 0, 0, 1, 2, 1, 4, 0, 1))
    x_fast <- c(p_fast, n_fast)
    names(x_fast) <- c("positive", "negative")
    out_fast <- forms_matrix_apply(x_fast)
    expect_equal(names(out_fast), "form")
    expect_true(terra::is.factor(out_fast))
    expect_equal(nrow(terra::levels(out_fast)[[1]]), 10L)
    expect_equal(nrow(terra::coltab(out_fast)[[1]]), 10L)
}

# One-row 3-d array input keeps its 1 x 5 shape and matches the equivalent matrix call
pos_row <- matrix(c(0, 1, 2, 0, 3), nrow = 1)
neg_row <- matrix(c(0, 0, 0, 1, 2), nrow = 1)
arr_row <- array(c(pos_row, neg_row), dim = c(1, 5, 2),
                 dimnames = list(NULL, NULL, c("positive", "negative")))
res_arr_row <- forms_matrix_apply(arr_row)
res_list_row <- forms_matrix_apply(list(positive = pos_row, negative = neg_row))
expect_equal(dim(res_arr_row), c(1L, 5L))
expect_equal(as.vector(res_arr_row), as.vector(res_list_row))
