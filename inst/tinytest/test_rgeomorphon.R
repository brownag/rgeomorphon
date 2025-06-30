library(tinytest)
library(rgeomorphon)

SEARCH <- 7        # outer search radius (cells)
SKIP <- 1          # inner skip radius (cells)
DIST <- 2          # flatness distance (cells)
FLAT <- 1          # flat angle threshold
MODE <- "anglev1"  # comparison mode

## classic volcano
data("volcano", package = "datasets")

# matrix interface
res0 <- geomorphons(volcano, xres = 10, simplify = TRUE)
expect_equivalent(dim(res0), dim(volcano))

# SpatRaster interface requires terra
if (requireNamespace("terra", quietly = TRUE)) {
    dem <- terra::rast(volcano)
    terra::crs(dem) <- terra::crs("EPSG:2193")
    terra::ext(dem) <- c(1756968, 1757578, 5917000, 5917870)
    names(dem) <- "elevation"

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

    expect_warning({
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
}
