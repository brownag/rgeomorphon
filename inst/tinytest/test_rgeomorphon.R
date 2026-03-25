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
}
