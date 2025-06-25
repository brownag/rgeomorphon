# Rcpp::sourceCpp("rgeomorphon") # for testing

#' Calculate Geomorphons
#'
#' Parallel implementation of the 'geomorphon' terrain classification algorithm
#' largely based on 'r.geomorphon' algorithm of Jasiewicz and Stepinski (2013)
#' from 'GRASS GIS'.
#'
#' The algorithm assumes planar distances and angles are calculated based on
#' cell resolutions, so it is strongly recommended that elevation data be in a
#' projected coordinate system.
#'
#' For reliable geomorphon classification, especially near study area
#' boundaries, it is recommended to use a raster that includes a buffer of at
#' least `search + 1` cells around the area of interest. This implementation
#' utilizes all available DEM data up to the specified search radius or the
#' boundary (whichever is encountered first).
#'
#' @param elevation matrix or SpatRaster object. Digital Elevation Model values.
#'   It is **STRONGLY** recommended to use a grid in a projected coordinate system.
#' @param search numeric. User input for search radius (default: `3`). Units
#'   depend on `use_meters`.
#' @param skip numeric. User input for skip radius (default: `0`). Units depend
#'   on `use_meters`.
#' @param flat_angle_deg numeric. Flatness angle threshold in **degrees**.
#'   Default: `1.0`.
#' @param dist numeric. Flatness distance (default: `0`). Units depend on
#'   `use_meters`.
#' @param comparison_mode Character. One of `"anglev1"`, `"anglev2"`,
#'   `"anglev2_distance"`. Default: `"anglev1"`.
#' @param tdist numeric. Terrain distance factor. When greater than 0, overrides
#'   Z tolerance from angular logic. Default: `0.0`.
#' @param use_meters Logical. Default: `FALSE` uses cell units. Set to `TRUE` to
#'   specify `search`, `skip`, and `dist` in units of meters.
#' @param forms Character. Number of geomorphon forms to identify. One of
#'   `"forms10` (default), `"forms6"`, `"forms5"`, or `"forms4`.
#' @param ternary Logical. Include "ternary" output? Default: `FALSE`
#' @param positive Logical. Include "positive" output? Default: `FALSE`
#' @param negative Logical. Include "negative" output? Default: `FALSE`
#' @param nodata_val numeric. NODATA value. Default: `NA_real_`.
#' @param xres numeric. X grid resolution (used only when `elevation` is a
#'   matrix). Default: `NULL`.
#' @param yres numeric. Y grid resolution (used only when `elevation` is a
#' matrix). Default: `xres`.
#' @param simplify Logical. If result is length `1` list, the first element is
#'   returned. Default: `FALSE`
#'
#' @return List of SpatRaster or matrix of geomorphon algorithm outputs. When
#'   more than one of `forms`, `ternary`, `positive`, `negative` are set the
#'   result is a list. For one result type, and default `simplify` argument, the
#'   result is the first (and only) element of the list.
#'
#' @seealso [geomorphon_theme()]
#'
#' @details
#'
#' This implementation achieves very high agreement with the classification
#' logic of GRASS GIS 'r.geomorphon' when using equivalent parameters and data
#' in a projected coordinate system.
#'
#' 'r.geomorphon' employs a row buffering strategy which can, for cells near the
#' edges of a raster, result in a truncated line-of-sight compared to the full
#' raster extent. This may lead GRASS to classify edge-region cells differently
#' or as NODATA where this implementation may produce a more 'valid' geomorphon
#' form given the available data.
#'
#' More information about the 'r.geomorphon' module can be found in the GRASS
#' GIS manual: \url{https://grass.osgeo.org/grass-stable/manuals/r.geomorphon.html}
#'
#' @references
#'
#' Stepinski, T., Jasiewicz, J., 2011, Geomorphons - a new approach to
#' classification of landform, in : Eds: Hengl, T., Evans, I.S., Wilson, J.P.,
#' and Gould, M., Proceedings of Geomorphometry 2011, Redlands, 109-112.
#' Available online:
#' \url{https://www.geomorphometry.org/uploads/pdf/pdf2011/StepinskiJasiewicz2011geomorphometry.pdf}
#'
#' Jasiewicz, J., Stepinski, T., 2013, Geomorphons - a pattern recognition
#' approach to classification and mapping of landforms, Geomorphology, vol. 182,
#' 147-156. (\doi{doi:10.1016/j.geomorph.2012.11.005})
#'
#' @export
#' @examplesIf requireNamespace("terra")
#' library(terra)
#' library(rgeomorphon)
#'
#' SEARCH = 7       # outer search radius (cells)
#' SKIP = 1         # inner skip radius (cells)
#' DIST = 0         # flatness distance (cells)
#' FLAT = 1         # flat angle threshold
#' MODE = "anglev1" # comparison mode
#'
#' ## classic volcano
#' data("volcano", package = "datasets")
#' dem <- terra::rast(volcano)
#' terra::crs(dem) <- terra::crs("EPSG:2193")
#' terra::ext(dem) <- c(1756968, 1757578, 5917000, 5917870)
#' names(dem) <- "elevation"
#'
#' system.time({
#'     rg <- geomorphons(
#'         dem,
#'         search = SEARCH,
#'         skip = SKIP,
#'         dist = DIST,
#'         flat = FLAT,
#'         comparison_mode = MODE
#'     )
#' })
#'
#' plot(c(dem, rg))
#'
geomorphons <- function(elevation,
                        search = 3,
                        skip = 0,
                        flat_angle_deg = 1.0,
                        dist = 0.0,
                        comparison_mode = "anglev1",
                        tdist = 0.0,
                        forms = TRUE,
                        ternary = FALSE,
                        positive = FALSE,
                        negative = FALSE,
                        use_meters = FALSE,
                        nodata_val = NA_real_,
                        xres = NULL,
                        yres = xres,
                        simplify = FALSE) {

    outputs <- character(0)
    forms_int <- -1

    if (!is.null(forms) && !isFALSE(forms)) {
        if (isTRUE(forms) || forms == "forms") {
            forms_int <- 10
        } else {
            fi <- as.integer(gsub("forms", "", forms))
            if (is.na(fi)) {
                if (!forms %in% c("forms10", "forms6", "forms5", "forms4")) {
                    message("Valid values for 'forms' include: `TRUE` ('forms10'; default), 'forms6', 'forms5', 'forms4'", call. = FALSE)
                }
                fi <- 10
            }
            forms_int <- fi
        }
        forms <- "forms"
    }

    if (is.null(ternary)) {
        ternary <- FALSE
    }

    if (is.null(positive)) {
        positive <- FALSE
    }

    if (is.null(negative)) {
        negative <- FALSE
    }

    outputs <- c(forms = forms,
                 ternary = if (isTRUE(ternary)) "ternary" else character(0),
                 positive = if (isTRUE(positive)) "positive" else character(0),
                 negative = if (isTRUE(negative)) "negative" else character(0))

    comparison_mode <- tolower(comparison_mode)

    if (!is.matrix(elevation)) {
        if (!requireNamespace("terra")) {
            stop("Package 'terra' is required for SpatRaster input")
        }

        if (!inherits(elevation, "SpatRaster")) {
            stop("Input 'elevation' must be a SpatRaster or matrix.")
        }

        if (terra::nlyr(elevation) > 1) {
            warning("Input raster has multiple layers. Using the first layer.")
            elevation <- elevation[[1]]
        }

        if (terra::is.lonlat(elevation, warn = FALSE)) {
            warning("Input SpatRaster appears to be in geographic coordinates (longitude/latitude).\n\n",
                    "This algorithm assumes planar distances and angles are calculated based on cell resolutions. ",
                    "Results will be inaccurate--especially for diagonal views and areas away from the equator.\n\n",
                    "It is STRONGLY recommended to project the DEM to a suitable Cartesian coordinate system before use.")
        }

        dem_values <- terra::as.matrix(elevation, wide = TRUE)

        x_res <- terra::xres(elevation)
        y_res <- terra::yres(elevation)

    } else {
        if (is.null(xres) || is.null(yres) ||
            is.na(xres) || is.na(yres) ||
            !is.numeric(xres) || !is.numeric(yres)) {
            stop("When `elevation` is a matrix `xres` and `yres` must be positive numeric values.")
        }

        dem_values <- elevation
        x_res <- xres
        y_res <- yres
    }

    nodata <- NA_real_

    if (is.na(nodata_val)) {
        nodata <- NA_real_
    } else {
        dem_values[is.na(dem_values)] <- nodata_val
        nodata <- as.double(nodata_val)
    }

    if ((x_res <= 0 || y_res <= 0)) {
        stop("Cell resolution (xres, yres) must be positive.")
    }

    valid_comp_modes <- c("anglev1", "anglev2", "anglev2_distance")
    if (!comparison_mode %in% valid_comp_modes) {
        stop(paste(
            "Invalid comparison_mode. Must be one of:",
            toString(valid_comp_modes)
        ))
    }

    forms_matrix_res <- geomorphons_cpp_worker(
        elevation = dem_values,
        search = as.double(search),
        skip = as.double(skip),
        flat_angle_deg = as.double(flat_angle_deg),
        dist = as.double(dist),
        comparison_mode = as.character(comparison_mode),
        tdist = as.double(tdist),
        use_meters = as.logical(use_meters),
        x_res_dem = as.double(x_res),
        y_res_dem = as.double(y_res),
        forms = as.integer(forms_int),
        ternary = as.integer(ternary),
        positive = as.integer(positive),
        negative = as.integer(negative),
        nodata = nodata
    )

    if (inherits(elevation, 'SpatRaster')) {
        forms_rast <- terra::rast(elevation)
        terra::values(forms_rast) <- forms_matrix_res[["forms"]]
        forms_rast <- geomorphon_theme(forms_rast, forms = forms)
        ternary_rast <- terra::rast(elevation)
        terra::values(ternary_rast) <- forms_matrix_res[["ternary"]]
        positive_rast <- terra::rast(elevation)
        terra::values(positive_rast) <- forms_matrix_res[["positive"]]
        negative_rast <- terra::rast(elevation)
        terra::values(negative_rast) <- forms_matrix_res[["negative"]]
    } else {
        forms_rast <- forms_matrix_res[["forms"]]
        ternary_rast <- forms_matrix_res[["ternary"]]
        positive_rast <- forms_matrix_res[["positive"]]
        negative_rast <- forms_matrix_res[["negative"]]
    }

    res <- list(
        forms = forms_rast,
        ternary = ternary_rast,
        positive = positive_rast,
        negative = negative_rast
    )[outputs]

    if (isTRUE(simplify) && length(res) == 1) {
        res <- res[[1]]
    }

    if (inherits(elevation, 'SpatRaster')) {
        terra::rast(res)
    } else {
        res
    }
}

#' @export
#' @rdname geomorphon_theme
geomorphon_categories <- function() {
    data.frame(
        value = 1:10,
        form = c(
            "Flat",
            "Peak",
            "Ridge",
            "Shoulder",
            "Spur",
            "Slope",
            "Hollow",
            "Footslope",
            "Valley",
            "Pit"
        )
    )
}

#' @export
#' @rdname geomorphon_theme
geomorphon_colors <- function() {
    data.frame(
        value = 1:10,
        col = c(
            "#DCDCDC",
            "#380000",
            "#C80000",
            "#FF5014",
            "#FAD23C",
            "#FFFF3C",
            "#B4E614",
            "#3CFA96",
            "#0000FF",
            "#000038"
        )
    )
}

#' Apply Geomorphon Theme to Result Object
#'
#' Applies standard class names and colors to a SpatRaster, or creates a factor
#' matrix. Input values should be integers between 1 and 10.
#'
#' @param x A SpatRaster or matrix object.
#' @param forms character. One of: `"forms10"` (default), `"forms6"`,
#'   `"forms5"`, `"forms4"`. These are themes corresponding to the built-in
#'   10-form, 6-form, 5-form, and 4-form `"forms`" outputs from [geomorphons()].
#'
#' @details
#'
#' When `x` is a matrix the result is a factor using `geomorphon_categories()`.
#' Values are integers 1 to 10 and labels are the geomorphon form names.
#'
#' @return A SpatRaster or matrix object with geomorphon class names (and colors
#'   for SpatRaster) applied.
#'
#' @export
#' @rdname geomorphon_theme
geomorphon_theme <- function(x, forms = "forms10") {
    if (inherits(x, 'SpatRaster')) {
        names(x) <- "geomorphon"
        levels(x) <- geomorphon_categories()
        terra::coltab(x) <- geomorphon_colors()
    } else {
        cls <- geomorphon_categories()
        x[] <- as.character(factor(x, levels = cls$value, labels = cls$form))
    }

    x
}
