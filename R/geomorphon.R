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
#'   It is **STRONGLY** recommended to use a grid in a projected coordinate
#'   system.
#' @param filename character. Output filename. Default `NULL` creates a
#'   temporary file.
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
#' @param forms character. Number of geomorphon forms to identify. One of
#'   `"forms10` (default), `"forms6"`, `"forms5"`, or `"forms4`.
#' @param ternary logical. Include "ternary" output? Default: `FALSE`
#' @param positive logical. Include "positive" output? Default: `FALSE`
#' @param negative logical. Include "negative" output? Default: `FALSE`
#' @param nodata_val numeric. NODATA value. Default: `NA_integer_`.
#' @param xres numeric. X grid resolution (used only when `elevation` is a
#'   matrix). Default: `NULL`.
#' @param yres numeric. Y grid resolution (used only when `elevation` is a
#'   matrix). Default: `xres`.
#' @param simplify logical. If result is length `1` list, the first element is
#'   returned. Default: `FALSE`
#' @param LAPPLY.FUN An [lapply()]-like function such as
#'   `future.apply::future_lapply()`. Default: `lapply()`.
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
                        filename = NULL,
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
                        nodata_val = NA_integer_,
                        xres = NULL,
                        yres = xres,
                        simplify = FALSE,
                        LAPPLY.FUN = lapply) {

    if (inherits(elevation, 'SpatRaster')) {
        if (!requireNamespace("terra")) {
            stop("Package 'terra' is required to process SpatRaster input.")
        }
        mi <- as.data.frame(t(terra::mem_info(elevation, print = FALSE)))
        nchunk <- as.integer(ceiling(nrow(elevation) / mi$chunk_rows))
        return(
            .geomorphons_tiled(
                elevation,
                y = nchunk,
                search = search,
                skip = skip,
                flat_angle_deg = flat_angle_deg,
                dist = dist,
                comparison_mode = comparison_mode,
                tdist = tdist,
                forms = forms,
                ternary = ternary,
                positive = positive,
                negative = negative,
                use_meters = use_meters,
                nodata_val = nodata_val,
                xres = xres,
                yres = yres,
                filename = filename,
                overwrite = TRUE,
                LAPPLY.FUN = LAPPLY.FUN,
                simplify = simplify
            )
        )
    } else {
        return(
            .geomorphons(
                elevation,
                search = search,
                skip = skip,
                flat_angle_deg = flat_angle_deg,
                dist = dist,
                comparison_mode = comparison_mode,
                tdist = tdist,
                forms = forms,
                ternary = ternary,
                positive = positive,
                negative = negative,
                use_meters = use_meters,
                nodata_val = nodata_val,
                xres = xres,
                yres = yres,
                filename = NULL,
                overwrite = TRUE,
                simplify = simplify
            )
        )
    }
}

.geomorphons <- function(elevation,
                         search,
                         skip,
                         flat_angle_deg,
                         dist,
                         comparison_mode,
                         tdist,
                         forms,
                         ternary,
                         positive,
                         negative,
                         use_meters,
                         nodata_val,
                         xres,
                         yres,
                         filename,
                         overwrite,
                         simplify) {
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

    outputs <- character(0)

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
        tmp <- terra::rast(elevation)
        forms_rast <- tmp
        terra::values(forms_rast) <- forms_matrix_res[["forms"]]
        forms_rast <- geomorphon_theme(forms_rast, forms = forms)
        ternary_rast <- tmp
        terra::values(ternary_rast) <- forms_matrix_res[["ternary"]]
        positive_rast <- tmp
        terra::values(positive_rast) <- forms_matrix_res[["positive"]]
        negative_rast <- tmp
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
        res <- terra::rast(res)
        if (!is.null(filename)) {
            res <- terra::writeRaster(
                res,
                filename = filename,
                overwrite = TRUE,
                datatype = ifelse("ternary" %in% outputs, "INT2U", "INT1U")
            )
        }
    }
    res
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

#' Calculate Geomorphons on a Large Raster via Tiling
#'
#' Processes a large Digital Elevation Model (DEM) that may not fit into memory
#' by dividing it into smaller tiles, calculating geomorphons on each tile
#' independently, and then mosaicking the results into a single, seamless output
#' raster. This function orchestrates the tiling, processing, and combining
#' steps.
#'
#' This chunked processing approach works by creating temporary, buffered tiles
#' on disk. A buffer region around each tile provides sufficient data for the
#' geomorphon calculation's search window, thus avoiding edge artifacts between
#' tiles. After processing, the buffer area is removed, and the final tiles are
#' merged.
#'
#' @param elevation `SpatRaster`. The DEM to be processed.
#' @param y integer. The number of tiles to create along each dimension (rows
#'   and columns). For example, `y = 2` creates a 2x2 grid of four tiles. A
#'   larger number creates more, smaller tiles, which can be useful for very
#'   large rasters or parallel processing, but may increase I/O overhead.
#'   Default: `2`.
#' @param filename character. The path and filename for the final, combined
#'   output SpatRaster file. This must be specified as the function writes its
#'   final result to disk.
#' @param overwrite logical. If `TRUE`, the final output file and any temporary
#'   tile files will be overwritten if they already exist. Default: `TRUE`.
#' @param ... Additional arguments passed directly to the core geomorphon
#'   processing function e.g [geomorphons()]. This includes parameters like
#'   `search`, `skip`, `flat_angle_deg`, etc.
#' @param LAPPLY.FUN An [lapply()]-like function such as
#'   `future.apply::future_lapply()`. Default: `lapply()`.
#'
#' @returns A `SpatRaster` object pointing to the final, combined geomorphon
#'   raster file specified by the `filename` argument.
#'
#' @seealso The unexported lower-level functions used by this wrapper:
#'   `.geomorphons_tile_extents()`, `.geomorphons_create_tiles()`,
#'   `.geomorphons_process_tiles()`, `.geomorphons_combine_tiles()`.
#'
#' @noRd
#' @examplesIf requireNamespace("terra")
#' # Load the 'salton' dataset
#' data("salton", package = "rgeomorphon")
#'
#' # Construct and georeference a SpatRaster object
#' dem <- terra::rast(salton)
#' terra::crs(dem) <- attr(salton, "crs")
#' terra::ext(dem) <- attr(salton, "extent")
#' names(dem) <- "Elevation (feet)"
#'
#' # Define a temporary output file
#' tf <- tempfile(fileext = ".tif")
#'
#' res <- .geomorphons_tiled(
#'     elevation = dem,
#'     filename = tf,
#'     overwrite = TRUE,
#'     search = 10,
#'     skip = 5,
#'     flat_angle_deg = 0.1,
#'     forms = "forms6",
#'     ternary = TRUE,
#'     positive = FALSE,
#'     negative = FALSE,
#'     comparison_mode = "anglev1",
#'     nodata_val = NA_real_,
#'     dist = 0,
#'     tdist = 0,
#'     use_meters = FALSE,
#'     simplify = FALSE
#' )
#'
#' res
#'
#' terra::plot(res, main = "Geomorphons from Tiled Processing")
#'
#' \dontshow{
#' # clean up temp files
#' related_files <- terra::sources(res, all=TRUE, relative=FALSE)
#' related_dirs <- unique(dirname(unlist(related_files)))
#' unlink(related_dirs, recursive = TRUE, force = TRUE)
#' }
.geomorphons_tiled <- function(elevation,
                               y = 2,
                               filename,
                               overwrite = TRUE,
                               ...,
                               LAPPLY.FUN = lapply) {

    if (missing(filename)) {
        stop("`filename` must be specified for output of combined tiled result")
    }

    args <- list(...)

    search <- 0
    if (!is.null(args[["search"]])) {
        search <- args[["search"]]
    }
    skip <- 0
    if (!is.null(args[["skip"]])) {
        skip <- args[["skip"]]
    }
    if (!is.null(args[["ternary"]])) {
        ternary <- isTRUE(args[["ternary"]])
    }
    dtype <- ifelse(ternary, "INT2U", "INT1U")
    BUFFER <- search + skip + 1
    if(length(BUFFER) == 0) {
        BUFFER <- 0
    }
    gte <- .geomorphons_tile_extents(elevation, y, cell_buffer = BUFFER)
    tiles <- .geomorphons_create_tiles(elevation, gte$buffered_polys, overwrite = overwrite)
    res <- .geomorphons_process_tiles(
        tiles,
        gte$tile_polys,
        cell_buffer = BUFFER,
        FUN = .geomorphons,
        filename = filename,
        datatype = dtype,
        ...,
        LAPPLY.FUN = LAPPLY.FUN
    )
    .geomorphons_combine_tiles(res, elevation, filename = filename, datatype = dtype, overwrite = overwrite)
}

#' Calculate Tiling Scheme and Buffered Extents for Chunked Processing
#'
#' This function defines the spatial layout for chunked raster processing. It
#' divides a large source SpatRaster into a grid of smaller tiles and then
#' calculates the extent of a buffered area around each tile. This buffered
#' extent is crucial for ensuring that neighborhood operations (like
#' geomorphons) have sufficient data around the edges of each tile to produce
#' correct results.
#'
#' @param x A `SpatRaster` object to be tiled.
#' @param y integer or `SpatVector`. The number of tiles to create along each
#'   dimension (rows and columns). For example, `y = 2` creates a 2x2 grid of
#'   tiles. When `y` is a `SpatVector`, the extent of each polygon is used to
#'   create the tiles.
#' @param cell_buffer integer. The number of cells to add as a buffer around
#'   each tile. This should be at least `search + skip + 1` from the main
#'   geomorphon function to ensure the search window is always filled with valid
#'   data.
#'
#' @return A `list` containing two `SpatVector` objects:
#'         \describe{
#'           \item{`tile_polys`}{`SpatVector` of polygons representing the unbuffered extents of each processing tile.}
#'           \item{`buffered_polys`}{`SpatVector` of polygons representing the buffered extents, which will be used to crop the source raster for processing.}
#'           \item{`cell_buffer`}{integer. Number of cells buffered around each tile.}
#'           \item{`buffer_distance`}{numeric. Distance, in map units, of buffer around each tike.}
#'         }
#'
#' @noRd
#'
#' @examplesIf requireNamespace(terra)
#' r <- terra::rast(volcano)
#' terra::ext(r) <- c(0, 610, 0, 870)
#' terra::crs(r) <- "local"
#'
#' # calculate a 2x2 tiling scheme with a 5-cell buffer
#' tile_scheme <- .geomorphons_tile_extents(r, y = 2, cell_buffer = 5)
#'
#' tile_scheme$tile_polys
#' tile_scheme$buffered_polys
#'
#' terra::plot(r)
#' terra::plot(
#'     tile_scheme$tile_polys,
#'     border = "blue",
#'     lwd = 2,
#'     add = TRUE
#' )
#' terra::plot(
#'     tile_scheme$buffered_polys,
#'     border = "red",
#'     lty = 2,
#'     add = TRUE
#' )
.geomorphons_tile_extents <- function(x, y, cell_buffer) {
    if (!inherits(x, "SpatRaster")) {
        stop("Input 'x' must be a SpatRaster object.")
    }
    if (cell_buffer < 0) {
        stop("'cell_buffer' must be a non-negative integer.")
    }

    if (is.numeric(y)) {
        tile_extents <- apply(terra::getTileExtents(x, ceiling(c(
            nrow(x) / y, ncol(x) / y
        ))), 1, FUN = terra::ext)

        tile_polygons <- do.call('rbind', lapply(
            tile_extents,
            terra::as.polygons,
            crs = terra::crs(x)
        ))
    } else {
        tile_polygons <- y
    }

    buffer_dist_map_units <- terra::res(x)[1] * cell_buffer
    buffered_polygons <- terra::buffer(tile_polygons, buffer_dist_map_units)

    return(list(
        tile_polys = tile_polygons,
        buffered_polys = buffered_polygons,
        cell_buffer = cell_buffer,
        buffer_distance = buffer_dist_map_units
    ))
}

#' Create Buffered Raster Tiles on Disk
#'
#' Using a set of buffered polygon extents, this function crops a source
#' SpatRaster into multiple smaller raster files (tiles) on disk. Each tile
#' contains the data needed for processing one chunk, including the necessary
#' buffer region. Storing tiles on disk is essential for processing DEMs that do
#' not fit into memory.
#'
#' @param x A `SpatRaster` object to be tiled.
#' @param y A `SpatVector` of polygons defining the extents of the buffered
#'   tiles to be created. This is typically the output from
#'   `geomorphons_tile_extents()`.
#' @param filename_template Character. A filename template for the output tile
#'   files. Should include a placeholder for the tile number. `terra::makeTiles`
#'   uses a `_x_y` suffix for tile indices if no placeholder is given. Example:
#'   `file.path(tempdir(), "dem_tile_.tif")`.
#' @param overwrite Logical. If `TRUE`, existing tile files will be overwritten.
#'
#' @return A `SpatRasterCollection` where each element is a SpatRaster object
#'   pointing to a tile file on disk.
#'
#' @noRd
#'
#' @examplesIf requireNamespace("terra")
#' r <- terra::rast(volcano)
#' terra::ext(r) <- c(0, 610, 0, 870)
#' terra::crs(r) <- "local"
#'
#' tile_scheme <- .geomorphons_tile_extents(r, y = 2, cell_buffer = 10)
#'
#' # Create the tile files in a temporary directory
#' tile_collection <- .geomorphons_create_tiles(
#'   x = terra::extend(r, tile_scheme$cell_buffer),
#'   y = tile_scheme$buffered_polys,
#'   filename_template = file.path(tempdir(), "tile_.tif"),
#'   overwrite = TRUE
#' )
#'
#' tile_collection
#'
#' terra::plot(tile_collection[1])
#' terra::lines(tile_scheme$tile_polys[1])
.geomorphons_create_tiles <- function(x,
                                      y,
                                      filename_template = file.path(tempdir(), "rgeomorphon_.tif"),
                                      overwrite = TRUE) {
    if (!inherits(x, "SpatRaster")) {
        stop("Input 'x' must be a SpatRaster object.")
    }
    if (!inherits(y, "SpatVector")) {
        stop("Input 'y' must be a SpatVector object.")
    }

    z <- terra::makeTiles(
        x = x,
        y = y,
        filename = filename_template,
        overwrite = overwrite
    )

    if (length(z) == 1) {
        z <- terra::rast(z)
    }

    terra::sprc(z)
}

#' Process Tiled Rasters to Calculate Geomorphons
#'
#' This function iterates through a collection of buffered raster tiles, applies
#' a specified geomorphon function to each, and saves the correctly cropped
#' results. The buffer is used in the calculation but removed from the final
#' output tile, ensuring seamless results when the tiles are later merged.
#'
#' @param x A `SpatRasterCollection` containing the buffered input tiles,
#'   typically from `geomorphons_create_tiles()`.
#' @param y A `SpatVector` of polygons defining the original, unbuffered extents
#'   of each tile. This is used for cropping the results. This is typically the
#'   output from `geomorphons_tile_extents()`.
#' @param FUN The geomorphon function to apply to each tile (e.g.,
#'   [geomorphons()]).
#' @param cell_buffer Integer. The number of cells to add as a buffer around
#'   each tile. Default `NULL` uses `search + skip + 1`.
#' @param ... Additional arguments to be passed to `FUN` (e.g., `search`,
#'   `skip`, `flat_angle_deg`, etc.).
#' @param LAPPLY.FUN An [lapply()]-like function such as
#'   `future.apply::future_lapply()`. Default: `lapply()`.
#'
#' @return A `SpatRasterCollection` where each element points to the processed,
#'   cropped result file for each corresponding input tile.
#' @noRd
.geomorphons_process_tiles <- function(x,
                                       y,
                                       cell_buffer = NULL,
                                       FUN = geomorphons,
                                       ...,
                                       datatype,
                                       LAPPLY.FUN = lapply) {
    if (!inherits(x, "SpatRasterCollection")) {
        stop("Input 'x' must be a SpatRasterCollection.")
    }
    if (!inherits(y, "SpatVector")) {
        stop("Input 'y' must be a SpatVector object.")
    }
    if (length(x) != length(y)) {
        stop("Number of tiles in 'x' must match number of polygons in 'y'.")
    }

    geomorphon_args <- list(...)

    if (is.null(cell_buffer)) {
        cell_buffer <- geomorphon_args$search + geomorphon_args$skip + 1
        if (is.null(cell_buffer)) {
            cell_buffer <- 0
        }
    }

    xs <- terra::sources(x)
    yw <- terra::geom(y, wkt = TRUE)
    ycrs <- terra::crs(y)

    processed_tiles_list <- LAPPLY.FUN(seq_len(length(x)), function(j) {

        input_source <- xs[j]
        input_tile <- terra::rast(input_source)

        processed_buffered_tile <- do.call(FUN, c(list(
            elevation = terra::extend(input_tile, cell_buffer)
        ), geomorphon_args))

        if (is.list(processed_buffered_tile)) {
            processed_buffered_tile <- processed_buffered_tile[[1]]
        }

        output_filename <- file.path(
            dirname(input_source),
            paste0("crop_", basename(input_source))
        )

        if( input_source == "") {
            output_filename <- NULL
        }

        cropped_result <- terra::crop(
            x = processed_buffered_tile,
            y = terra::vect(yw[j], crs = ycrs),
            filename = output_filename,
            overwrite = TRUE,
            datatype = datatype
        )
        terra::varnames(cropped_result) <- ""
        return(cropped_result)
    })

    terra::sprc(processed_tiles_list)
}

#' Combine Processed Tiles into a Single Output Raster
#'
#' Mosaics a collection of processed raster tiles into a single, seamless
#' SpatRaster. The final raster is cropped to the exact extent of the original
#' source DEM to ensure a perfect match.
#'
#' @param x `SpatRasterCollection` of the final result tiles, typically from
#'   [geomorphons_process_tiles()].
#' @param y `SpatExtent` object from the original, full DEM. Used to ensure
#'   the final output has the exact same extent. Default `NULL` uses full extent
#'   of `x`.
#' @param filename character. The path and filename for the final output raster.
#' @param overwrite logical. If `TRUE`, the final output file will be
#'   overwritten.
#'
#' @return A `SpatRaster` object pointing to the final, merged file on disk.
#' @noRd
.geomorphons_combine_tiles <- function(x, y = NULL, filename, datatype, overwrite = TRUE) {
    if (!inherits(x, "SpatRasterCollection")) {
        stop("Input 'x' must be a SpatRasterCollection.")
    }

    if (is.null(y)) {
        y <- terra::ext(x)
    }

    terra::merge(
        terra::crop(x, y),
        datatype = datatype,
        na.rm = TRUE,
        filename = filename,
        overwrite = overwrite
    )

}
