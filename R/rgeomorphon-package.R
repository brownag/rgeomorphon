#' @keywords internal
#' @useDynLib rgeomorphon
#' @import RcppParallel
#' @importFrom Rcpp getRcppVersion
#' @seealso [geomorphons()]
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' Bathymetric Information on California's Salton Sea
#'
#' @description
#'
#' Matrix derived from one foot contours of the Salton Sea floor. These data were
#' created with the vertical datum NGVD29 and NAD83 California Teale Albers
#' (`EPSG:3110`) projection. Each value in the matrix represents the elevation, in
#' feet, of a 300 m x 300 m cell. Cell values are interpolated using a thin
#' plate spline fit to an exhaustive sample of contour line vertices.
#'
#' @format matrix, with cells representing X, Y grid locations, and attributes
#'   `"crs"` (containing WKT2019 string with coordinate reference system
#'   information) and `"extent"` (named numeric of length 4, containing xmin,
#'   xmax, ymin, ymax)
#'
#' @source California Division of Fish and Wildlife. 2007. Bathymetric Contours
#'   (1 foot) - Salton Sea (ds426). Available online:
#'   <https://map.dfg.ca.gov/metadata/ds0426.html>
#'
#' @examplesIf requireNamespace("terra")
#'
#' str(salton)
#'
#' # construct and georeference a SpatRaster object
#' dem <- terra::rast(salton)
#' terra::crs(dem) <- attr(salton, "crs")
#' terra::ext(dem) <- attr(salton, "extent")
#' names(dem) <- "Elevation (feet)"
#'
#' dem
"salton"

#' Elevation and SSURGO Soil Map Units for an Erosional Remnant in Calaveras County, California
#'
#' @description
#'
#' Co-registered 10 m Digital Elevation Model (DEM) and Soil Survey Geographic
#' Database (SSURGO) soil map unit grid for an old erosional remnant landform
#' that resembles a scorpion in Calaveras County, California, projected in
#' NAD83 / Conus Albers (`EPSG:5070`). Elevation values (in meters) are
#' aggregated from USGS 3DEP 1 m LiDAR data, and soil map unit keys (`mukey`)
#' are rasterized from SSURGO map unit polygons onto the same 235 x 339 grid.
#'
#' @format A named `list` with 3 elements and spatial attributes `"crs"`
#'   (WKT2019 string for `EPSG:5070`) and `"extent"` (numeric vector of length 4:
#'   `xmin`, `xmax`, `ymin`, `ymax`):
#'   \describe{
#'     \item{`dem`}{Numeric matrix (235 x 339) of 10 m elevation values in meters.}
#'     \item{`mukey`}{Integer matrix (235 x 339) of SSURGO map unit keys (`mukey`).}
#'     \item{`mapunit`}{Data frame of SSURGO `mapunit` and `muaggatt` attributes
#'       with columns `mukey` (map unit key, used as the raster category ID),
#'       `musym` (map unit symbol), `muname` (map unit name), `mukind` (map unit
#'       kind), `farmlndcl` (farmland classification), `slopegraddcp`
#'       (representative slope gradient, dominant component, %), `slopegradwta`
#'       (weighted average slope gradient, %), `drclassdcd` (dominant drainage
#'       class), `hydgrpdcd` (dominant hydrologic soil group), and `aws0100wta`
#'       (available water storage 0-100 cm, cm). Suitable for assigning directly
#'       to a `SpatRaster` via `levels(r) <- scorpion$mapunit`.}
#'   }
#'
#' @source
#' U.S. Geological Survey (USGS) 3D Elevation Program (3DEP) 1 m Digital
#' Elevation Model, accessed via Google Earth Engine (`USGS/3DEP/1m`).
#'
#' USDA Natural Resources Conservation Service (NRCS) Soil Survey Geographic
#' Database (SSURGO), accessed via Soil Data Access (`soilDB::SDA_spatialQuery()`
#' and `soilDB::SDA_query()`).
#'
#' @examplesIf requireNamespace("terra")
#'
#' data("scorpion", package = "rgeomorphon")
#' str(scorpion)
#'
#' # reconstruct elevation and soil map unit SpatRaster objects
#' dem <- terra::rast(scorpion$dem)
#' terra::crs(dem) <- attr(scorpion, "crs")
#' terra::ext(dem) <- attr(scorpion, "extent")
#' names(dem) <- "elevation"
#'
#' soil <- terra::rast(scorpion$mukey)
#' terra::crs(soil) <- attr(scorpion, "crs")
#' terra::ext(soil) <- attr(scorpion, "extent")
#' levels(soil) <- scorpion$mapunit
#'
#' c(dem, soil)
"scorpion"
