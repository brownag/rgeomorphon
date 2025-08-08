#' @keywords internal
#' @useDynLib rgeomorphon
#' @import RcppParallel
#' @importFrom Rcpp getRcppVersion
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' Bathymetric Information on California's Salton Sea
#'
#' @description
#'
#' Matrix derived from one foot contours of the Salton Sea floor. This data was
#' created with the vertical datum NGVD29 and NAD83 California Teale Albers
#' (EPSG:3110) projection. Each value in the matrix represents the elevation, in
#' meters, of a 300 m x 300 m cell. Cell values are interpolated using a thin
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
