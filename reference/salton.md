# Bathymetric Information on California's Salton Sea

Matrix derived from one foot contours of the Salton Sea floor. This data
was created with the vertical datum NGVD29 and NAD83 California Teale
Albers (EPSG:3110) projection. Each value in the matrix represents the
elevation, in meters, of a 300 m x 300 m cell. Cell values are
interpolated using a thin plate spline fit to an exhaustive sample of
contour line vertices.

## Usage

``` r
salton
```

## Format

matrix, with cells representing X, Y grid locations, and attributes
`"crs"` (containing WKT2019 string with coordinate reference system
information) and `"extent"` (named numeric of length 4, containing xmin,
xmax, ymin, ymax)

## Source

California Division of Fish and Wildlife. 2007. Bathymetric Contours (1
foot) - Salton Sea (ds426). Available online:
<https://map.dfg.ca.gov/metadata/ds0426.html>

## Examples

``` r
str(salton)
#>  num [1:161, 1:165] NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ...
#>  - attr(*, "crs")= chr "PROJCRS[\"NAD83 / California Albers\",\n    BASEGEOGCRS[\"NAD83\",\n        DATUM[\"North American Datum 1983\""| __truncated__
#>  - attr(*, "extent")= Named num [1:4] 363300 412800 -538100 -489800
#>   ..- attr(*, "names")= chr [1:4] "xmin" "xmax" "ymin" "ymax"

# construct and georeference a SpatRaster object
dem <- terra::rast(salton)
terra::crs(dem) <- attr(salton, "crs")
terra::ext(dem) <- attr(salton, "extent")
names(dem) <- "Elevation (feet)"

dem
#> class       : SpatRaster 
#> size        : 161, 165, 1  (nrow, ncol, nlyr)
#> resolution  : 300, 300  (x, y)
#> extent      : 363300, 412800, -538100, -489800  (xmin, xmax, ymin, ymax)
#> coord. ref. : NAD83 / California Albers (EPSG:3310) 
#> source(s)   : memory
#> name        : Elevation (feet) 
#> min value   :        -278.1904 
#> max value   :        -216.2552 
```
