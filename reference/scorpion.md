# Elevation and SSURGO Soil Map Units for an Erosional Remnant in Calaveras County, California

Co-registered 10 m Digital Elevation Model (DEM) and Soil Survey
Geographic Database (SSURGO) soil map unit grid for an old erosional
remnant landform that resembles a scorpion in Calaveras County,
California, projected in NAD83 / Conus Albers (`EPSG:5070`). Elevation
values (in meters) are aggregated from USGS 3DEP 1 m LiDAR data, and
soil map unit keys (`mukey`) are rasterized from SSURGO map unit
polygons onto the same 235 x 339 grid.

## Usage

``` r
scorpion
```

## Format

A named `list` with 3 elements and spatial attributes `"crs"` (WKT2019
string for `EPSG:5070`) and `"extent"` (numeric vector of length 4:
`xmin`, `xmax`, `ymin`, `ymax`):

- `dem`:

  Numeric matrix (235 x 339) of 10 m elevation values in meters.

- `mukey`:

  Integer matrix (235 x 339) of SSURGO map unit keys (`mukey`).

- `mapunit`:

  Data frame of SSURGO `mapunit` and `muaggatt` attributes with columns
  `mukey` (map unit key, used as the raster category ID), `musym` (map
  unit symbol), `muname` (map unit name), `mukind` (map unit kind),
  `farmlndcl` (farmland classification), `slopegraddcp` (representative
  slope gradient, dominant component, %), `slopegradwta` (weighted
  average slope gradient, %), `drclassdcd` (dominant drainage class),
  `hydgrpdcd` (dominant hydrologic soil group), and `aws0100wta`
  (available water storage 0-100 cm, cm). Suitable for assigning
  directly to a `SpatRaster` via `levels(r) <- scorpion$mapunit`.

## Source

U.S. Geological Survey (USGS) 3D Elevation Program (3DEP) 1 m Digital
Elevation Model, accessed via Google Earth Engine (`USGS/3DEP/1m`).

USDA Natural Resources Conservation Service (NRCS) Soil Survey
Geographic Database (SSURGO), accessed via Soil Data Access
(`soilDB::SDA_spatialQuery()` and `soilDB::SDA_query()`).

## Examples

``` r

data("scorpion", package = "rgeomorphon")
str(scorpion)
#> List of 3
#>  $ dem    : num [1:235, 1:339] 94.3 89.6 84.9 81.6 79.8 ...
#>  $ mukey  : int [1:235, 1:339] 2600529 2600529 2600529 2600529 2600529 2766838 2766838 2766838 2766838 2766838 ...
#>  $ mapunit:'data.frame': 7 obs. of  10 variables:
#>   ..$ mukey       : int [1:7] 2462630 2600529 2766836 2766838 2766850 2924880 2924908
#>   ..$ musym       : chr [1:7] "W" "5013" "5201" "5012" ...
#>   ..$ muname      : chr [1:7] "Water" "Miltonhills-Amador complex, 15 to 45 percent slopes" "Pardee-Amador complex, 1 to 15 percent slopes" "Amador sandy loam, 2 to 15 percent slopes" ...
#>   ..$ mukind      : chr [1:7] "Consociation" "Complex" "Consociation" "Consociation" ...
#>   ..$ farmlndcl   : chr [1:7] "Not prime farmland" "Not prime farmland" "Not prime farmland" "Not prime farmland" ...
#>   ..$ slopegraddcp: num [1:7] 0 40 4 7 2 33 1
#>   ..$ slopegradwta: num [1:7] 0 35.6 5.4 6.1 2.8 29.7 1.1
#>   ..$ drclassdcd  : chr [1:7] NA "Well drained" "Well drained" "Well drained" ...
#>   ..$ hydgrpdcd   : chr [1:7] "D" "C" "D" "D" ...
#>   ..$ aws0100wta  : num [1:7] NA 7.14 5.04 4.66 3.42 ...
#>  - attr(*, "crs")= chr "PROJCRS[\"NAD83 / Conus Albers\",\n    BASEGEOGCRS[\"NAD83\",\n        DATUM[\"North American Datum 1983\",\n  "| __truncated__
#>  - attr(*, "extent")= Named num [1:4] -2141632 -2138242 1962844 1965194
#>   ..- attr(*, "names")= chr [1:4] "xmin" "xmax" "ymin" "ymax"

# reconstruct elevation and soil map unit SpatRaster objects
dem <- terra::rast(scorpion$dem)
terra::crs(dem) <- attr(scorpion, "crs")
terra::ext(dem) <- attr(scorpion, "extent")
names(dem) <- "elevation"

soil <- terra::rast(scorpion$mukey)
terra::crs(soil) <- attr(scorpion, "crs")
terra::ext(soil) <- attr(scorpion, "extent")
levels(soil) <- scorpion$mapunit

c(dem, soil)
#> class       : SpatRaster
#> size        : 235, 339, 2  (nrow, ncol, nlyr)
#> resolution  : 10, 10  (x, y)
#> extent      : -2141632, -2138242, 1962844, 1965194  (xmin, xmax, ymin, ymax)
#> coord. ref. : NAD83 / Conus Albers (EPSG:5070)
#> source(s)   : memory
#> names       : elevation, musym
#> min values  :     65.21,     W
#> max values  :    119.51,   301
```
