# Reference Attribute Table Target Layer Generator

Derives aligned categorical raster layers from an integer ID raster and
its Reference Attribute Table (RAT), with support for categorical
columns and binned numeric columns. Attaches coverage statistics
describing valid cell proportions and class counts.

## Usage

``` r
geomorphon_targets(
  x,
  attrs = NULL,
  id_col = 1,
  vars = NULL,
  breaks = list(),
  labels = list(),
  right = TRUE
)
```

## Arguments

- x:

  A single-layer categorical or integer `SpatRaster`.

- attrs:

  Optional `data.frame` containing the attribute table. If `NULL`
  (default), the table is extracted from `terra::cats(x)[[1]]`.

- id_col:

  Column index or column name in `attrs` containing the keys that match
  raster values in `x`. Default is `1`.

- vars:

  Character vector of column names in `attrs` to convert into target
  raster layers. Default is all columns except `id_col`.

- breaks:

  Named list of numeric break vectors for binning continuous attribute
  columns using [`base::cut()`](https://rdrr.io/r/base/cut.html).

- labels:

  Named list of character vectors of category labels corresponding to
  `breaks`.

- right:

  Logical. Indicating if the intervals should be closed on the right
  (and open on the left) or vice versa. Default is `TRUE`. Passed to
  [`base::cut()`](https://rdrr.io/r/base/cut.html).

## Value

A multi-layer `SpatRaster` where each layer corresponds to a derived
target variable with categories set via
[`terra::levels()`](https://rspatial.github.io/terra/reference/factors.html).
The returned raster has an attached attribute `"coverage"`, a tidy
`data.frame` containing:

- target:

  Target layer name.

- n_total:

  Total number of cells.

- n_valid:

  Number of valid non-NA cells.

- n_na:

  Number of missing (NA) cells.

- prop_valid:

  Proportion of valid cells (`n_valid / n_total`).

- n_classes:

  Number of factor levels / distinct categories.

## Examples

``` r
library(terra)
r <- rast(matrix(c(1, 2, 1, 3), nrow = 2, ncol = 2))
rat <- data.frame(id = 1:3, landuse = c("Forest", "Water", "Forest"), slope = c(12, 1, 25))
targets <- geomorphon_targets(r, attrs = rat, id_col = "id", vars = "landuse")
attr(targets, "coverage")
#>    target n_total n_valid n_na prop_valid n_classes
#> 1 landuse       4       4    0          1         2
```
