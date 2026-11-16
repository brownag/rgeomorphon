# Parameter Sweep for Geomorphon Terrain Classification

Evaluates a grid of parameter combinations over a Digital Elevation
Model (DEM), taking advantage of line-of-sight ray-tracing count reuse
across multiple form taxonomies in each evaluation. Records execution
runtime and parameter metadata in a tidy summary table.

## Usage

``` r
geomorphon_sweep(
  elevation,
  param_grid,
  forms = c("forms10", "forms6", "forms5", "forms4"),
  ...
)
```

## Arguments

- elevation:

  A `SpatRaster` or numeric `matrix` containing elevation data.

- param_grid:

  A `data.frame` or `list` containing parameter combinations. May
  contain columns or list elements for `search`, `skip`,
  `flat_angle_deg` (or `flat`), `dist`, `comparison_mode`, and `tdist`.
  If a list of vectors is supplied, it is expanded via
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html). Any
  combinations where `skip >= search` are pruned with a warning.

- forms:

  character or numeric vector. Landform taxonomies to compute for each
  parameter combination. Default is
  `c("forms10", "forms6", "forms5", "forms4")`.

- ...:

  Additional arguments passed to
  [`geomorphons()`](http://humus.rocks/rgeomorphon/reference/geomorphons.md).

## Value

An S3 object of class `geomorphon_sweep` containing:

- summary:

  A tidy `data.frame` with columns `run_id`, parameter values,
  `resolution`, `elapsed_sec` (runtime in seconds), and `n_cells`.

- rasters:

  A named list of classification outputs (one per run_id).

## Examples

``` r
if (FALSE) { # \dontrun{
library(terra)
data("salton", package = "rgeomorphon")
dem <- rast(salton)
crs(dem) <- attr(salton, "crs")
ext(dem) <- attr(salton, "extent")

grid <- list(
  search = c(3, 5),
  skip = c(0, 1),
  flat_angle_deg = 1.0
)

sw <- geomorphon_sweep(dem, param_grid = grid, forms = c("forms10", "forms6"))
print(sw)
summary(sw)
} # }
```
