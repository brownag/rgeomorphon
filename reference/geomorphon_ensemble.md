# Ensemble Aggregation of Geomorphon Classifications

Computes cell-wise ensemble summaries across parameter sweep runs,
supporting deterministic modal classification and per-class proportion
layers.

## Usage

``` r
geomorphon_ensemble(
  sweep_obj,
  method = c("modal", "proportions"),
  form = NULL,
  classes = NULL,
  ties = "lowest",
  ...
)
```

## Arguments

- sweep_obj:

  An object of class `geomorphon_sweep`, a list of `SpatRaster` objects,
  or a multi-layer `SpatRaster` representing runs of a parameter sweep.

- method:

  character. Aggregation method: `"modal"` (default) computes the modal
  landform class per cell; `"proportions"` computes per-class
  proportions across non-NA runs.

- form:

  character or numeric. Which form taxonomy to aggregate when multiple
  exist in `sweep_obj` (e.g. `"forms10"`, `"forms6"`, `"forms5"`,
  `"forms4"`). Defaults to the first or only form available.

- classes:

  integer or character vector. Which landform classes to compute
  proportions for when `method = "proportions"`. Defaults to `NULL`,
  which evaluates all reachable classes for `form` (10 for `forms10`, 6
  for `forms6`, 5 for `forms5`, 4 for `forms4`), matching
  [`geomorphon_categories()`](https://humus.rocks/rgeomorphon/reference/geomorphon_theme.md).

- ties:

  character. Tie-breaking rule for `method = "modal"`. Passed to
  [`terra::modal()`](https://rspatial.github.io/terra/reference/modal.html),
  default is `"lowest"` for deterministic tie-breaking.

- ...:

  Additional arguments passed to methods.

## Value

A `SpatRaster` (or `matrix` / list of matrices) containing the
aggregated ensemble. For `method = "modal"`, returns a single layer with
geomorphon factor levels and color table applied. For
`method = "proportions"`, returns multi-layer proportions named by
class.

## Examples

``` r
if (FALSE) { # \dontrun{
library(terra)
data("salton", package = "rgeomorphon")
dem <- rast(salton)
crs(dem) <- attr(salton, "crs")
ext(dem) <- attr(salton, "extent")

grid <- data.frame(search = c(3, 5), skip = c(0, 1))
sw <- geomorphon_sweep(dem, param_grid = grid, forms = "forms10")

ens_mode <- geomorphon_ensemble(sw, method = "modal")
ens_prop <- geomorphon_ensemble(sw, method = "proportions")
} # }
```
