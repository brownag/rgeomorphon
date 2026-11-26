# Cell-Wise Divergence Between Landform Proportion Ensembles

Computes cell-wise divergence and summary statistics between two aligned
multi-layer class-proportion ensembles (e.g. from
[`geomorphon_ensemble()`](https://humus.rocks/rgeomorphon/reference/geomorphon_ensemble.md)
with `method = "proportions"`), or two lists of proportion matrices /
numeric matrices.

## Usage

``` r
geomorphon_divergence(p, q, method = "jensen-shannon", unit = "log2")
```

## Arguments

- p:

  A multi-layer class-proportion `SpatRaster` or numeric matrix / matrix
  list.

- q:

  A multi-layer class-proportion `SpatRaster` or numeric matrix / matrix
  list aligned with `p`.

- method:

  Character distance metric. Default is `"jensen-shannon"`.

- unit:

  Logarithm unit. Default is `"log2"`.

## Value

A list containing:

- raster:

  A single-layer `SpatRaster` named `"divergence"`, returned for raster
  input.

- matrix:

  The cell-vector divergence as a numeric matrix, returned for matrix or
  list input.

- mean_divergence:

  Arithmetic mean divergence across valid cells.

- median_divergence:

  Median divergence across valid cells.

- q90_divergence:

  90th percentile divergence across valid cells.

- n_cells:

  Count of valid cells compared.

## Examples

``` r
library(terra)
p <- rast(array(c(0.8, 0.2, 0.5, 0.5), dim = c(2, 2, 2)))
q <- rast(array(c(0.2, 0.8, 0.5, 0.5), dim = c(2, 2, 2)))
names(p) <- names(q) <- c("Flat", "Slope")
div <- geomorphon_divergence(p, q)
div$mean_divergence
#> [1] 0
```
