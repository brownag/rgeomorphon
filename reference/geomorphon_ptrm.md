# Calculate Perceived Terrain Realism Metric (PTRM)

Calculates the Perceived Terrain Realism Metric based on the
distribution of geomorphon landform features. The PTRM provides a
quantitative estimate of human-perceived terrain realism by relating the
proportions of terrain features to their perceived realism.

## Usage

``` r
geomorphon_ptrm(
  x,
  mode = c("dem", "class"),
  search = 3,
  skip = 0,
  flat_angle_deg = 1,
  dist = 0,
  forms = "forms10",
  coefficients = NULL,
  intercept = -38.02,
  normalizer = 69.96,
  ...
)
```

## Arguments

- x:

  `SpatRaster`. Either a Digital Elevation Model (when `mode = "dem"`)
  or a geomorphon classification raster (when `mode = "class"`).

- mode:

  Character. One of `"dem"` or `"class"`. When `"dem"`, geomorphons are
  calculated from the input DEM. When `"class"`, the input is assumed to
  be a pre-computed geomorphon classification. Default: `"dem"`.

- search:

  Numeric. Search radius for geomorphon calculation (only used when
  `mode = "dem"`). Units depend on the CRS of the input. Default: `3`.
  This is a convenience default, not a value calibrated against the PTRM
  paper, so it should not be read as a validated choice.

- skip:

  Numeric. Skip radius for geomorphon calculation (only used when
  `mode = "dem"`). Default: `0`.

- flat_angle_deg:

  Numeric. Flatness angle threshold in degrees (only used when
  `mode = "dem"`). Default: `1.0`.

- dist:

  Numeric. Flatness distance (only used when `mode = "dem"`). Default:
  `0`.

- forms:

  Character. Form system to use: `"forms10"`, `"forms6"`, `"forms5"`, or
  `"forms4"`. Default: `"forms10"`.

- coefficients:

  Named numeric vector. Custom PTRM coefficients for each geomorphon
  class. If `NULL`, uses default values from Rajasekaran et al. (2022).
  Names should match geomorphon form names. A ten-form vector is reduced
  to the requested `forms` system. Default: `NULL`.

- intercept:

  Numeric. PTRM intercept. Default: `-38.02`.

- normalizer:

  Numeric. PTRM divisor. Must be finite and non-zero. Default: `69.96`.

- ...:

  Additional arguments passed to
  [`geomorphons()`](http://humus.rocks/rgeomorphon/reference/geomorphons.md)
  when `mode = "dem"`. Ignored when `mode = "class"`.

## Value

A list containing:

- `ptrm_score`:

  Numeric. PTRM value obtained by dividing the raw score by the
  normalizer. Scores are raw and may fall outside `[0, 1]`; higher
  values indicate more realistic terrain.

- `raw_score`:

  Numeric. Unnormalized score before division.

- `proportions`:

  Named numeric vector. Proportion of each observed geomorphon form in
  the input, named with the canonical system labels.

- `n_cells`:

  Integer. Total number of valid (non-NA) cells analyzed.

## Distance Calculation and Coordinate Reference Systems

The geomorphon classification algorithm assumes planar distances, so it
is strongly recommended that elevation data be in a projected coordinate
system.

## Tiled Processing

When `mode = "dem"`, geomorphons are calculated internally. For large
DEMs, this may trigger the automatic tiled processing workflow. See
[`geomorphons()`](http://humus.rocks/rgeomorphon/reference/geomorphons.md)
for details.

The PTRM is calculated using the formula:

\$\$PTRM = \frac{c_0 + \sum\_{i=1}^{n} c_i \times p_i}{D}\$\$

where \\c_0\\ is the intercept (-38.02 by default), \\c_i\\ are the
coefficients for each geomorphon form, \\p_i\\ are the proportions of
each form in the terrain, and \\D\\ is the divisor (69.96 by default).
The intercept, normalizer, and coefficients can all be overridden.

The returned score is not clamped or rescaled, so it may fall outside
`[0, 1]` depending on the input distribution and the model parameters.
Higher PTRM scores indicate more realistic terrain from a human
perception perspective. Rajasekaran et al. (2022) found that real
terrain from Arizona achieved a PTRM of 0.76, ranking in the 78th
percentile for perceived realism.

## References

Rajasekaran, S.D., Kang, H., Čadík, M., Galin, E., Guérin, E., Peytavie,
A., Slavík, P., & Mech, R. (2022). PTRM: Perceived Terrain Realism
Metric. ACM Transactions on Applied Perception, 19(2), Article 5.
[doi:10.1145/3514244](https://doi.org/10.1145/3514244)

Jasiewicz, J., Stepinski, T. (2013). Geomorphons - a pattern recognition
approach to classification and mapping of landforms. Geomorphology, 182,
147-156.
[doi:10.1016/j.geomorph.2012.11.005](https://doi.org/10.1016/j.geomorph.2012.11.005)

## See also

[`geomorphons()`](http://humus.rocks/rgeomorphon/reference/geomorphons.md),
[`geomorphon_categories()`](http://humus.rocks/rgeomorphon/reference/geomorphon_theme.md)

## Examples

``` r
library(terra)
library(rgeomorphon)

# Load example DEM
data("volcano", package = "datasets")
dem <- rast(volcano)
crs(dem) <- crs("EPSG:2193")
ext(dem) <- c(1756968, 1757578, 5917000, 5917870)

# Calculate PTRM from DEM
ptrm1 <- geomorphon_ptrm(dem, mode = "dem", search = 7)
print(ptrm1$ptrm_score)
#> [1] -0.2654868

# Or calculate geomorphons first, then PTRM
geo <- geomorphons(dem, search = 7, simplify = TRUE)
ptrm2 <- geomorphon_ptrm(geo, mode = "class")

# Results should be identical
all.equal(ptrm1$ptrm_score, ptrm2$ptrm_score)
#> [1] TRUE
```
