# Calculate Focal (Moving Window) PTRM

Calculates Perceived Terrain Realism Metric (PTRM) for each pixel using
a circular moving window. Produces a spatially-varying measure of local
terrain complexity and realism.

## Usage

``` r
geomorphon_ptrm_focal(
  x,
  radius = 3,
  mode = c("dem", "class"),
  coefficients = NULL,
  forms = "forms10",
  intercept = -38.02,
  normalizer = 69.96,
  search = 3,
  skip = 0,
  flat_angle_deg = 1,
  dist = 0,
  filename = NULL,
  overwrite = FALSE,
  ...
)
```

## Arguments

- x:

  `SpatRaster`. Digital Elevation Model (when `mode = "dem"`) or
  pre-computed geomorphon classification (when `mode = "class"`).

- radius:

  Numeric. Integer window radius in cells (not map units). Default: 3.
  Larger values compute PTRM from larger neighborhoods.

- mode:

  Character. Input type: `"dem"` to calculate geomorphons internally, or
  `"class"` if input is pre-computed geomorphon classification. Default:
  `"dem"`.

- coefficients:

  Named numeric vector. Custom PTRM coefficients, either named for the
  requested form system (length k) or a named ten-form vector, which is
  reduced to the requested system. If `NULL`, uses the paper weights
  from Rajasekaran et al. (2022), derived for reduced form systems.
  Default: `NULL`.

- forms:

  Character. Form system to use: `"forms10"`, `"forms6"`, `"forms5"`, or
  `"forms4"`. Default: `"forms10"`.

- intercept:

  Numeric. PTRM intercept. Default: `-38.02`.

- normalizer:

  Numeric. PTRM normalizer; must be finite and non-zero. Default:
  `69.96`.

- search:

  Numeric. Search radius for geomorphon calculation (only used when
  `mode = "dem"`). Default: 3.

- skip:

  Numeric. Skip radius for geomorphon calculation (only used when
  `mode = "dem"`). Default: 0.

- flat_angle_deg:

  Numeric. Flatness angle threshold in degrees (only used when
  `mode = "dem"`). Default: 1.0.

- dist:

  Numeric. Flatness distance parameter (only used when `mode = "dem"`).
  Default: 0.

- filename:

  Character. Output filename. If provided, the result is written to disk
  instead of returned in memory. Default: `NULL`.

- overwrite:

  Logical. Overwrite existing file? Default: `FALSE`.

- ...:

  Additional arguments passed to
  [`geomorphons()`](http://humus.rocks/rgeomorphon/reference/geomorphons.md)
  when `mode = "dem"`.

## Value

A `SpatRaster` with raw focal PTRM scores. Values are not normalized to
`[0, 1]` and may fall outside that range. Border cells are NA.

## Details

The PTRM score measures local terrain realism on a continuous scale.
Higher values indicate more realistic terrain patterns according to the
perception model of Rajasekaran et al. (2022). The focal version applies
this metric locally to map spatial variation in terrain complexity.

**Computational Complexity:**

- Geomorphons: O(N × 8R) where R = search radius

- Focal PTRM: O(N × W) where W ≈ π × radius²

- Total: O(N × (8R + W))

For a 1000×1000 DEM with radius=5, search=7:

- Expected runtime: 3-6 seconds on modern multi-core system

- Memory: ~20 MB

## Algorithm

Uses a two-pass approach for efficiency:

1.  **Calculate geomorphons globally** (if `mode = "dem"`)

    - Applies the full geomorphon classification algorithm once

    - Produces classification raster with classes 1-10

2.  **Apply focal PTRM calculation** to each pixel's neighborhood

    - For each interior pixel, extracts geomorphon classes within
      circular window

    - Calculates proportions of each geomorphon form

    - Applies PTRM formula using empirical coefficients

    - Returns spatially-varying PTRM scores

This two-pass approach is much more efficient than recalculating
geomorphons for each window (~38x speedup for typical parameters).

## Window Type

Uses a circular window based on Euclidean distance. For an integer
`radius = r`, includes all cells where `sqrt(dx^2 + dy^2) <= r`.

## Border Handling

Cells within `radius` distance of raster edges are set to NA, as their
neighborhoods are incomplete. This prevents edge artifacts in the
output.

## Memory Requirements

For small/medium rasters, uses the fast C++ implementation requiring
approximately 4 bytes per input cell plus 8 bytes per output cell (about
12 bytes per cell) in RAM. For large rasters that exceed half of the
available memory, automatically falls back to block-processed
[`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html),
which is slower but avoids full-raster R matrices. Without `filename`,
terra-backed temporary storage is used for intermediate results.

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

[`geomorphon_ptrm()`](http://humus.rocks/rgeomorphon/reference/geomorphon_ptrm.md),
[`geomorphons()`](http://humus.rocks/rgeomorphon/reference/geomorphons.md),
[`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)

## Examples

``` r
library(terra)
library(rgeomorphon)

# Load example DEM
data("volcano", package = "datasets")
dem <- rast(volcano)
crs(dem) <- crs("EPSG:2193")
ext(dem) <- c(1756968, 1757578, 5917000, 5917870)

# Calculate focal PTRM with 5-cell radius
ptrm_focal <- geomorphon_ptrm_focal(dem, radius = 5, mode = "dem")
#> Warning: [subst] all 'from' values are missing, returning a copy of 'x'

# Compare different window sizes
ptrm_r3 <- geomorphon_ptrm_focal(dem, radius = 3)
#> Warning: [subst] all 'from' values are missing, returning a copy of 'x'
ptrm_r7 <- geomorphon_ptrm_focal(dem, radius = 7)
#> Warning: [subst] all 'from' values are missing, returning a copy of 'x'
```
