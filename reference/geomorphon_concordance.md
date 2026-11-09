# Calculate Structural Concordance Between Geomorphon and Categorical Maps

Measures the degree of structural agreement between a geomorphon
classification and another categorical map (e.g., soil types, vegetation
classes) using Normalized Mutual Information and boundary alignment
metrics.

## Usage

``` r
geomorphon_concordance(
  geo_rast,
  cat_rast,
  method = c("both", "nmi", "boundary"),
  threshold_cells = 2,
  spatial_null = TRUE
)
```

## Arguments

- geo_rast:

  `SpatRaster`. Geomorphon classification raster (typically output from
  [`geomorphons()`](http://humus.rocks/rgeomorphon/reference/geomorphons.md)).

- cat_rast:

  `SpatRaster`. Categorical comparison raster (e.g., soil map,
  vegetation map). Must have matching extent and resolution to
  `geo_rast`.

- method:

  Character. Which concordance metrics to calculate. One of `"nmi"`
  (Normalized Mutual Information only), `"boundary"` (boundary alignment
  only), or `"both"` (all metrics). Default: `"both"`.

- threshold_cells:

  Integer. Buffer distance in cells for boundary concordance
  calculation, converted to map units with the cell size. Categorical
  boundary cells within this distance of the nearest geomorphon boundary
  cell are counted as concordant. Must be a finite, non-negative whole
  number when `method` includes boundary metrics. Geographic
  (longitude/latitude) CRS inputs are rejected with an error when
  `method` includes boundary metrics because boundary distances require
  a planar Cartesian CRS. Default: `2`.

- spatial_null:

  Logical. Whether to calculate the spatial-null ensemble. If `FALSE`,
  skip it and return `NA_real_` null means/deltas and zero shift counts.
  Default: `TRUE`.

## Value

A list with the following components (depending on `method`):

- `nmi`:

  Numeric. Normalized Mutual Information (0-1). Higher values indicate
  stronger statistical association between classifications. Existing NMI
  behavior is unchanged. Only included when method includes "nmi".

- `homogeneity`:

  Numeric. `I(geo; cat) / H(cat)`, measuring how well the geomorphon
  classes determine the categorical reference. `NA_real_` when
  `H(cat) = 0`; otherwise zero for a constant geomorphon map against a
  varying categorical map. Included only for methods "nmi" and "both".

- `completeness`:

  Numeric. `I(geo; cat) / H(geo)`, measuring how completely the
  categorical classes capture geomorphon information. `NA_real_` when
  `H(geo) = 0`; otherwise zero for a constant categorical map against a
  varying geomorphon map. Included only for methods "nmi" and "both".

- `edge_offset`:

  Numeric. Mean distance (in map units) from each categorical boundary
  cell to the nearest geomorphon boundary cell. Lower values indicate
  better boundary alignment. Only included when method includes
  "boundary".

- `edge_concordance`:

  Numeric. Proportion (0-1) of categorical boundary cells within
  `threshold_cells` of the nearest geomorphon boundary cell. Higher
  values indicate better alignment. Only included when method includes
  "boundary".

- `edge_f1`:

  Numeric. Symmetric threshold-based boundary F1 (0-1). `NA_real_` if
  either shared-domain boundary set is empty. Included when method
  includes "boundary".

- `edge_mhd`:

  Numeric. Maximum of the two directed mean nearest-boundary distances,
  in map units; not classical maximum Hausdorff distance. `NA_real_` if
  either shared-domain boundary set is empty. Included when method
  includes "boundary".

- `edge_chamfer`:

  Numeric. Mean of the two directed mean nearest-boundary distances, in
  map units. `NA_real_` if either shared-domain boundary set is empty.
  Included when method includes "boundary".

- `n_cells`:

  Integer. Total number of valid cells analyzed.

- `n_edges_geo`:

  Integer. Number of geomorphon boundary cells. Only included when
  method includes "boundary".

- `n_edges_cat`:

  Integer. Number of categorical boundary cells. Only included when
  method includes "boundary".

- `contingency_table`:

  Matrix. Cross-tabulation of geomorphon vs. categorical classes. Only
  included when method includes "nmi".

- `mean_null_nmi`:

  Numeric. Arithmetic mean NMI across eligible deterministic spatial
  translations. `NA_real_` if no shifts have paired-valid cells.
  Included only for methods "nmi" and "both".

- `delta_nmi`:

  Numeric. Observed `nmi` minus `mean_null_nmi`; may be negative and is
  `NA_real_` if either score is `NA`. Included only for methods "nmi"
  and "both".

- `nmi_null_shifts`:

  Integer. Number of evaluated shifts eligible for the NMI null mean
  (shifts with paired-valid cells). Included only for methods "nmi" and
  "both".

- `mean_null_edge_f1`:

  Numeric. Arithmetic mean edge F1 across shifts with nonempty edge
  sets. `NA_real_` if no shifts are eligible. Included only for methods
  "boundary" and "both".

- `delta_edge_f1`:

  Numeric. Observed `edge_f1` minus `mean_null_edge_f1`; may be negative
  and is `NA_real_` if either score is `NA`. Included only for methods
  "boundary" and "both".

- `edge_f1_null_shifts`:

  Integer. Number of evaluated shifts eligible for the edge-F1 null mean
  (both edge sets nonempty). Included only for methods "boundary" and
  "both".

- `null_shifts_evaluated`:

  Integer. Number of non-identity shifts evaluated for the null
  ensemble, shared by the requested metrics.

## Note

**Memory Requirements:** This function performs edge detection on both
input rasters. For small to medium rasters (\< 1GB each), a fast C++
implementation is used. For larger rasters, the function automatically
falls back to
[`terra::boundaries()`](https://rspatial.github.io/terra/reference/boundaries.html)
which is slower but more memory-efficient.

Since both rasters are processed simultaneously, the total memory
requirement can approach 4x the size of a single raster in the worst
case. For rasters larger than available RAM:

- Process in tiles using
  [`terra::makeTiles()`](https://rspatial.github.io/terra/reference/makeTiles.html)

- Increase available system memory

- Consider splitting the analysis into smaller regions

## Information-Theoretic Concordance

Normalized Mutual Information (NMI) measures the mutual dependence
between two categorical variables, independent of the specific class
mappings. An NMI of 1 indicates perfect association, while 0 indicates
independence. When NMI is requested, the result also includes
directional `homogeneity` and `completeness`. With `cat_rast` as the
reference, `homogeneity = I(geo; cat) / H(cat)` and
`completeness = I(geo; cat) / H(geo)`. Swapping the rasters swaps these
two values but leaves NMI unchanged. When both entropies are positive,
their harmonic mean equals NMI. A component with a zero entropy
denominator is `NA_real_`; the other component is zero when defined. If
both maps are constant, both components are `NA_real_` while existing
NMI remains 0. These additive components do not change the existing NMI
calculation or values.

## Boundary Alignment

Boundary alignment measures physical proximity between geomorphon and
categorical boundaries. The `edge_offset` is the mean distance from each
categorical boundary cell to the nearest geomorphon boundary cell, in
map units of the input CRS. The `edge_concordance` is the proportion of
categorical boundary cells within a threshold distance of the nearest
geomorphon boundary cell. These legacy summaries retain their
directional definitions and per-input valid masks.

The additive `edge_f1`, `edge_mhd`, and `edge_chamfer` summaries are
symmetric. They detect boundaries for both maps only within their shared
non-`NA` cells. `edge_f1` is threshold-based: precision and recall are
the proportions of boundary cells in each map within the threshold of a
boundary in the other map, and F1 is their harmonic mean (range 0-1). It
is not a kernel F1 or a one-to-one matched-boundary score. `edge_mhd` is
the maximum of the two directed mean nearest-boundary distances; despite
its name, it is not the classical maximum Hausdorff distance.
`edge_chamfer` is the mean of those directed means. The distance
summaries are non-negative and use input CRS map units. Each new score
is `NA_real_` if either shared-domain boundary set is empty.

Shared-domain edges use the existing eight-neighbor rule and exclude
raster borders. Mask holes are not boundaries and do not block
nearest-boundary distances. `threshold_cells` is multiplied by the
x-resolution, including for rectangular cells; distance calculations use
both cell dimensions. When input NA masks differ, legacy counts and
directional summaries still use their existing per-input masks, so they
need not describe the edges used by the new shared-domain scores.

For rasters to be comparable, they must have matching row and column
counts, extents, and resolutions. Use
[`terra::crop()`](https://rspatial.github.io/terra/reference/crop.html),
[`terra::resample()`](https://rspatial.github.io/terra/reference/resample.html),
or
[`terra::disagg()`](https://rspatial.github.io/terra/reference/disaggregate.html)
to align rasters before calling this function.

The NMI is symmetric and normalized to the range 0-1, with higher values
indicating stronger association regardless of how classes map to each
other. Homogeneity and completeness are directional: `cat_rast` is the
reference for their definitions, so swapping inputs exchanges the
scores. If both marginal entropies are positive, their harmonic mean
equals NMI. A zero denominator yields `NA_real_` for that component; for
two constant maps both components are `NA_real_` and NMI remains 0.
These fields are returned only by `method = "nmi"` or `"both"`, not by
`"boundary"`. The boundary metrics are useful for validating whether
geomorphological transitions correspond to ecological or pedological
boundaries.

The additive spatial-null summaries compare the observed alignment with
a deterministic sample of up to 128 non-identity circular translations
of `cat_rast`; all candidate shifts are used when fewer than 128 exist.
The schedule is defined for a raster with `nr` rows and `nc` columns by
`M = nr * nc - 1`, `S = min(128, M)`, and
`k_i = floor((i + 0.5) * M / S)` for `i = 0, ..., S - 1`. For each
index, `q = k_i + 1` and the row/column shift is
`(dr, dc) = (q %/% nc, q %% nc)`. Using zero-based coordinates, the
source category at destination `(r, c)` is
`[(r - dr) mod nr, (c - dc) mod nc]`; the geomorphon map stays fixed.
Outputs contain aggregate summaries only, not per-shift scores. The
category values and their `NA` mask move together, and paired-valid
cells are recomputed for each shift. NMI uses base-2 entropy and is zero
when both marginal entropies are zero. Null edge F1 uses shared-valid
eight-neighbor edges (excluding image borders) and toroidal Euclidean
distances in cell coordinates with `threshold_cells`; a shift with
either edge set empty is ineligible for its mean. Null summaries are
arithmetic means over eligible shifts; deltas are observed minus mean
and are not clamped. With no eligible shifts, the corresponding summary
and delta are `NA_real_` and its eligible count is zero. These are
descriptive contrasts against this translation ensemble, not adjusted
mutual information, calibrated scores, p-values, or evidence of
statistical significance. NMI-null fields are returned only for
`method = "nmi"` or `"both"`; edge-F1-null fields only for
`method = "boundary"` or `"both"`. `null_shifts_evaluated` is always
returned. When `spatial_null = FALSE`, the ensemble is skipped, null
means and deltas are `NA_real_`, and shift counts are `0L`.

## See also

[`geomorphons()`](http://humus.rocks/rgeomorphon/reference/geomorphons.md),
[`terra::boundaries()`](https://rspatial.github.io/terra/reference/boundaries.html)

## Examples

``` r
library(terra)
library(rgeomorphon)

# Load example DEM
data("volcano", package = "datasets")
dem <- rast(volcano)
crs(dem) <- crs("EPSG:2193")
ext(dem) <- c(1756968, 1757578, 5917000, 5917870)

# Calculate geomorphons
geo <- geomorphons(dem, search = 7, simplify = TRUE)

# Create synthetic categorical map for demonstration
cat_map <- dem
values(cat_map) <- cut(values(dem), breaks = 5, labels = FALSE)

# Calculate concordance
conc <- geomorphon_concordance(geo, cat_map)
print(conc$nmi)
#> [1] 0.123542
print(conc$edge_concordance)
#> [1] 0.9707773
print(conc[c("edge_f1", "edge_mhd", "edge_chamfer")])
#> $edge_f1
#> [1] 0.7801114
#> 
#> $edge_mhd
#> [1] 26.19394
#> 
#> $edge_chamfer
#> [1] 14.6822
#> 
```
