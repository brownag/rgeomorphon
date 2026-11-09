---
id: "0007"
title: "Structural concordance conventions"
status: "Accepted"
date: 2026-11-01
authors:
  - "Andrew Gene Brown"
supersedes: null
superseded_by: null
source_context: "EXECUTION_PLAN_3, EXECUTION_PLAN_4, EXECUTION_PLAN_5"
---

# ADR 0007: Structural concordance conventions

## Context and Problem Statement

`rgeomorphon` provides `geomorphon_concordance()` to evaluate structural
agreement between a geomorphon classification (`geo_rast`) and an aligned
categorical raster (`cat_rast`, such as soil map units or ecological classes).
Evaluating categorical map correspondence requires complementary perspectives:
areal class association, physical boundary alignment, and baseline contrast
against spatial autocorrelation and boundary-density saturation.

Without an explicit record, entropy normalization, reference orientation for
directional components, edge-detection masks, map-unit distance scaling, and
spatial-null translation rules could drift between the R and C++
implementations. This ADR defines those conventions.

## Decision

### Input alignment

Require `geo_rast` and `cat_rast` to have compatible coordinate reference
systems (either both undefined or matching via `terra::same.crs()`), identical
extents, and matching resolutions within tolerance. Validate `threshold_cells`
as a finite, non-negative whole number and reject geographic (lon/lat) CRS
inputs with an error when boundary metrics are requested (`method = "boundary"`
or `"both"`).

### Information-theoretic concordance (`method = "nmi"` or `"both"`)

Compute Shannon entropy in base 2 over cells that are valid (non-`NA`) in both
rasters:

- `nmi = 2 * I(geo; cat) / (H(geo) + H(cat))`, clamped to `[0, 1]` for numerical
  precision, and `0` when `H(geo) + H(cat) == 0`.
- With `cat_rast` as the reference classification, report directional
  `homogeneity = I(geo; cat) / H(cat)` and
  `completeness = I(geo; cat) / H(geo)`.
- When both marginal entropies are positive, the harmonic mean of `homogeneity`
  and `completeness` equals `nmi`. A component with a zero entropy denominator
  is `NA_real_` (so two constant maps yield `nmi = 0` and `NA_real_` for both
  directional components).

### Boundary alignment (`method = "boundary"` or `"both"`)

Detect 8-connected interior boundary cells (where a valid cell differs from at
least one valid 8-connected neighbor; outer raster borders and `NA`/`NaN` cells
are set to `NA`) using the parallel C++ kernel `focal_edge_cpp()` when the
working set fits in half of available RAM, or `terra::boundaries()` with border
masking as a fallback.

Report two groups of boundary metrics in planar map units of the projected
input CRS (using `threshold_mapunits = threshold_cells * terra::res(geo_rast)[1]`):

1. **Directional (per-input valid masks):**
   - `n_edges_geo` and `n_edges_cat`: total edge cells detected on each input's
     own valid domain.
   - `edge_offset`: mean distance from each categorical boundary cell to the
     nearest geomorphon boundary cell.
   - `edge_concordance`: proportion of categorical boundary cells within
     `threshold_mapunits` of the nearest geomorphon boundary cell.
2. **Symmetric (shared non-`NA` domain):**
   - Detect edges for both rasters only within their shared non-`NA` mask so
     that asymmetric missing-data holes do not create one-sided boundaries.
   - `edge_f1`: harmonic mean of thresholded precision (categorical-to-geomorphon
     hit rate) and recall (geomorphon-to-categorical hit rate), or `0` when
     both hit rates are zero.
   - `edge_mhd`: maximum of the two directed mean nearest-boundary distances
     (modified Hausdorff summary, not classical maximum Hausdorff distance).
   - `edge_chamfer`: arithmetic mean of the two directed mean nearest-boundary
     distances.
   - Return `NA_real_` for `edge_f1`, `edge_mhd`, and `edge_chamfer` if either
     shared-domain boundary set is empty.

### Toroidal spatial-null ensemble (`spatial_null`)

When `spatial_null = TRUE` (the default), compare observed concordance against a
deterministic sample of at most 128 non-identity circular translations of
`cat_rast` evaluated in `spatial_null_cpp()`:

- For `M = nr * nc - 1` candidate shifts and `S = min(128, M)`, select 0-based
  candidate rank `k_i = floor((i + 0.5) * M / S)` using exact integer
  quotient/remainder arithmetic, map `q = k_i + 1` to
  `(dr, dc) = (q %/% nc, q %% nc)`, and translate `cat_rast` values and their
  `NA` mask together on the cell-index torus while keeping `geo_rast` fixed.
- When `M` and `S` are zero (such as for a one-cell raster), return zero shift
  counts and `NA_real_` null means without calculating candidate ranks or
  dividing by `S`.
- Recompute the paired-valid domain for each shift. Exclude shifts with zero
  paired-valid cells from `mean_null_nmi`, and shifts with an empty shared-domain
  edge set from `mean_null_edge_f1`.
- Compute null edge F1 using shared-valid 8-neighbor edges (excluding outer
  raster borders) and toroidal Euclidean distances in cell coordinates against
  `threshold_cells`.
- Report arithmetic means (`mean_null_nmi`, `mean_null_edge_f1`), un-clamped
  observed-minus-null deltas (`delta_nmi`, `delta_edge_f1`), and shift counts
  (`nmi_null_shifts`, `edge_f1_null_shifts`, `null_shifts_evaluated`) matching
  the requested `method`.
- When `spatial_null = FALSE`, skip the ensemble and return `NA_real_` for null
  means/deltas and `0L` for shift counts.

## Consequences

- `geomorphon_concordance()` is an independent structural-comparison tool for
  categorical maps; it is not part of the PTRM perceptual model (ADR 0006).
- Directional boundary summaries (`edge_offset`, `edge_concordance`) can
  saturate when one map has a much higher boundary density than the other;
  symmetric metrics (`edge_f1`, `edge_mhd`, `edge_chamfer`) and spatial-null
  deltas (`delta_nmi`, `delta_edge_f1`) provide diagnostic counterweights.
- Spatial-null summaries are descriptive contrasts against the fixed 128-shift
  toroidal translation ensemble, not p-values or formal significance tests.
- Changing the entropy base, reference orientation, edge connectivity, border
  exclusion, or spatial-null shift schedule requires updating this ADR and the
  concordance test suite.
