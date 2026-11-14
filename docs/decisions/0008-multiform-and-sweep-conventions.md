---
id: "0008"
title: "Multi-form count reuse, fast classification, and parameter sweep conventions"
status: "Accepted"
date: 2026-11-08
authors:
  - "Andrew Gene Brown"
supersedes: null
superseded_by: null
source_context: "EXECUTION_PLAN_9, RESEARCH_6_FORM_ENSEMBLES"
---

# ADR 0008: Multi-form count reuse, fast classification, and parameter sweep conventions

## Context and Problem Statement

Geomorphon classification requires ray-tracing line-of-sight profiles across
8 principal directions from every valid cell up to an outer search radius
`search` (`L`). In `rgeomorphon`, this ray-tracing stage represents >99% of CPU
execution time.

The standard geomorphon algorithm classifies terrain into 10 fundamental
landform classes (`forms10`; Jasiewicz & Stepinski, 2013). Reduced form systems
(`forms6`, `forms5`, `forms4`; Masetti et al., 2018) map the exact same
`(negative, positive)` overlook count tuple `(nn, np)` onto simplified
geomorphic taxonomies via fixed 9x9 lookup tables.

Previously:

1. Evaluating multiple form taxonomies required calling `geomorphons()` repeatedly
   with different `forms` arguments, repeating identical ray-tracing calculations
   from scratch for each taxonomy at $O(K \cdot N)$ compute cost.
2. Reclassifying exported `positive` and `negative` overlook rasters via
   `forms_matrix_apply()` relied on `terra::subst` and table lookups in R,
   creating a substantial bottleneck for large datasets and bathymetric workflows.
3. Multi-scale terrain modeling requires sweeping parameter grids (varying
   `search`, `skip`, and `flat_angle_deg`) and aggregating runs into multi-scale
   ensembles (modal classifications and class proportions per cell; Issue #6).
4. Direct aggregation with `terra::modal()` drops categorical factor levels and
   color tables, exhibits version-dependent tie-breaking, and returns `NaN` for
   tied or missing cells.

This ADR establishes architectural conventions for single-pass multi-form count
reuse in C++, fast C++ overlook reclassification dispatch, vector `forms` API
contracts, and parameter sweep / ensemble aggregation schemas.

## Decision

### Multi-form count reuse architecture

Perform line-of-sight ray tracing exactly once per cell geometry. During ray
tracing in `GeomorphonWorker` (`src/geomorphon.cpp`), compute the negative and
positive overlook counts `(nn, np)` ($0 \le nn, np \le 8$).

When multiple form taxonomies are requested, populate each requested output array
(`forms10`, `forms6`, `forms5`, `forms4`) directly from the computed `(nn, np)`
counts within the same parallel cell iteration by referencing the constant 9x9
lookup tables (`forms_table10`, `forms_table6`, `forms_table5`, `forms_table4`).

Evaluating 4 taxonomies in a single C++ ray-tracing pass incurs little additional
kernel overhead ($< 1.4\times$ of a single-form kernel call, enforced by the
gated benchmark in `test_sweep.R`). End to end, a 4-form `geomorphons()` pass
measured about 45% less time than 4 separate single-form passes on `salton`
(search 30 and 40; multi-form pass about $0.55\times$ the separate-pass time).
This is well short of the ~75% estimate in the original goal; the earlier figure
is retired. Pixel-for-pixel identity across all four taxonomies (`forms10`, `forms6`,
`forms5`, `forms4`) to individual single-form evaluations is verified across both
`SpatRaster` and `matrix` inputs (closing D4).

### Vector `forms` API contract and backward compatibility

Extend `geomorphons()` and `.geomorphons()` to accept scalar or vector `forms`:

- Accept character tokens (`"forms10"`, `"forms6"`, `"forms5"`, `"forms4"`),
  integer codes (`10`, `6`, `5`, `4`), `TRUE` (canonical 10-form), `FALSE`,
  and `NULL`.
- Normalize requests to unique integer arrays (`10L`, `6L`, `5L`, `4L`).
- **Single-form backward compatibility:**
  - When `forms = TRUE`, `forms = "forms10"`, or any single form is requested,
    maintain existing return structures: return a 1-layer `SpatRaster` named
    `"form"` styled with `geomorphon_theme()`.
  - For numeric matrix inputs, return a single character factor matrix (when
    `simplify = TRUE`) or a named list `list(forms = ...)` (when `simplify = FALSE`).
- **Multi-form return structure:**
  - For `SpatRaster` inputs with multiple requested forms, return a multi-layer
    `SpatRaster` with canonical layer names (`forms10`, `forms6`, `forms5`,
    `forms4`), applying `geomorphon_theme()` across all layers.
  - For numeric matrix inputs with multiple requested forms, return a named list
    containing each requested taxonomy (`forms10`, `forms6`, etc.) plus a
    backward-compatible `"forms"` element pointing to the primary requested form.
  - Additional layers (`ternary`, `positive`, `negative`) are appended cleanly
    as extra layers.
- **Tiled execution:**
  - Ensure `.geomorphons_tiled()` and `.geomorphons_combine_tiles()` propagate
    multi-layer form stacks across spatial chunks, restoring layer names,
    factor categories, and color tables on the mosaicked raster.

### Fast reclassification dispatch in `forms_matrix_apply()`

Refactor `forms_matrix_apply()` in `R/forms_matrix.R`:

- Inspect the reclassification matrix `rcl`. When `rcl` matches any of the
  standard 9x9 forms matrices (`10`, `6`, `5`, `4`), bypass `terra::subst` and
  dispatch directly to the parallel C++ kernel `forms_matrix_apply_cpp()`.
  The kernel stays internal (not exported): users reach it through
  `forms_matrix_apply()`. Export it only if a user needs standalone
  `(positive, negative)` reclassification.
- Handle missing values (`NA` / `NaN`) by propagating `NA` cleanly without
  generating warnings.
- Support both `SpatRaster` inputs (converting to matrix buffers, calling C++,
  and reconstructing themed `SpatRaster` objects) and `matrix` / `list` inputs.
- When `rcl` is a custom non-standard matrix, fall back to the legacy
  expansion formula, preserving full backward compatibility for user-defined
  matrices.

### Parameter sweep runner (`geomorphon_sweep()`)

Implement `geomorphon_sweep()` in `R/sweep.R`:

- Accept `param_grid` as a `data.frame` or `list` (expanded via [expand.grid()]).
- Supply standard defaults: `search = 3`, `skip = 0`, `flat_angle_deg = 1.0`,
  `dist = 0.0`, `comparison_mode = "anglev1"`, `tdist = 0.0`.
- Prune invalid combinations where `skip >= search` with an informative warning;
  raise an error if no valid combinations remain.
- For each unique parameter tuple, execute a single ray-tracing evaluation across
  all requested `forms` simultaneously.
- Record execution runtime (`elapsed_sec`), raster resolution, and cell count.
- Return an S3 object of class `geomorphon_sweep` containing:
  - `summary`: tidy `data.frame` (`run_id`, parameters, `resolution`, `elapsed_sec`, `n_cells`).
  - `rasters`: named list of outputs (`run_1`, `run_2`, ...).
- Provide concise S3 `print()` and `summary()` methods.

### Ensemble aggregation conventions (`geomorphon_ensemble()`)

Implement `geomorphon_ensemble()` in `R/sweep.R` (resolving Issue #6):

- Support `geomorphon_sweep` objects, lists of rasters or matrices, and
  multi-layer `SpatRaster` stacks.
- **`method = "modal"`:**
  - Compute the cell-wise mode across ensemble runs using
    `terra::modal(..., ties = ties, na.rm = TRUE)`.
  - Default to `ties = "lowest"` for deterministic, version-independent
    tie-breaking across terra releases.
  - Mask cells with zero valid runs explicitly to `NA`, eliminating `NaN`
    artifacts.
  - Apply `geomorphon_theme()`, restoring categorical factor levels and
    color tables.
- **`method = "proportions"`:**
  - Compute per-class proportion layers $p_k = n_k / n_{\text{valid}}$ across
    valid non-NA runs.
  - Mask cells with zero valid runs to `NA_real_`.
  - Name layers with standard category labels (`Flat`, `Peak`, etc.) matching
    [geomorphon_categories()].
  - Guarantee valid cells sum to $1.0$.
  - Default reachable class subsets for reduced forms:
    - `forms10`: 1 to 10 (10 layers).
    - `forms6`: 1, 3, 4, 6, 8, 9 (6 layers: Flat, Ridge, Shoulder, Slope, Footslope, Valley).
    - `forms5`: 1, 2, 3, 6, 9 (5 layers: Flat, Peak, Ridge, Slope, Valley).
    - `forms4`: 1, 3, 6, 9 (4 layers: Flat, Ridge, Slope, Valley).
- Support full matrix workflow when `terra` is unavailable.

## Consequences

- Multi-taxonomy terrain analysis is faster when evaluating all 4 taxonomies in
  one pass: about 45% less end-to-end time than four separate passes on `salton`
  (measured; see Decision above).
- Legacy workflows and tests continue to work without modification; single-form
  calls return identical 1-layer rasters named `"form"`.
- Reclassification of `(positive, negative)` overlook rasters is 10–50x faster
  via C++ dispatch.
- Multi-run sweeps and ensemble aggregations provide reproducible, deterministic
  multi-scale mapping with explicit tie resolution and valid proportion surfaces.
- Downstream modules (e.g., PTRM perceptual modeling and structural concordance)
  can consume sweep and ensemble outputs directly.
