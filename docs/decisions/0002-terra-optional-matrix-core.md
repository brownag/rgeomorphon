---
id: "0002"
title: "Keep Runtime Dependencies Minimal and Make terra Optional"
status: "Accepted"
date: 2025-06-03
authors:
  - "Andrew Gene Brown"
supersedes: null
superseded_by: null
source_context: "8b58bc5; DESCRIPTION; R/geomorphon.R; NAMESPACE; 7021c14"
---

# ADR 0002: Keep Dependencies Small and terra Optional

## Context and Problem Statement

`rgeomorphon` is a low-dependency alternative to GRASS GIS. Many R users pass
georeferenced `terra::SpatRaster` objects, but requiring `terra` would add a
large geospatial stack to every installation, including matrix-only use and
builds that need only the classification kernel.

The package must set the boundary between its core and spatial data handling,
including coordinates, raster metadata, and file output.

## Decision Drivers

* Keep installed dependencies small.
* Run the core algorithm on a numeric matrix.
* Support `terra` without requiring `raster`, `sf`, or GRASS.
* Keep optional spatial features from blocking installation or checks.

## Considered Options

1. Depend unconditionally on `terra` and accept only `SpatRaster` input.
2. Support the older `raster` package as well as `terra`.
3. Require callers to supply matrices only and leave all geospatial I/O to
   other packages.
4. Accept matrices in the core and provide an optional `terra` adapter through
   `Suggests` plus `requireNamespace()`.

## Decision Outcome

Chosen option: "Accept matrices in the core and provide an optional `terra`
adapter through `Suggests` plus `requireNamespace()`." This keeps the core small
and supports the common R raster API.

Evidence from the initial commit (`8b58bc5`):

* `DESCRIPTION` lists `Imports: RcppParallel` and `Suggests: terra, tinytest`;
  `terra` is not imported.
* `R/geomorphon.R` documents `elevation` as "matrix or SpatRaster object" and
  branches on `!is.matrix(elevation)`.
* The `SpatRaster` branch checks `requireNamespace("terra")` and stops with
  "Package 'terra' is required for SpatRaster input" when it is absent, then
  converts the raster to a matrix with `terra::as.matrix()`.
* Matrix input uses explicit `xres`/`yres` arguments when no raster metadata is
  available.

Tiled processing kept the same boundary. It runs only for `SpatRaster` input;
matrices use the in-memory kernel (`7021c14`). Optional parallelism follows the
same rule: `future` and `future.apply` are in `Suggests` and `Enhances`, not
`Imports`.

### Consequences

* **Positive**: `geomorphons()` runs on a numeric matrix with only compiled
  dependencies installed.
* **Positive**: `terra` users get raster input and output without making it
  mandatory for everyone.
* **Positive**: `requireNamespace()` guards optional features without adding
  them to `Imports`.
* **Negative**: Matrix users must supply cell resolution and handle CRS
  interpretation themselves.
* **Negative**: Tiling, focal processing, and raster output need `terra` and must
  report clearly when it is missing.
* **Negative**: The matrix and raster paths must stay consistent, including at
  borders and in output names.

## Validation & Compliance

Examples that need `terra` use `@examplesIf requireNamespace("terra")`. Tests
cover matrix and `SpatRaster` input when the optional package is available.
