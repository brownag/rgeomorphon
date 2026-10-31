---
id: "0005"
title: "Use GRASS-Compatible Integer Class Codes with a Shared Display Theme"
status: "Accepted"
date: 2025-06-03
authors:
  - "Andrew Gene Brown"
supersedes: null
superseded_by: null
source_context: "8b58bc5; a37d6af; R/geomorphon.R; src/geomorphon.cpp"
---

# ADR 0005: Use GRASS Class Codes and a Shared Display Theme

## Context and Problem Statement

The package classifies terrain into named forms. It needs stable, compact codes
for raster storage, plots, and comparison with GRASS GIS. Character labels use
more space and are harder to exchange. Labels and colors can also drift when
each example defines its own.

The package needs one public set of class names and a shared way to apply them.

## Decision Drivers

* Match `r.geomorphon` integer codes so results are comparable with GRASS GIS.
* Store classes as small integers, with labels kept separately.
* Use the same form names and colors in examples and plots.
* Support `terra::SpatRaster` and matrix output.

## Considered Options

1. Return character labels from the classification kernel.
2. Define package-specific integer codes unrelated to GRASS.
3. Reuse GRASS integer codes for the ten forms and expose separate category and
   color tables plus a helper that applies them.
4. Leave all labeling and coloring to users.

## Decision Outcome

Chosen option: "Reuse GRASS integer codes for the ten forms and expose separate
category and color tables plus a helper that applies them." This keeps results
compact and comparable, with display rules in one place.

The initial commit (`8b58bc5`) establishes the vocabulary:

* `src/geomorphon.cpp` defines `enum FORMS_GRASS` with `G_NONE = 0` and
  `G_FL = 1` through `G_PT = 10`, matching the `r.geomorphon` class order.
* `R/geomorphon.R` exports `geomorphon_categories()` (values `1:10`, names
  Flat, Peak, Ridge, Shoulder, Spur, Slope, Hollow, Footslope, Valley, Pit) and
  `geomorphon_colors()`.
* `geomorphon_theme()` applies the names and colors to a `SpatRaster` through
  `terra::levels()` and `terra::coltab()`, or converts a matrix to a factor
  with the same category labels.
* `NAMESPACE` exports all three helpers alongside `geomorphons()`.

Commit `a37d6af` aligned the palette with `r.geomorphon` using hex values
instead of generic R color names. It changed the display table, not the class
codes.

### Consequences

* **Positive**: Integer rasters have a stable, documented mapping to forms.
* **Positive**: GRASS results use the same class values.
* **Positive**: Examples and plots share labels and colors.
* **Positive**: The theme helper works with rasters and matrices.
* **Negative**: Class codes are public. Renaming or reordering forms breaks
  stored rasters and downstream code.
* **Negative**: Without `geomorphon_theme()`, results contain only integers.
* **Negative**: Color changes can alter published figures and comparisons.
* **Negative**: Positive, negative, and ternary counts are not form classes and
  must not use the ten-class theme.

## Validation & Compliance

Tests check that themed rasters retain their values and use the ten labels from
`geomorphon_categories()`. Treat changes to class order or labels as breaks to
the public data contract.
