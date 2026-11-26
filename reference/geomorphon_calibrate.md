# Calibration of Geomorphon Parameters Against Reference Targets

Evaluates spatial concordance between multiple runs of a
[`geomorphon_sweep()`](https://humus.rocks/rgeomorphon/reference/geomorphon_sweep.md)
and reference categorical target raster layers. Combines sweep parameter
metadata with concordance metrics, boundary-density penalties, and
Pareto-optimal frontier identification.

## Usage

``` r
geomorphon_calibrate(
  sweep_obj,
  targets,
  forms = NULL,
  method = c("both", "nmi", "boundary"),
  threshold_cells = 2,
  spatial_null = TRUE
)
```

## Arguments

- sweep_obj:

  An object of class `geomorphon_sweep`.

- targets:

  A single- or multi-layer categorical `SpatRaster` representing
  reference classification targets (e.g. from
  [`geomorphon_targets()`](https://humus.rocks/rgeomorphon/reference/geomorphon_targets.md)).

- forms:

  Character vector of form taxonomies to evaluate (e.g.
  `c("forms10", "forms6")`). Default is `NULL`, which evaluates all
  taxonomies found in `sweep_obj`.

- method:

  Character concordance method: `"both"` (default), `"nmi"`, or
  `"boundary"`. Passed to
  [`geomorphon_concordance()`](https://humus.rocks/rgeomorphon/reference/geomorphon_concordance.md).

- threshold_cells:

  Integer or numeric cell distance tolerance for boundary metrics.
  Default is `2`.

- spatial_null:

  Logical. Whether to evaluate spatial-null circular shift ensembles in
  [`geomorphon_concordance()`](https://humus.rocks/rgeomorphon/reference/geomorphon_concordance.md).
  Default is `TRUE`.

## Value

A tidy `data.frame` combining sweep parameter columns (`run_id`,
`search`, `skip`, `flat_angle_deg`, `dist`, `comparison_mode`, `tdist`,
`resolution`, `elapsed_sec`, `n_cells`), `forms`, `target`, concordance
metrics, `edge_ratio`, `density_penalty`, `is_pareto`, and
`pareto_score`.

## Examples

``` r
library(terra)
data(salton, package = "rgeomorphon")
dem <- rast(salton)
crs(dem) <- attr(salton, "crs")
ext(dem) <- attr(salton, "extent")
sw <- geomorphon_sweep(dem, param_grid = list(search = c(3, 4)), forms = "forms6")
cat_rast <- (dem > 50) + 1
levels(cat_rast) <- data.frame(id = 1:2, category = c("Low", "High"))
names(cat_rast) <- "elevation_class"
cal <- geomorphon_calibrate(sw, targets = cat_rast, spatial_null = FALSE)
print(cal)
#>   run_id search skip flat_angle_deg dist comparison_mode tdist resolution
#> 1      1      3    0              1    0         anglev1     0        300
#> 2      2      4    0              1    0         anglev1     0        300
#>   elapsed_sec n_cells  forms          target nmi homogeneity completeness
#> 1       0.028   26565 forms6 elevation_class   0          NA            0
#> 2       0.026   26565 forms6 elevation_class   0          NA            0
#>   n_edges_geo n_edges_cat edge_ratio density_penalty edge_offset
#> 1         858           0         NA              NA          NA
#> 2         903           0         NA              NA          NA
#>   edge_concordance edge_f1 edge_mhd edge_chamfer mean_null_nmi delta_nmi
#> 1               NA      NA       NA           NA            NA        NA
#> 2               NA      NA       NA           NA            NA        NA
#>   nmi_null_shifts mean_null_edge_f1 delta_edge_f1 edge_f1_null_shifts
#> 1               0                NA            NA                   0
#> 2               0                NA            NA                   0
#>   null_shifts_evaluated is_pareto pareto_score
#> 1                     0      TRUE            1
#> 2                     0      TRUE            1
```
