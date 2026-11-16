# Multi-Objective Pareto Frontier and Composite Calibration Score

Identifies non-dominated Pareto-optimal parameter sets and computes
composite harmonic-mean scores over multi-objective evaluation tables.

## Usage

``` r
geomorphon_pareto(
  cal_df,
  objectives = NULL,
  by = "target",
  min_edge_ratio = 0,
  max_edge_ratio = Inf
)
```

## Arguments

- cal_df:

  A `data.frame` containing calibration evaluation metrics (such as the
  output of
  [`geomorphon_calibrate()`](http://humus.rocks/rgeomorphon/reference/geomorphon_calibrate.md)).

- objectives:

  Named character vector defining optimization directions for each
  objective (e.g.
  `c(delta_nmi = "max", delta_edge_f1 = "max", edge_mhd = "min")`).
  Values must be `"max"` or `"min"`. If `NULL`, default objectives are
  selected based on available columns.

- by:

  Character vector of column names defining independent Pareto grouping
  domains. Default is `"target"`.

- min_edge_ratio:

  Minimum permissible `edge_ratio` (`n_edges_geo / n_edges_cat`).
  Default is `0`.

- max_edge_ratio:

  Maximum permissible `edge_ratio`. Default is `Inf`.

## Value

The input `cal_df` augmented with:

- is_pareto:

  Logical. `TRUE` if the candidate belongs to the non-dominated Pareto
  frontier within its group.

- pareto_score:

  Numeric in \\\[0, 1\]\\. Normalized composite score computed as the
  harmonic mean across normalized objectives.

## Examples

``` r
df <- data.frame(
  target = "musym",
  run_id = 1:3,
  edge_ratio = c(1.1, 1.2, 0.9),
  delta_nmi = c(0.12, 0.15, 0.10),
  delta_edge_f1 = c(0.20, 0.18, 0.22),
  edge_mhd = c(45, 40, 50)
)
geomorphon_pareto(df)
#>   target run_id edge_ratio delta_nmi delta_edge_f1 edge_mhd is_pareto
#> 1  musym      1        1.1      0.12          0.20       45      TRUE
#> 2  musym      2        1.2      0.15          0.18       40      TRUE
#> 3  musym      3        0.9      0.10          0.22       50      TRUE
#>   pareto_score
#> 1    0.4615385
#> 2    0.0000000
#> 3    0.0000000
```
