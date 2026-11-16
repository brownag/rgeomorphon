---
id: "0009"
title: "Target aggregation, Jensen-Shannon divergence, and multi-objective Pareto calibration"
status: "Accepted"
date: 2026-11-09
authors:
  - "Andrew Gene Brown"
supersedes: null
superseded_by: null
source_context: "EXECUTION_PLAN_10"
---

# ADR 0009: Target aggregation, Jensen-Shannon divergence, and multi-objective Pareto calibration

## Context and Problem Statement

Evaluating geomorphon landform classifications against ecological, pedological,
and geological mapping themes requires comparing parameter combinations across
multiple objectives. Calibration requires:

1. Deriving aligned categorical reference layers from integer ID rasters with
   Reference Attribute Tables (RAT) without enforcing domain-specific schemas
   (such as SSURGO).
2. Diagnosing whether specific landform classes partition target themes
   cohesively or redundantly using distribution dissimilarity metrics.
3. Comparing continuous spatial distributions (such as class-proportion ensemble
   surfaces) without external heavy dependency requirements while retaining the
   ability to dispatch to established ecology/information-theory packages.
4. Ranking parameter sweep configurations against reference maps without collapsing
   conflicting criteria (information-theoretic overlap, physical boundary
   focal alignment, and boundary saturation) into a brittle arbitrary scalar sum.

## Decision

### Domain-agnostic Reference Attribute Table (RAT) target derivation (`geomorphon_targets`)

Implement `geomorphon_targets(x, attrs = NULL, id_col = 1, vars = NULL, breaks = list(), labels = list(), right = TRUE)`
operating on any integer ID `SpatRaster` and its RAT (`terra::cats(x)[[1]]` or an
explicit user-provided `attrs` `data.frame`):

- Support categorical variables as well as continuous numeric variables partitioned
  by user-specified `breaks` via `cut()`.
- Attach deterministic factor levels (`levels(lyr)`) to each derived `SpatRaster` layer.
- Propagate missing attribute values deterministically without distorting the raster
  coordinate reference system or resolution.
- Attach an audit data frame `attr(res, "coverage")` recording:
  `target`, `n_total`, `n_valid`, `n_na`, `prop_valid`, and `n_classes`.

### Hybrid internal base-2 Jensen-Shannon divergence with optional `philentropy` dispatch

Provide `geomorphon_profile_dist()` and `geomorphon_divergence()` with an exact,
zero-dependency internal base-2 Jensen-Shannon divergence ($D_{\text{JS}} \in [0, 1]$):

- Calculate $D_{\text{JS}}(P \parallel Q) = \frac{1}{2} D_{\text{KL}}(P \parallel M) + \frac{1}{2} D_{\text{KL}}(Q \parallel M)$
  in base 2 where $M = \frac{1}{2}(P + Q)$.
- Implement the information-theoretic convention $0 \log_2 0 = 0$ strictly by
  restricting summations to elements with positive probability, guaranteeing exact
  $D_{\text{JS}} = 1.0$ for distributions with disjoint support without floating-point
  regularization noise ($\epsilon$).
- In `geomorphon_profile_dist()`, compute pairwise dissimilarity between conditional
  probability profiles from a contingency table or `geomorphon_concordance()` result
  along either margin ($P(\text{Category} \mid \text{Geomorphon})$ for `margin = 1`
  or $P(\text{Geomorphon} \mid \text{Category})$ for `margin = 2`).
- In `geomorphon_divergence()`, compute per-cell Jensen-Shannon divergence between
  two class-proportion ensemble stacks (`geomorphon_ensemble(..., method = "proportions")`),
  reconciling layer names across taxonomies (e.g., `forms10` vs. `forms6` on the
  10-class union with unreached classes assigned proportion 0 on valid cells).
- When distance methods other than `"jensen-shannon"` or units other than `"log2"`
  are requested, dispatch dynamically via `requireNamespace("philentropy", quietly = TRUE)`
  to `philentropy::distance()` or `philentropy::dist_one_one()`, keeping `philentropy`
  under `Suggests:` to preserve minimal runtime dependencies (ADR 0002).

### Multi-objective Pareto calibration and boundary diagnostics

Provide `geomorphon_calibrate()` and `geomorphon_pareto()` for parameter sweeps:

- Cross-tabulate and evaluate `geomorphon_concordance()` across each run and requested
  form taxonomy in a `geomorphon_sweep()` object against each layer in `targets`.
- Record physical boundary density diagnostics:
  - `edge_ratio = n_edges_geo / n_edges_cat`
  - `density_penalty = exp(-abs(log(edge_ratio)))`, which penalizes over-fragmentation
    (`edge_ratio >> 1`) and over-generalization (`edge_ratio << 1`) symmetrically.
- Evaluate non-dominated parameter sets via `geomorphon_pareto(cal_df, objectives = NULL, by = "target", min_edge_ratio = 0, max_edge_ratio = Inf)`:
  - Default to 3-objective Pareto evaluation: `c(delta_nmi = "max", delta_edge_f1 = "max", edge_mhd = "min")`
    when spatial null contrasts are available, or `c(nmi = "max", edge_f1 = "max", edge_mhd = "min")`
    when `spatial_null = FALSE`.
  - Filter candidate parameter sets by boundary-density bounds (`[min_edge_ratio, max_edge_ratio]`).
  - Flag non-dominated solutions with `is_pareto = TRUE`.
  - Compute a composite `pareto_score` via the harmonic mean of min-max normalized
    objective values within each evaluation group.

## Consequences

- Thematic calibration and sensitivity analysis operate seamlessly on any RAT-enabled
  raster data source.
- Information divergence calculations require no external libraries by default,
  maintaining CRAN lightweight compliance while supporting extensive distance families
  when `philentropy` is available.
- Parameter calibration balances information-theoretic agreement against physical
  boundary alignment and edge-density balance across landform taxonomies.
