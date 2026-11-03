---
id: "0006"
title: "PTRM conventions"
status: "Accepted"
date: 2026-10-31
authors:
  - "Andrew Gene Brown"
supersedes: null
superseded_by: null
source_context: "EXECUTION_PLAN_1"
---

# ADR 0006: PTRM conventions

## Context and Problem Statement

`rgeomorphon` exposes the Perceived Terrain Realism Metric (PTRM) of
Rajasekaran et al. (2022) through `geomorphon_ptrm()` and the focal variant
`geomorphon_ptrm_focal()`. The metric combines the proportions of geomorphon
forms with fitted weights. The package also supports the reduced form systems
`forms6`, `forms5`, and `forms4`, and keeps `geomorphon_concordance()` as a
boundary-alignment companion.

Without a record, the formula, coefficient source, reduced-system derivation,
model constants, and focal-window rules can change differently in the global,
focal, C++, and R code. This ADR fixes them.

## Decision

### PTRM formula and coefficients

PTRM is computed as `PTRM = (c0 + sum_i(c_i * p_i)) / D`. Here `c0` is the
intercept, `c_i` is the coefficient for form `i`, `p_i` is the proportion of
form `i`, and `D` is the normalizer.

Use the ten-form coefficients from Rajasekaran et al. (2022), Table 1, as the
only source for the PTRM model:

- flat 25.12
- peak 1.75
- ridge 7.59
- shoulder 7.31
- spur 9.02
- slope 28.95
- hollow 6.71
- footslope 7.63
- valley 9.61
- depression 3.55

### Reduced form systems

Reduced systems have no published coefficients. Their default weights are the
composition-weighted mean of the ten-form weights over the multiset of reachable
(negative, positive) count pairs for each reduced class. Count a repeated pair
each time it occurs.

A user coefficient vector must be named for the requested form system and cover
exactly its reachable classes. A complete named ten-form vector is also accepted
and reduced by the same rule. Reject coefficients for another form system. Do
not reuse them.

### Model constants and output scale

`intercept` and `normalizer` are tunable. Their defaults are `-38.02` and
`69.96`, which match the paper's model. The normalizer must be finite and
non-zero.

Return PTRM without clamping. Scores may fall outside `[0, 1]` with different
form distributions or parameters. Do not rescale or truncate them.

### Focal window

The focal window uses cells. Include a cell at offset `(dx, dy)` when
`sqrt(dx^2 + dy^2) <= radius`. The `radius` is an integer cell count, not a
distance in map units.

Set cells within `radius` of a raster edge to `NA` because their neighborhoods
are incomplete. The fast C++ path and the memory-efficient `terra::focal()` path
must use the same mask, border handling, and formula. Keep them equal.

### Concordance companion

Keep `geomorphon_concordance()` as an experimental companion to the PTRM
functions. It measures structural agreement with normalized mutual information
and boundary alignment. It does not produce a PTRM score.

Report boundary distances in map units of the input CRS. `edge_offset` is the
mean distance from each categorical boundary cell to the nearest geomorphon
boundary cell. `edge_concordance` is the proportion of categorical boundary
cells within the requested threshold of a geomorphon boundary.

## Consequences

- The `forms4`, `forms5`, and `forms6` coefficients are derived approximations.
  The PTRM paper does not validate them. Do not present them as paper-published
  values.
- PTRM is a perceptual metric from the cited study. It is not a validated
  measure of terrain quality.
- This change does not benchmark against the paper's reported `0.76` score.
  Absolute scores depend on the study context, not package calibration.
- Testing whether PTRM is useful for this package's applications is out of
  scope.
- Changing the coefficients, constants, reduced-system rule, focal window, or
  concordance definitions requires updating this ADR and the tests.
- Keep testing the C++ and `terra::focal()` paths for matching results,
  including border behavior.

## References

Rajasekaran, S. D., Kang, H., Čadík, M., Galin, E., Guérin, E., Peytavie, A.,
Slavík, P., and Mech, R. (2022). PTRM (Perceived Terrain Realism Metric). ACM
Transactions on Applied Perception, 19(2), Article 5.
<https://doi.org/10.1145/3514244>
