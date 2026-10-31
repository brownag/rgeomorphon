---
id: "0003"
title: "Represent Reduced Form Systems as Composable Classification Matrices"
status: "Accepted"
date: 2025-06-25
authors:
  - "Andrew Gene Brown"
supersedes: null
superseded_by: null
source_context: "3c531e3; ee5a2d0; src/geomorphon.cpp; R/forms_matrix.R"
---

# ADR 0003: Use Classification Matrices for Reduced Form Systems

## Context and Problem Statement

The original geomorphon output has ten classes based on positive and negative
overlook counts. Users also need the 4-, 5-, and 6-form systems from Masetti et
al. (2018), as well as custom groupings of count cells. A separate pass for each
system would duplicate the cell loop and rule out custom systems.

One representation must support the published systems, custom classes, and the
ten-form default without changing the compiled kernel for each new scheme.

## Decision Drivers

* Support the 10-, 6-, 5-, and 4-form systems with one engine.
* Let users group positive/negative count pairs.
* Keep ten forms as the default and preserve their class IDs.
* Change classification without repeating the neighbor search.

## Considered Options

1. Add a separate branch in the C++ worker for each published form system.
2. Always compute the ten-form result and reclassify it in R by form label.
3. Compute positive and negative overlook counts in the worker, then map each
   `(negative, positive)` count pair through a 9x9 classification matrix.
4. Require users to post-process a ternary pattern outside the package.

## Decision Outcome

Chosen option: "Compute positive and negative overlook counts in the worker,
then map each `(negative, positive)` pair through a 9x9 classification matrix."
This keeps the neighbor search in one place and stores class mappings as data.

The count outputs and the reduced-system idea first appear in `3c531e3`, which
added positive, negative, and ternary outputs to `GeomorphonWorker` and recorded
the 4-form table as a commented alternative. `ee5a2d0` then implemented the
general mechanism:

* `src/geomorphon.cpp` holds `forms_table10`, `forms_table6`, `forms_table5`,
  and `forms_table4`, all 9x9 arrays indexed by negative and positive counts.
* `form_from_counts(num_neg, num_pos, num_forms, none_val)` selects the table,
  and `get_forms_matrix_cpp()` returns the chosen table to R.
* `get_forms_grass_enum()` exposes the stable `G_*` integer identifiers.
* `R/forms_matrix.R` defines the `forms_matrix` S3 class, which pairs a 9x9
  matrix with named levels, plus `forms_matrix_get()` for the standard systems
  and `forms_matrix_apply()` for reclassifying positive/negative raster layers
  with `terra::subst()`.

### Consequences

* **Positive**: Published and custom systems use one engine and matrix format.
* **Positive**: Reduced forms can reuse positive and negative layers without
  repeating the geomorphon search.
* **Positive**: Users can inspect and edit the 9x9 matrix to define custom
  classes without compiled code.
* **Negative**: To use `forms_matrix_apply()` later, callers must request
  `positive = TRUE` and `negative = TRUE`.
* **Negative**: The C++ tables and R matrices must match. Changes affect the
  public classification contract.
* **Negative**: Reduced systems use ten-form pixel assignments, not a separate
  reduced-form algorithm.
* **Negative**: PTRM coefficients use the same count mapping and must change
  when the tables do.

## Validation & Compliance

Tests cover the standard systems in `forms_matrix_get()` and check
`forms_matrix_apply()` against `geomorphons()` positive and negative outputs.
Table changes must also update the tests and PTRM conventions record.
