---
id: "0001"
title: "Use a C++ RcppParallel Worker for Geomorphon Classification"
status: "Accepted"
date: 2025-06-03
authors:
  - "Andrew Gene Brown"
supersedes: null
superseded_by: null
source_context: "8b58bc5; src/geomorphon.cpp; src/Makevars; DESCRIPTION"
---

# ADR 0001: Use an RcppParallel Worker for Classification

## Context and Problem Statement

`rgeomorphon` reimplements the `r.geomorphon` algorithm as an R package. It
visits every raster cell and compares nearby elevations. An R implementation
would spend too much time in the interpreter on realistic DEMs. The package
needs a compiled kernel that can use multiple CPU cores without requiring users
to set up a cluster.

The initial package description promises a "fast and low-dependency"
implementation, with a parallel C++ algorithm that has no external dependencies
(`8b58bc5`, `DESCRIPTION`). The threading method must fit that limit.

## Decision Drivers

* Avoid interpreter overhead in the per-cell loop.
* Use multiple cores on one machine by default.
* Keep runtime and build dependencies small and portable across R platforms.
* Keep the R interface simple. Users should not manage native pointers or
  threads.

## Considered Options

1. Implement the whole algorithm in R.
2. Shell out to GRASS GIS `r.geomorphon`.
3. Write a serial C++ kernel and rely on R-level parallelism only.
4. Write a C++ kernel as an `RcppParallel::Worker` and dispatch it with
   `RcppParallel::parallelFor()`.
5. Use raw OpenMP loops directly from Rcpp.

## Decision Outcome

Chosen option: "Write a C++ kernel as an `RcppParallel::Worker` and dispatch it
with `RcppParallel::parallelFor()`." This uses shared-memory parallelism while
keeping the public interface in R and avoiding a system GRASS dependency.

The initial commit (`8b58bc5`) sets this pattern:

* `DESCRIPTION` imports `RcppParallel` and links to `Rcpp`, `RcppArmadillo`,
  and `RcppParallel`.
* `src/geomorphon.cpp` defines `GeomorphonWorker`, which reads a shared
  `RcppParallel::RMatrix<double>` DEM and writes preallocated output matrices.
  Its `operator()` handles a contiguous range of rows.
* The exported `geomorphons_cpp_worker()` calls
  `RcppParallel::parallelFor(0, elevation.nrow(), geomorphon_worker)`.
* `src/Makevars` and `src/Makevars.win` add the platform OpenMP flags supplied
  by R (`$(SHLIB_OPENMP_CXXFLAGS)`).

Later, `3c531e3` added ternary, positive, and negative output matrices to the
same worker instead of adding another execution model.

### Consequences

* **Positive**: The compiled loop can spread rows across threads without
  changing the R API.
* **Positive**: R's build tools handle the worker. No separate scheduler or
  geospatial program is needed.
* **Positive**: R allocates output matrices, and the worker writes to them.
* **Negative**: `RcppParallel` is needed even for matrix-only use.
* **Negative**: Each output type adds fields and arguments to the worker.
* **Negative**: Builds use R's OpenMP flags where available. The package must
  still work without them.

## Validation & Compliance

The tinytest suite covers the worker and its exported wrapper. Kernel changes
must preserve the row-range contract and write only to per-cell output matrices.
