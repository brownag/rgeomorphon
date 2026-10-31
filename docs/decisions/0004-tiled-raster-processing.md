---
id: "0004"
title: "Process Large Rasters as Buffered Tiles with Two-Level Parallelism"
status: "Accepted"
date: 2025-06-30
authors:
  - "Andrew Gene Brown"
supersedes: null
superseded_by: null
source_context: "7021c14; 54bbe95; 119db3b; bfd01cd; bd4631d; R/geomorphon.R; R/utils.R"
---

# ADR 0004: Process Large Rasters as Buffered Tiles

## Context and Problem Statement

A DEM can exceed available memory. Each output cell also needs nearby elevation
data out to the search radius. Splitting a raster without overlap leaves gaps in
that data at tile edges and can change the results.

The package already parallelizes rows in one C++ call. It needs a way to process
larger rasters without changing classifications at tile edges.

## Decision Drivers

* Process rasters larger than memory without requiring users to split them.
* Match single-pass results at tile edges.
* Reuse the C++ kernel instead of writing a streaming algorithm.
* Let users override the memory heuristic.
* Let users distribute tiles across R sessions or machines with their own
  parallel backend.

## Considered Options

1. Require the full raster to fit in memory.
2. Split the raster into unbuffered tiles and accept edge artifacts.
3. Rely on a generic `terra` block-processing loop without geomorphon-specific
   buffering.
4. Split into tiles, add a buffer of `search + skip + 1` cells, process each
   buffered tile with the existing kernel, crop the buffer, and mosaic.
5. Stream cell-by-cell and maintain state across the raster.

## Decision Outcome

Chosen option: "Split into tiles, add a buffer of `search + skip + 1` cells,
process each buffered tile with the existing kernel, crop the buffer, and
mosaic." This limits memory use while keeping each cell's full neighborhood.

Commit `7021c14` introduced tiled processing. Later commits refined it:

* `7021c14` adds `.geomorphons_tiled()` and helper functions for tile extents,
  tile creation, processing, and combination. It computes a cell buffer of
  `search + skip + 1` and accepts an `LAPPLY.FUN` argument so callers can pass
  `future.apply::future_lapply()`.
* `54bbe95` adds a memory-based tile-count heuristic, later exported as
  `geomorphon_chunks_needed()` and documented in `119db3b`. The heuristic uses
  `terra::mem_info()` plus configurable scaling factors and worker counts.
* `bfd01cd` adds the `nchunk` argument so users can bypass the heuristic.
* `bd4631d` adds internal support for a custom `cell_buffer`.
* `geomorphons()` applies the tiled path to `SpatRaster` input, with a special
  extend-then-crop path when `nchunk == 1` (`R/geomorphon.R`).

This gives two levels of parallelism, as described in
`vignettes/geomorphons-tiled.Rmd`. `RcppParallel` works within each tile, and
`LAPPLY.FUN` can distribute tiles across R sessions.

### Consequences

* **Positive**: Large rasters use the existing kernel.
* **Positive**: Buffers prevent edge errors from splitting the raster.
* **Positive**: Users can set the tile count or use the memory heuristic.
* **Positive**: The same `geomorphons()` call supports sequential, multisession,
  or cluster backends.
* **Negative**: Temporary tiles add file I/O, so this is slower than one
  in-memory pass.
* **Negative**: Buffers must cover the search and skip windows.
* **Negative**: Tiling requires `terra`; matrices use the in-memory kernel.
* **Negative**: The memory heuristic cannot guarantee a tile will fit every
  system.

## Validation & Compliance

Chunked tests compare tiled results with the in-memory path. The vignette
describes both parallelism layers. Check buffer and tile-combination changes
against single-pass results at tile edges.
