# Estimate Tile Processing Needs

`geomorphon_chunks_needed()` is a heuristic for number of tiles needed
to calculate geomorphons on larger-than-memory rasters. Allows for
scaling by number of parallel workers, a multiplicative factor for the
memory needs, and a multiplicative factor for worker needs.

## Usage

``` r
geomorphon_chunks_needed(
  x,
  workers = Sys.getenv("R_RGEOMORPHON_N_WORKERS", unset = 1),
  scl_need = Sys.getenv("R_RGEOMORPHON_MEM_SCALE_NEED", unset = 10),
  scl_workers = Sys.getenv("R_RGEOMORPHON_MEM_SCALE_WORKERS", unset = 1),
  pow_total = Sys.getenv("R_RGEOMORPHON_MEM_POWER", unset = 0.5)
)
```

## Arguments

- x:

  A *SpatRaster* object.

- workers:

  *integer*. Number of parallel workers. Default uses value of
  environment variable `R_RGEOMORPHON_N_WORKERS`. If unset, `1`

- scl_need:

  *numeric*. Scaling factor for memory needs. Default uses value of
  environment variable `R_RGEOMORPHON_MEM_SCALE_NEED`. If unset, `10`.

- scl_workers:

  *numeric*. Scaling factor for each worker. Default uses value of
  environment variable `R_RGEOMORPHON_MEM_SCALE_WORKERS`. If unset, `1`.

- pow_total:

  *numeric*. Exponent for scaling total number of chunks. Default uses
  value of environment variable `R_RGEOMORPHON_MEM_POWER`. If unset,
  `1`.

## Value

*integer*. Number of tile chunks to divide `x` into.

## Examples

``` r
data("salton", package = "rgeomorphon")

x <- terra::rast(salton)
terra::ext(x) <- attr(salton, "extent")
terra::crs(x) <- attr(salton, "crs")

geomorphon_chunks_needed(x)
#> [1] 1
```
