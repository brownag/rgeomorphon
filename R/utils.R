#' Estimate Tile Processing Needs
#'
#' `geomorphon_chunks_needed()` is a heuristic for number of tiles needed to
#' calculate geomorphons on larger-than-memory rasters. Allows for scaling by
#' number of parallel workers, a multiplicative factor for the memory needs, and
#' a multiplicative factor for worker needs.
#'
#' @param x A _SpatRaster_ object.
#' @param workers _integer_. Number of parallel workers. Default uses value of
#'   environment variable `R_RGEOMORPHON_N_WORKERS`. If unset, `1`
#' @param scl_need _numeric_. Scaling factor for memory needs. Default uses
#'   value of environment variable `R_RGEOMORPHON_MEM_SCALE_NEED`. If unset,
#'   `10`.
#' @param scl_workers _numeric_. Scaling factor for each worker. Default uses
#'   value of environment variable `R_RGEOMORPHON_MEM_SCALE_WORKERS`. If unset,
#'   `1`.
#' @param pow_total _numeric_. Exponent for scaling total number of chunks.
#'   Default uses value of environment variable `R_RGEOMORPHON_MEM_POWER`. If
#'   unset, `1`.
#'
#' @returns _integer_. Number of tile chunks to divide `x` into.
#' @export
geomorphon_chunks_needed <- function(x,
                                     workers = Sys.getenv("R_RGEOMORPHON_N_WORKERS", unset = 1),
                                     scl_need = Sys.getenv("R_RGEOMORPHON_MEM_SCALE_NEED", unset = 10),
                                     scl_workers = Sys.getenv("R_RGEOMORPHON_MEM_SCALE_WORKERS", unset = 1),
                                     pow_total = Sys.getenv("R_RGEOMORPHON_MEM_POWER", unset = 0.5)) {
    mi <- as.data.frame(t(terra::mem_info(x, print = FALSE)))
    ceiling((mi$needed * as.numeric(scl_need) /
                (mi$available * (mi$memfrac / (as.numeric(workers) * as.numeric(scl_workers)))))^as.numeric(pow_total))
}
