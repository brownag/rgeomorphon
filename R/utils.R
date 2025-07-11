#' Estimate Tile Processing Needs
#'
#' `.terra_mem_chunks_needed()` is a heuristic for number of tiles needed to
#' calculate geomorphons on larger-than-memory rasters. Allows for scaling by
#' number of parallel workers, a multiplicative factor for the memory needs, and
#' a multiplicative factor for worker needs.
#'
#' @param x A _SpatRaster_ object.
#' @param workers _integer_. Number of parallel workers. Default: `1`
#' @param scl_need _numeric_. Scaling factor for memory needs. Default: `10`
#' @param scl_workers _numeric_ .Scaling factor for each worker. Default: `1`
#'
#' @returns _integer_. Number of tile chunks to divide `x` into.
#' @noRd
.terra_mem_chunks_needed <- function(x, workers = 1, scl_need = 10, scl_workers = 1) {
    mi <- as.data.frame(t(terra::mem_info(x, print = FALSE)))
    ceiling(mi$needed * scl_need / (mi$available * (mi$memfrac / (workers * scl_workers)))) ^ 2
}
