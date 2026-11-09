library(tinytest)
library(rgeomorphon)

# ==============================================================================
# 1. geomorphon_sweep on matrix interface (volcano)
# ==============================================================================
data("volcano", package = "datasets")

grid_mat <- list(
    search = c(3, 5),
    skip = c(0, 1),
    flat_angle_deg = 1.0
)

sw_mat <- geomorphon_sweep(
    volcano,
    param_grid = grid_mat,
    xres = 10,
    forms = c("forms10", "forms6")
)

expect_true(inherits(sw_mat, "geomorphon_sweep"))
expect_true(is.data.frame(sw_mat$summary))
expect_equal(nrow(sw_mat$summary), 4L)
expected_cols <- c("run_id", "search", "skip", "flat_angle_deg", "dist",
                   "comparison_mode", "tdist", "resolution", "elapsed_sec", "n_cells")
expect_true(all(expected_cols %in% names(sw_mat$summary)))
expect_true(all(sw_mat$summary$elapsed_sec >= 0))
expect_equal(sw_mat$summary$n_cells[1], length(volcano))
expect_equal(length(sw_mat$rasters), 4L)
expect_equal(names(sw_mat$rasters), paste0("run_", 1:4))
expect_true(all(c("forms10", "forms6") %in% names(sw_mat$rasters[[1]])))

# Provenance records what produced the sweep
expect_true(all(c("package_version", "created", "forms", "grid_columns", "extra_args", "n_runs")
                %in% names(sw_mat$provenance)))
expect_equal(sw_mat$provenance$forms, c("forms10", "forms6"))
expect_equal(sw_mat$provenance$n_runs, 4L)
expect_equal(sw_mat$provenance$package_version, as.character(packageVersion("rgeomorphon")))

# S3 methods for geomorphon_sweep
out_p <- capture.output(print(sw_mat))
expect_true(any(grepl("<geomorphon_sweep>", out_p)))
expect_true(any(grepl("Runs: 4", out_p)))
expect_equal(summary(sw_mat), sw_mat$summary)

# ==============================================================================
# 2. geomorphon_sweep edge cases (invalid params, pruning, errors)
# ==============================================================================
# Missing param_grid
expect_error(geomorphon_sweep(volcano))

# Empty param_grid
expect_error(geomorphon_sweep(volcano, param_grid = data.frame()))

# Pruning invalid combinations (skip >= search)
grid_prune <- data.frame(
    search = c(5, 3),
    skip = c(1, 3), # row 2 has skip 3 >= search 3
    flat_angle_deg = c(1.0, 1.0)
)
expect_warning({
    sw_prune <- geomorphon_sweep(volcano, param_grid = grid_prune, xres = 10, forms = "forms10")
})
expect_equal(nrow(sw_prune$summary), 1L)
expect_equal(sw_prune$summary$search, 5)
expect_equal(sw_prune$summary$skip, 1)

# All combinations invalid
grid_all_bad <- data.frame(search = c(3, 4), skip = c(3, 5))
expect_error({
    suppressWarnings(geomorphon_sweep(volcano, param_grid = grid_all_bad, xres = 10))
})

# ==============================================================================
# 3. geomorphon_ensemble on matrix interface
# ==============================================================================
m1 <- matrix(c(1, 3, 1, NA), nrow = 2, ncol = 2)
m2 <- matrix(c(1, 6, 3, NA), nrow = 2, ncol = 2)
m3 <- matrix(c(3, 6, 1, NA), nrow = 2, ncol = 2)

ens_m_mat <- geomorphon_ensemble(list(m1, m2, m3), method = "modal", ties = "lowest")
expect_true(is.matrix(ens_m_mat))
expect_equal(ens_m_mat[1, 1], "Flat")
expect_equal(ens_m_mat[2, 1], "Slope")
expect_true(is.na(ens_m_mat[2, 2]))

ens_p_mat <- geomorphon_ensemble(list(m1, m2, m3), method = "proportions")
expect_true(is.list(ens_p_mat))
expect_equal(ens_p_mat[["Flat"]][1, 1], 2 / 3)
expect_equal(ens_p_mat[["Ridge"]][1, 1], 1 / 3)
expect_true(is.na(ens_p_mat[["Flat"]][2, 2]))
# Proportions default to the reachable classes of the selected form (forms6 = 6 classes)
ens_p6_mat <- geomorphon_ensemble(list(m1, m2, m3), method = "proportions", form = "forms6")
expect_equal(names(ens_p6_mat), c("Flat", "Ridge", "Shoulder", "Slope", "Footslope", "Valley"))

# ==============================================================================
# 4. SpatRaster tests: geomorphon_sweep, geomorphon_ensemble, and benchmarks
# ==============================================================================
if (requireNamespace("terra", quietly = TRUE)) {
    data("salton", package = "rgeomorphon")
    dem <- terra::rast(salton)
    terra::crs(dem) <- attr(salton, "crs")
    terra::ext(dem) <- attr(salton, "extent")

    # Sweep on SpatRaster
    grid_rast <- list(
        search = c(3, 5),
        skip = c(0, 1),
        flat_angle_deg = 1.0
    )
    sw_rast <- geomorphon_sweep(
        dem,
        param_grid = grid_rast,
        forms = c("forms10", "forms6", "forms5", "forms4")
    )
    expect_equal(nrow(sw_rast$summary), 4L)
    expect_true(inherits(sw_rast$rasters[[1]], "SpatRaster"))
    expect_equal(terra::nlyr(sw_rast$rasters[[1]]), 4L)
    expect_equal(names(sw_rast$rasters[[1]]), c("forms10", "forms6", "forms5", "forms4"))
    expect_equal(terra::cats(sw_rast$rasters[[1]][["forms10"]])[[1]]$form, geomorphon_categories()$form)

    # --------------------------------------------------------------------------
    # Hand-calculated synthetic ensemble fixture
    # --------------------------------------------------------------------------
    # 5 cells:
    # Cell 1: 1, 1, 3  -> mode=1 (Flat), prop Flat=2/3, prop Ridge=1/3
    # Cell 2: 3, 6, 6  -> mode=6 (Slope), prop Slope=2/3, prop Ridge=1/3
    # Cell 3: 1, 3, NA -> tie between 1 and 3 over 2 valid runs
    # Cell 4: NA, NA, NA -> all NA
    # Cell 5: 1, NA, 1 -> 2 valid runs of 1 -> mode=1, prop Flat=1.0
    syn1 <- terra::rast(nrows = 1, ncols = 5, vals = c(1, 3, 1, NA, 1), crs = "EPSG:4326", ext = c(0, 5, 0, 1))
    syn2 <- terra::rast(nrows = 1, ncols = 5, vals = c(1, 6, 3, NA, NA), crs = "EPSG:4326", ext = c(0, 5, 0, 1))
    syn3 <- terra::rast(nrows = 1, ncols = 5, vals = c(3, 6, NA, NA, 1), crs = "EPSG:4326", ext = c(0, 5, 0, 1))
    names(syn1) <- "forms10"
    names(syn2) <- "forms10"
    names(syn3) <- "forms10"
    syn_list <- list(syn1, syn2, syn3)

    # Modal aggregation with deterministic tie-breaking
    ens_syn_low <- geomorphon_ensemble(syn_list, method = "modal", ties = "lowest")
    ens_syn_high <- geomorphon_ensemble(syn_list, method = "modal", ties = "highest")

    vals_low <- as.vector(terra::values(ens_syn_low))
    vals_high <- as.vector(terra::values(ens_syn_high))

    expect_equal(vals_low[1], 1)      # Flat
    expect_equal(vals_low[2], 6)      # Slope
    expect_equal(vals_low[3], 1)      # Tie resolved to lowest: 1
    expect_equal(vals_high[3], 3)     # Tie resolved to highest: 3
    expect_true(is.na(vals_low[4]))   # All NA -> NA
    expect_equal(vals_low[5], 1)      # Partial NA -> 1

    # Preserves factor categories and color table
    expect_equal(terra::cats(ens_syn_low)[[1]]$form, geomorphon_categories()$form)
    expect_true(!is.null(terra::coltab(ens_syn_low)[[1]]))

    # Proportions aggregation
    ens_syn_prop <- geomorphon_ensemble(syn_list, method = "proportions")
    expect_equal(terra::nlyr(ens_syn_prop), 10L)
    expect_equal(names(ens_syn_prop), geomorphon_categories()$form)

    prop_flat <- as.vector(terra::values(ens_syn_prop[["Flat"]]))
    prop_ridge <- as.vector(terra::values(ens_syn_prop[["Ridge"]]))
    prop_slope <- as.vector(terra::values(ens_syn_prop[["Slope"]]))

    expect_equal(prop_flat[1], 2 / 3, tolerance = 1e-6)
    expect_equal(prop_ridge[1], 1 / 3, tolerance = 1e-6)
    expect_equal(prop_slope[1], 0.0, tolerance = 1e-6)

    expect_equal(prop_slope[2], 2 / 3, tolerance = 1e-6)
    expect_equal(prop_ridge[2], 1 / 3, tolerance = 1e-6)

    expect_equal(prop_flat[3], 0.5, tolerance = 1e-6)
    expect_equal(prop_ridge[3], 0.5, tolerance = 1e-6)

    expect_true(is.na(prop_flat[4]))
    expect_true(is.na(prop_ridge[4]))

    expect_equal(prop_flat[5], 1.0, tolerance = 1e-6)

    # Sum of proportions equals 1 on all valid cells
    sum_prop <- as.vector(terra::values(sum(ens_syn_prop)))
    expect_equal(sum_prop[c(1, 2, 3, 5)], rep(1.0, 4), tolerance = 1e-6)
    expect_true(is.na(sum_prop[4]))

    # --------------------------------------------------------------------------
    # Real DEM ensemble (salton)
    # --------------------------------------------------------------------------
    ens_mode_salton <- geomorphon_ensemble(sw_rast, method = "modal")
    expect_true(inherits(ens_mode_salton, "SpatRaster"))
    expect_equal(terra::nlyr(ens_mode_salton), 1L)
    expect_equal(names(ens_mode_salton), "modal")
    expect_equal(terra::cats(ens_mode_salton)[[1]]$form, geomorphon_categories()$form)

    ens_prop_salton <- geomorphon_ensemble(sw_rast, method = "proportions")
    expect_equal(terra::nlyr(ens_prop_salton), 10L)
    expect_equal(names(ens_prop_salton), geomorphon_categories()$form)
    sum_prop_salton <- as.vector(terra::values(sum(ens_prop_salton)))
    expect_true(all(abs(sum_prop_salton[!is.na(sum_prop_salton)] - 1.0) < 1e-5))

    # Reduced form ensemble: forms6 reachable classes (6 layers)
    ens_prop6 <- geomorphon_ensemble(sw_rast, method = "proportions", form = "forms6")
    expect_equal(terra::nlyr(ens_prop6), 6L)
    expect_equal(names(ens_prop6), c("Flat", "Ridge", "Shoulder", "Slope", "Footslope", "Valley"))

    # Per-run filenames: each run writes its own file, so no run overwrites another
    fn_out <- file.path(tempdir(), "sweep_fn.tif")
    sw_fn <- geomorphon_sweep(dem, param_grid = list(search = c(3, 5), skip = 0, flat_angle_deg = 1),
        forms = "forms10", filename = fn_out)
    expect_true(file.exists(file.path(tempdir(), "sweep_fn_run_1.tif")))
    expect_true(file.exists(file.path(tempdir(), "sweep_fn_run_2.tif")))
    ref_s3 <- geomorphons(dem, search = 3, skip = 0, flat_angle_deg = 1, forms = "forms10")
    ref_s5 <- geomorphons(dem, search = 5, skip = 0, flat_angle_deg = 1, forms = "forms10")
    expect_equal(as.vector(terra::values(sw_fn$rasters[[1]])), as.vector(terra::values(ref_s3)))
    expect_equal(as.vector(terra::values(sw_fn$rasters[[2]])), as.vector(terra::values(ref_s5)))

    # Direct SpatRaster input: a single run and a multi-layer stack of runs
    ens_one <- geomorphon_ensemble(sw_fn$rasters[[1]], method = "modal", ties = "lowest")
    expect_true(inherits(ens_one, "SpatRaster"))
    stk_runs <- c(sw_fn$rasters[[1]], sw_fn$rasters[[2]])
    ens_stk <- geomorphon_ensemble(stk_runs, method = "modal", ties = "lowest")
    ens_lst <- geomorphon_ensemble(list(sw_fn$rasters[[1]], sw_fn$rasters[[2]]), method = "modal", ties = "lowest")
    expect_equal(as.vector(terra::values(ens_stk)), as.vector(terra::values(ens_lst)))

    # Invalid ensemble inputs
    expect_error(geomorphon_ensemble(list()))
    expect_error(geomorphon_ensemble("invalid_input"))

    # --------------------------------------------------------------------------
    # Runtime benchmark assertions
    # Wall-clock timing is machine-dependent: run only when RGEOMORPHON_BENCH is set.
    # --------------------------------------------------------------------------
    if (nzchar(Sys.getenv("RGEOMORPHON_BENCH"))) {
        # Warmup
        invisible(rgeomorphon:::geomorphons_cpp_worker(
            salton, search = 40, skip = 0, flat_angle_deg = 1, dist = 0,
            comparison_mode = "anglev1", tdist = 0, use_meters = FALSE,
            x_res_dem = 30, y_res_dem = 30, forms = 10L, ternary = 0L,
            positive = 0L, negative = 0L, nodata = NA_integer_
        ))

        t0 <- proc.time()[["elapsed"]]
        for (i in 1:5) {
            invisible(rgeomorphon:::geomorphons_cpp_worker(
                salton, search = 40, skip = 0, flat_angle_deg = 1, dist = 0,
                comparison_mode = "anglev1", tdist = 0, use_meters = FALSE,
                x_res_dem = 30, y_res_dem = 30, forms = 10L, ternary = 0L,
                positive = 0L, negative = 0L, nodata = NA_integer_
            ))
        }
        t_cpp_1 <- (proc.time()[["elapsed"]] - t0) / 5

        t0 <- proc.time()[["elapsed"]]
        for (i in 1:5) {
            invisible(rgeomorphon:::geomorphons_cpp_worker(
                salton, search = 40, skip = 0, flat_angle_deg = 1, dist = 0,
                comparison_mode = "anglev1", tdist = 0, use_meters = FALSE,
                x_res_dem = 30, y_res_dem = 30, forms = c(10L, 6L, 5L, 4L), ternary = 0L,
                positive = 0L, negative = 0L, nodata = NA_integer_
            ))
        }
        t_cpp_4 <- (proc.time()[["elapsed"]] - t0) / 5

        ratio_cpp <- t_cpp_4 / max(t_cpp_1, 0.001)
        # C++ multi-form ray-tracing count reuse: 4 forms evaluates in < 1.4x of 1 form
        expect_true(ratio_cpp < 1.4)

        # 2. End-to-end evaluation benchmark
        # Single multi-form pass vs 4 individual calls
        invisible(geomorphons(dem, search = 30, forms = "forms10"))
        invisible(geomorphons(dem, search = 30, forms = c("forms10", "forms6", "forms5", "forms4")))

        t0 <- proc.time()[["elapsed"]]
        for (i in 1:3) {
            invisible(geomorphons(dem, search = 30, forms = c("forms10", "forms6", "forms5", "forms4")))
        }
        t_multi_pass <- (proc.time()[["elapsed"]] - t0) / 3

        t0 <- proc.time()[["elapsed"]]
        for (i in 1:3) {
            invisible(geomorphons(dem, search = 30, forms = "forms10"))
            invisible(geomorphons(dem, search = 30, forms = "forms6"))
            invisible(geomorphons(dem, search = 30, forms = "forms5"))
            invisible(geomorphons(dem, search = 30, forms = "forms4"))
        }
        t_separate_passes <- (proc.time()[["elapsed"]] - t0) / 3

        # Multi-form pass takes significantly less than 4 separate passes
        # (at least 20% faster than 4 separate passes)
        expect_true(t_multi_pass < 0.80 * t_separate_passes)
    }
}
