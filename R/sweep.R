#' Parameter Sweep for Geomorphon Terrain Classification
#'
#' Evaluates a grid of parameter combinations over a Digital Elevation Model
#' (DEM), taking advantage of line-of-sight ray-tracing count reuse across
#' multiple form taxonomies in each evaluation. Records execution runtime and
#' parameter metadata in a tidy summary table.
#'
#' @param elevation A `SpatRaster` or numeric `matrix` containing elevation data.
#' @param param_grid A `data.frame` or `list` containing parameter combinations.
#'   May contain columns or list elements for `search`, `skip`, `flat_angle_deg`
#'   (or `flat`), `dist`, `comparison_mode`, and `tdist`. If a list of vectors
#'   is supplied, it is expanded via [expand.grid()]. Any combinations where
#'   `skip >= search` are pruned with a warning.
#' @param forms character or numeric vector. Landform taxonomies to compute
#'   for each parameter combination. Default is
#'   `c("forms10", "forms6", "forms5", "forms4")`.
#' @param ... Additional arguments passed to [geomorphons()].
#'
#' @return An S3 object of class `geomorphon_sweep` containing:
#'   \item{summary}{A tidy `data.frame` with columns `run_id`, parameter values,
#'     `resolution`, `elapsed_sec` (runtime in seconds), and `n_cells`.}
#'   \item{rasters}{A named list of classification outputs (one per run_id).}
#'
#' @export
#' @rdname geomorphon_sweep
#'
#' @examples
#' \dontrun{
#' library(terra)
#' data("salton", package = "rgeomorphon")
#' dem <- rast(salton)
#' crs(dem) <- attr(salton, "crs")
#' ext(dem) <- attr(salton, "extent")
#'
#' grid <- list(
#'   search = c(3, 5),
#'   skip = c(0, 1),
#'   flat_angle_deg = 1.0
#' )
#'
#' sw <- geomorphon_sweep(dem, param_grid = grid, forms = c("forms10", "forms6"))
#' print(sw)
#' summary(sw)
#' }
geomorphon_sweep <- function(elevation,
                             param_grid,
                             forms = c("forms10", "forms6", "forms5", "forms4"),
                             ...) {
    if (missing(param_grid) || is.null(param_grid)) {
        stop("'param_grid' must be provided as a data.frame or list of parameters.", call. = FALSE)
    }

    # Normalize param_grid to data.frame
    if (is.data.frame(param_grid)) {
        grid_df <- as.data.frame(param_grid)
    } else if (is.list(param_grid)) {
        if (all(vapply(param_grid, is.atomic, logical(1)))) {
            grid_df <- expand.grid(param_grid, stringsAsFactors = FALSE)
        } else {
            grid_df <- do.call(rbind, lapply(param_grid, as.data.frame))
        }
    } else {
        stop("'param_grid' must be a data.frame or a list.", call. = FALSE)
    }

    if (nrow(grid_df) == 0) {
        stop("'param_grid' has 0 rows.", call. = FALSE)
    }

    # Parameter normalization and defaults
    col_names <- names(grid_df)
    if ("flat" %in% col_names && !"flat_angle_deg" %in% col_names) {
        grid_df$flat_angle_deg <- grid_df$flat
        grid_df$flat <- NULL
    }

    if (!"search" %in% names(grid_df)) {
        grid_df$search <- 3
    }
    if (!"skip" %in% names(grid_df)) {
        grid_df$skip <- 0
    }
    if (!"flat_angle_deg" %in% names(grid_df)) {
        grid_df$flat_angle_deg <- 1.0
    }
    if (!"dist" %in% names(grid_df)) {
        grid_df$dist <- 0.0
    }
    if (!"comparison_mode" %in% names(grid_df)) {
        grid_df$comparison_mode <- "anglev1"
    }
    if (!"tdist" %in% names(grid_df)) {
        grid_df$tdist <- 0.0
    }

    # Validate and prune invalid parameter combinations where skip >= search
    invalid <- grid_df$skip >= grid_df$search
    if (any(invalid)) {
        warning(
            sprintf(
                "Pruning %d invalid parameter combination(s) where skip >= search.",
                sum(invalid)
            ),
            call. = FALSE
        )
        grid_df <- grid_df[!invalid, , drop = FALSE]
    }

    if (nrow(grid_df) == 0) {
        stop("No valid parameter combinations remaining (all had skip >= search).", call. = FALSE)
    }

    rownames(grid_df) <- NULL
    n_runs <- nrow(grid_df)
    run_id <- seq_len(n_runs)

    # Resolution and cell count
    dots <- list(...)
    if (inherits(elevation, "SpatRaster")) {
        res_val <- terra::res(elevation)[1]
        n_cells_val <- terra::ncell(elevation)
        xres_arg <- NULL
        yres_arg <- NULL
    } else if (is.matrix(elevation)) {
        res_val <- if (!is.null(dots$xres)) as.numeric(dots$xres) else 1.0
        n_cells_val <- nrow(elevation) * ncol(elevation)
        xres_arg <- res_val
        yres_arg <- if (!is.null(dots$yres)) as.numeric(dots$yres) else res_val
    } else {
        res_val <- 1.0
        n_cells_val <- length(elevation)
        xres_arg <- NULL
        yres_arg <- NULL
    }

    rasters <- vector("list", n_runs)
    elapsed_sec <- numeric(n_runs)

    for (i in seq_len(n_runs)) {
        s  <- grid_df$search[i]
        sk <- grid_df$skip[i]
        fl <- grid_df$flat_angle_deg[i]
        di <- grid_df$dist[i]
        cm <- as.character(grid_df$comparison_mode[i])
        td <- grid_df$tdist[i]

        call_args <- list(
            elevation = elevation,
            search = s,
            skip = sk,
            flat_angle_deg = fl,
            dist = di,
            comparison_mode = cm,
            tdist = td,
            forms = forms
        )
        if (!is.null(xres_arg)) {
            call_args$xres <- xres_arg
            call_args$yres <- yres_arg
        }
        extra_args <- dots[setdiff(names(dots), names(call_args))]
        if (length(extra_args) > 0) {
            call_args <- c(call_args, extra_args)
        }
        if (!is.null(call_args$filename)) {
            # Give each run its own file so later runs do not overwrite earlier ones
            fn <- as.character(call_args$filename)
            m <- regexpr("\\.[^./]+$", fn)
            base <- if (m > 0) substr(fn, 1L, m - 1L) else fn
            ext <- if (m > 0) substr(fn, m, nchar(fn)) else ""
            call_args$filename <- sprintf("%s_run_%d%s", base, i, ext)
        }

        t0 <- proc.time()[["elapsed"]]
        res_i <- do.call(geomorphons, call_args)
        t1 <- proc.time()[["elapsed"]]

        elapsed_sec[i] <- t1 - t0
        rasters[[i]] <- res_i
    }

    names(rasters) <- paste0("run_", run_id)

    # Reorder columns: run_id, search, skip, flat_angle_deg, dist, comparison_mode, tdist, other columns..., resolution, elapsed_sec, n_cells
    std_params <- c("search", "skip", "flat_angle_deg", "dist", "comparison_mode", "tdist")
    other_params <- setdiff(names(grid_df), std_params)
    summary_df <- data.frame(
        run_id = run_id,
        grid_df[, c(std_params, other_params), drop = FALSE],
        resolution = res_val,
        elapsed_sec = elapsed_sec,
        n_cells = n_cells_val,
        stringsAsFactors = FALSE
    )

    structure(
        list(
            summary = summary_df,
            rasters = rasters
        ),
        class = "geomorphon_sweep"
    )
}

#' @export
print.geomorphon_sweep <- function(x, ...) {
    cat("<geomorphon_sweep>\n")
    n_runs <- nrow(x$summary)
    cat("Runs:", n_runs, "\n")
    if (n_runs > 0 && length(x$rasters) > 0) {
        first_r <- x$rasters[[1]]
        form_names <- if (inherits(first_r, "SpatRaster")) {
            names(first_r)
        } else if (is.list(first_r)) {
            names(first_r)
        } else {
            "forms"
        }
        cat("Forms:", paste(form_names, collapse = ", "), "\n")
        cat("Total elapsed:", round(sum(x$summary$elapsed_sec), 3), "sec\n")
        cat("Mean elapsed per run:", round(mean(x$summary$elapsed_sec), 3), "sec\n")
    }
    cat("\nSummary:\n")
    print(x$summary, ...)
    invisible(x)
}

#' @export
summary.geomorphon_sweep <- function(object, ...) {
    object$summary
}

#' Ensemble Aggregation of Geomorphon Classifications
#'
#' Computes cell-wise ensemble summaries across parameter sweep runs,
#' supporting deterministic modal classification and per-class proportion layers.
#'
#' @param sweep_obj An object of class `geomorphon_sweep`, a list of `SpatRaster`
#'   objects, or a multi-layer `SpatRaster` representing runs of a parameter sweep.
#' @param method character. Aggregation method: `"modal"` (default) computes the
#'   modal landform class per cell; `"proportions"` computes per-class proportions
#'   across non-NA runs.
#' @param form character or numeric. Which form taxonomy to aggregate when multiple
#'   exist in `sweep_obj` (e.g. `"forms10"`, `"forms6"`, `"forms5"`, `"forms4"`).
#'   Defaults to the first or only form available.
#' @param classes integer or character vector. Which landform classes to compute
#'   proportions for when `method = "proportions"`. Defaults to `NULL`, which evaluates
#'   all reachable classes for `form` (10 for `forms10`, 6 for `forms6`, 5 for
#'   `forms5`, 4 for `forms4`), matching [geomorphon_categories()].
#' @param ties character. Tie-breaking rule for `method = "modal"`. Passed to
#'   [terra::modal()], default is `"lowest"` for deterministic tie-breaking.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `SpatRaster` (or `matrix` / list of matrices) containing the aggregated
#'   ensemble. For `method = "modal"`, returns a single layer with geomorphon factor
#'   levels and color table applied. For `method = "proportions"`, returns multi-layer
#'   proportions named by class.
#'
#' @export
#' @rdname geomorphon_ensemble
#'
#' @examples
#' \dontrun{
#' library(terra)
#' data("salton", package = "rgeomorphon")
#' dem <- rast(salton)
#' crs(dem) <- attr(salton, "crs")
#' ext(dem) <- attr(salton, "extent")
#'
#' grid <- data.frame(search = c(3, 5), skip = c(0, 1))
#' sw <- geomorphon_sweep(dem, param_grid = grid, forms = "forms10")
#'
#' ens_mode <- geomorphon_ensemble(sw, method = "modal")
#' ens_prop <- geomorphon_ensemble(sw, method = "proportions")
#' }
geomorphon_ensemble <- function(sweep_obj,
                                method = c("modal", "proportions"),
                                form = NULL,
                                classes = NULL,
                                ties = "lowest",
                                ...) {
    method <- match.arg(method, c("modal", "proportions", "proportion"))
    if (method == "proportion") {
        method <- "proportions"
    }

    if (inherits(sweep_obj, "geomorphon_sweep")) {
        runs_list <- sweep_obj$rasters
    } else if (inherits(sweep_obj, "SpatRaster")) {
        runs_list <- list(sweep_obj)
    } else if (is.list(sweep_obj)) {
        runs_list <- sweep_obj
    } else {
        stop(
            "'sweep_obj' must be a 'geomorphon_sweep' object, a list of rasters/matrices, or a SpatRaster.",
            call. = FALSE
        )
    }

    if (length(runs_list) == 0) {
        stop("'sweep_obj' contains no runs.", call. = FALSE)
    }

    is_raster <- inherits(runs_list[[1]], "SpatRaster") || inherits(sweep_obj, "SpatRaster")

    if (is_raster) {
        if (inherits(sweep_obj, "SpatRaster") && length(runs_list) == 1 && is.null(sweep_obj$rasters)) {
            stack_runs <- sweep_obj
            form_name <- if (!is.null(form)) as.character(form) else "forms10"
        } else {
            first_r <- runs_list[[1]]
            avail <- names(first_r)
            if (is.null(form)) {
                if ("forms10" %in% avail) {
                    form_key <- "forms10"
                } else if ("form" %in% avail) {
                    form_key <- "form"
                } else {
                    form_key <- avail[1]
                }
            } else {
                if (is.numeric(form)) {
                    form_str <- paste0("forms", form)
                } else {
                    form_str <- as.character(form)
                }
                if (form_str %in% avail) {
                    form_key <- form_str
                } else if (as.character(form) %in% avail) {
                    form_key <- as.character(form)
                } else if (length(avail) == 1 && ("form" %in% avail || "geomorphon" %in% avail)) {
                    form_key <- avail[1]
                } else {
                    form_key <- form_str
                }
            }
            form_name <- form_key

            stack_runs <- terra::rast(lapply(runs_list, function(r) {
                if (form_key %in% names(r)) {
                    r[[form_key]]
                } else if (terra::nlyr(r) == 1) {
                    r[[1]]
                } else {
                    r[[form_key]]
                }
            }))
        }

        runs_num <- stack_runs
        levels(runs_num) <- NULL

        if (method == "modal") {
            m <- terra::modal(runs_num, ties = ties, na.rm = TRUE)
            n_valid <- sum(!is.na(runs_num))
            m[n_valid == 0] <- NA
            m[is.nan(m)] <- NA
            m <- geomorphon_theme(m, forms = form_name)
            names(m) <- "modal"
            return(m)
        } else {
            # method == "proportions"
            cats <- geomorphon_categories()
            if (is.null(classes)) {
                f_clean <- gsub("^forms", "", tolower(trimws(form_name)))
                if (f_clean == "6") {
                    classes_to_calc <- c(1L, 3L, 4L, 6L, 8L, 9L)
                } else if (f_clean == "5") {
                    classes_to_calc <- c(1L, 2L, 3L, 6L, 9L)
                } else if (f_clean == "4") {
                    classes_to_calc <- c(1L, 3L, 6L, 9L)
                } else {
                    classes_to_calc <- 1:10
                }
            } else if (is.numeric(classes)) {
                classes_to_calc <- as.integer(classes)
            } else if (is.character(classes)) {
                idx <- match(tolower(trimws(classes)), tolower(cats$form))
                classes_to_calc <- cats$value[idx[!is.na(idx)]]
            } else {
                classes_to_calc <- 1:10
            }

            n_valid <- sum(!is.na(runs_num))
            counts <- terra::rast(lapply(classes_to_calc, function(k) {
                sum(runs_num == k, na.rm = TRUE)
            }))
            prop <- counts / n_valid
            prop[n_valid == 0] <- NA
            names(prop) <- cats$form[match(classes_to_calc, cats$value)]
            return(prop)
        }
    } else {
        # Matrix workflow
        cats <- geomorphon_categories()
        first_r <- runs_list[[1]]
        if (is.list(first_r) && !is.data.frame(first_r)) {
            avail <- names(first_r)
            if (is.null(form)) {
                if ("forms10" %in% avail) {
                    form_key <- "forms10"
                } else if ("forms" %in% avail) {
                    form_key <- "forms"
                } else {
                    form_key <- avail[1]
                }
            } else {
                if (is.numeric(form)) {
                    form_str <- paste0("forms", form)
                } else {
                    form_str <- as.character(form)
                }
                if (form_str %in% avail) {
                    form_key <- form_str
                } else if (as.character(form) %in% avail) {
                    form_key <- as.character(form)
                } else if (length(avail) == 1) {
                    form_key <- avail[1]
                } else {
                    form_key <- form_str
                }
            }
            to_code_matrix <- function(val) {
                nr <- nrow(val)
                nc <- ncol(val)
                if (is.character(val)) {
                    matrix(match(val, cats$form), nrow = nr, ncol = nc)
                } else if (is.factor(val)) {
                    matrix(as.integer(val), nrow = nr, ncol = nc)
                } else {
                    matrix(as.integer(val), nrow = nr, ncol = nc)
                }
            }

            mats <- lapply(runs_list, function(r) {
                val <- if (form_key %in% names(r)) r[[form_key]] else r[[1]]
                to_code_matrix(val)
            })
        } else if (is.matrix(first_r)) {
            to_code_matrix <- function(val) {
                nr <- nrow(val)
                nc <- ncol(val)
                if (is.character(val)) {
                    matrix(match(val, cats$form), nrow = nr, ncol = nc)
                } else if (is.factor(val)) {
                    matrix(as.integer(val), nrow = nr, ncol = nc)
                } else {
                    matrix(as.integer(val), nrow = nr, ncol = nc)
                }
            }

            mats <- lapply(runs_list, function(r) {
                to_code_matrix(r)
            })
        } else {
            stop("Unsupported sweep_obj raster/matrix format.", call. = FALSE)
        }

        nr <- nrow(mats[[1]])
        nc <- ncol(mats[[1]])
        mat_runs <- do.call(cbind, lapply(mats, as.vector))

        if (method == "modal") {
            mode_vec <- apply(mat_runs, 1, function(vals) {
                v <- vals[!is.na(vals)]
                if (length(v) == 0) return(NA_integer_)
                tbl <- tabulate(v, nbins = 10)
                max_cnt <- max(tbl)
                if (max_cnt == 0) return(NA_integer_)
                which(tbl == max_cnt)[1] # lowest code
            })
            res_mat <- matrix(mode_vec, nrow = nr, ncol = nc)
            return(geomorphon_theme(res_mat))
        } else {
            cats <- geomorphon_categories()
            if (is.null(classes)) {
                classes_to_calc <- 1:10
            } else if (is.numeric(classes)) {
                classes_to_calc <- as.integer(classes)
            } else if (is.character(classes)) {
                idx <- match(tolower(trimws(classes)), tolower(cats$form))
                classes_to_calc <- cats$value[idx[!is.na(idx)]]
            } else {
                classes_to_calc <- 1:10
            }

            n_valid <- rowSums(!is.na(mat_runs))
            res_prop_list <- lapply(classes_to_calc, function(k) {
                cnt <- rowSums(mat_runs == k, na.rm = TRUE)
                p <- cnt / n_valid
                p[n_valid == 0] <- NA_real_
                matrix(p, nrow = nr, ncol = nc)
            })
            names(res_prop_list) <- cats$form[match(classes_to_calc, cats$value)]
            return(res_prop_list)
        }
    }
}
