# Internal base-2 Jensen-Shannon Divergence and Distance Helpers

.jsd_vec <- function(p, q) {
  if (length(p) != length(q) || length(p) == 0L) {
    return(NA_real_)
  }
  if (anyNA(p) || anyNA(q)) {
    return(NA_real_)
  }
  if (any(p < 0) || any(q < 0)) {
    return(NA_real_)
  }
  sum_p <- sum(p)
  sum_q <- sum(q)
  if (sum_p <= 0 || sum_q <= 0) {
    return(NA_real_)
  }
  p <- p / sum_p
  q <- q / sum_q
  m <- 0.5 * (p + q)
  
  kl_pm <- 0.0
  idx_p <- which(p > 0)
  if (length(idx_p) > 0L) {
    kl_pm <- sum(p[idx_p] * log2(p[idx_p] / m[idx_p]))
  }
  
  kl_qm <- 0.0
  idx_q <- which(q > 0)
  if (length(idx_q) > 0L) {
    kl_qm <- sum(q[idx_q] * log2(q[idx_q] / m[idx_q]))
  }
  
  jsd <- 0.5 * (kl_pm + kl_qm)
  # Clamp to [0, 1] against numerical floating-point inaccuracies
  if (jsd < 0) jsd <- 0.0
  if (jsd > 1) jsd <- 1.0
  jsd
}

.dist_pair <- function(p, q, method = "jensen-shannon", unit = "log2") {
  if (identical(tolower(method), "jensen-shannon") && identical(tolower(unit), "log2")) {
    return(.jsd_vec(p, q))
  }
  if (!requireNamespace("philentropy", quietly = TRUE)) {
    stop(sprintf("Package 'philentropy' is required for method = '%s' / unit = '%s'.", method, unit), call. = FALSE)
  }
  # philentropy::dist_one_one
  philentropy::dist_one_one(as.numeric(p), as.numeric(q), method = method, unit = unit)
}

.dist_matrix <- function(mat, method = "jensen-shannon", unit = "log2") {
  nr <- nrow(mat)
  if (is.null(nr) || nr == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }
  rnames <- rownames(mat)
  if (identical(tolower(method), "jensen-shannon") && identical(tolower(unit), "log2")) {
    dmat <- matrix(0, nrow = nr, ncol = nr)
    if (nr > 1L) {
      for (i in 1:(nr - 1L)) {
        pi <- mat[i, ]
        for (j in (i + 1L):nr) {
          val <- .jsd_vec(pi, mat[j, ])
          dmat[i, j] <- val
          dmat[j, i] <- val
        }
      }
    }
    rownames(dmat) <- rnames
    colnames(dmat) <- rnames
    return(dmat)
  }
  
  if (!requireNamespace("philentropy", quietly = TRUE)) {
    stop(sprintf("Package 'philentropy' is required for method = '%s' / unit = '%s'.", method, unit), call. = FALSE)
  }
  
  dmat <- philentropy::dist_many_many(mat, mat, method = method, unit = unit)
  rownames(dmat) <- rnames
  colnames(dmat) <- rnames
  dmat
}

#' Reference Attribute Table Target Layer Generator
#'
#' Derives aligned categorical raster layers from an integer ID raster and its
#' Reference Attribute Table (RAT), with support for categorical columns and
#' binned numeric columns. Attaches coverage statistics describing valid cell
#' proportions and class counts.
#'
#' @param x A single-layer categorical or integer `SpatRaster`.
#' @param attrs Optional `data.frame` containing the attribute table. If `NULL`
#'   (default), the table is extracted from `terra::cats(x)[[1]]`.
#' @param id_col Column index or column name in `attrs` containing the keys that
#'   match raster values in `x`. Default is `1`.
#' @param vars Character vector of column names in `attrs` to convert into target
#'   raster layers. Default is all columns except `id_col`.
#' @param breaks Named list of numeric break vectors for binning continuous
#'   attribute columns using [base::cut()].
#' @param labels Named list of character vectors of category labels corresponding
#'   to `breaks`.
#' @param right Logical. Indicating if the intervals should be closed on the right
#'   (and open on the left) or vice versa. Default is `TRUE`. Passed to [base::cut()].
#'
#' @return A multi-layer `SpatRaster` where each layer corresponds to a derived
#'   target variable with categories set via [terra::levels()]. The returned raster
#'   has an attached attribute `"coverage"`, a tidy `data.frame` containing:
#'   \item{target}{Target layer name.}
#'   \item{n_total}{Total number of cells.}
#'   \item{n_valid}{Number of valid non-NA cells.}
#'   \item{n_na}{Number of missing (NA) cells.}
#'   \item{prop_valid}{Proportion of valid cells (`n_valid / n_total`).}
#'   \item{n_classes}{Number of factor levels / distinct categories.}
#'
#' @export
#' @examplesIf requireNamespace("terra", quietly = TRUE)
#' library(terra)
#' r <- rast(matrix(c(1, 2, 1, 3), nrow = 2, ncol = 2))
#' rat <- data.frame(id = 1:3, landuse = c("Forest", "Water", "Forest"), slope = c(12, 1, 25))
#' targets <- geomorphon_targets(r, attrs = rat, id_col = "id", vars = "landuse")
#' attr(targets, "coverage")
geomorphon_targets <- function(x,
                               attrs = NULL,
                               id_col = 1,
                               vars = NULL,
                               breaks = list(),
                               labels = list(),
                               right = TRUE) {
  if (!inherits(x, "SpatRaster")) {
    stop("'x' must be a SpatRaster.", call. = FALSE)
  }
  if (terra::nlyr(x) != 1L) {
    stop("'x' must be a single-layer SpatRaster.", call. = FALSE)
  }
  
  if (is.null(attrs)) {
    c_list <- terra::cats(x)
    if (length(c_list) == 0L || is.null(c_list[[1]]) || nrow(c_list[[1]]) == 0L) {
      stop("No attribute table found in 'x' and 'attrs' was not provided.", call. = FALSE)
    }
    attrs <- c_list[[1]]
  } else {
    attrs <- as.data.frame(attrs)
  }
  
  if (nrow(attrs) == 0L) {
    stop("'attrs' must contain at least one row.", call. = FALSE)
  }
  
  if (is.numeric(id_col)) {
    if (id_col < 1 || id_col > ncol(attrs)) {
      stop("Column index 'id_col' is out of bounds.", call. = FALSE)
    }
    id_col_name <- names(attrs)[id_col]
  } else if (is.character(id_col)) {
    if (!id_col %in% names(attrs)) {
      stop(sprintf("Column '%s' specified by 'id_col' not found in 'attrs'.", id_col), call. = FALSE)
    }
    id_col_name <- id_col
  } else {
    stop("'id_col' must be an integer column index or column name.", call. = FALSE)
  }
  
  if (is.null(vars)) {
    vars <- setdiff(names(attrs), id_col_name)
    if (length(vars) == 0L) {
      stop("No target variables found in 'attrs' besides 'id_col'.", call. = FALSE)
    }
  } else {
    missing_vars <- setdiff(vars, names(attrs))
    if (length(missing_vars) > 0L) {
      stop(sprintf("Requested variable(s) not found in 'attrs': %s", paste(missing_vars, collapse = ", ")), call. = FALSE)
    }
  }
  
  id_vals <- terra::values(x, mat = FALSE)
  n_total <- length(id_vals)
  row_idx <- match(id_vals, attrs[[id_col_name]])
  
  layers <- vector("list", length(vars))
  cov_target <- character(length(vars))
  cov_n_total <- integer(length(vars))
  cov_n_valid <- integer(length(vars))
  cov_n_na <- integer(length(vars))
  cov_prop_valid <- numeric(length(vars))
  cov_n_classes <- integer(length(vars))
  
  for (i in seq_along(vars)) {
    v <- vars[i]
    col_v <- attrs[[v]]
    
    if (v %in% names(breaks)) {
      b <- breaks[[v]]
      lbl <- labels[[v]]
      col_v <- cut(col_v, breaks = b, labels = lbl, right = right, include.lowest = TRUE)
    }
    
    non_na_vals <- col_v[!is.na(col_v)]
    if (is.factor(col_v)) {
      present_levels <- levels(col_v)[levels(col_v) %in% as.character(non_na_vals)]
      lvls <- present_levels
    } else {
      lvls <- sort(unique(as.character(non_na_vals)))
    }
    
    code_map <- match(as.character(col_v), lvls)
    mapped_vals <- code_map[row_idx]
    
    lyr <- terra::rast(x)
    terra::values(lyr) <- mapped_vals
    names(lyr) <- v
    
    if (length(lvls) >= 1L) {
      levels(lyr) <- data.frame(id = seq_along(lvls), category = lvls, stringsAsFactors = FALSE)
    }
    
    n_val <- sum(!is.na(mapped_vals))
    n_miss <- sum(is.na(mapped_vals))
    
    layers[[i]] <- lyr
    cov_target[i] <- v
    cov_n_total[i] <- n_total
    cov_n_valid[i] <- n_val
    cov_n_na[i] <- n_miss
    cov_prop_valid[i] <- if (n_total > 0L) n_val / n_total else NA_real_
    cov_n_classes[i] <- length(lvls)
  }
  
  out <- if (length(layers) == 1L) layers[[1]] else do.call(c, layers)
  coverage_df <- data.frame(
    target = cov_target,
    n_total = cov_n_total,
    n_valid = cov_n_valid,
    n_na = cov_n_na,
    prop_valid = cov_prop_valid,
    n_classes = cov_n_classes,
    stringsAsFactors = FALSE
  )
  attr(out, "coverage") <- coverage_df
  out
}

#' Pairwise Dissimilarity of Concordance Profiles
#'
#' Computes the pairwise dissimilarity matrix between conditional probability
#' profiles from a contingency table or [geomorphon_concordance()] output using
#' exact base-2 Jensen-Shannon divergence (default) or any method supported by
#' `philentropy`.
#'
#' @param x A `geomorphon_concordance` result list (containing `contingency_table`)
#'   or a 2D contingency `table` or numeric `matrix`.
#' @param margin Integer indicating the profile direction: `1` for row profiles
#'   (e.g. \eqn{P(\text{Category} \mid \text{Geomorphon})}), or `2` for column
#'   profiles (e.g. \eqn{P(\text{Geomorphon} \mid \text{Category})}). Default is `1`.
#' @param method Character string specifying distance metric. Default is `"jensen-shannon"`.
#' @param unit Logarithm unit (`"log2"`, `"log"`, `"log10"`). Default is `"log2"`.
#'
#' @return A symmetric numeric distance matrix with row and column names matching
#'   the profile labels.
#'
#' @export
#' @examples
#' ct <- matrix(c(10, 0, 5, 5, 0, 10), nrow = 3, byrow = TRUE)
#' rownames(ct) <- c("Flat", "Slope", "Ridge")
#' colnames(ct) <- c("CatA", "CatB")
#' geomorphon_profile_dist(ct, margin = 1)
geomorphon_profile_dist <- function(x,
                                    margin = 1,
                                    method = "jensen-shannon",
                                    unit = "log2") {
  if (is.list(x) && !is.null(x$contingency_table)) {
    x <- x$contingency_table
  }
  if (!is.matrix(x) && !is.table(x)) {
    stop("'x' must be a matrix, table, or geomorphon_concordance result containing 'contingency_table'.", call. = FALSE)
  }
  if (!margin %in% c(1L, 2L)) {
    stop("'margin' must be 1 (rows) or 2 (columns).", call. = FALSE)
  }
  
  mat <- as.matrix(x)
  if (margin == 2L) {
    mat <- t(mat)
  }
  
  rsums <- rowSums(mat, na.rm = TRUE)
  # Normalize rows to probabilities
  prob_mat <- mat
  for (i in seq_len(nrow(prob_mat))) {
    if (rsums[i] > 0) {
      prob_mat[i, ] <- prob_mat[i, ] / rsums[i]
    } else {
      prob_mat[i, ] <- NA_real_
    }
  }
  
  .dist_matrix(prob_mat, method = method, unit = unit)
}

#' Cell-Wise Divergence Between Landform Proportion Ensembles
#'
#' Computes cell-wise divergence and summary statistics between two aligned
#' multi-layer class-proportion ensembles (e.g. from [geomorphon_ensemble()] with
#' `method = "proportions"`), or two lists of proportion matrices / numeric matrices.
#'
#' @param p A multi-layer class-proportion `SpatRaster` or numeric matrix / matrix list.
#' @param q A multi-layer class-proportion `SpatRaster` or numeric matrix / matrix list
#'   aligned with `p`.
#' @param method Character distance metric. Default is `"jensen-shannon"`.
#' @param unit Logarithm unit. Default is `"log2"`.
#'
#' @return A list containing:
#'   \item{raster}{A single-layer `SpatRaster` (or numeric matrix) named `"divergence"`.}
#'   \item{mean_divergence}{Arithmetic mean divergence across valid cells.}
#'   \item{median_divergence}{Median divergence across valid cells.}
#'   \item{q90_divergence}{90th percentile divergence across valid cells.}
#'   \item{n_cells}{Count of valid cells compared.}
#'
#' @export
#' @examplesIf requireNamespace("terra", quietly = TRUE)
#' library(terra)
#' p <- rast(array(c(0.8, 0.2, 0.5, 0.5), dim = c(2, 2, 2)))
#' q <- rast(array(c(0.2, 0.8, 0.5, 0.5), dim = c(2, 2, 2)))
#' names(p) <- names(q) <- c("Flat", "Slope")
#' div <- geomorphon_divergence(p, q)
#' div$mean_divergence
geomorphon_divergence <- function(p,
                                  q,
                                  method = "jensen-shannon",
                                  unit = "log2") {
  is_rast <- inherits(p, "SpatRaster") && inherits(q, "SpatRaster")
  
  if (is_rast) {
    if (!terra::compareGeom(p, q, stopOnError = FALSE, messages = FALSE)) {
      stop("Geometry of 'p' and 'q' does not match.", call. = FALSE)
    }
    
    p_names <- names(p)
    q_names <- names(q)
    cats_all <- geomorphon_categories()$form
    
    has_known_p <- any(p_names %in% cats_all)
    has_known_q <- any(q_names %in% cats_all)
    
    if ((has_known_p && has_known_q) || (length(intersect(p_names, q_names)) > 0L)) {
      union_names <- union(p_names, q_names)
      if (has_known_p || has_known_q) {
        union_names <- cats_all[cats_all %in% union_names]
        rem <- setdiff(union(p_names, q_names), union_names)
        if (length(rem) > 0L) union_names <- c(union_names, rem)
      }
      
      align_stack <- function(r, u_names) {
        r_names <- names(r)
        missing_names <- setdiff(u_names, r_names)
        if (length(missing_names) == 0L) {
          return(r[[u_names]])
        }
        valid_mask <- !is.na(r[[1]])
        zero_lyr <- terra::rast(r[[1]])
        terra::values(zero_lyr) <- ifelse(terra::values(valid_mask, mat = FALSE), 0, NA_real_)
        
        all_lyrs <- vector("list", length(u_names))
        for (idx in seq_along(u_names)) {
          nm <- u_names[idx]
          if (nm %in% r_names) {
            all_lyrs[[idx]] <- r[[nm]]
          } else {
            all_lyrs[[idx]] <- zero_lyr
            names(all_lyrs[[idx]]) <- nm
          }
        }
        do.call(c, all_lyrs)
      }
      
      p_aligned <- align_stack(p, union_names)
      q_aligned <- align_stack(q, union_names)
    } else {
      if (terra::nlyr(p) != terra::nlyr(q)) {
        stop("Number of layers in 'p' and 'q' must match when layer names are disjoint or unrecognised.", call. = FALSE)
      }
      p_aligned <- p
      q_aligned <- q
    }
    
    mat_p <- terra::values(p_aligned, mat = TRUE)
    mat_q <- terra::values(q_aligned, mat = TRUE)
    n_cells_total <- nrow(mat_p)
    
    res_vec <- numeric(n_cells_total)
    
    if (identical(tolower(method), "jensen-shannon") && identical(tolower(unit), "log2")) {
      for (k in seq_len(n_cells_total)) {
        res_vec[k] <- .jsd_vec(mat_p[k, ], mat_q[k, ])
      }
    } else {
      for (k in seq_len(n_cells_total)) {
        pk <- mat_p[k, ]
        qk <- mat_q[k, ]
        if (anyNA(pk) || anyNA(qk) || sum(pk) <= 0 || sum(qk) <= 0) {
          res_vec[k] <- NA_real_
        } else {
          res_vec[k] <- .dist_pair(pk, qk, method = method, unit = unit)
        }
      }
    }
    
    out_r <- terra::rast(p[[1]])
    terra::values(out_r) <- res_vec
    names(out_r) <- "divergence"
    
    valid_vals <- res_vec[!is.na(res_vec)]
    n_valid <- length(valid_vals)
    
    return(list(
      raster = out_r,
      mean_divergence = if (n_valid > 0L) mean(valid_vals) else NA_real_,
      median_divergence = if (n_valid > 0L) stats::median(valid_vals) else NA_real_,
      q90_divergence = if (n_valid > 0L) stats::quantile(valid_vals, probs = 0.90, names = FALSE) else NA_real_,
      n_cells = n_valid
    ))
  }
  
  # Matrix input fallback
  if (is.list(p)) {
    p <- do.call(cbind, lapply(p, as.vector))
  }
  if (is.list(q)) {
    q <- do.call(cbind, lapply(q, as.vector))
  }
  
  mat_p <- as.matrix(p)
  mat_q <- as.matrix(q)
  if (nrow(mat_p) != nrow(mat_q) || ncol(mat_p) != ncol(mat_q)) {
    stop("Dimensions of 'p' and 'q' must match.", call. = FALSE)
  }
  
  n_rows <- nrow(mat_p)
  res_vec <- numeric(n_rows)
  if (identical(tolower(method), "jensen-shannon") && identical(tolower(unit), "log2")) {
    for (k in seq_len(n_rows)) {
      res_vec[k] <- .jsd_vec(mat_p[k, ], mat_q[k, ])
    }
  } else {
    for (k in seq_len(n_rows)) {
      pk <- mat_p[k, ]
      qk <- mat_q[k, ]
      if (anyNA(pk) || anyNA(qk) || sum(pk) <= 0 || sum(qk) <= 0) {
        res_vec[k] <- NA_real_
      } else {
        res_vec[k] <- .dist_pair(pk, qk, method = method, unit = unit)
      }
    }
  }
  
  valid_vals <- res_vec[!is.na(res_vec)]
  n_valid <- length(valid_vals)
  
  list(
    matrix = res_vec,
    mean_divergence = if (n_valid > 0L) mean(valid_vals) else NA_real_,
    median_divergence = if (n_valid > 0L) stats::median(valid_vals) else NA_real_,
    q90_divergence = if (n_valid > 0L) stats::quantile(valid_vals, probs = 0.90, names = FALSE) else NA_real_,
    n_cells = n_valid
  )
}

#' Multi-Objective Pareto Frontier and Composite Calibration Score
#'
#' Identifies non-dominated Pareto-optimal parameter sets and computes composite
#' harmonic-mean scores over multi-objective evaluation tables.
#'
#' @param cal_df A `data.frame` containing calibration evaluation metrics (such as
#'   the output of [geomorphon_calibrate()]).
#' @param objectives Named character vector defining optimization directions for
#'   each objective (e.g. `c(delta_nmi = "max", delta_edge_f1 = "max", edge_mhd = "min")`).
#'   Values must be `"max"` or `"min"`. If `NULL`, default objectives are selected
#'   based on available columns.
#' @param by Character vector of column names defining independent Pareto grouping
#'   domains. Default is `"target"`.
#' @param min_edge_ratio Minimum permissible `edge_ratio` (`n_edges_geo / n_edges_cat`).
#'   Default is `0`.
#' @param max_edge_ratio Maximum permissible `edge_ratio`. Default is `Inf`.
#'
#' @return The input `cal_df` augmented with:
#'   \item{is_pareto}{Logical. `TRUE` if the candidate belongs to the non-dominated
#'     Pareto frontier within its group.}
#'   \item{pareto_score}{Numeric in \eqn{[0, 1]}. Normalized composite score
#'     computed as the harmonic mean across normalized objectives.}
#'
#' @export
#' @examples
#' df <- data.frame(
#'   target = "musym",
#'   run_id = 1:3,
#'   edge_ratio = c(1.1, 1.2, 0.9),
#'   delta_nmi = c(0.12, 0.15, 0.10),
#'   delta_edge_f1 = c(0.20, 0.18, 0.22),
#'   edge_mhd = c(45, 40, 50)
#' )
#' geomorphon_pareto(df)
geomorphon_pareto <- function(cal_df,
                              objectives = NULL,
                              by = "target",
                              min_edge_ratio = 0,
                              max_edge_ratio = Inf) {
  if (!is.data.frame(cal_df)) {
    stop("'cal_df' must be a data.frame.", call. = FALSE)
  }
  if (nrow(cal_df) == 0L) {
    cal_df$is_pareto <- logical(0)
    cal_df$pareto_score <- numeric(0)
    return(cal_df)
  }
  
  if (is.null(objectives)) {
    has_finite <- function(col) {
      col %in% names(cal_df) && any(is.finite(cal_df[[col]]))
    }
    if (has_finite("delta_nmi") && has_finite("delta_edge_f1") && has_finite("edge_mhd")) {
      objectives <- c(delta_nmi = "max", delta_edge_f1 = "max", edge_mhd = "min")
    } else if (has_finite("delta_nmi") && has_finite("delta_edge_f1")) {
      objectives <- c(delta_nmi = "max", delta_edge_f1 = "max")
    } else if (has_finite("nmi") && has_finite("edge_f1") && has_finite("edge_mhd")) {
      objectives <- c(nmi = "max", edge_f1 = "max", edge_mhd = "min")
    } else if (has_finite("nmi") && has_finite("edge_f1")) {
      objectives <- c(nmi = "max", edge_f1 = "max")
    } else if (has_finite("delta_nmi")) {
      objectives <- c(delta_nmi = "max")
    } else if (has_finite("nmi")) {
      objectives <- c(nmi = "max")
    } else if (has_finite("delta_edge_f1")) {
      objectives <- c(delta_edge_f1 = "max")
    } else if (has_finite("edge_f1")) {
      objectives <- c(edge_f1 = "max")
    } else {
      stop("Could not determine default objectives from 'cal_df'. Please provide 'objectives'.", call. = FALSE)
    }
  }
  
  obj_names <- names(objectives)
  if (is.null(obj_names) || any(obj_names == "")) {
    stop("'objectives' must be a named vector, e.g. c(delta_nmi = 'max', edge_mhd = 'min').", call. = FALSE)
  }
  missing_objs <- setdiff(obj_names, names(cal_df))
  if (length(missing_objs) > 0L) {
    stop(sprintf("Objectives not found in 'cal_df': %s", paste(missing_objs, collapse = ", ")), call. = FALSE)
  }
  
  directions <- tolower(trimws(as.character(objectives)))
  invalid_dirs <- setdiff(directions, c("max", "min"))
  if (length(invalid_dirs) > 0L) {
    stop("Objective directions must be 'max' or 'min'.", call. = FALSE)
  }
  
  by <- intersect(by, names(cal_df))
  if (length(by) == 0L) {
    group_ids <- rep(1L, nrow(cal_df))
  } else {
    group_ids <- interaction(cal_df[by], drop = TRUE)
  }
  
  is_pareto <- rep(FALSE, nrow(cal_df))
  pareto_score <- rep(NA_real_, nrow(cal_df))
  
  # Check eligibility
  has_edge_ratio <- "edge_ratio" %in% names(cal_df)
  ratio_ok <- rep(TRUE, nrow(cal_df))
  if (has_edge_ratio) {
    er <- cal_df$edge_ratio
    ratio_ok <- !is.na(er) & (er >= min_edge_ratio) & (er <= max_edge_ratio)
  }
  
  objs_finite <- rep(TRUE, nrow(cal_df))
  for (nm in obj_names) {
    objs_finite <- objs_finite & is.finite(cal_df[[nm]])
  }
  
  eligible <- ratio_ok & objs_finite
  
  for (grp in unique(group_ids)) {
    idx <- which(group_ids == grp)
    idx_elig <- idx[eligible[idx]]
    
    if (length(idx_elig) == 0L) {
      next
    }
    
    # Extract objective matrix for eligible candidates
    m <- as.matrix(cal_df[idx_elig, obj_names, drop = FALSE])
    
    # Transform to all "max"
    for (j in seq_along(obj_names)) {
      if (directions[j] == "min") {
        m[, j] <- -m[, j]
      }
    }
    
    n_elig <- length(idx_elig)
    pareto_elig <- rep(TRUE, n_elig)
    
    if (n_elig > 1L) {
      for (i in seq_len(n_elig)) {
        for (j in seq_len(n_elig)) {
          if (i == j) next
          # Does candidate j dominate candidate i?
          # j is at least as good on all objectives, and strictly better on >= 1
          if (all(m[j, ] >= m[i, ]) && any(m[j, ] > m[i, ])) {
            pareto_elig[i] <- FALSE
            break
          }
        }
      }
    }
    is_pareto[idx_elig] <- pareto_elig
    
    # Compute pareto_score: harmonic mean of normalized objectives
    # Normalization: for each objective, map values in this group to [0, 1]
    norm_m <- matrix(0, nrow = length(idx), ncol = length(obj_names))
    for (j in seq_along(obj_names)) {
      raw_vals <- cal_df[[obj_names[j]]][idx]
      vals_elig <- cal_df[[obj_names[j]]][idx_elig]
      min_v <- min(vals_elig)
      max_v <- max(vals_elig)
      
      if (max_v == min_v) {
        norm_v <- rep(1, length(idx))
      } else {
        if (directions[j] == "max") {
          norm_v <- (raw_vals - min_v) / (max_v - min_v)
        } else {
          norm_v <- (max_v - raw_vals) / (max_v - min_v)
        }
      }
      norm_m[, j] <- norm_v
    }
    
    # Harmonic mean across objectives for each row in group
    # H = K / sum(1 / (norm_v + eps))
    # If any normalized value <= 0, harmonic mean is 0 (or close to 0)
    eps <- 1e-6
    h_scores <- numeric(length(idx))
    for (r in seq_along(idx)) {
      row_idx_global <- idx[r]
      if (!eligible[row_idx_global]) {
        h_scores[r] <- NA_real_
      } else {
        nv <- norm_m[r, ]
        if (any(nv < 0)) nv[nv < 0] <- 0
        if (any(nv == 0)) {
          h_scores[r] <- 0.0
        } else {
          h_scores[r] <- length(nv) / sum(1.0 / nv)
        }
      }
    }
    pareto_score[idx] <- h_scores
  }
  
  cal_df$is_pareto <- is_pareto
  cal_df$pareto_score <- pareto_score
  cal_df
}

#' Calibration of Geomorphon Parameters Against Reference Targets
#'
#' Evaluates spatial concordance between multiple runs of a [geomorphon_sweep()]
#' and reference categorical target raster layers. Combines sweep parameter
#' metadata with concordance metrics, boundary-density penalties, and
#' Pareto-optimal frontier identification.
#'
#' @param sweep_obj An object of class `geomorphon_sweep`.
#' @param targets A single- or multi-layer categorical `SpatRaster` representing
#'   reference classification targets (e.g. from [geomorphon_targets()]).
#' @param forms Character vector of form taxonomies to evaluate (e.g.
#'   `c("forms10", "forms6")`). Default is `NULL`, which evaluates all taxonomies
#'   found in `sweep_obj`.
#' @param method Character concordance method: `"both"` (default), `"nmi"`, or
#'   `"boundary"`. Passed to [geomorphon_concordance()].
#' @param threshold_cells Integer or numeric cell distance tolerance for boundary
#'   metrics. Default is `2`.
#' @param spatial_null Logical. Whether to evaluate spatial-null circular shift
#'   ensembles in [geomorphon_concordance()]. Default is `TRUE`.
#'
#' @return A tidy `data.frame` combining sweep parameter columns (`run_id`, `search`,
#'   `skip`, `flat_angle_deg`, `dist`, `comparison_mode`, `tdist`, `resolution`,
#'   `elapsed_sec`, `n_cells`), `forms`, `target`, concordance metrics,
#'   `edge_ratio`, `density_penalty`, `is_pareto`, and `pareto_score`.
#'
#' @export
#' @examplesIf requireNamespace("terra", quietly = TRUE)
#' library(terra)
#' data(salton, package = "rgeomorphon")
#' dem <- rast(salton)
#' crs(dem) <- attr(salton, "crs")
#' ext(dem) <- attr(salton, "extent")
#' sw <- geomorphon_sweep(dem, param_grid = list(search = c(3, 4)), forms = "forms6")
#' cat_rast <- (dem > 50) + 1
#' levels(cat_rast) <- data.frame(id = 1:2, category = c("Low", "High"))
#' names(cat_rast) <- "elevation_class"
#' cal <- geomorphon_calibrate(sw, targets = cat_rast, spatial_null = FALSE)
#' print(cal)
geomorphon_calibrate <- function(sweep_obj,
                                 targets,
                                 forms = NULL,
                                 method = c("both", "nmi", "boundary"),
                                 threshold_cells = 2,
                                 spatial_null = TRUE) {
  if (!inherits(sweep_obj, "geomorphon_sweep")) {
    stop("'sweep_obj' must be an object of class 'geomorphon_sweep'.", call. = FALSE)
  }
  if (!inherits(targets, "SpatRaster")) {
    stop("'targets' must be a SpatRaster.", call. = FALSE)
  }
  method <- match.arg(method)
  
  summary_df <- sweep_obj$summary
  rasters_list <- sweep_obj$rasters
  n_runs <- length(rasters_list)
  if (n_runs == 0L) {
    stop("'sweep_obj' contains no rasters.", call. = FALSE)
  }
  
  target_names <- names(targets)
  if (is.null(target_names) || length(target_names) == 0L) {
    target_names <- paste0("target_", seq_len(terra::nlyr(targets)))
  }
  
  # Determine forms to evaluate
  prov_forms <- sweep_obj$provenance$forms
  avail_taxonomies <- c("forms10", "forms6", "forms5", "forms4")
  
  if (!is.null(forms)) {
    target_forms <- as.character(forms)
  } else if (!is.null(prov_forms)) {
    target_forms <- as.character(prov_forms)
  } else {
    first_r <- rasters_list[[1]]
    first_names <- names(first_r)
    matched_forms <- first_names[first_names %in% avail_taxonomies]
    if (length(matched_forms) > 0L) {
      target_forms <- matched_forms
    } else {
      target_forms <- "forms10"
    }
  }
  
  res_list <- list()
  
  for (i in seq_len(n_runs)) {
    run_r <- rasters_list[[i]]
    run_names <- names(run_r)
    param_row <- summary_df[i, , drop = FALSE]
    
    for (f in target_forms) {
      # Resolve layer key for form f in run_r
      form_clean <- if (is.numeric(f)) paste0("forms", f) else as.character(f)
      
      if (form_clean %in% run_names) {
        g_lyr <- run_r[[form_clean]]
      } else if (f %in% run_names) {
        g_lyr <- run_r[[f]]
      } else if (length(run_names) == 1L && run_names %in% c("form", "forms", "geomorphon")) {
        g_lyr <- run_r[[1]]
      } else if (length(run_names) == 1L) {
        g_lyr <- run_r[[1]]
      } else {
        next
      }
      
      for (t_idx in seq_along(target_names)) {
        t_nm <- target_names[t_idx]
        t_lyr <- targets[[t_idx]]
        
        conc <- geomorphon_concordance(
          g_lyr,
          t_lyr,
          method = method,
          threshold_cells = threshold_cells,
          spatial_null = spatial_null
        )
        
        n_edges_g <- if (!is.null(conc$n_edges_geo)) conc$n_edges_geo else NA_integer_
        n_edges_c <- if (!is.null(conc$n_edges_cat)) conc$n_edges_cat else NA_integer_
        
        e_ratio <- if (!is.na(n_edges_c) && n_edges_c > 0) {
          n_edges_g / n_edges_c
        } else {
          NA_real_
        }
        
        d_penalty <- if (!is.na(e_ratio) && e_ratio > 0) {
          exp(-abs(log(e_ratio)))
        } else {
          NA_real_
        }
        
        row_df <- data.frame(
          run_id = if (!is.null(param_row$run_id)) param_row$run_id else i,
          search = if (!is.null(param_row$search)) param_row$search else NA_real_,
          skip = if (!is.null(param_row$skip)) param_row$skip else NA_real_,
          flat_angle_deg = if (!is.null(param_row$flat_angle_deg)) param_row$flat_angle_deg else NA_real_,
          dist = if (!is.null(param_row$dist)) param_row$dist else NA_real_,
          comparison_mode = if (!is.null(param_row$comparison_mode)) as.character(param_row$comparison_mode) else NA_character_,
          tdist = if (!is.null(param_row$tdist)) param_row$tdist else NA_real_,
          resolution = if (!is.null(param_row$resolution)) param_row$resolution else NA_real_,
          elapsed_sec = if (!is.null(param_row$elapsed_sec)) param_row$elapsed_sec else NA_real_,
          n_cells = if (!is.null(param_row$n_cells)) param_row$n_cells else NA_integer_,
          forms = form_clean,
          target = t_nm,
          nmi = if (!is.null(conc$nmi)) conc$nmi else NA_real_,
          homogeneity = if (!is.null(conc$homogeneity)) conc$homogeneity else NA_real_,
          completeness = if (!is.null(conc$completeness)) conc$completeness else NA_real_,
          n_edges_geo = n_edges_g,
          n_edges_cat = n_edges_c,
          edge_ratio = e_ratio,
          density_penalty = d_penalty,
          edge_offset = if (!is.null(conc$edge_offset)) conc$edge_offset else NA_real_,
          edge_concordance = if (!is.null(conc$edge_concordance)) conc$edge_concordance else NA_real_,
          edge_f1 = if (!is.null(conc$edge_f1)) conc$edge_f1 else NA_real_,
          edge_mhd = if (!is.null(conc$edge_mhd)) conc$edge_mhd else NA_real_,
          edge_chamfer = if (!is.null(conc$edge_chamfer)) conc$edge_chamfer else NA_real_,
          mean_null_nmi = if (!is.null(conc$mean_null_nmi)) conc$mean_null_nmi else NA_real_,
          delta_nmi = if (!is.null(conc$delta_nmi)) conc$delta_nmi else NA_real_,
          nmi_null_shifts = if (!is.null(conc$nmi_null_shifts)) conc$nmi_null_shifts else NA_integer_,
          mean_null_edge_f1 = if (!is.null(conc$mean_null_edge_f1)) conc$mean_null_edge_f1 else NA_real_,
          delta_edge_f1 = if (!is.null(conc$delta_edge_f1)) conc$delta_edge_f1 else NA_real_,
          edge_f1_null_shifts = if (!is.null(conc$edge_f1_null_shifts)) conc$edge_f1_null_shifts else NA_integer_,
          null_shifts_evaluated = if (!is.null(conc$null_shifts_evaluated)) conc$null_shifts_evaluated else NA_integer_,
          stringsAsFactors = FALSE
        )
        res_list[[length(res_list) + 1L]] <- row_df
      }
    }
  }
  
  if (length(res_list) == 0L) {
    stop("No calibration evaluations could be performed.", call. = FALSE)
  }
  
  cal_df <- do.call(rbind, res_list)
  rownames(cal_df) <- NULL
  
  geomorphon_pareto(cal_df, by = "target")
}
