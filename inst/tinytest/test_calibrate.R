# inst/tinytest/test_calibrate.R
# Unit tests for Reference Attribute Table (RAT) targets, base-2 JSD,
# divergence, and multi-objective Pareto calibration.

library(tinytest)
library(rgeomorphon)

if (requireNamespace("terra", quietly = TRUE)) {

# -----------------------------------------------------------------------------
# 1. geomorphon_targets()
# -----------------------------------------------------------------------------

# Load packaged scorpion fixture
data("scorpion", package = "rgeomorphon")

soil <- terra::rast(scorpion$mukey)
terra::crs(soil) <- attr(scorpion, "crs")
terra::ext(soil) <- attr(scorpion, "extent")
levels(soil) <- scorpion$mapunit

# Derive targets for categorical columns + binned numeric column
targets <- geomorphon_targets(
  soil,
  vars = c("musym", "drclassdcd", "mukind", "hydgrpdcd", "slopegraddcp"),
  breaks = list(slopegraddcp = c(-Inf, 2, 5, 15, 30, Inf))
)

expect_inherits(targets, "SpatRaster")
expect_equal(terra::nlyr(targets), 5L)

cov <- attr(targets, "coverage")
expect_true(is.data.frame(cov))
expect_equal(nrow(cov), 5L)
expect_equal(cov$target, c("musym", "drclassdcd", "mukind", "hydgrpdcd", "slopegraddcp"))
expect_true(all(cov$n_total == 79665L))

# Deterministic NA propagation:
# In scorpion$mapunit, mukey = 2462630 (Water) has drclassdcd = NA and aws0100wta = NA.
# Therefore, drclassdcd must have n_valid < n_total and n_na > 0.
expect_true(cov$n_valid[cov$target == "drclassdcd"] < 79665L)
expect_true(cov$n_na[cov$target == "drclassdcd"] > 0L)
expect_equal(cov$n_na[cov$target == "drclassdcd"], 1135L)
expect_equal(cov$n_classes[cov$target == "drclassdcd"], 2L)

# musym covers all 7 mapunits (Water has musym = "W")
expect_equal(cov$n_classes[cov$target == "musym"], 7L)

# slopegraddcp binned into 4 present USDA slope categories
expect_equal(cov$n_classes[cov$target == "slopegraddcp"], 4L)

# Factor levels are attached to layers
cats_list <- terra::cats(targets)
expect_equal(length(cats_list), 5L)
expect_true("Well drained" %in% cats_list[[2]]$category)
expect_true("Moderately well drained" %in% cats_list[[2]]$category)

# Error handling: missing attribute table
r_no_rat <- terra::rast(matrix(1:9, 3, 3))
expect_error(geomorphon_targets(r_no_rat), "No attribute table")

# Error handling: missing requested variable in attrs
expect_error(
  geomorphon_targets(soil, vars = c("nonexistent_column")),
  "Requested variable"
)

# Explicit attrs data.frame override
custom_rat <- data.frame(
  id = c(2462630, 2600529),
  zone = c("wetland", "upland")
)
tgt_custom <- geomorphon_targets(soil, attrs = custom_rat, id_col = "id", vars = "zone")
expect_equal(terra::nlyr(tgt_custom), 1L)
cov_custom <- attr(tgt_custom, "coverage")
expect_equal(cov_custom$n_classes, 2L)


# -----------------------------------------------------------------------------
# 2. geomorphon_profile_dist()
# -----------------------------------------------------------------------------

# Hand-calculated 3-row contingency matrix:
# Row 1: (1, 1, 0) -> normalized: (0.5, 0.5, 0.0)
# Row 2: (2, 2, 0) -> normalized: (0.5, 0.5, 0.0) -> identical to Row 1 (D_JS = 0.0)
# Row 3: (0, 0, 5) -> normalized: (0.0, 0.0, 1.0) -> disjoint support from Row 1 & 2 (D_JS = 1.0)
m3 <- matrix(c(
  1, 1, 0,
  2, 2, 0,
  0, 0, 5
), nrow = 3, byrow = TRUE, dimnames = list(c("r1", "r2", "r3"), c("c1", "c2", "c3")))

d3 <- geomorphon_profile_dist(m3, margin = 1)
expect_equal(nrow(d3), 3L)
expect_equal(ncol(d3), 3L)
expect_equal(rownames(d3), c("r1", "r2", "r3"))
expect_equal(colnames(d3), c("r1", "r2", "r3"))

# Identical conditional distributions: D_JS == 0
expect_equal(d3["r1", "r2"], 0.0)
expect_equal(d3["r2", "r1"], 0.0)
expect_equal(diag(d3), c(r1 = 0.0, r2 = 0.0, r3 = 0.0))

# Disjoint support: D_JS == 1.0 exactly (verifying 0 * log2(0) = 0 without epsilon drift)
expect_equal(d3["r1", "r3"], 1.0)
expect_equal(d3["r2", "r3"], 1.0)

# Intermediate analytical pair:
# Row A: (1, 1, 0) -> p_A = (0.5, 0.5, 0)
# Row B: (0, 2, 0) -> p_B = (0.0, 1.0, 0)
# Midpoint m = (0.25, 0.75, 0)
# KL(p_A || m) = 0.5 * log2(0.5 / 0.25) + 0.5 * log2(0.5 / 0.75) = 0.5 - 0.5 * log2(1.5)
# KL(p_B || m) = 1.0 * log2(1.0 / 0.75) = log2(4/3)
# Analytical base-2 JSD = 1.5 - 0.75 * log2(3) ≈ 0.3112781244591328
m_inter <- matrix(c(
  1, 1, 0,
  2, 2, 0,
  0, 2, 0
), nrow = 3, byrow = TRUE)
d_inter <- geomorphon_profile_dist(m_inter, margin = 1)
jsd_analytical <- 1.5 - 0.75 * log2(3)
expect_equal(d_inter[1, 2], 0.0)
expect_true(abs(d_inter[1, 3] - jsd_analytical) < 1e-12)
expect_true(abs(d_inter[2, 3] - jsd_analytical) < 1e-12)

# Margin = 2 (columns normalized across rows)
d_col <- geomorphon_profile_dist(m3, margin = 2)
expect_equal(nrow(d_col), 3L)
expect_equal(colnames(d_col), c("c1", "c2", "c3"))
# Col 1 (1, 2, 0) and Col 2 (1, 2, 0) have identical distributions across rows -> JSD = 0
expect_equal(d_col["c1", "c2"], 0.0)
# Col 1 and Col 3 have disjoint support -> JSD = 1.0
expect_equal(d_col["c1", "c3"], 1.0)

# Direct acceptance of geomorphon_concordance() result object
synth_g <- terra::rast(matrix(c(1, 1, 2, 2), 2, 2))
synth_c <- terra::rast(matrix(c(1, 2, 1, 2), 2, 2))
conc_res <- geomorphon_concordance(synth_g, synth_c, method = "nmi")
d_from_conc <- geomorphon_profile_dist(conc_res, margin = 1)
expect_true(is.matrix(d_from_conc))
expect_equal(nrow(d_from_conc), nrow(conc_res$contingency_table))

# Optional philentropy dispatch
if (requireNamespace("philentropy", quietly = TRUE)) {
  d_hellinger <- geomorphon_profile_dist(m3, method = "hellinger")
  expect_equal(nrow(d_hellinger), 3L)
  expect_equal(d_hellinger["r1", "r2"], 0.0)
  expect_true(d_hellinger["r1", "r3"] > 0)

  d_manhattan <- geomorphon_profile_dist(m3, method = "manhattan")
  expect_equal(nrow(d_manhattan), 3L)
  expect_equal(d_manhattan["r1", "r2"], 0.0)
  expect_equal(d_manhattan["r1", "r3"], 2.0)
}


# -----------------------------------------------------------------------------
# 3. geomorphon_divergence()
# -----------------------------------------------------------------------------

# Synthetic proportion SpatRasters:
# Stack 1: Identical proportion layers
p_mat1 <- matrix(c(0.7, 0.3, 0.5, 0.5), 2, 2)
p_mat2 <- matrix(c(0.3, 0.7, 0.5, 0.5), 2, 2)
r_p1 <- terra::rast(p_mat1)
r_p2 <- terra::rast(p_mat2)
stack_a <- c(r_p1, r_p2)
names(stack_a) <- c("class1", "class2")

# Identical stacks yield 0 divergence everywhere
div_zero <- geomorphon_divergence(stack_a, stack_a)
expect_inherits(div_zero$raster, "SpatRaster")
expect_equal(div_zero$mean_divergence, 0.0)
expect_equal(div_zero$median_divergence, 0.0)
expect_equal(div_zero$q90_divergence, 0.0)
expect_equal(as.vector(terra::values(div_zero$raster)), rep(0.0, 4L))

# Disjoint one-hot stacks yield 1.0 divergence everywhere
stack_onehot1 <- c(terra::rast(matrix(1, 2, 2)), terra::rast(matrix(0, 2, 2)))
names(stack_onehot1) <- c("class1", "class2")
stack_onehot2 <- c(terra::rast(matrix(0, 2, 2)), terra::rast(matrix(1, 2, 2)))
names(stack_onehot2) <- c("class1", "class2")

div_one <- geomorphon_divergence(stack_onehot1, stack_onehot2)
expect_equal(div_one$mean_divergence, 1.0)
expect_equal(div_one$median_divergence, 1.0)
expect_equal(as.vector(terra::values(div_one$raster)), rep(1.0, 4L))

# Matrix input support
mat_p <- matrix(c(0.5, 0.5, 0.0, 1.0), nrow = 2, byrow = TRUE)
mat_q <- matrix(c(0.5, 0.5, 1.0, 0.0), nrow = 2, byrow = TRUE)
div_mat <- geomorphon_divergence(mat_p, mat_q)
expect_true(is.numeric(div_mat$matrix) || is.matrix(div_mat$matrix))
expect_equal(div_mat$n_cells, 2L)
expect_equal(div_mat$matrix[1], 0.0)
expect_equal(div_mat$matrix[2], 1.0)

# forms10 vs forms6 proportion stacks align on 10-class union cleanly
dem_full <- terra::rast(scorpion$dem)
terra::crs(dem_full) <- attr(scorpion, "crs")
terra::ext(dem_full) <- attr(scorpion, "extent")
sub_ext <- terra::ext(-2141500, -2140500, 1963500, 1964500)
sub_dem <- terra::crop(dem_full, sub_ext)

swp_multi <- geomorphon_sweep(
  sub_dem,
  param_grid = data.frame(search = c(5, 8), flat_angle_deg = c(1, 1)),
  forms = c("forms10", "forms6")
)
ens10 <- geomorphon_ensemble(swp_multi, form = "forms10", method = "proportions")
ens6 <- geomorphon_ensemble(swp_multi, form = "forms6", method = "proportions")

expect_equal(terra::nlyr(ens10), 10L)
expect_equal(terra::nlyr(ens6), 6L)
div_ens <- geomorphon_divergence(ens10, ens6)
expect_inherits(div_ens$raster, "SpatRaster")
expect_true(div_ens$mean_divergence >= 0.0 && div_ens$mean_divergence <= 1.0)
expect_true(div_ens$n_cells > 0L)


# -----------------------------------------------------------------------------
# 4. geomorphon_calibrate() and geomorphon_pareto()
# -----------------------------------------------------------------------------

# Compact 2-run sweep on cropped scorpion with target layers
sub_soil <- terra::crop(soil, sub_ext)
targets_sub <- geomorphon_targets(sub_soil, vars = c("musym", "drclassdcd"))
names(targets_sub) <- c("musym", "drclassdcd")

swp <- geomorphon_sweep(
  sub_dem,
  param_grid = data.frame(search = c(5, 8), flat_angle_deg = c(1, 1)),
  forms = c("forms10", "forms6")
)

cal <- geomorphon_calibrate(
  swp,
  targets_sub,
  spatial_null = FALSE
)

expect_true(is.data.frame(cal))

# A run holding only a forms6 layer must not be scored as forms10
fake_run <- geomorphons(sub_dem, search = 5, skip = 0, flat_angle_deg = 1, forms = "forms6")
names(fake_run) <- "forms6"  # one layer that names a taxonomy
fake_sweep <- structure(
  list(summary = swp$summary[1, , drop = FALSE], rasters = list(run_1 = fake_run),
       provenance = swp$provenance),
  class = "geomorphon_sweep"
)
cal_fake <- geomorphon_calibrate(fake_sweep, targets_sub, forms = c("forms10", "forms6"),
                                 spatial_null = FALSE)
expect_true(all(cal_fake$forms == "forms6"))
expect_equal(nrow(cal_fake), 2L)  # 1 run * 1 layer * 2 targets

# 2 runs * 2 forms * 2 targets = 8 rows
expect_equal(nrow(cal), 8L)

required_cols <- c(
  "run_id", "search", "flat_angle_deg", "forms", "target",
  "nmi", "homogeneity", "completeness",
  "n_edges_geo", "n_edges_cat", "edge_ratio", "density_penalty",
  "edge_f1", "edge_mhd", "is_pareto", "pareto_score"
)
expect_true(all(required_cols %in% colnames(cal)))
expect_true(is.logical(cal$is_pareto))
expect_true(is.numeric(cal$pareto_score))
expect_true(all(cal$edge_ratio > 0, na.rm = TRUE))
expect_true(all(cal$density_penalty > 0 & cal$density_penalty <= 1.0, na.rm = TRUE))

# Synthetic Pareto domination testing:
# Objective defaults: c(delta_nmi = "max", delta_edge_f1 = "max", edge_mhd = "min")
syn_cal <- data.frame(
  run_id = 1:4,
  target = "soil_type",
  forms = "forms10",
  delta_nmi = c(0.40, 0.50, 0.30, 0.45),
  delta_edge_f1 = c(0.30, 0.40, 0.20, 0.35),
  edge_mhd = c(15.0, 10.0, 25.0, 12.0),
  edge_ratio = c(1.5, 2.0, 1.2, 8.0)
)
# Note:
# Row 2 (0.50, 0.40, 10.0) strictly dominates Row 1 (0.40, 0.30, 15.0) and Row 3 (0.30, 0.20, 25.0).
# Row 4 (0.45, 0.35, 12.0) is dominated by Row 2.
# So only Row 2 is non-dominated.

par_res <- geomorphon_pareto(syn_cal)
expect_equal(par_res$is_pareto, c(FALSE, TRUE, FALSE, FALSE))
expect_equal(par_res$pareto_score[2], 1.0)
expect_true(par_res$pareto_score[3] < par_res$pareto_score[1])

# max_edge_ratio filtering:
# If max_edge_ratio = 1.8, Row 2 (edge_ratio = 2.0) and Row 4 (8.0) are filtered out!
# Candidates become Row 1 and Row 3. Row 1 dominates Row 3.
par_filtered <- geomorphon_pareto(syn_cal, max_edge_ratio = 1.8)
expect_equal(par_filtered$is_pareto, c(TRUE, FALSE, FALSE, FALSE))

# NA edge ratios (no boundaries, e.g. NMI-only): with no bounds set they stay eligible
syn_na <- syn_cal
syn_na$edge_ratio <- NA_real_
par_na <- geomorphon_pareto(syn_na)
expect_equal(par_na$is_pareto, c(FALSE, TRUE, FALSE, FALSE))
expect_true(all(!is.na(par_na$pareto_score)))
# With bounds requested, NA edge ratios are not eligible
par_na_bounded <- geomorphon_pareto(syn_na, max_edge_ratio = 1.8)
expect_equal(par_na_bounded$is_pareto, c(FALSE, FALSE, FALSE, FALSE))

}

