// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <vector>
#include <string>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

const double GEOM_PI = 3.14159265358979323846;

const double EPSILON = 1e-9;

enum FORMS_GRASS {
    G_NONE = 0,
    G_FL = 1,
    G_PK = 2,
    G_RI = 3,
    G_SH = 4,
    G_SP = 5,
    G_SL = 6,
    G_HL = 7,
    G_FS = 8,
    G_VL = 9,
    G_PT = 10
};

FORMS_GRASS form_from_counts(int num_neg, int num_pos) {
    if (num_neg < 0 || num_neg > 8 || num_pos < 0 || num_pos > 8) {
        return G_FL; // TODO: G_NONE?
    }

    const FORMS_GRASS forms_table[9][9] = {
        /* 0 */ {G_FL, G_FL, G_FL, G_FS, G_FS, G_VL, G_VL, G_VL, G_PT},
        /* 1 */ {G_FL, G_FL, G_FS, G_FS, G_FS, G_VL, G_VL, G_VL, G_NONE},
        /* 2 */ {G_FL, G_SH, G_SL, G_SL, G_HL, G_HL, G_VL, G_NONE, G_NONE},
        /* 3 */ {G_SH, G_SH, G_SL, G_SL, G_SL, G_HL, G_NONE, G_NONE, G_NONE},
        /* 4 */ {G_SH, G_SH, G_SP, G_SL, G_SL, G_NONE, G_NONE, G_NONE, G_NONE},
        /* 5 */ {G_RI, G_RI, G_SP, G_SP, G_NONE, G_NONE, G_NONE, G_NONE, G_NONE},
        /* 6 */ {G_RI, G_RI, G_RI, G_NONE, G_NONE, G_NONE, G_NONE, G_NONE, G_NONE},
        /* 7 */ {G_RI, G_RI, G_NONE, G_NONE, G_NONE, G_NONE, G_NONE, G_NONE, G_NONE},
        /* 8 */ {G_PK, G_NONE, G_NONE, G_NONE, G_NONE, G_NONE, G_NONE, G_NONE, G_NONE}
    };

    FORMS_GRASS result = forms_table[num_neg][num_pos];

    if (result == G_NONE)
        return G_FL;

    return result;
}

enum COMPARISON_MODE {
    ANGLEV1,
    ANGLEV2,
    ANGLEV2_DISTANCE
};

struct GeomorphonWorker : public RcppParallel::Worker {
    const RcppParallel::RMatrix<double> dem_access_mat_w;
    RcppParallel::RMatrix<double> forms_output_mat_w;

    int n_rows_w, n_cols_w;
    double x_res_w, y_res_w;
    int search_cells_w;
    double search_dist_w;
    int skip_cells_w;
    double flat_angle_rad_w;
    double planim_flat_dist_w;
    double flat_thresh_height_w;
    double t_dist_factor_w;
    COMPARISON_MODE comp_mode_w;
    double nodata_w;

    const int grass_nextr[8] = {-1, -1, -1,  0,  1,  1,  1,  0};
    const int grass_nextc[8] = { 1,  0, -1, -1, -1,  0,  1,  1};

    GeomorphonWorker(
        const Rcpp::NumericMatrix& dem_r_in,
        Rcpp::NumericMatrix& forms_r_out,
        double dem_x_res,
        double dem_y_res,
        int search_r_cells,
        double search_dist,
        int skip_r_cells,
        double los_flat_angle_rad,
        double planim_flat_dist,
        double flat_thresh_h,
        double t_dist_f,
        COMPARISON_MODE comp_mode,
        double nodata_v
    ) : dem_access_mat_w(dem_r_in),
        forms_output_mat_w(forms_r_out),
        n_rows_w(dem_r_in.nrow()),
        n_cols_w(dem_r_in.ncol()),
        x_res_w(dem_x_res),
        y_res_w(dem_y_res),
        search_cells_w(search_r_cells),
        search_dist_w(search_dist),
        skip_cells_w(skip_r_cells),
        flat_angle_rad_w(los_flat_angle_rad),
        planim_flat_dist_w(planim_flat_dist),
        flat_thresh_height_w(flat_thresh_h),
        t_dist_factor_w(t_dist_f),
        comp_mode_w(comp_mode),
        nodata_w(nodata_v)
    { }

    inline double get_dem(int r, int c) const {
        if (r < 0 || r >= n_rows_w || c < 0 || c >= n_cols_w) {
            return R_NaReal;
        }
        return dem_access_mat_w(r, c);
    }

    inline bool check_is_nodata(double val) const {
        if (std::isnan(nodata_w)) {
            return std::isnan(val);
        }
        if (!std::isnan(val) && std::abs(val - nodata_w) < EPSILON) {
            return true;
        }
        return false;
    }

    double planimetric_distance(int r_start, int c_start, int r_end, int c_end) const {
        double dr = static_cast<double>(r_end - r_start) * y_res_w;
        double dc = static_cast<double>(c_end - c_start) * x_res_w;
        return std::sqrt(dr * dr + dc * dc);
    }

    void operator()(std::size_t begin_row, std::size_t end_row) {
        std::vector<int> tp_grass_order(8);

        for (std::size_t r_idx = begin_row; r_idx < end_row; ++r_idx) {
            int r_center_cell = static_cast<int>(r_idx);
            for (int c_center_cell = 0; c_center_cell < n_cols_w; ++c_center_cell) {

                double center_height = get_dem(r_center_cell, c_center_cell);
                if (check_is_nodata(center_height)) {
                    forms_output_mat_w(r_center_cell, c_center_cell) = nodata_w;
                    continue;
                }

                if (r_center_cell <= skip_cells_w || r_center_cell >= n_rows_w - (skip_cells_w + 1) ||
                    c_center_cell <= skip_cells_w || c_center_cell >= n_cols_w - (skip_cells_w + 1) ) {
                    forms_output_mat_w(r_center_cell, c_center_cell) = nodata_w;
                    continue;
                }

                for (int i_dir_grass = 0; i_dir_grass < 8; ++i_dir_grass) {
                    double cur_zenith_angle = -M_PI_2;
                    double cur_nadir_angle = M_PI_2;
                    double cur_zenith_height_diff = 0.0;
                    double cur_nadir_height_diff = 0.0;
                    double cur_zenith_dist = 0.0;
                    double cur_nadir_dist = 0.0;

                    bool line_of_sight_hasid_points = false;
                    tp_grass_order[i_dir_grass] = 0;

                    int r_first_adj = r_center_cell + grass_nextr[i_dir_grass];
                    int c_first_adj = c_center_cell + grass_nextc[i_dir_grass];

                    if (r_first_adj < 0 || r_first_adj >= n_rows_w ||
                        c_first_adj < 0 || c_first_adj >= n_cols_w ||
                        check_is_nodata(get_dem(r_first_adj, c_first_adj))) {
                        continue;
                    }

                    for (int j_step = skip_cells_w + 1; j_step <= search_cells_w; ++j_step) {
                        int r_target_cell = r_center_cell + j_step * grass_nextr[i_dir_grass];
                        int c_target_cell = c_center_cell + j_step * grass_nextc[i_dir_grass];

                        double map_dist_to_target = planimetric_distance(r_center_cell, c_center_cell, r_target_cell, c_target_cell);

                        // while cur_distance < search_distance
                        if (map_dist_to_target >= search_dist_w - EPSILON) {
                            break;
                        }

                        if (r_target_cell < 0 || r_target_cell >= n_rows_w ||
                            c_target_cell < 0 || c_target_cell >= n_cols_w) {
                            break;
                        }

                        double target_height = get_dem(r_target_cell, c_target_cell);
                        if (check_is_nodata(target_height)) {
                            break;
                        }

                        line_of_sight_hasid_points = true;
                        double height_diff_from_center = target_height - center_height;
                        if (map_dist_to_target < EPSILON) {
                            continue;
                        }

                        double angle_to_target = std::atan2(height_diff_from_center, map_dist_to_target);

                        if (angle_to_target > cur_zenith_angle) {
                            cur_zenith_angle = angle_to_target;
                            cur_zenith_height_diff = height_diff_from_center;
                            cur_zenith_dist = map_dist_to_target;
                        }

                        if (angle_to_target < cur_nadir_angle) {
                            cur_nadir_angle = angle_to_target;
                            cur_nadir_height_diff = height_diff_from_center;
                            cur_nadir_dist = map_dist_to_target;
                        }
                    }

                    if (!line_of_sight_hasid_points) {
                        continue;
                    }

                    int tp_for_dir = 0;
                    double dom_h_diff_for_tdist = 0.0;
                    double dom_m_dist_for_tdist = search_dist_w;
                    double zenith_thresh_angle = flat_angle_rad_w;

                    if (planim_flat_dist_w > EPSILON && cur_zenith_dist > EPSILON &&
                        planim_flat_dist_w < cur_zenith_dist) {
                        zenith_thresh_angle = std::atan2(flat_thresh_height_w, cur_zenith_dist);
                    }

                    double nadir_thresh_angle = flat_angle_rad_w;
                    if (planim_flat_dist_w > EPSILON && cur_nadir_dist > EPSILON &&
                        planim_flat_dist_w < cur_nadir_dist) {
                        nadir_thresh_angle = std::atan2(flat_thresh_height_w, cur_nadir_dist);
                    }

                    bool zenith_signif = (cur_zenith_angle > zenith_thresh_angle + EPSILON);
                    bool nadir_signif = (cur_nadir_angle < -nadir_thresh_angle - EPSILON);

                    if (comp_mode_w == ANGLEV1) {
                        if (zenith_signif || nadir_signif) {
                            if (std::abs(cur_zenith_angle) > std::abs(cur_nadir_angle) + EPSILON) {
                                tp_for_dir = 1;
                                dom_h_diff_for_tdist = cur_zenith_height_diff;
                                dom_m_dist_for_tdist = cur_zenith_dist;
                            } else if (std::abs(cur_nadir_angle) > std::abs(cur_zenith_angle) + EPSILON) {
                                tp_for_dir = -1;
                                dom_h_diff_for_tdist = cur_nadir_height_diff;
                                dom_m_dist_for_tdist = cur_nadir_dist;
                            }
                        }
                    } else {
                        if (!zenith_signif && !nadir_signif) {
                            tp_for_dir = 0;
                        } else if (zenith_signif && !nadir_signif) {
                            tp_for_dir = 1;
                            dom_h_diff_for_tdist = cur_zenith_height_diff;
                            dom_m_dist_for_tdist = cur_zenith_dist;
                        } else if (!zenith_signif && nadir_signif) {
                            tp_for_dir = -1;
                            dom_h_diff_for_tdist = cur_nadir_height_diff;
                            dom_m_dist_for_tdist = cur_nadir_dist;
                        } else {
                            if (std::abs(cur_zenith_angle) > std::abs(cur_nadir_angle) + EPSILON) {
                                tp_for_dir = 1;
                                dom_h_diff_for_tdist = cur_zenith_height_diff;
                                dom_m_dist_for_tdist = cur_zenith_dist;
                            } else if (std::abs(cur_nadir_angle) > std::abs(cur_zenith_angle) + EPSILON) {
                                tp_for_dir = -1;
                                dom_h_diff_for_tdist = cur_nadir_height_diff;
                                dom_m_dist_for_tdist = cur_nadir_dist;
                            } else {
                                if (comp_mode_w == ANGLEV2_DISTANCE) {
                                    if (cur_zenith_dist > cur_nadir_dist + EPSILON) {
                                        tp_for_dir = 1;
                                        dom_h_diff_for_tdist = cur_zenith_height_diff;
                                        dom_m_dist_for_tdist = cur_zenith_dist;
                                    } else if (cur_nadir_dist > cur_zenith_dist + EPSILON) {
                                        tp_for_dir = -1;
                                        dom_h_diff_for_tdist = cur_nadir_height_diff;
                                        dom_m_dist_for_tdist = cur_nadir_dist;
                                    } else {
                                        tp_for_dir = 1;
                                        dom_h_diff_for_tdist = cur_zenith_height_diff;
                                        dom_m_dist_for_tdist = cur_zenith_dist;
                                    }
                                } else { // ANGLEV2
                                    tp_for_dir = 1;
                                    dom_h_diff_for_tdist = cur_zenith_height_diff;
                                    dom_m_dist_for_tdist = cur_zenith_dist;
                                }
                            }
                        }
                    }

                    if (t_dist_factor_w > EPSILON && dom_m_dist_for_tdist > EPSILON) {

                        double representative_res_w = std::min(x_res_w, y_res_w);
                        if (representative_res_w < EPSILON) {
                            representative_res_w = 1.0;
                        }

                        double dom_fvp_pixel_dist = dom_m_dist_for_tdist / representative_res_w;
                        if (dom_fvp_pixel_dist < 1.0) {
                            dom_fvp_pixel_dist = 1.0;
                        }

                        double tdist_vertical_thresh = t_dist_factor_w * dom_fvp_pixel_dist;

                        if (dom_h_diff_for_tdist > tdist_vertical_thresh + EPSILON) {
                            tp_for_dir = 1;
                        } else if (dom_h_diff_for_tdist < -tdist_vertical_thresh - EPSILON) {
                            tp_for_dir = -1;
                        } else {
                            tp_for_dir = 0;
                        }
                    }

                    tp_grass_order[i_dir_grass] = tp_for_dir;
                }

                int np = 0;
                int nn = 0;

                for (int k_tp = 0; k_tp < 8; ++k_tp) {
                    if (tp_grass_order[k_tp] == 1) {
                        np++;
                    } else if (tp_grass_order[k_tp] == -1) {
                        nn++;
                    }
                }

                FORMS_GRASS final_form_code = form_from_counts(nn, np);
                forms_output_mat_w(r_center_cell, c_center_cell) = static_cast<double>(final_form_code);
            }
        }
    }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix geomorphons_cpp_worker(
        Rcpp::NumericMatrix elevation,
        double search,
        double skip,
        double flat_angle_deg,
        double dist,
        std::string comparison_mode,
        double tdist,
        bool use_meters,
        double x_res_dem,
        double y_res_dem,
        double nodata
    ) {

    double ns_res = y_res_dem;
    double max_res = std::max(x_res_dem, y_res_dem);

    if (x_res_dem <= EPSILON || y_res_dem <= EPSILON) {
        Rcpp::stop("Cell resolutions (x_res, y_res) must be positive");
    }

    int search_cells;
    int skip_cells;
    double search_dist;
    double skip_dist;

    if (use_meters) {
        search_dist = search;
        skip_dist = skip;

        if (max_res <= EPSILON) {
            Rcpp::stop("Max resolution is near zero for map unit to cell conversion");
        }

        search_cells = static_cast<int>(search_dist / max_res);
        skip_cells = static_cast<int>(skip_dist / max_res);
    } else {
        search_cells = static_cast<int>(search);
        skip_cells = static_cast<int>(skip);

        search_dist = static_cast<double>(search_cells) * ns_res;
        skip_dist = static_cast<double>(skip_cells) * ns_res;
    }

    // sanity checks
    if (search_cells < 1) {
        Rcpp::stop("Search radius must be at least 1 cell");
    }

    if (skip_cells >= search_cells && search_cells > 0) { // Check search_cells > 0
        Rcpp::stop("Skip radius (cells) must be less than search radius (cells)");
    }

    if (skip_cells < 0) {
        Rcpp::stop("Skip radius (cells) must be >= 0");
    }

    if (flat_angle_deg <= EPSILON) {
        Rcpp::stop("Flatness threshold (degrees) must be > 0");
    }

    double flat_angle_rad = flat_angle_deg * GEOM_PI / 180.0;

    double planim_flat_dist = static_cast<double>(dist);
    if (!use_meters) {
        planim_flat_dist *= ns_res;
    }

    if (planim_flat_dist > EPSILON) {
        if ((planim_flat_dist <= skip_dist + EPSILON) ||
            (search_dist > EPSILON && planim_flat_dist >= search_dist - EPSILON)) {
            Rcpp::Rcout << "Flatness distance ('dist') is not between skip distance and search distance. Planimetric distance set to 0."
                        << std::endl;
            planim_flat_dist = 0.0;
        }
    }

    double flat_thresh_height = std::tan(flat_angle_rad) * planim_flat_dist;

    COMPARISON_MODE cmode = ANGLEV1;
    if (comparison_mode == "anglev2") {
        cmode = ANGLEV2;
    } else if (comparison_mode == "anglev2_distance") {
        cmode = ANGLEV2_DISTANCE;
    }

    Rcpp::NumericMatrix forms_out_mat(elevation.nrow(), elevation.ncol());
    std::fill(forms_out_mat.begin(), forms_out_mat.end(), nodata);

    GeomorphonWorker geomorphon_worker(elevation,
                                       forms_out_mat,
                                       x_res_dem,
                                       y_res_dem,
                                       search_cells,
                                       search_dist,
                                       skip_cells,
                                       flat_angle_rad,
                                       planim_flat_dist,
                                       flat_thresh_height,
                                       tdist,
                                       cmode,
                                       nodata);

    RcppParallel::parallelFor(0, elevation.nrow(), geomorphon_worker);
    return forms_out_mat;
}
