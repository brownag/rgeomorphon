// Focal PTRM - Moving window PTRM calculation using RcppParallel
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>

// Circular mask for window extraction.
// Canonical rule: include cell (dr, dc) when sqrt(dr^2 + dc^2) <= radius.
// With radius = 0 this yields only the centre cell (0, 0).
struct CircularMask {
    std::vector<std::pair<int, int>> offsets;  // Pre-computed (row_offset, col_offset) pairs
    int radius;

    CircularMask(int r) : radius(r) {
        for (int dr = -r; dr <= r; ++dr) {
            for (int dc = -r; dc <= r; ++dc) {
                double dist = std::sqrt(dr * dr + dc * dc);
                if (dist <= r) {
                    offsets.push_back({dr, dc});
                }
            }
        }
    }
};

// Worker class for parallel focal PTRM calculation
struct PTRMFocalWorker : public RcppParallel::Worker {
    // Input: geomorphon classification matrix (class codes 1..k)
    const RcppParallel::RMatrix<int> geo_mat;

    // Output: PTRM scores
    RcppParallel::RMatrix<double> output_mat;

    // Pre-computed circular mask
    const CircularMask& mask;

    // PTRM coefficients (length k, indexed 0..k-1 for classes 1..k)
    std::vector<double> coeffs;

    // PTRM formula constants
    double intercept;   // -38.02
    double normalizer;  // 69.96
    double na_value;

    // Constructor
    PTRMFocalWorker(
        const Rcpp::IntegerMatrix& geo,
        Rcpp::NumericMatrix& output,
        const CircularMask& m,
        const Rcpp::NumericVector& c,
        double interc,
        double norm,
        double na_val
    ) : geo_mat(geo), output_mat(output), mask(m),
        intercept(interc), normalizer(norm), na_value(na_val) {
        // Convert Rcpp vector to std::vector
        coeffs.assign(c.begin(), c.end());
    }

    // Process rows from begin_row to end_row in parallel
    void operator()(std::size_t begin_row, std::size_t end_row) {
        int nrows = geo_mat.nrow();
        int ncols = geo_mat.ncol();
        int k = static_cast<int>(coeffs.size());

        for (std::size_t r = begin_row; r < end_row; ++r) {
            for (int c = 0; c < ncols; ++c) {

                // Handle border cells - set to NA
                if (r < (std::size_t)mask.radius ||
                    r >= (std::size_t)(nrows - mask.radius) ||
                    c < mask.radius ||
                    c >= ncols - mask.radius) {
                    output_mat(r, c) = na_value;
                    continue;
                }

                // Count geomorphon classes in circular window
                // class codes 1..k are valid (indices 0..k-1)
                std::vector<int> counts(k, 0);
                int total_valid = 0;

                // Iterate through offset pairs in circular mask
                for (const auto& offset : mask.offsets) {
                    int nr = r + offset.first;
                    int nc = c + offset.second;

                    // Bounds check (should always pass due to border handling above)
                    if (nr >= 0 && nr < nrows && nc >= 0 && nc < ncols) {
                        int geo_class = geo_mat(nr, nc);

                        // Check if value is not NA (NA_INTEGER in R)
                        // R uses INT_MIN as NA
                        if (geo_class != NA_INTEGER && geo_class >= 1 && geo_class <= k) {
                            counts[geo_class - 1]++;
                            total_valid++;
                        }
                    }
                }

                // If no valid cells in window, set output to NA
                if (total_valid == 0) {
                    output_mat(r, c) = na_value;
                    continue;
                }

                // Calculate PTRM score
                double weighted_sum = 0.0;
                for (int i = 0; i < k; ++i) {
                    double proportion = counts[i] / (double)total_valid;
                    weighted_sum += proportion * coeffs[i];
                }

                // Apply PTRM formula: (intercept + weighted_sum) / normalizer
                double raw_score = intercept + weighted_sum;
                double ptrm_score = raw_score / normalizer;

                output_mat(r, c) = ptrm_score;
            }
        }
    }
};

// Focal PTRM kernel.
//
// Window: canonical circular mask, i.e. include (dr, dc) when
//   sqrt(dr^2 + dc^2) <= radius.
// Classes: generic k = coefficients.size() classes. Valid input codes are
// 1..k; NA_INTEGER and codes outside 1..k are ignored. The R wrapper is
// responsible for remapping reduced-system codes to 1..k before calling
// this kernel.
//
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix geomorphon_ptrm_focal_cpp(
    Rcpp::IntegerMatrix geomorphon_class,
    int radius,
    Rcpp::NumericVector coefficients,  // Length k, indexed 0..k-1 for classes 1..k
    double intercept,                   // -38.02
    double normalizer,                  // 69.96
    double na_value                     // NA_REAL
) {
    // Validate inputs
    if (radius < 0) {
        Rcpp::stop("radius must be non-negative");
    }

    if (coefficients.size() < 1) {
        Rcpp::stop("coefficients must have length >= 1");
    }

    // Create circular mask
    CircularMask mask(radius);

    // Create output matrix
    int nrows = geomorphon_class.nrow();
    int ncols = geomorphon_class.ncol();
    Rcpp::NumericMatrix result(nrows, ncols);

    // Initialize with NA
    std::fill(result.begin(), result.end(), na_value);

    // Create worker
    PTRMFocalWorker worker(
        geomorphon_class,
        result,
        mask,
        coefficients,
        intercept,
        normalizer,
        na_value
    );

    // Run in parallel
    RcppParallel::parallelFor(0, nrows, worker);

    return result;
}
