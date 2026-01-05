// Focal edge detection using RcppParallel for high performance
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

#include <RcppArmadillo.h>
#include <RcppParallel.h>

// Worker class for parallel edge detection
struct EdgeDetectionWorker : public RcppParallel::Worker {
  // Input matrix (read-only, thread-safe access)
  const RcppParallel::RMatrix<double> input_mat;

  // Output matrix (write-only, no race conditions since each thread writes different rows)
  RcppParallel::RMatrix<double> output_mat;

  // NA value for comparison
  double na_value;

  // Constructor
  EdgeDetectionWorker(
      const Rcpp::NumericMatrix& input,
      Rcpp::NumericMatrix& output,
      double na_val
  ) : input_mat(input), output_mat(output), na_value(na_val) { }

  // Operator: process rows from begin to end in parallel
  void operator()(std::size_t begin_row, std::size_t end_row) {
    int nrows = input_mat.nrow();
    int ncols = input_mat.ncol();

    // Process each row
    for (std::size_t r = begin_row; r < end_row; ++r) {
      // Handle border rows (set to NA)
      if (r == 0 || r == (std::size_t)(nrows - 1)) {
        for (int c = 0; c < ncols; ++c) {
          output_mat(r, c) = na_value;
        }
        continue;
      }

      // Process interior cells
      for (int c = 0; c < ncols; ++c) {
        // Handle border columns (set to NA)
        if (c == 0 || c == ncols - 1) {
          output_mat(r, c) = na_value;
          continue;
        }

        double center = input_mat(r, c);

        // Skip if center is NA
        if (R_IsNA(center)) {
          output_mat(r, c) = na_value;
          continue;
        }

        // Check 8-connected neighbors for difference from center
        bool is_edge = false;

        for (int dr = -1; dr <= 1 && !is_edge; ++dr) {
          for (int dc = -1; dc <= 1 && !is_edge; ++dc) {
            // Skip center cell
            if (dr == 0 && dc == 0) continue;

            double neighbor = input_mat(r + dr, c + dc);

            // If neighbor is not NA and different from center, it's an edge
            if (!R_IsNA(neighbor) && neighbor != center) {
              is_edge = true;
            }
          }
        }

        // Store result: 1 for edge, 0 for non-edge
        output_mat(r, c) = is_edge ? 1.0 : 0.0;
      }
    }
  }
};

//' Detect Edges in a Raster (C++ Parallel Implementation)
//'
//' Fast edge detection for a raster using 8-connected neighborhood.
//' A cell is considered an edge if it differs from any of its 8 neighbors.
//' Border cells and NA values are set to NA in the output.
//'
//' @param x Numeric matrix representing raster values
//' @param na_value Value to use for NA cells (default: NA_REAL)
//'
//' @return Numeric matrix with same dimensions as input.
//'   Values: 1 = edge, 0 = non-edge, NA = border or NA input cell
//'
//' @noRd
//'
//' @examples
//' \dontrun{
//' # Create sample raster
//' m <- matrix(c(1,1,2,1,1,2,2,2,3), nrow=3)
//' edges <- focal_edge_cpp(m)
//' }
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix focal_edge_cpp(
    Rcpp::NumericMatrix x,
    double na_value = NA_REAL
) {
  int nrows = x.nrow();
  int ncols = x.ncol();

  // Create output matrix initialized to na_value
  Rcpp::NumericMatrix result(nrows, ncols);
  std::fill(result.begin(), result.end(), na_value);

  // Create worker
  EdgeDetectionWorker worker(x, result, na_value);

  // Process in parallel
  RcppParallel::parallelFor(0, nrows, worker);

  return result;
}
