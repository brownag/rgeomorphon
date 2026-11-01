#include <Rcpp.h>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

using namespace Rcpp;

namespace {
double nmi_score(const std::vector<double>& g, const std::vector<double>& c,
                 const std::vector<unsigned char>& valid, int n) {
  std::map<double, int> fg, fc;
  std::map<std::pair<double, double>, int> fj;
  int count = 0;
  for (int i = 0; i < n; ++i) if (valid[i]) {
    ++fg[g[i]]; ++fc[c[i]]; ++fj[std::make_pair(g[i], c[i])]; ++count;
  }
  if (!count) return NA_REAL;
  const auto entropy = [count](const auto& f) {
    double h = 0;
    for (const auto& x : f) { double p = double(x.second) / count; h -= p * std::log2(p); }
    return h;
  };
  double hg = entropy(fg), hc = entropy(fc);
  if (hg + hc == 0) return 0;
  double hj = entropy(fj), mi = std::max(0.0, hg + hc - hj);
  return std::min(1.0, std::max(0.0, 2 * mi / (hg + hc)));
}

std::vector<unsigned char> edges(const std::vector<double>& x,
                                 const std::vector<unsigned char>& valid,
                                 int nr, int nc) {
  std::vector<unsigned char> out(nr * nc, 0);
  for (int r = 1; r + 1 < nr; ++r) for (int col = 1; col + 1 < nc; ++col) {
    int i = r * nc + col;
    if (!valid[i]) continue;
    for (int dr = -1; dr <= 1 && !out[i]; ++dr) for (int dc = -1; dc <= 1; ++dc) {
      if (!dr && !dc) continue;
      int j = (r + dr) * nc + col + dc;
      if (valid[j] && x[i] != x[j]) { out[i] = 1; break; }
    }
  }
  return out;
}

double f1_score(const std::vector<unsigned char>& a,
                const std::vector<unsigned char>& b,
                int nr, int nc, double threshold) {
  std::vector<int> ea, eb;
  for (int i = 0; i < nr * nc; ++i) { if (a[i]) ea.push_back(i); if (b[i]) eb.push_back(i); }
  if (ea.empty() || eb.empty()) return NA_REAL;
  const auto hit = [nr, nc, threshold](const std::vector<int>& from,
                                        const std::vector<int>& to) {
    int hits = 0;
    const int max_radius = std::max(nr / 2, nc / 2);
    const int radius = threshold >= static_cast<double>(max_radius)
      ? max_radius : static_cast<int>(std::floor(threshold));
    const double threshold2 = threshold * threshold;
    for (int i : from) {
      const int r = i / nc, c = i % nc;
      bool found = false;
      for (int dr = -radius; dr <= radius && !found; ++dr) {
        const int row_offset = std::abs(dr) % nr;
        const int rr = std::min(row_offset, nr - row_offset);
        for (int dc = -radius; dc <= radius; ++dc) {
          const int col_offset = std::abs(dc) % nc;
          const int cc = std::min(col_offset, nc - col_offset);
          if (double(rr * rr + cc * cc) > threshold2) continue;
          const int row = (r + dr % nr + nr) % nr;
          const int col = (c + dc % nc + nc) % nc;
          if (to[row * nc + col]) { found = true; break; }
        }
      }
      if (found) ++hits;
    }
    return double(hits) / from.size();
  };
  double precision = hit(eb, ea), recall = hit(ea, eb);
  return precision + recall == 0 ? 0 : 2 * precision * recall / (precision + recall);
}
}

// [[Rcpp::export]]
List spatial_null_cpp(NumericVector geo, NumericVector cat, int nr, int nc,
                      double threshold, bool do_nmi, bool do_edge) {
  int n = nr * nc;
  std::vector<double> g(geo.begin(), geo.end()), c(cat.begin(), cat.end());
  std::vector<unsigned char> gv(n), cv(n);
  for (int i = 0; i < n; ++i) { gv[i] = !NumericVector::is_na(g[i]); cv[i] = !NumericVector::is_na(c[i]); }
  NumericVector nmi_vals, edge_vals;
  int evaluated = 0, nmi_eligible = 0, edge_eligible = 0;
  double nmi_sum = 0, edge_sum = 0;
  const int shifts = std::min(128, n - 1);
  for (int i = 0; i < shifts; ++i) {
    const int rank = static_cast<int>(std::floor((i + 0.5) * (n - 1.0) / shifts));
    const int q = rank + 1;
    const int dr = q / nc, dc = q % nc;
    ++evaluated;
    std::vector<double> shifted(n);
    std::vector<unsigned char> valid(n);
    for (int r = 0; r < nr; ++r) for (int col = 0; col < nc; ++col) {
      int dst = r * nc + col;
      int sr = (r - dr + nr) % nr, sc = (col - dc + nc) % nc, src = sr * nc + sc;
      shifted[dst] = c[src]; valid[dst] = gv[dst] && cv[src];
    }
    if (do_nmi) { double v = nmi_score(g, shifted, valid, n); if (!NumericVector::is_na(v)) { nmi_sum += v; ++nmi_eligible; } }
    if (do_edge) {
      auto ge = edges(g, valid, nr, nc), ce = edges(shifted, valid, nr, nc);
      double v = f1_score(ge, ce, nr, nc, threshold);
      if (!NumericVector::is_na(v)) { edge_sum += v; ++edge_eligible; }
    }
  }
  return List::create(_["evaluated"] = evaluated,
    _["nmi_mean"] = do_nmi && nmi_eligible ? nmi_sum / nmi_eligible : NA_REAL,
    _["nmi_eligible"] = nmi_eligible,
    _["edge_mean"] = do_edge && edge_eligible ? edge_sum / edge_eligible : NA_REAL,
    _["edge_eligible"] = edge_eligible);
}
