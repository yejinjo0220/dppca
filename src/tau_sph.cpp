#include <Rcpp.h>
#include <cfloat>
#include <cmath>
#include <vector>

// [[Rcpp::export]]
Rcpp::NumericMatrix tau_sph_cpp(
    const Rcpp::NumericMatrix& X
) {
  const int n = X.nrow();
  const int d = X.ncol();

  Rcpp::NumericMatrix out(d, d);
  std::vector<double> a(d);

  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double norm_sq = 0.0;

      for (int l = 0; l < d; ++l) {
        a[l] = X(j, l) - X(i, l);
        norm_sq += a[l] * a[l];
      }

      if (!std::isfinite(norm_sq) || norm_sq <= DBL_EPSILON) {
        Rcpp::stop(
          "Cannot normalize a vector with zero or non-finite norm. "
          "The spherical Kendall mechanism is undefined for "
          "duplicate observations."
        );
      }

      for (int l = 0; l < d; ++l) {
        a[l] /= std::sqrt(norm_sq);
      }

      for (int l = 0; l < d; ++l) {
        for (int m = 0; m <= l; ++m) {
          out(l, m) += a[l] * a[m];
        }
      }
    }
  }

  const double scale =
    2.0 / (static_cast<double>(n) * (n - 1));

  for (int l = 0; l < d; ++l) {
    for (int m = 0; m <= l; ++m) {
      out(l, m) *= scale;
      out(m, l) = out(l, m);
    }
  }

  return out;
}
