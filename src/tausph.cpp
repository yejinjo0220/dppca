#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Calculate multivariate Kendall's tau matrix (Han and Liu, 2017).
// 
// INPUT
// X: n by p data matrix. (Each column is an observed p-dim data)
// m: number of PC directions to select
// OUTPUT
// hK: p by p multivariate Kendall's tau matrix


// [[Rcpp::export]]
arma::mat tau_sph_cpp(arma::mat X) {
  
  int p = X.n_cols;
  int n = X.n_rows;
  mat hK(p, p, fill::zeros);
  rowvec a;
  
  for (int i=0; i < n-1; i++) {
    for (int j=i+1; j<n; j++) {
      a = X.row(j) - X.row(i);
      a = a / norm(a);
      hK = hK + a.t() * a;
    }
  }
  
  hK = hK * (2.0 / (n * (n-1.0))); // we need to convert n to double type
  
  return hK;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
//
// /*** R
// X <- as.matrix(USArrests)
// multiKendallcpp(X)
// */
