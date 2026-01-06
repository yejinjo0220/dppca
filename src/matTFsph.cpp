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
arma::mat matTFsph(arma::mat X) {
  
  int n = X.n_rows;

  for (int i=0; i<n; i++){
    X.row(i) = X.row(i) / norm(X.row(i));
  }
  
  return X;
}
