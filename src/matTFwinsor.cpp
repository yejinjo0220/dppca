#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Winsorize each row of given matrix X (Han and Liu, 2017).
// 
// INPUT
// X: n by p data matrix. (Each column is an observed p-dim data)
// R: positive scalar. winsorization radius.
// OUTPUT
// hK: p by p multivariate Kendall's tau matrix

// [[Rcpp::export]]
arma::mat matTFwinsor(arma::mat X, double R) {
  
  int n = X.n_rows;

  for (int i=0; i<n; i++){
    X.row(i) = X.row(i) * (std::min(R/norm(X.row(i)), 1.0));
  }
  
  return X;
}


