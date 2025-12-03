# helper functions for DP-PCA ----------------------------------------------------
normalize <- function(x){x/norm2(x)}

norm2 <- function(x){sqrt(sum(x^2))}
vec2mat <- function(v, p){
  # vector to matrix as in Schwartzman (2016)
  # p and q should be comparable
  q <- length(v)
  if(q != p*(p+1)/2){
    stop('dimension is not compatible')
  }

  tmp_mat <- matrix(0, p, p)
  diag(tmp_mat) <- v[1:p]
  tmp_mat[lower.tri(tmp_mat)] <- v[(p+1):q] / sqrt(2)
  tmp_mat[upper.tri(tmp_mat)] <- (t(tmp_mat))[upper.tri(tmp_mat)]

  return(tmp_mat)
}

tau_sph <- function(X, cpp.option = FALSE) {

  if (cpp.option) {
    stop("cpp.option = TRUE is not yet supported in the dppca package. ",
         "Please use cpp.option = FALSE.")
  }

  n <- nrow(X)
  d <- ncol(X)
  hK <- matrix(0, d, d)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      a <- X[j, ] - X[i, ]
      a <- normalize(a)
      hK <- hK + a %*% t(a)
    }
  }
  hK * (2 / (n * (n - 1)))
}



mech_tau_sph <- function(X, sig, cpp.option=F){
  # simple additive mechanism; just add a Gaussian error to sample Kendall's multivariate tau matrix
  # For (eps, delta)-DP: sig <- 4 * sqrt(2 * log(1.25 / delta)) / (n * eps)
  # For rho-zCDP: sig <- (4 / n) / sqrt(2 * rho)

  n <- nrow(X)
  d <- ncol(X)
  q <- d*(d+1)/2

  # generate a sample multivariate Kendall's tau matrix
  hK <- tau_sph(X, cpp.option)

  # generate a random noise for DP
  zeta <- stats::rnorm(q, 0, sig)
  zeta_mat <- vec2mat(zeta, d)

  # add a generated noise to hat_K
  hK_dp <- hK + zeta_mat
  return(hK_dp)
}











