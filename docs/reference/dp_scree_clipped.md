# Differentially private scree estimation via clipped mean

Internal helper that estimates the leading scree values using a
clipped-mean differentially private estimator applied to the centered
squared scores of each principal component.

The procedure is as follows:

1.  preprocess the data matrix for PCA,

2.  compute a non-private or private direction matrix,

3.  form the projected score matrix,

4.  for each component \\\ell\\, compute centered squared scores
    \\w\_{i\ell} = (y\_{i\ell} - \bar y\_\ell)^2\\,

5.  clip each \\w\_{i\ell}\\ at `C_clip`,

6.  release the clipped mean with Gaussian noise,

7.  apply the \\n/(n-1)\\ correction to obtain the final scree estimate.

If `mono = TRUE`, the final scree vector is post-processed so that it is
nonnegative and nonincreasing.

## Usage

``` r
dp_scree_clipped(
  X,
  k,
  eps_total,
  delta_total,
  center = TRUE,
  standardize = TRUE,
  C_clip,
  g_dppca = FALSE,
  cpp.option = FALSE,
  mono = TRUE
)
```

## Arguments

- X:

  Numeric data matrix with observations in rows.

- k:

  Integer number of leading principal components.

- eps_total:

  Total privacy epsilon allocated to the scree routine.

- delta_total:

  Total privacy delta allocated to the scree routine.

- center:

  Logical; whether to center columns before PCA.

- standardize:

  Logical; whether to standardize columns before PCA.

- C_clip:

  Positive clipping threshold applied to centered squared scores.

- g_dppca:

  Logical; whether to privatize the PCA direction matrix.

- cpp.option:

  Logical passed to
  [`mech_tau_sph()`](https://yejinjo0220.github.io/dppca/reference/mech_tau_sph.md)
  when private directions are computed.

- mono:

  Logical; whether to post-process the final scree vector so that it is
  nonnegative and nonincreasing.

## Value

A list with components:

- `scree`:

  Private scree estimates.

- `scree_np`:

  Non-private clipped scree estimates.

- `evr`:

  Explained variance ratios computed from the private scree estimates.
