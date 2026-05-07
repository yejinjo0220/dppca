# Compute principal component directions

Returns the principal component directions actually used in downstream
routines. By default, it returns the usual non-private principal
component directions. When `g_dppca = TRUE`, it returns differentially
private principal component directions based on the spherical Kendall
mechanism.

The input data `X` are preprocessed internally by
[`prep_matrix_for_pca()`](https://yejinjo0220.github.io/dppca/reference/prep_matrix_for_pca.md)
using the `center` and `standardize` options.

## Usage

``` r
dp_pc_dir(
  X,
  k,
  center = TRUE,
  standardize = FALSE,
  g_dppca = FALSE,
  eps_dir = NULL,
  delta_dir = NULL,
  cpp.option = FALSE
)
```

## Arguments

- X:

  Numeric matrix-like object.

- k:

  Number of leading components.

- center:

  Logical; whether to center columns before PCA.

- standardize:

  Logical; whether to scale columns by their standard deviations after
  optional centering.

- g_dppca:

  Logical; whether to privatize the principal component directions.

- eps_dir:

  Privacy epsilon for releasing principal component directions.

- delta_dir:

  Privacy delta for releasing principal component directions.

- cpp.option:

  Logical passed to
  [`mech_tau_sph()`](https://yejinjo0220.github.io/dppca/reference/mech_tau_sph.md).

## Value

A direction matrix whose columns are the principal component directions
used.
