# Compute principal component directions

Internal helper that returns the principal component directions actually
used in downstream routines. By default, it returns the usual
non-private principal component directions. When `g_dppca = TRUE`, it
returns differentially private principal component directions based on
the spherical Kendall mechanism.

## Usage

``` r
compute_pc_dir(
  X_proc,
  k,
  g_dppca = FALSE,
  eps_dir = NULL,
  delta_dir = NULL,
  cpp.option = FALSE
)
```

## Arguments

- X_proc:

  Preprocessed numeric matrix.

- k:

  Number of leading components.

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
