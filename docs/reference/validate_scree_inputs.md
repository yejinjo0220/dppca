# Validate common inputs for scree-related routines

Internal helper used at the beginning of scree-related routines to
validate the most common user inputs. In particular, this function
checks that:

- `X` can be coerced to a numeric matrix,

- the sample size and ambient dimension are valid,

- `k` is a valid number of leading components,

- `eps_total` is a positive scalar, and

- `delta_total` is a scalar in `(0, 1)`.

This helper does not return a transformed object for downstream use; its
role is simply to fail early with informative error messages before
private scree estimation begins.

## Usage

``` r
validate_scree_inputs(X, k, eps_total, delta_total)
```

## Arguments

- X:

  A matrix-like object containing the data, with observations in rows
  and variables in columns.

- k:

  Integer number of leading principal components for which scree values
  will be estimated.

- eps_total:

  Total privacy budget \\\varepsilon\\ allocated to the scree-estimation
  routine.

- delta_total:

  Total privacy budget \\\delta\\ allocated to the scree-estimation
  routine.

## Value

Invisibly returns `TRUE` if all checks pass.
