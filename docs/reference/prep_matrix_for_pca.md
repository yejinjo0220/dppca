# Preprocess matrix for PCA

Internal helper that preprocesses a numeric matrix before PCA by
applying optional centering and standardization.

## Usage

``` r
prep_matrix_for_pca(X, center = TRUE, standardize = FALSE)
```

## Arguments

- X:

  Numeric matrix-like object.

- center:

  Logical; whether to center columns.

- standardize:

  Logical; whether to scale columns by their standard deviations.

## Value

A numeric matrix after preprocessing.
