# Five Gaussian clusters

A simulated 20-dimensional Gaussian cluster dataset used as an example
for principal component analysis and differentially private PCA
visualization.

## Usage

``` r
gau
```

## Format

A data frame with 5,000 rows and 20 columns:

- V1:

  Simulated numerical variable.

- V2:

  Simulated numerical variable.

- V3:

  Simulated numerical variable.

- V4:

  Simulated numerical variable.

- V5:

  Simulated numerical variable.

- V6:

  Simulated numerical variable.

- V7:

  Simulated numerical variable.

- V8:

  Simulated numerical variable.

- V9:

  Simulated numerical variable.

- V10:

  Simulated numerical variable.

- V11:

  Simulated numerical variable.

- V12:

  Simulated numerical variable.

- V13:

  Simulated numerical variable.

- V14:

  Simulated numerical variable.

- V15:

  Simulated numerical variable.

- V16:

  Simulated numerical variable.

- V17:

  Simulated numerical variable.

- V18:

  Simulated numerical variable.

- V19:

  Simulated numerical variable.

- V20:

  Simulated numerical variable.

## Source

Simulated by the authors of the package. The data-generating code is
available in `data-raw/make_datasets.R` in the package source
repository.

## Details

The dataset contains 5,000 observations in 20 dimensions. It consists of
five groups, with 1,000 observations in each group. Each group is
generated from a multivariate normal distribution. The group mean
vectors and covariance matrices are chosen differently across groups so
that the data contain both separated and partially overlapping cluster
structures.
