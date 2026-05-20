# Adult numeric data

A numerical subset of the Adult dataset from the UCI Machine Learning
Repository. The original Adult dataset is based on data extracted from
the 1994 United States Census database.

## Usage

``` r
adult
```

## Format

A data frame with 32,561 rows and 5 columns:

- age:

  Age of the individual.

- education_num:

  Number of years of education.

- capital_gain:

  Capital gain.

- capital_loss:

  Capital loss.

- hours_per_week:

  Number of working hours per week.

## Source

UCI Machine Learning Repository, Adult dataset. The Adult dataset is
licensed under the Creative Commons Attribution 4.0 International (CC BY
4.0) license.

## Details

This package dataset retains five numerical variables: `age`,
`education_num`, `capital_gain`, `capital_loss`, and `hours_per_week`.
These selected numerical variables contain no missing values in the
original data file used here. The resulting dataset contains 32,561
observations.

Since the variables have substantially different units and scales,
scaling is recommended before applying PCA-based methods.

## References

Becker, B. and Kohavi, R. (1996). Adult dataset. UCI Machine Learning
Repository.
