# Vowel Recognition Dataset

Acoustic measurements of eleven British English vowels spoken by fifteen
speakers.

## Usage

``` r
vowel
```

## Format

A data frame with 990 observations and 14 variables:

- Split:

  Training or test set (Train, Test).

- Speaker:

  Speaker identifier (15 levels).

- Sex:

  Speaker sex (Male, Female).

- F0:

  Acoustic feature 0.

- F1:

  Acoustic feature 1.

- F2:

  Acoustic feature 2.

- F3:

  Acoustic feature 3.

- F4:

  Acoustic feature 4.

- F5:

  Acoustic feature 5.

- F6:

  Acoustic feature 6.

- F7:

  Acoustic feature 7.

- F8:

  Acoustic feature 8.

- F9:

  Acoustic feature 9.

- Class:

  Vowel category (11 levels).

## Source

UCI Machine Learning Repository. <https://doi.org/10.24432/C58P4S>

## Details

The dataset contains 990 observations and 14 variables, including ten
numeric acoustic features and four categorical variables.

This dataset was obtained from the UCI file `vowel-context.data`.
Categorical codes were converted to factors; numeric acoustic
measurements were preserved.

Original dataset license: CC BY 4.0.
<https://creativecommons.org/licenses/by/4.0/>

## References

Deterding, D., Niranjan, M., and Robinson, T. (1988). Connectionist
Bench (Vowel Recognition - Deterding Data). UCI Machine Learning
Repository.

Turney, P. (1993). Robust Classification With Context-Sensitive
Features.
