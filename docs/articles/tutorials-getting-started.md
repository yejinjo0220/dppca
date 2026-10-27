# Getting Started with dppca

This guide covers the essentials for getting started with `dppca`:
installing and loading the package, understanding its example datasets,
preparing your own data, and allocating privacy budgets with
[`privacy_control()`](https://yejinjo0220.github.io/dppca/reference/privacy_control.md).
If you are looking for a particular visualization, use the guide below
to choose a tutorial and its corresponding theory article.

## Choose a starting point

The **tutorials** explain how to use `dppca` in R, while the **theory
articles** describe the statistical concepts and methods. For an
introduction to ordinary PCA, see [PCA
Background](https://yejinjo0220.github.io/dppca/articles/theory-pca-background.md).

| What you want to do | Tutorial | Theory and Methods |
|:---|:---|:---|
| Explore private PC directions and their coefficients | [DP Loading Plot](https://yejinjo0220.github.io/dppca/articles/tutorials-loading-plots.md) | [Theory of DP Loading](https://yejinjo0220.github.io/dppca/articles/theory-dp-loading.md) |
| Examine scree values and proportions of variance explained (PVE) | [DP Scree Plot](https://yejinjo0220.github.io/dppca/articles/tutorials-scree-plots.md) | [Theory of DP Scree](https://yejinjo0220.github.io/dppca/articles/theory-dp-scree.md) |
| Visualize private score distributions and compare groups | [DP Score Plot](https://yejinjo0220.github.io/dppca/articles/tutorials-score-plots.md) | [Theory of DP Score](https://yejinjo0220.github.io/dppca/articles/theory-dp-score.md) |
| Explore loading, scree, score, variable, and biplot views in one workflow | [DP PCA Visualization](https://yejinjo0220.github.io/dppca/articles/tutorials-dp-visualization.md) | [Theory of DP PCA Visualizations](https://yejinjo0220.github.io/dppca/articles/theory-dp-visualization.md) |

For complete argument descriptions and return values, see
[Functions](https://yejinjo0220.github.io/dppca/reference/index.md).

## 1. Install and load the package

Install the released version of `dppca` from CRAN:

``` r
install.packages("dppca")
```

Alternatively, install the development version from GitHub:

``` r
install.packages("remotes")  # Only needed once
remotes::install_github("yejinjo0220/dppca")
```

Once installed, load the package:

``` r
library(dppca)
```

## 2. Prepare your input data

In `dppca`, rows represent observations and columns represent variables.
Numeric features are used for PCA. For functions supporting `group`
(such as
[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
and
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)),
a data frame may also contain categorical columns; the selected grouping
column is excluded from PCA features. Other non-numeric columns are not
included in the PCA calculation.

The private PCA method implemented in `dppca` is based on pairwise
differences of observations and therefore does not require explicit
centering. Also, if standardization is needed, it should be performed
before applying `dppca`. When standardization involves estimating
parameters from sensitive data, the preprocessing must satisfy
differential privacy, and its privacy budget must be included in the
total privacy budget.

## 3. Explore the Example Datasets

The `dppca` package provides several example datasets for demonstrating
its functions. These include numeric datasets and datasets containing
categorical variables for grouped visualizations.

| Dataset | Rows | Columns | What it contains |
|:---|---:|---:|:---|
| `gau` | 1,000 | 20 | Numeric variables `V1`–`V20` from five synthetic Gaussian clusters, with 200 observations per cluster. |
| `gau_g` | 1,000 | 21 | The same 20 numeric variables plus the character label `group`, with values `group1`–`group5`. |
| `vowel` | 990 | 14 | Ten acoustic measurements and four categorical columns describing British English vowel observations. |

The Gaussian clusters differ in means and covariance scales. `gau_g`
adds grouping labels to the observations in `gau`; it is not a separate
sample.

### 3.1. Gaussian data

``` r
data("gau", package = "dppca")
data("gau_g", package = "dppca")

dim(gau)
#> [1] 1000   20
dim(gau_g)
#> [1] 1000   21
head(gau[, 1:4, drop = FALSE])
#>           V1         V2          V3         V4
#> 1  0.7019527 -0.2136231 -0.36868434 -0.8423263
#> 2  0.3361815  1.1978761  0.97822807  0.1018881
#> 3  0.7498257  0.2318031 -0.30707361 -0.8979258
#> 4 -0.8008823 -0.5028415 -0.05840928  1.3939254
#> 5 -0.1227414  0.6304571  0.35253375 -2.4865239
#> 6  0.6642886  0.9572975 -0.18232763  0.4012941
table(gau_g$group)
#> 
#> group1 group2 group3 group4 group5 
#>    200    200    200    200    200
```

The variables `V1`–`V20` are numeric features. The `group` column in
`gau_g` identifies five groups and can be used as a categorical label in
analyses that support grouping.

For more details about the datasets, see the documentation for
[`gau`](https://yejinjo0220.github.io/dppca/reference/gau.md) and
[`gau_g`](https://yejinjo0220.github.io/dppca/reference/gau_g.md).

### 3.2. Vowel data

The `vowel` dataset contains 990 observations of spoken British English
vowels from 15 speakers, covering 11 vowel categories. Each observation
consists of ten numeric acoustic measurements (`F0`–`F9`) and four
categorical variables (`Split`, `Speaker`, `Sex`, and `Class`).

``` r
data("vowel", package = "dppca")

numeric_columns <- vapply(vowel, is.numeric, logical(1))
dim(vowel)
#> [1] 990  14
names(vowel)[numeric_columns]
#>  [1] "F0" "F1" "F2" "F3" "F4" "F5" "F6" "F7" "F8" "F9"
names(vowel)[!numeric_columns]
#> [1] "Split"   "Speaker" "Sex"     "Class"
table(vowel$Class)
#> 
#>  0  1  2  3  4  5  6  7  8  9 10 
#> 90 90 90 90 90 90 90 90 90 90 90
```

The variables in the dataset are summarized below.

| Columns | Role |
|:---|:---|
| `F0`–`F9` | Ten numeric acoustic features used for PCA. |
| `Split` | Training/test-set label. |
| `Speaker` | Speaker identifier, with 15 levels. |
| `Sex` | Speaker descriptor. |
| `Class` | Vowel category, with 11 levels. This variable can be used for grouped score visualizations. |

For more details about the dataset, see the documentation for
[`vowel`](https://yejinjo0220.github.io/dppca/reference/vowel.md).

## 4. Set and Allocate the Privacy Budget

The
[`privacy_control()`](https://yejinjo0220.github.io/dppca/reference/privacy_control.md)
function specifies the privacy budget for PCA estimation in `dppca`. The
privacy parameters, `eps` and `delta`, are allocated to three
components:

| Component | Purpose                                               |
|:----------|:------------------------------------------------------|
| `loading` | Estimation of private principal component directions. |
| `scree`   | Estimation of private scree values and PVE.           |
| `score`   | Construction of private score visualizations.         |

The `dppca` package provides three functions for individual analyses:
[`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md),
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md),
and
[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md).
The
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
function combines these analyses into a single workflow. Since these
functions involve different estimation steps, the required privacy
budget components depend on the function used.

### 4.1. Default privacy budget allocation

The total privacy budget can be specified using `eps` and `delta`:

``` r
privacy <- privacy_control(
  eps = 6,
  delta = 1e-4
)
```

When `split` is omitted, the fitting function automatically divides the
total budget equally among its required components.

| Fitting function | Required components | Default allocation |
|:---|:---|:---|
| [`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md) | `loading` | 100% |
| [`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md) | `loading`, `scree` | 50% each |
| [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md) | `loading`, `score` | 50% each |
| [`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md) | `loading`, `scree`, `score` | 1/3 each |

For example, when the specified total budget is `eps = 6` and
`delta = 1e-4`:

- [`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
  allocates the entire budget to direction estimation.
- [`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
  and
  [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
  allocate epsilon `3` and delta `5e-5` to each of their two components.
- [`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
  allocates epsilon `2` and delta `1e-4 / 3` to each of its three
  components.

The total budget and its allocations can be inspected using:

``` r
privacy$total
#>   eps delta 
#> 6e+00 1e-04
privacy$budgets
#> $loading
#>          eps        delta 
#> 2.000000e+00 3.333333e-05 
#> 
#> $scree
#>          eps        delta 
#> 2.000000e+00 3.333333e-05 
#> 
#> $score
#>          eps        delta 
#> 2.000000e+00 3.333333e-05
```

The `total` field contains the overall budget, and `budgets` contains
its component allocations.

### 4.2. Customize the budget allocation

To change the default allocation, use the `split` argument to specify
the proportion of the total budget assigned to each component. The
proportions must be positive, named, and sum to one. The same
proportions are applied to both epsilon and delta.

For example, to allocate more budget to scree and score estimation in
[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md):

``` r
pca_privacy_custom <- privacy_control(
  eps = 6,
  delta = 1e-4,
  split = c(loading = 0.2, scree = 0.4, score = 0.4)
)

pca_privacy_custom$budgets
#> $loading
#>     eps   delta 
#> 1.20000 0.00002 
#> 
#> $scree
#>     eps   delta 
#> 2.40000 0.00004 
#> 
#> $score
#>     eps   delta 
#> 2.40000 0.00004
```

Similarly, for
[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md),
a different allocation can be specified as follows:

``` r
score_privacy_custom <- privacy_control(
  eps = 6,
  delta = 1e-4,
  split = c(loading = 0.4, score = 0.6)
)

score_privacy_custom$budgets
#> $loading
#>     eps   delta 
#> 2.40000 0.00004 
#> 
#> $score
#>     eps   delta 
#> 3.60000 0.00006
```

### 4.3. Specify individual component budgets

Alternatively, epsilon and delta can be specified separately for each
component, without using a total budget or `split`.

``` r
pca_privacy_absolute <- privacy_control(
  loading = c(eps = 2, delta = 2e-5),
  scree   = c(eps = 3, delta = 3e-5),
  score   = c(eps = 5, delta = 5e-5)
)

pca_privacy_absolute$total
#>   eps delta 
#> 1e+01 1e-04
```

In this example, the total budget is epsilon `10` and delta `1e-4`.
Unlike the proportional allocation approach, this allows epsilon and
delta to have different allocation proportions.

The total-budget and individual-component approaches should not be
combined in a single call.

## 5. Next steps

With the package installed, the input data prepared, and the privacy
budget specified, you are ready to explore differentially private PCA
visualizations. Start with the [DP Loading Plot
Tutorial](https://yejinjo0220.github.io/dppca/articles/tutorials-loading-plots.md)
to estimate and visualize private PC directions, or choose another
tutorial from the overview.

------------------------------------------------------------------------

**Next:** [DP Loading
Plot](https://yejinjo0220.github.io/dppca/articles/tutorials-loading-plots.md)  
**Tutorial overview:** [Choose a starting
point](#choose-a-starting-point)
