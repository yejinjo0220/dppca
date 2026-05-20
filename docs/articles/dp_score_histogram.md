# Differentially Private PCA Score Histograms

## Overview

A PCA score plot is a standard visualization tool for understanding the
low-dimensional structure of high-dimensional data. In ordinary PCA, one
often plots the first two score coordinates directly as a scatter plot.

However, directly releasing all score points can reveal individual-level
information. The goal of a differentially private PCA score
visualization is to display the empirical structure of the projected
data while protecting the privacy of individual observations.

The `dppca` package represents the private PCA score plot through a
two-dimensional differentially private histogram. Instead of releasing
all score points, the package releases a privatized approximation of the
empirical score distribution.

This article describes the main steps:

1.  compute PCA scores,
2.  choose two score coordinates,
3.  construct a private plotting frame,
4.  divide the frame into bins,
5.  privatize the bin counts, and
6.  optionally repeat the procedure by group.

## PCA score matrix

Let
``` math

X_c \in \mathbb{R}^{n \times p}
```
be the preprocessed data matrix. Typically, this means that the data
have been centered, and possibly standardized.

Let
``` math

V_k = [v_1,\ldots,v_k] \in \mathbb{R}^{p \times k}
```
be the principal component direction matrix, where each column
$`v_\ell`$ is a principal component direction.

The $`k`$-dimensional PCA score for the $`i`$-th observation is
``` math

z_i
=
V_k^\top x_{c,i}
\in
\mathbb{R}^k,
```
where $`x_{c,i}`$ denotes the $`i`$-th row of $`X_c`$, viewed as a
column vector.

Equivalently, the score matrix is
``` math

Z = X_c V_k
\in
\mathbb{R}^{n \times k}.
```

The $`i`$-th row of $`Z`$ contains the PCA score coordinates of
observation $`i`$.

## Two-dimensional score coordinates

In visualization, we usually select two components rather than plotting
all $`k`$ score dimensions.

Let
``` math

(a,b)
```
be the two selected component indices. A common choice is
``` math

(a,b) = (1,2).
```

For each observation, define the two-dimensional score point
``` math

s_i
=
(z_{i,a}, z_{i,b})
\in
\mathbb{R}^2,
\qquad
i=1,\ldots,n.
```

The collection
``` math

S = \{s_i\}_{i=1}^n
```
is the two-dimensional score point cloud.

A non-private score plot would draw all points
``` math

s_1,\ldots,s_n
```
directly. In the private setting, we instead approximate the empirical
distribution of $`S`$ through a private histogram.

## Goal of private score visualization

The goal is to produce a visualization that approximates the empirical
distribution of
``` math

S = \{s_i\}_{i=1}^n
\subset
\mathbb{R}^2
```
without directly releasing the individual score points.

In `dppca`, the private score visualization is based on the following
idea:

1.  Construct a private plotting frame $`F \subset \mathbb{R}^2`$.
2.  Divide $`F`$ into rectangular bins.
3.  Count how many score points fall into each bin.
4.  Add privacy noise to the bin counts.
5.  Normalize and visualize the noisy bin frequencies.

Thus, the released object is a privatized score histogram rather than
the raw score point cloud.

## Private plotting frame

Before building a two-dimensional histogram, we need to determine the
plotting region. This region is called the plotting frame.

If the frame is too narrow, many points may be excluded. If it is too
wide, the histogram may become too sparse. Therefore, the frame should
cover the relevant range of score coordinates while remaining stable
under privacy constraints.

The `dppca` package constructs a square plotting frame using private
quantile estimation.

## Stacking score coordinates

Let
``` math

S \in \mathbb{R}^{n \times 2}
```
be the matrix of two-dimensional scores, where the $`i`$-th row is
``` math

s_i^\top = (z_{i,a}, z_{i,b}).
```

To construct one common range for both axes, we stack the two score
coordinates into a single vector:
``` math

z
=
(z_1,\ldots,z_m)^\top
\in
\mathbb{R}^m,
\qquad
m = 2n.
```

Here, the vector $`z`$ contains all values from both selected score
coordinates:
``` math

\{z_{i,a}\}_{i=1}^n
\cup
\{z_{i,b}\}_{i=1}^n.
```

This converts the problem of choosing a two-dimensional plotting range
into a one-dimensional private quantile estimation problem.

## Private lower and upper quantiles

To determine the lower and upper boundaries of the plotting frame,
define
``` math

q_{\min} = \frac{1}{m},
\qquad
q_{\max} = \frac{m-1}{m}.
```

These levels correspond approximately to the minimum and maximum of the
stacked score coordinates.

Let
``` math

\widetilde z_{\min}
\quad\text{and}\quad
\widetilde z_{\max}
```
be private estimates of the $`q_{\min}`$- and $`q_{\max}`$-quantiles of
the stacked vector $`z`$.

In `dppca`, these private quantiles may be computed using a smooth
sensitivity based DP quantile estimator.

## Square plotting frame

Given the private lower and upper quantile estimates, define the center
``` math

c
=
\frac{\widetilde z_{\min}+\widetilde z_{\max}}{2}
```
and the half-length
``` math

L
=
\frac{\widetilde z_{\max}-\widetilde z_{\min}}{2}.
```

To add a visual margin and reduce boundary effects, introduce an inflate
parameter
``` math

\alpha > 0.
```

The inflated half-length is
``` math

L'
=
(1+\alpha)L.
```

The final square plotting frame is
``` math

F
=
[c-L',c+L']
\times
[c-L',c+L'].
```

Using the same range for the horizontal and vertical axes simplifies
interpretation and makes group-wise comparisons easier.

## Functions for the plotting frame

The plotting frame construction can be summarized by two conceptual
functions.

| Function | Role |
|----|----|
| [`dp_quantile_ss()`](https://yejinjo0220.github.io/dppca/reference/dp_quantile_ss.md) | Computes a private quantile using smooth sensitivity |
| [`dp_frame()`](https://yejinjo0220.github.io/dppca/reference/dp_frame.md) | Constructs a square plotting frame from private lower and upper quantiles |

The exact function names may depend on the package implementation.
Conceptually, the frame step takes score coordinates and privacy
parameters as input, and returns a private square plotting region.

## Choosing the number of bins

Once the plotting frame $`F`$ has been determined, it must be divided
into histogram bins.

Let
``` math

m_{\mathrm{axis}}
```
be the number of bins per axis. Then the total number of two-dimensional
bins is
``` math

m = m_{\mathrm{axis}}^2.
```

The number of bins can be chosen by the user or selected using a
heuristic rule. For a $`d`$-dimensional histogram, common rules have the
form
``` math

m_{\mathrm{axis}}
\asymp
n^{1/(2+d)}
```
or
``` math

m_{\mathrm{axis}}
\asymp
\left(
\frac{n}{\log n}
\right)^{1/(1+d)}.
```

Since a score histogram is two-dimensional, we typically take
``` math

d = 2.
```

These rules should be interpreted as starting points rather than strict
optimal choices. In practice, the best bin resolution depends on the
sample size, privacy budget, and the amount of structure in the data.

## Two-dimensional histogram

Let the private plotting frame be divided into bins
``` math

B_1,\ldots,B_m.
```

For the score point set
``` math

S = \{s_i\}_{i=1}^n,
```
the non-private count in bin $`B_k`$ is
``` math

c_k
=
\sum_{i=1}^n
\mathbf{1}\{s_i \in B_k\},
\qquad
k=1,\ldots,m.
```

The count vector is
``` math

c = (c_1,\ldots,c_m)
\in
\mathbb{N}^m.
```

If the histogram is normalized to a frequency vector, then
``` math

q_k
=
\frac{c_k}{\sum_{j=1}^m c_j}
=
\frac{c_k}{n},
\qquad
k=1,\ldots,m,
```
assuming all points are assigned to bins in the frame.

The final private score visualization displays a noisy version of the
frequency vector
``` math

q = (q_1,\ldots,q_m).
```

## Sensitivity of histogram counts

Under row-level adjacency, two neighboring datasets differ in one
observation. Changing one observation can move one score point from one
bin to another.

Therefore, the count vector can change by at most $`+1`$ in one bin and
$`-1`$ in another bin. Hence,
``` math

\Delta_1(c)
\leq
2
```
and
``` math

\Delta_2(c)
\leq
\sqrt{2}.
```

These sensitivity bounds are used to calibrate privacy noise for the
histogram mechanism.

## Additive DP histogram

A simple DP histogram mechanism adds independent noise to each bin count
and then post-processes the resulting vector to obtain a valid
histogram-like output.

### Algorithm: Additive DP histogram

**Input:** score points $`S=\{s_i\}_{i=1}^n`$, bins $`\{B_k\}_{k=1}^m`$,
privacy parameters $`(\epsilon,\delta)`$.

**Output:** a private normalized histogram
``` math

\hat q = (\hat q_1,\ldots,\hat q_m).
```

``` text
1. For each bin k = 1,...,m:
   a. Compute the count
      c_k = sum_i 1{s_i in B_k}.
   b. Add Gaussian noise:
      \tilde c_k = c_k + eta_k.
   c. Truncate negative values:
      \hat c_k = max(\tilde c_k, 0).

2. Let Z = sum_j \hat c_j.

3. Normalize:
   \hat q_k = \hat c_k / Z.

4. Return \hat q.
```

Using Gaussian noise, the noisy count can be written as
``` math

\widetilde c_k
=
c_k + \eta_k,
\qquad
\eta_k \sim N(0,\sigma^2).
```

With $`\ell_2`$-sensitivity $`\sqrt{2}`$, one possible Gaussian noise
scale is
``` math

\sigma
=
\sqrt{2}
\cdot
\frac{\sqrt{2\log(1.25/\delta)}}{\epsilon}.
```

After adding noise, negative counts are truncated to zero and the vector
is renormalized. These steps are post-processing and therefore do not
consume additional privacy budget.

### Interpretation

The additive mechanism is straightforward and works well when many bins
have moderate counts. However, when the grid is fine, many bins are
empty, and adding noise to every bin can create artificial mass in empty
regions.

This motivates sparse DP histogram mechanisms.

## Sparse DP histogram

When many bins are empty, adding noise to every bin can dominate the
visualization. A sparse histogram aims to report only bins whose counts
are large enough to be distinguishable from noise.

The `dppca` package may use a count-based sparse histogram mechanism
based on thresholding.

### Algorithm: Count-based sparse DP histogram

**Input:** score points $`S=\{s_i\}_{i=1}^n`$, bins $`\{B_k\}_{k=1}^m`$,
privacy parameters $`(\epsilon,\delta)`$.

**Output:** a private normalized histogram
``` math

\hat q = (\hat q_1,\ldots,\hat q_m).
```

``` text
1. Set the threshold
   t = 2 log(2/delta) / epsilon + 1.

2. For each bin k = 1,...,m:
   a. Compute the count
      c_k = sum_i 1{s_i in B_k}.

   b. If c_k = 0, set \hat c_k = 0.

   c. If c_k > 0:
      i.   Add Laplace noise:
           \tilde c_k = c_k + eta_k,
           eta_k ~ Laplace(0, 2/epsilon).
      ii.  If \tilde c_k < t, set \tilde c_k = 0.
      iii. Set \hat c_k = max(\tilde c_k, 0).

3. Let Z = sum_j \hat c_j.

4. Normalize:
   \hat q_k = \hat c_k / Z.

5. Return \hat q.
```

The threshold is
``` math

t
=
\frac{2\log(2/\delta)}{\epsilon}
+
1.
```

The Laplace noise scale is
``` math

\frac{2}{\epsilon}.
```

After thresholding and truncation, the remaining noisy counts are
normalized.

### Interpretation

The sparse mechanism is especially useful when the histogram has many
bins but only a small number of occupied bins.

Compared with the additive mechanism, the sparse mechanism can reduce
visual noise in empty regions. However, small but real clusters may be
suppressed if their counts are below the threshold.

## Privacy accounting

The full private score histogram procedure has two main
privacy-consuming steps:

1.  private quantile estimation for constructing the plotting frame,
2.  private histogram release.

Let the frame construction budget be
``` math

(\epsilon_{\mathrm{frame}},\delta_{\mathrm{frame}})
```
and the histogram budget be
``` math

(\epsilon_{\mathrm{hist}},\delta_{\mathrm{hist}}).
```

By basic composition, the combined procedure satisfies
``` math

(\epsilon_{\mathrm{frame}}+\epsilon_{\mathrm{hist}},
\delta_{\mathrm{frame}}+\delta_{\mathrm{hist}})
\text{-DP}.
```

If the PCA directions are also estimated privately, then the
direction-estimation budget must be included as well. For example, if
the private direction step uses
``` math

(\epsilon_{\mathrm{dir}},\delta_{\mathrm{dir}}),
```
then the full procedure satisfies
``` math

(\epsilon_{\mathrm{dir}}
+
\epsilon_{\mathrm{frame}}
+
\epsilon_{\mathrm{hist}},
\delta_{\mathrm{dir}}
+
\delta_{\mathrm{frame}}
+
\delta_{\mathrm{hist}})
\text{-DP}.
```

## Group-wise DP score histograms

In many applications, one may want to compare PCA score distributions
across groups. For example, observations may have group labels
``` math

g_i \in \mathcal{G}.
```

The data are then represented as
``` math

\{(s_i,g_i)\}_{i=1}^n,
```
where $`s_i \in \mathbb{R}^2`$ is the two-dimensional PCA score and
$`g_i`$ is the group label.

A group-wise DP score histogram releases one private histogram for each
group:
``` math

\{\hat q^{(g)}\}_{g\in\mathcal{G}}.
```

## Group-wise additive DP histogram

For each group $`g\in\mathcal{G}`$, define the group-specific bin count
``` math

c_k^{(g)}
=
\sum_{i=1}^n
\mathbf{1}\{s_i \in B_k,\; g_i = g\}.
```

The additive group-wise mechanism adds Gaussian noise to each group-bin
count:
``` math

\widetilde c_k^{(g)}
=
c_k^{(g)}
+
\eta_k^{(g)}.
```

Then the noisy counts are truncated and normalized within each group:
``` math

\hat q_k^{(g)}
=
\frac{\hat c_k^{(g)}}{\sum_{j=1}^m \hat c_j^{(g)}}.
```

### Algorithm: Additive DP histogram with groups

``` text
1. Initialize c_k^(g) = 0 for every group g and bin k.

2. For each observation i:
   a. Find the bin k(i) such that s_i in B_{k(i)}.
   b. Add one count to c_{k(i)}^(g_i).

3. For each group g:
   a. For each bin k:
      i.   Add Gaussian noise:
           \tilde c_k^(g) = c_k^(g) + eta_k^(g).
      ii.  Set \hat c_k^(g) = max(\tilde c_k^(g), 0).
   b. Normalize:
      \hat q_k^(g) = \hat c_k^(g) / sum_j \hat c_j^(g).

4. Return {\hat q^(g)}_{g in G}.
```

## Group-wise sparse DP histogram

A group-wise sparse histogram applies thresholding separately to each
group and bin.

### Algorithm: Count-based sparse DP histogram with groups

``` text
1. Initialize c_k^(g) = 0 for every group g and bin k.

2. For each observation i:
   a. Find the bin k(i) such that s_i in B_{k(i)}.
   b. Add one count to c_{k(i)}^(g_i).

3. For each group g:
   a. For each bin k:
      i.   If c_k^(g) = 0, set \hat c_k^(g) = 0.
      ii.  If c_k^(g) > 0:
           - Add Laplace noise:
             \tilde c_k^(g) = c_k^(g) + eta_k^(g),
             eta_k^(g) ~ Laplace(0, 2/epsilon).
           - Set threshold
             t = 2 log(2/delta) / epsilon + 1.
           - If \tilde c_k^(g) < t, set \tilde c_k^(g) = 0.
           - Set \hat c_k^(g) = max(\tilde c_k^(g), 0).
   b. Normalize:
      \hat q_k^(g) = \hat c_k^(g) / sum_j \hat c_j^(g).

4. Return {\hat q^(g)}_{g in G}.
```

Group-wise histograms are useful for comparing low-dimensional structure
across categories while still releasing only privatized aggregate
information.

## Practical interpretation

The private score histogram should be interpreted as a noisy
approximation of the empirical score distribution.

- Regions with larger noisy mass indicate areas where many projected
  observations are concentrated.
- Empty or low-count regions may be affected by privacy noise or
  thresholding.
- Smaller privacy budgets produce stronger privacy but noisier
  histograms.
- Larger numbers of bins give finer resolution but increase sparsity.
- Sparse mechanisms are useful when most bins are empty.

Therefore, the choice of bin number, privacy budget split, and histogram
mechanism should be made jointly.

## Example usage

The exact function arguments may depend on the installed version of
`dppca`. A typical workflow is as follows.

``` r
library(dppca)

# x: numeric data matrix or data frame
# epsilon, delta: privacy parameters

out <- dp_score_plot(
  x,
  k = 2,
  epsilon = 1,
  delta = 1e-6,
  hist_method = "additive"
)

out
```

A sparse histogram version may look like:

``` r
out_sparse <- dp_score_plot(
  x,
  k = 2,
  epsilon = 1,
  delta = 1e-6,
  hist_method = "sparse"
)
```

For group-wise score visualization:

``` r
out_group <- group_dp_score_plot(
  x,
  group = group_label,
  k = 2,
  epsilon = 1,
  delta = 1e-6,
  hist_method = "sparse"
)
```

These examples are not evaluated in this vignette so that the article
can be built without requiring a specific dataset.

## Summary

The private PCA score histogram workflow can be summarized as follows.

1.  Compute or obtain principal component directions $`V_k`$.
2.  Project the data:
    ``` math

    Z = X_c V_k.
    ```
3.  Select two score coordinates $`(a,b)`$.
4.  Construct a private square plotting frame using private quantiles.
5.  Divide the frame into bins.
6.  Compute bin counts.
7.  Privatize the histogram using an additive or sparse mechanism.
8.  Normalize the noisy counts.
9.  Visualize the resulting private score distribution.

The key idea is that the package releases a private aggregate
representation of the score distribution rather than the raw score
points.

## References

Karwa, V. and Vadhan, S. P. (2017). Finite sample differentially private
confidence intervals. *CoRR*, abs/1711.03908.

Nissim, K., Raskhodnikova, S., and Smith, A. (2007). Smooth sensitivity
and sampling in private data analysis. In *Proceedings of the
Thirty-Ninth Annual ACM Symposium on Theory of Computing*, 75–84.

Wasserman, L. and Zhou, S. (2010). A statistical framework for
differential privacy. *Journal of the American Statistical Association*,
105, 375–389.

Kim, M. and Jung, S. (2025). *Robust and Differentially Private
Principal Component Analysis*.
