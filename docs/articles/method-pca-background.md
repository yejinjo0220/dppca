# Principal Component Analysis

## 1. Introduction

Principal component analysis (PCA) is a method that represents
multivariate observations in a lower-dimensional space while preserving
as much variation as possible. It constructs mutually orthogonal axes,
called principal component (PC) directions, and projects observations
onto these axes.

Suppose we observe $`n`$ vectors $`x_1,\ldots,x_n\in\mathbb{R}^p`$,
collected as the rows of an $`n\times p`$ data matrix

``` math

X=
\begin{bmatrix}
x_1^\top\\
\vdots\\
x_n^\top
\end{bmatrix}.
```

Each row represents an observation, and each column represents a
variable.

## 2. Population PCA

Let $`\xi\in\mathbb{R}^p`$ be a random vector with mean $`\mu=E[\xi]`$
and covariance matrix $`\Sigma
=
\operatorname{Cov}(\xi),`$ assuming finite second moments.

For a unit vector $`v\in\mathbb{R}^p`$, the variance of the projection
$`v^\top\xi`$ is $`\operatorname{Var}(v^\top\xi)=v^\top\Sigma v`$.

The first PC direction maximizes this projected variance:

``` math

v_1\in\arg\max_{\|v\|_2=1}v^\top\Sigma v.
```
Each subsequent direction maximizes the same quantity subject to being
orthogonal to the preceding directions.

Since $`\Sigma`$ is symmetric and positive semidefinite, it has an
eigen-decomposition

``` math

\Sigma=V\Lambda V^\top,
\quad
V=[v_1,\ldots,v_p],
\quad
\Lambda=\operatorname{diag}(\lambda_1,\ldots,\lambda_p),
```

where the $`v_\ell`$ are orthonormal eigenvectors and
$`\lambda_1\ge\cdots\ge\lambda_p\ge0`$. In particular,

``` math

\Sigma v_\ell=\lambda_\ell v_\ell,
\quad
\operatorname{Var}(v_\ell^\top\xi)=\lambda_\ell.
```

Thus, **eigenvectors determine the PC directions**, while **eigenvalues
quantify the variance along those directions**. The first few principal
components account for the largest amounts of variation along mutually
orthogonal directions.

## 3. Sample PCA

In practice, the population covariance matrix $`\Sigma`$ is generally
unknown, so PCA is performed using the sample covariance matrix

``` math

\widehat \Sigma
=
\frac{1}{n-1}
\sum_{i=1}^n
(x_i-\bar x)(x_i-\bar x)^\top,
\qquad
\bar x=\frac1n\sum_{i=1}^n x_i.
```

Equivalently, letting $`X_c=X-\mathbf{1}_n\bar x^\top`$ denote the
centered data matrix, where $`\mathbf{1}_n`$ is the $`n`$-dimensional
vector of ones, we can write

``` math

\widehat \Sigma =\frac{1}{n-1}X_c^\top X_c.
```

### 3.1 Directions

The **PC directions** are obtained from the eigen-decomposition of the
sample covariance matrix:
``` math

\widehat \Sigma =\widehat V\widehat\Lambda\widehat V^\top,
\quad
\widehat V=[\widehat v_1,\ldots,\widehat v_p],
\quad
\widehat\Lambda
=\operatorname{diag}(\widehat\lambda_1,\ldots,\widehat\lambda_p),
```

where $`\widehat\lambda_1\ge\cdots\ge\widehat\lambda_p\ge0`$.

The columns of $`\widehat V`$ are orthonormal eigenvectors satisfying

``` math

\widehat \Sigma \widehat v_\ell
=
\widehat\lambda_\ell\widehat v_\ell,
\qquad \ell=1,\ldots,p.
```

Each direction $`\widehat v_\ell`$ defines a PC axis, with the first
direction corresponding to the largest sample variance.

### 3.2 Eigenvalues

The **eigenvalues** of $`\widehat \Sigma`$ represent the sample
variances along the corresponding PC directions. In particular,
``` math

\widehat\lambda_\ell
=
\widehat v_\ell^\top \widehat \Sigma\widehat v_\ell,
\qquad \ell=1,\ldots,p.
```

The **proportion of variance explained (PVE)** by the $`\ell`$th
principal component is defined as

``` math

\operatorname{PVE}_\ell
=
\frac{\widehat\lambda_\ell}
{\sum_{j=1}^p\widehat\lambda_j},
```

where the denominator equals $`\operatorname{tr}(\widehat \Sigma)`$, the
total sample variance, assumed to be positive.

### 3.3 Scores

The **PC scores** are obtained by projecting the centered observations
onto the PC directions. The score matrix is defined as
``` math

S=X_c\widehat V,
```

where the $`\ell`$th PC score of observation $`i`$ is $`s_{i\ell}
=
\widehat v_\ell^\top(x_i-\bar x)`$.

For dimensionality reduction, we may retain only the first $`k`$
principal components. Letting
$`\widehat V_k=[\widehat v_1,\ldots,\widehat v_k]`$, the corresponding
reduced score matrix is $`S_k=X_c\widehat V_k`$.

## 4. PCA visualizations

PCA provides several visualization methods for exploring different
aspects of multivariate data. The five commonly used displays are
loading plots, scree plots, score plots, variable plots, and biplots.

For two-dimensional visualizations, let $`a`$ and $`b`$ denote two
selected principal component indices, with $`1\le a<b\le p`$. Although
the first two components are commonly used, any two components can be
selected depending on the purpose of the analysis.

### 4.1 Loading plot

A **loading plot** displays the coefficients of the original variables
in the PC directions. Specifically, let $`\widehat v_{j\ell}`$ denote
the coefficient of variable $`j`$ in PC direction $`\ell`$. For two
selected components $`a`$ and $`b`$, the loading plot represents
variable $`j`$ by

``` math

L_j=(\widehat v_{ja},\widehat v_{jb}).
```

The magnitude of each coefficient indicates the contribution of the
corresponding variable to a PC direction. Larger absolute coefficients
indicate greater contributions, while the signs distinguish
contributions in opposite directions.

### 4.2 Scree plot

A **scree plot** displays the eigenvalues $`\widehat\lambda_\ell`$
against the component indices $`\ell=1,\ldots,p`$. Since each eigenvalue
represents the sample variance along its corresponding PC direction, the
plot illustrates how the total variation is distributed across the
principal components.

Alternatively, the proportion of variance explained (PVE) can be plotted
against the component index. These displays help identify the components
that account for most of the variation and guide the choice of how many
components to retain.

### 4.3 Score plot

A **score plot** displays the projections of observations onto selected
PC directions. For two selected components $`a`$ and $`b`$, the score
plot represents observation $`i`$ by

``` math

(s_{ia},s_{ib}),
```

where $`s_{i\ell}=\widehat v_\ell^\top(x_i-\bar x)`$ is the score of
observation $`i`$ on the $`\ell`$th principal component.

Score plots are useful for exploring patterns, group structures, and
potential outliers in the data.

### 4.4 Variable plot (correlation circle)

A **variable plot** displays the correlations between the original
variables and the PC scores. Let
$`\widehat \sigma_j=\sqrt{\widehat \Sigma_{jj}}`$ denote the sample
standard deviation of variable $`j`$. For $`\widehat \sigma>0`$ and
$`\widehat\lambda_\ell>0`$, the sample correlation between variable
$`j`$ and the $`\ell`$th PC scores is
``` math

\operatorname{Cor}(X_{\cdot j},S_{\cdot\ell})
=
\frac{\sqrt{\widehat\lambda_\ell}\,\widehat v_{j\ell}}{\widehat \sigma_j},
```

where $`X_{\cdot j}`$ and $`S_{\cdot\ell}`$ denote the $`j`$th and
$`\ell`$th columns of $`X`$ and $`S`$, respectively.

For two selected components $`a`$ and $`b`$, the variable plot
represents variable $`j`$ by

``` math

R_j=
\left(
\operatorname{Cor}(X_{\cdot j},S_{\cdot a}),
\operatorname{Cor}(X_{\cdot j},S_{\cdot b})
\right).
```

Since $`\|R_j\|_2^2
=
\operatorname{Cor}(X_{\cdot j},S_{\cdot a})^2
+
\operatorname{Cor}(X_{\cdot j},S_{\cdot b})^2
\le 1,`$ the variable coordinates lie inside or on the unit circle.

Variables near the unit circle are well represented by the two selected
components, whereas variables near the origin have relatively little of
their variation represented in the displayed plane. Also, the angles
between variable vectors can provide information about associations
between variables, particularly when both variables are well represented
by the selected components.

### 4.5 Biplot

A **biplot** displays observations and variables simultaneously in a
low-dimensional representation of the data. To describe its
construction, consider the reduced singular value decomposition (SVD) of
the centered data matrix:

``` math

X_c=UDV^\top,
```

where $`r=\operatorname{rank}(X_c)`$, $`U\in\mathbb R^{n\times r}`$ and
$`V\in\mathbb R^{p\times r}`$ have orthonormal columns, and

``` math

D=\operatorname{diag}(d_1,\ldots,d_r),
\qquad d_1\ge\cdots\ge d_r>0.
```

For a two-dimensional biplot, let $`I=\{a,b\}`$ denote the selected
component indices, where $`1\le a<b\le r`$. Define

``` math

U_I=[u_a,u_b],
\qquad
D_I=\operatorname{diag}(d_a,d_b),
\qquad
V_I=[v_a,v_b].
```

The corresponding rank-two approximation of the centered data matrix is

``` math

X_c\approx U_ID_IV_I^\top.
```

A general family of biplot scalings is defined by

``` math

G_\alpha=U_ID_I^{1-\alpha},
\qquad
H_\alpha=V_ID_I^\alpha,
\qquad 0\le\alpha\le1.
```

Recall that the PC score matrix for the selected components is
$`S_I=X_cV_I=U_ID_I`$. Therefore,
$`G_\alpha=U_ID_I^{1-\alpha}=S_ID_I^{-\alpha}`$ showing that the
observation coordinates in a biplot are obtained by scaling the ordinary
PC scores according to $`\alpha`$.

These matrices satisfy $`G_\alpha H_\alpha^\top\approx X_c`$. The rows
of $`G_\alpha`$ represent observations, while the rows of $`H_\alpha`$
represent variables. The parameter $`\alpha`$ determines how the
singular-value scaling is distributed between them:

- $`\alpha=0`$: Observation coordinates are the ordinary PC scores on
  the selected components, and variable coordinates are the
  corresponding PC direction coefficients.
- $`\alpha=1/2`$: Singular-value scaling is shared symmetrically between
  observations and variables.
- $`\alpha=1`$: Variable coordinates are scaled by the singular values,
  emphasizing the variance and covariance structure of the original
  variables.

**Note.** Since $`\widehat\Sigma=\frac{1}{n-1}VD^2V^\top`$, the columns
of $`V`$ are eigenvectors of $`\widehat\Sigma`$, with corresponding
positive eigenvalues $`\widehat\lambda_\ell=\frac{d_\ell^2}{n-1}`$,
$`\ell=1,\ldots,r`$. Thus, the right singular vectors of $`X_c`$ provide
the sample PC directions, while its singular values determine the
corresponding eigenvalues. Hence, PCA can equivalently be performed
using the SVD of $`X_c`$.

### Summary of PCA visualizations

The five visualization methods provide complementary information about
the data:

| Display | Main purpose |
|----|----|
| Loading plot | Examine the contributions of variables to PC directions |
| Scree plot | Assess the variance explained by each principal component |
| Score plot | Explore the distribution and patterns of observations in PC space |
| Variable plot | Examine correlations between original variables and PC scores |
| Biplot | Explore observations and variables in a common display |

## References

Abdi, H. and Williams, L. J. (2010). Principal component analysis.
*WIREs Computational Statistics*, 2(4), 433–459.
<https://doi.org/10.1002/wics.101>

Gabriel, K. R. (1971). The biplot graphic display of matrices with
application to principal component analysis. *Biometrika*, 58(3),
453–467. <https://doi.org/10.1093/biomet/58.3.453>

Greenacre, M., Groenen, P. J. F., Hastie, T., Iodice D’Enza, A., Markos,
A., and Tuzhilina, E. (2022). Principal component analysis. *Nature
Reviews Methods Primers*, 2, 100.
<https://doi.org/10.1038/s43586-022-00184-w>

Jolliffe, I. T. (2002). *Principal Component Analysis* (2nd ed.).
Springer, New York. <https://doi.org/10.1007/b98835>

Jolliffe, I. T. and Cadima, J. (2016). Principal component analysis: a
review and recent developments. *Philosophical Transactions of the Royal
Society A*, 374(2065), 20150202.
<https://doi.org/10.1098/rsta.2015.0202>
