# PC Directions in dppca

In ordinary PCA, the principal component directions are obtained from
the eigenvectors of the sample covariance matrix. In `dppca`, these
directions can be computed in two different ways.

1.  **Non-private PC directions**: eigenvectors of the sample covariance
    matrix.
2.  **Differentially private PC directions**: private principal
    component directions obtained through the g-DPPCA procedure.

## Notation

Let

``` math

X =
\begin{bmatrix}
X_1^\top \\
X_2^\top \\
\vdots \\
X_n^\top
\end{bmatrix}
\in \mathbb{R}^{n \times p}
```

be the data matrix used for PCA, where $`X_i \in \mathbb{R}^p`$is the
$`i`$-th observation. We assume that $`X`$ has been centered, and
optionally standardized.

The principal component direction matrix is denoted by

``` math

V_k = [v_1,\ldots,v_k] \in \mathbb{R}^{p \times k},
```

where each column $`v_\ell`$ is a unit vector representing the
$`\ell`$-th pc direction.

The corresponding score matrix is $`Z = X V_k`$.

## 1. Non-private PC directions

The classical sample covariance matrix is

``` math

\hat\Sigma
=
\frac{1}{n-1}X^\top X.
```

The non-private PCA directions are obtained from the eigenvalue
decomposition

``` math

\hat\Sigma
=
\hat V \hat\Lambda \hat V^\top,
```

where

``` math

\hat V = [\hat v_1,\ldots,\hat v_p],
\quad
\hat\Lambda
=
\operatorname{diag}(\hat\lambda_1,\ldots,\hat\lambda_p)
\quad \text{with} \quad
\hat\lambda_1 \geq \hat\lambda_2 \geq \cdots \geq \hat\lambda_p \geq 0.
```

The $`\ell`$-th sample principal component direction is $`\hat v_\ell`$.

Equivalently,

``` math

\hat v_\ell
=
\arg\max_{\|v\|_2 = 1}
v^\top \hat\Sigma v
\quad
\text{subject to}
\quad
v^\top \hat v_j = 0,
\qquad j = 1,\ldots,\ell-1.
```

In the non-private option of `dppca`, the direction matrix used for
projection is

``` math

\hat V_k = [\hat v_1,\ldots,\hat v_k].
```

## 2. DP PC directions

[Kim and Jung (2025)](#ref-Kim2025) proposed `g-DPPCA` by adding matrix
Gaussian mechanism on the generalized multivariate Kendall’s tau matrix
which based on the robust data transformation called generalized spatial
sign proposed by [Raymakers and Rousseeuw (2019)](#ref-Raymaekers2019).

For a positive valued scale function \$ : (0, ) (0, ) \$, consider a map
\$ g\_: ^d ^d \$ defined as

``` math

g_\xi(t) = \xi(\|t\|_2)\cdot \frac{t}{\|t\|_2}.
```

$`g_{\xi}`$ is called as a *generalized spatial sign* with respect to
$`\xi`$.

The *generalized multivariate Kendall’s tau* matrix with respect to \$
g\_\$ is defined as

``` math

K_{g_\xi} = \mathbb{E}_{X, \widetilde X}\left[ g_\xi\left( \frac{X - \widetilde X}{\sqrt{2}}\right) 
            g_\xi\left( \frac{X - \widetilde X}{\sqrt{2}}\right)^\top ~ \right],
```

where $`\widetilde X`$ is an independent copy of $`X`$. Importantly, if
$`X`$ follows an elliptical distribution (which including Gaussian and
multivariate $`t`$-distributions), $`K_{g_\xi}`$ shares the same
eigenvectors with same order to the $`\mbox{cov}(X)`$. So, one can
conduct a PCA by estimating $`K_{g_\xi}`$ and then get eigenvectors of
it.

For a convenience, we write $`g`$ as the given sign function. For a
random sample $`S = (X_1, \dots, X_n)`$, the second order \$ U
\$-statistic of \$ K_g \$ is could be written as

``` math

\widehat{K}_g(S) = \frac{2}{n(n-1)} \sum_{i < j} g\left(\frac{X_j - X_i}{\sqrt{2}}\right)
     g\left(\frac{X_j - X_i}{\sqrt{2}}\right)^\top.
```

Note that the sensitivity of $`\widehat K_g`$ with respec to the
Frobenius norm can be upper bounded by

``` math

\Delta_F(\widehat{K}_g) 
= \sup_{S \sim S'} \|\widehat{K}_g(S) - \widehat{K}_g(S')\|_F
\le \frac{4\|g\|_\infty^2}{n}.
```

So, for a dataset $`S = (x_1, \dots, x_n)`$ the randomized mechanism
$`\widetilde{K}_g`$ defined as

``` math

\widetilde K_g(S) :=   
\frac{2}{n(n-1)} \sum_{i < j} g\left(\frac{x_j-x_i}{\sqrt{2}}\right)g\left(\frac{x_j-x_i}{\sqrt{2}}\right)^\top + \mbox{vecd}^{-1}(\xi),
```
where \$N\_{d(d+1)/2}(0, *{, }^2 I*{d(d+1)/2}) \$ and \$ \_{, } = \$,
satisfies $`(\varepsilon, \delta)`$-DP.

Define \$ V\_{g, m}(S) (d, m)\$ as the matrix of the first $`m`$
eigenvectors of $`\widetilde K_g(S)`$. Then, \$ V\_{g, m}(S) \$
satisfies $`(\varepsilon, \delta)`$-DP due to the post-processing
property, and it can be served as a DP principal components. Kim and
Jung (2025) calls these process as a `g-DPPCA`.

In the implementation of the function `dp_pc_dir` with option
`g_dppca=TRUE`, we use the spherical transformation
$`g_{sph}(t) = t/\|t\|_2`$ to output differentially private PC
directions $`\widetilde{V}_{sph,m}`$. In this case, it holds that
$`\|g_{sph}\|_{\infty} = 1`$, and thus the variance of additive Gaussian
noise is set as \$ \_{, } = \$.

## Summary

The principal component direction step in `dppca` can be summarized as
follows.

1.  Start with a preprocessed data matrix $`X`$.
2.  Choose a direction estimation method.
3.  Obtain a direction matrix $`V_k`$.
4.  Compute projected scores $`Y = X V_k`$.
5.  Use the scores for private scree estimation or private score
    visualization.

The main distinction is whether $`V_k`$ is obtained from the ordinary
sample covariance matrix or from a differentially private robust PC
direction estimator.

## References

Minwoo Kim and Sungkyu Jung (2025), “Robust and differentially private
principal component analysis,” *Statistical Analysis and Data Mining*,
18(6), <https://doi.org/10.1002/sam.70053>

Jakob Raymaekers and Peter Rousseeuw (2019), “A generalized spatial sign
covariance matrix,” *Journal of Multivariate Analysis*, 171:94–111,
<https://doi.org/10.1016/j.jmva.2018.11.010>
