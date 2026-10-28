# Differentially Private PC Directions

In the `dppca` package, a first step in DP-PCA is to estimate the PC
directions privately. For this purpose, [Kim and Jung
(2025)](#ref-Kim2025) proposed **g-DPPCA**, which applies a matrix
Gaussian mechanism to the generalized multivariate Kendall’s tau matrix
based on the generalized spatial-sign transformations introduced by
[Raymaekers and Rousseeuw (2019)](#ref-Raymaekers2019). The `dppca`
package implements g-DPPCA using the spherical transformation. The
resulting DP PC directions are used in subsequent PCA visualizations and
can be directly visualized through loading plots.

For the definitions of ordinary PCA, see [PCA
Background](https://yejinjo0220.github.io/dppca/articles/method-pca-background.md).
For R examples, arguments, and plot customization, see the [Tutorial of
DP Loading
Plot](https://yejinjo0220.github.io/dppca/articles/tutorials-loading-plots.md).

## 1. Spherical Kendall matrix and PC directions

Let X\in\mathbb R^{n\times p} contain observations
x_1,\ldots,x_n\in\mathbb R^p in its rows. As described in the [PCA
Background](https://yejinjo0220.github.io/dppca/articles/method-pca-background.md),
ordinary sample PC directions are the eigenvectors of the sample
covariance matrix \widehat\Sigma =\frac{1}{n-1}\sum\_{i=1}^n (x_i-\bar
x)(x_i-\bar x)^\top, \qquad \bar x=\frac1n\sum\_{i=1}^n x_i.

Equivalently, the sample covariance can be expressed using pairwise
differences:

\widehat\Sigma =\frac{1}{n(n-1)}\sum\_{i\<j} (x_i-x_j)(x_i-x_j)^\top.

The **sample spherical Kendall matrix** similarly uses pairwise
differences, but normalizes each difference to unit length:
\displaystyle \widehat K\_{\mathrm{sph}} =\frac{2}{n(n-1)}\sum\_{i\<j}
\frac{(x_i-x_j)(x_i-x_j)^\top} {\\x_i-x_j\\\_2^2}.

This expression is defined when the pairwise differences are nonzero.
For numerical stability, the implementation requires the stronger
condition \\x_i-x_j\\\_2\>\sqrt{\texttt{.Machine\\double.eps}} for every
pair, and rejects nonfinite norms. Thus, distinct observations can still
be rejected when their separation is too small in the chosen measurement
units.

The matrices \widehat K\_{\mathrm{sph}} and \widehat\Sigma generally
have different sample eigenvectors. However, under suitable elliptical
distributional assumptions, their population counterparts share
principal eigendirections. Spherical Kendall PCA therefore provides an
alternative way to estimate PC directions, rather than an exact
reconstruction of the sample covariance eigenvectors.

More generally, a generalized spatial sign can be written as

g_h(z)=h(\\z\\\_2)\frac{z}{\\z\\\_2}, \qquad z\ne0,

where h is a radial weight function. For a bounded transformation g_h,
the corresponding **generalized multivariate Kendall’s tau matrix** is
estimated by

\widehat K\_{g_h} =\frac{2}{n(n-1)}\sum\_{i\<j}
g_h\\\left(\frac{x_i-x_j}{\sqrt2}\right)
g_h\\\left(\frac{x_i-x_j}{\sqrt2}\right)^\top.

The spherical case corresponds to h(r)=1, giving g_h(z)=z/\\z\\\_2 for
z\ne0. The
[`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
function implements this spherical case.

### Translation invariance

The spherical Kendall matrix is invariant to a common translation of the
observations. For any c\in\mathbb R^p, (x_i-c)-(x_j-c)=x_i-x_j, \qquad
\widehat K\_{\mathrm{sph}}(X-\mathbf1_n c^\top) =\widehat
K\_{\mathrm{sph}}(X).

Consequently, estimating spherical Kendall directions does not require
computing or releasing the sample mean.

## 2. DP PC directions

We consider fixed-size replacement adjacency, in which two neighboring
datasets have the same number of observations and differ in one row.
Since changing one observation affects at most n-1 terms in the pairwise
sum, the Frobenius sensitivity of the spherical Kendall matrix satisfies
\Delta_F =\sup\_{X\sim X'} \bigl\\\widehat K\_{\mathrm{sph}}(X)
-\widehat K\_{\mathrm{sph}}(X')\bigr\\\_F \le\frac4n.

The mechanism adds symmetric Gaussian noise to the spherical Kendall
matrix. For privacy parameters \varepsilon_L\>0 and 0\<\delta_L\<1,
define \sigma_K =\frac{4\sqrt{2\log(1.25/\delta_L)}}{n\varepsilon_L}.

This classical Gaussian calibration has a standard sufficient privacy
guarantee for 0\<\varepsilon_L\<1. The implementation also accepts
larger positive values, but a separate calibration check is required to
establish their privacy guarantee; the function does not perform that
check.

A symmetric noise matrix E\in\mathbb R^{p\times p} is generated with
independent entries on and below the diagonal, such that

E\_{jj}\sim N(0,\sigma_K^2), \qquad E\_{j\ell}\sim
N(0,\sigma_K^2/2)\quad(j\<\ell), \qquad E\_{\ell j}=E\_{j\ell}.

The perturbed matrix is \widetilde{K}=\widehat{K}\_{\mathrm{sph}}+E. The
**private PC directions** are obtained from the eigen-decomposition

\widetilde{K} =\widetilde{V} \widetilde{\Theta}\widetilde{V}^\top,
\qquad \widetilde{V}=\[\widetilde{v}\_1,\ldots,\widetilde{v}\_p\],

where the columns of \widetilde{V} are orthonormal eigenvectors ordered
by decreasing eigenvalues of \widetilde{K}.

Importantly, the eigenvalues of \widehat{K} are not same as the
covariance eigenvalues \widehat\lambda\_\ell defined in the non-private
ordinary PCA. The [DP
Scree](https://yejinjo0220.github.io/dppca/articles/method-dp-scree.md)
procedure estimates variance along private PC directions separately.

## 3. DP loading plots in `dppca`

A **DP loading plot** displays the coefficients of the private PC
direction vectors. Let \widetilde{v}\_{j\ell} denote the coefficient of
variable j in private PC direction \ell. For two selected components a
and b, with distinct a,b\in\\1,\ldots,p\\, variable j is represented by

L_j^{\mathrm{DP}} =(\widetilde{v}\_{ja},\widetilde{v}\_{jb}).

Large absolute coefficients indicate that the corresponding variables
have larger weights in a selected PC direction. The signs distinguish
opposing contributions along an axis, but eigenvector signs are
arbitrary: replacing \widetilde{v}\_\ell with -\widetilde{v}\_\ell does
not change the underlying direction. For comparison, the non-private
ordinary PCA loading coordinates are L_j^{\mathrm{NP}}=(\widehat
v\_{ja},\widehat v\_{jb}) where \widehat V consists of eigenvectors of
\widehat\Sigma.

The package supports **vector**, **point**, and **heatmap** displays of
loading coefficients. For a selected pair of PCs, vector and point plots
encode the same coordinates L_j^{\mathrm{DP}}. A **vector plot** draws
an arrow from (0,0) to L_j^{\mathrm{DP}}; a **point plot** places a
point at that coordinate. A **heatmap** displays the signed entries
\widetilde{v}\_{j\ell} for selected variables and PCs.

For examples and plotting options, see the [DP Loading Plot
Tutorial](https://yejinjo0220.github.io/dppca/articles/tutorials-loading-plots.md).

## References

- Kim, M. and Jung, S. (2025). *Robust and differentially private
  principal component analysis*. Statistical Analysis and Data Mining,
  18(6), e70053. <https://doi.org/10.1002/sam.70053>

- Raymaekers, J. and Rousseeuw, P. J. (2019). *A generalized spatial
  sign covariance matrix*. Journal of Multivariate Analysis, 171,
  94–111. <https://doi.org/10.1016/j.jmva.2018.11.010>

- Dwork, C. and Roth, A. (2014). *The Algorithmic Foundations of
  Differential Privacy*. Foundations and Trends in Theoretical Computer
  Science, 9(3–4), 211–407. <https://doi.org/10.1561/0400000042>
