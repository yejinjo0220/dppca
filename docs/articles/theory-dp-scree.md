# Differentially Private Eigenvalues

After obtaining differentially private (DP) principal component (PC)
directions, the next step is to privately estimate the variation along
these directions. In ordinary PCA, the eigenvalues of the sample
covariance matrix represent the variance explained by the corresponding
PC directions. For this purpose, the `dppca` package uses a
paired-difference construction and three DP mean estimators: **clipped
mean**, **private modified winsorized mean (PMWM)**, and **Huber mean**.
The resulting estimates are post-processed to obtain nonnegative scree
values and the corresponding proportions of variance explained (PVE).

This article describes the mathematical formulations, tuning parameters,
and privacy accounting of these methods. For R implementations and
visualization examples, see the [DP Scree Plot
Tutorial](https://yejinjo0220.github.io/dppca/articles/tutorials-scree-plots.md).
For background, see [PCA
Background](https://yejinjo0220.github.io/dppca/articles/theory-pca-background.md)
and [Theory of DP
Loading](https://yejinjo0220.github.io/dppca/articles/theory-dp-loading.md).

## 1. DP eigenvalue estimation

Let X_1 and X_2 be independent random vectors in \mathbb R^p, with
covariance matrix \Sigma. For any fixed unit vector v,

\mathbb E\left\[\frac12\\v^\top(X_1-X_2)\\^2\right\] = v^\top\Sigma v.

When v is a population covariance eigenvector, this quantity is its
eigenvalue.

Motivated by this identity, the `dppca` package estimates the variances
along privately estimated PC directions using squared differences of
projected observations.

For observations x_1,\ldots,x_n, let
\widetilde{V}\_k=(\widetilde{v}\_1,\ldots,\widetilde{v}\_k) contain k
privately estimated orthonormal directions. Define the projected
coordinates

Y=X\widetilde{V}\_k,\qquad y\_{i\ell}=x_i^\top\widetilde{v}\_\ell.

Randomly partition the observations into m=\lfloor n/2\rfloor disjoint
pairs (a_r,b_r), r=1,\ldots,m. One observation is unused when n is odd.
For each pair and direction, define

w\_{r\ell}=\frac12(y\_{a_r,\ell}-y\_{b_r,\ell})^2
=\frac12\\\widetilde{v}\_\ell^\top(x\_{a_r}-x\_{b_r})\\^2.

The average \bar w\_\ell=m^{-1}\sum\_{r=1}^m w\_{r\ell} is a natural
estimator of projected variation. For fixed observations and directions,

\mathbb E\_{\mathrm{pair}}\[\bar w\_\ell\mid X,\widetilde{V}\_k\]
=\widetilde{v}\_\ell^\top\widehat\Sigma\widetilde{v}\_\ell, \qquad
\widehat\Sigma=\frac1{n-1}\sum\_{i=1}^n (x_i-\bar x)(x_i-\bar x)^\top.

**Derivation of the pairing identity**

Conditional on the data and private PC directions, each of the
\binom{n}{2} distinct pairs is equally likely to appear at any given
pair position.

\begin{aligned} \mathbb E\_{\mathrm{pair}} \[\bar w\_\ell\mid
X,\widetilde{V}\_k\] &= \frac1m\sum\_{r=1}^m \mathbb E\_{\mathrm{pair}}
\[w\_{r\ell}\mid X,\widetilde{V}\_k\]\\\[6pt\] &= \frac1m\sum\_{r=1}^m
\frac1{\binom n2} \sum\_{i\<j} \frac12(y\_{i\ell}-y\_{j\ell})^2\\\[6pt\]
&= \frac{1}{n(n-1)} \sum\_{i\<j}(y\_{i\ell}-y\_{j\ell})^2\\\[6pt\] &=
\frac{1}{n-1} \sum\_{i=1}^n (y\_{i\ell}-\bar y\_\ell)^2 \quad (\because
\sum\_{i\<j}(y\_{i\ell}-y\_{j\ell})^2 =n\sum_i(y\_{i\ell}-\bar
y\_\ell)^2) \\\[6pt\] &= \frac{1}{n-1} \sum\_{i=1}^n \left\\
\widetilde{v}\_\ell^\top(x_i-\bar x) \right\\^2 \quad (\because
y\_{i\ell}=x_i^\top\widetilde{v}\_\ell \quad \\ \quad \bar y\_\ell=\bar
x^\top\widetilde{v}\_\ell)\\\[6pt\] &= \widetilde{v}\_\ell^\top
\underbrace{ \left\\ \frac{1}{n-1} \sum\_{i=1}^n (x_i-\bar x)(x_i-\bar
x)^\top \right\\ }\_{=\widehat\Sigma} \widetilde{v}\_\ell\\\[6pt\] &=
\widetilde{v}\_\ell^\top \widehat\Sigma \widetilde{v}\_\ell.
\end{aligned}

Therefore, the paired-difference mean is unbiased for the sample
variance along each private PC direction with respect to random pairing.
This property holds before applying any private mean estimation
procedure. —

To ensure differential privacy, the `dppca` package applies a private
mean estimation procedure to the columns of the paired-difference matrix
W=(w\_{r\ell})\in\mathbb R^{m\times k}:

(\widetilde{\lambda}\_1,\ldots,\widetilde{\lambda}\_k) =
\operatorname{DP\text{-}Mean}(W).

Here, \operatorname{DP\text{-}Mean} denotes one of three private mean
estimation methods: **clipped mean, private modified winsorized mean
(PMWM), or Huber mean**.

## 2. Clipped mean

The clipped mean bounds each nonnegative paired squared difference at a
publicly chosen threshold C\_{\mathrm{clip}}\>0:

w\_{r\ell}^{C}=\min\\w\_{r\ell},C\_{\mathrm{clip}}\\,\qquad \bar
w\_\ell^{C}=\frac1m\sum\_{r=1}^m w\_{r\ell}^{C}.

Let \bar w^{C}=(\bar w_1^{C},\ldots,\bar w_k^{C})^\top denote the vector
of clipped means. Since the observations are divided into disjoint
pairs, changing one observation affects at most one row of W,
conditional on the private PC directions and a fixed pairing. As each
clipped value lies in \[0,C\_{\mathrm{clip}}\], the \ell_2-sensitivity
of the mean vector is bounded by \Delta_2(\bar w^{C})\le
\frac{\sqrt{k}\\C\_{\mathrm{clip}}}{m}.

To ensure DP, Gaussian noise is added to the vector of clipped means:

\widetilde{\lambda}^{C,\mathrm{raw}}=\bar w^{C}+Z,\qquad Z\sim
N_k(0,\sigma_C^2I_k),\qquad
\sigma_C=\frac{\sqrt{k}\\C\_{\mathrm{clip}}}{m\mu_C}.

Here, \mu_C is the Gaussian-DP parameter calibrated from the scree
privacy allocation
(\varepsilon\_{\mathrm{scree}},\delta\_{\mathrm{scree}}).

**Tuning parameter.** The clipped estimator has one required tuning
parameter, C\_{\mathrm{clip}}. A smaller threshold reduces sensitivity
and the noise scale, but can increase downward clipping bias. A larger
threshold retains more paired values but requires more noise. The
threshold must be chosen using public information or a separately
accounted private procedure; there is no universally suitable value.

## 3. Private modified winsorized mean (PMWM)

PMWM limits the influence of extreme paired values using privately
estimated winsorization cutoffs. Following [Ramsay and Spicker
(2025)](#ref-ramsay-spicker), the method first estimates private lower
and upper quantiles and then releases a noisy mean of the winsorized
data. Private quantiles use the one-sided unbounded mechanism of [Durfee
(2023)](#ref-durfee).

### 3.1. Allocation of pair rows

The PMWM procedure consists of two stages: quantile estimation and
winsorized mean estimation. The m pair rows of W can either be reused in
both stages or divided into disjoint subsets.

Let I_q denote the indices used for quantile estimation and I_M those
used for mean estimation, with m_q=\|I_q\| and m_M=\|I_M\|. The two
stages may reuse all m pair rows, in which case I_q=I_M=\\1,\ldots,m\\,
or use disjoint subsets, with m_q=\lfloor m/2\rfloor and m_M=m-m_q.

### 3.2. Private quantiles and cutoff bounds

For each component, define the trimming proportion as

\zeta = \min\left\\ \max\left( \frac{c\_{\mathrm{trim}}}{m_Q},
\eta\_{\mathrm{trim}} \right), 0.49 \right\\,

where c\_{\mathrm{trim}}\>0 is a baseline trimming constant and
0\le\eta\_{\mathrm{trim}}\<0.5 is a lower-bound tuning value.

The lower and upper cutoffs are estimated using the private quantiles at
levels \zeta and 1-\zeta, respectively:

L\_\ell^{\mathrm{priv}} = \mathcal Q\_{\varepsilon_Q,\beta}
(W\_{I_Q,\ell},\zeta), \qquad U\_\ell^{\mathrm{priv}} = \mathcal
Q\_{\varepsilon_Q,\beta} (W\_{I_Q,\ell},1-\zeta).

Here, \mathcal Q\_{\varepsilon_Q,\beta} denotes the one-sided unbounded
private quantile mechanism of [Durfee (2023)](#ref-durfee), with privacy
parameter \varepsilon_Q and geometric grid parameter \beta\>1. Since the
squared pair differences are nonnegative, the quantile mechanism uses
zero as a public lower bound for its search. The parameter \beta
controls the grid resolution, with values closer to one producing a
finer grid.

After quantile estimation, the private cutoffs are restricted to a
prespecified public interval \[a,b\], where a\<b:

L\_\ell = \min\\b,\max\\a,L\_\ell^{\mathrm{priv}}\\\\, \qquad U\_\ell =
\min\\b,\max\\a,U\_\ell^{\mathrm{priv}}\\\\.

If the resulting upper cutoff is below the lower cutoff, the interval is
reset to \[a,b\]. These are bounds on the estimated cutoffs; the
original paired values do not have to lie inside \[a,b\].

### 3.3. Winsorized mean and Gaussian perturbation

Define the winsorization rule and component means by

\phi\_{L\_\ell,U\_\ell}(w)=\min\\U\_\ell,\max\\L\_\ell,w\\\\, \qquad
\bar w\_\ell^{\mathrm{PMWM}}= \frac1{m_M}\sum\_{r\in
I_M}\phi\_{L\_\ell,U\_\ell}(w\_{r\ell}).

Conditional on the privately selected cutoffs, the joint mean
sensitivity is bounded by

\Delta_2\le \frac1{m_M}\sqrt{\sum\_{\ell=1}^k
(U\_\ell-L\_\ell)^2}\le\frac{\sqrt{k}(b-a)}{m_M}.

To ensure DP, Gaussian noise is added to the vector of winsorized means:

\widetilde{\lambda}^{\mathrm{PMWM,raw}}= \bar w^{\mathrm{PMWM}}+Z,\qquad
Z\_\ell\sim N(0,\sigma\_\ell^2),\qquad
\sigma\_\ell=\frac{\sqrt{k}(U\_\ell-L\_\ell)}{m_M\mu_M},

Here, \mu_M is the Gaussian-DP parameter calibrated from the
mean-estimation privacy allocation (\varepsilon_M,\delta_M).

### 3.4. Privacy allocation and tuning choices

The scree epsilon budget is divided equally between private quantile
estimation and noisy winsorized mean estimation. For each of the k
components, two private quantiles are released to estimate the lower and
upper cutoffs. Thus, each of the 2k quantile releases uses \varepsilon_Q
=\frac{\varepsilon\_{\mathrm{scree}}}{4k}, with no delta allocation.

The k-dimensional winsorized mean vector is then released using the
Gaussian mechanism with privacy budget

(\varepsilon_M,\delta_M) = \left(
\frac{\varepsilon\_{\mathrm{scree}}}{2}, \delta\_{\mathrm{scree}}
\right).

By sequential composition, the two stages together satisfy
(\varepsilon\_{\mathrm{scree}},\delta\_{\mathrm{scree}})-DP, provided
that the respective mechanism guarantees hold.

The PMWM tuning parameters are:

| Parameter | Meaning and constraint |
|:---|:---|
| a,b | Public limits for the private cutoff estimates; a\<b. |
| c\_{\mathrm{trim}} | Positive baseline trimming constant in \zeta. |
| \eta\_{\mathrm{trim}} | Tail-trimming lower bound, 0\le\eta\_{\mathrm{trim}}\<0.5. |
| \beta | Geometric-grid base for private quantiles, \beta\>1. |
| I_Q,I_M | Choice between pair-row reuse and disjoint quantile/mean subsets. |

## 4. Huber mean

The Huber procedure estimates the means of squared pair differences
using private scale estimation and noisy gradient descent, drawing on
[Yu, Ren, and Zhou (2024)](#ref-yu-ren-zhou). Unlike clipping or PMWM,
the mean is found iteratively using a bounded Huber score.

### 4.1. Huber loss and score

For a threshold \tau\>0, the Huber loss function and its derivative are
defined as

\rho\_\tau(u)= \begin{cases} \frac12u^2,& \|u\|\le\tau,\\
\tau\|u\|-\frac12\tau^2,& \|u\|\>\tau, \end{cases} \qquad
\psi\_\tau(u)=\rho\_\tau'(u) =\min\\\tau,\max\\-\tau,u\\\\.

For each component \ell=1,\ldots,k, the corresponding non-private Huber
mean estimator is defined by

\widehat\mu\_{\tau\_\ell,\ell} \in \underset{\theta\in\mathbb
R}{\arg\min} \frac1m \sum\_{r=1}^m
\rho\_{\tau\_\ell}(w\_{r\ell}-\theta).

Since \|\psi\_{\tau\_\ell}(u)\|\le\tau\_\ell, extreme residuals have
bounded influence on the iterative updates.

### 4.2. Private scale and Huber thresholds

Before estimating the Huber mean, we privately estimate the variability
of the paired squared differences for each component. Then, the
resulting **private scale proxy** \widetilde{m}\_{2,\ell} is used to
determine the Huber threshold \tau\_\ell.

To obtain this scale proxy, the m paired squared differences are paired
again to form additional half squared differences. Their block medians
are then summarized using a private histogram over dyadic scales 2^r.
The public range r\_{\min}\le r\le r\_{\max} determines the candidate
scales, while M controls the number of blocks. The detailed procedure is
provided in the [private robust estimator for
m_2](https://yejinjo0220.github.io/dppca/articles/algorithms.html#alg-m2).

Let \rho\in(0,1) denote the proportion of the scree epsilon budget
allocated to private scale estimation. For k components, the privacy
allocations are

\varepsilon\_{\mathrm{scale},\ell} =
\frac{\rho\varepsilon\_{\mathrm{scree}}}{k}, \qquad
\varepsilon\_{\mathrm{GD}} = (1-\rho)\varepsilon\_{\mathrm{scree}}.

Using the private scale proxy and the privacy budget allocated to
gradient descent, the Huber threshold for each component is defined as

\tau\_\ell = \sqrt{\widetilde{m}\_{2,\ell}} \sqrt{ \frac{
(\varepsilon\_{\mathrm{GD}}/k)m }{ \sqrt{(1+\log m)\log m} } }.

### 4.3. Noisy gradient descent

Let \theta\_\ell^{(0)}=\mu_0 be a public initial value, \eta_0\>0 the
step size, and T the number of iterations. The private iterates are
projected to the nonnegative half-line:

\theta\_\ell^{(t+1)}= \max\left\\0,\\\theta\_\ell^{(t)}+
\frac{\eta_0}{m}\sum\_{r=1}^m
\psi\_{\tau\_\ell}(w\_{r\ell}-\theta\_\ell^{(t)})+
Z\_\ell^{(t)}\right\\.

The Gaussian perturbations are drawn independently across coordinates
and iterations with

Z\_\ell^{(t)}\sim N(0,\sigma\_{\ell,H}^2),\qquad \sigma\_{\ell,H}=
\frac{2\eta_0\tau\_\ell\sqrt{kT}}{m\mu\_{\mathrm{GD}}},

where \mu\_{\mathrm{GD}} is the Gaussian-DP parameter calibrated from
the gradient-descent privacy allocation
(\varepsilon\_{\mathrm{GD}},\delta\_{\mathrm{scree}}). After T steps,
\widetilde{\lambda}\_\ell^{\mathrm{Huber,raw}}=\theta\_\ell^{(T)}.

The scale step consumes \rho\varepsilon\_{\mathrm{scree}} and no delta;
the gradient step consumes
((1-\rho)\varepsilon\_{\mathrm{scree}},\delta\_{\mathrm{scree}}). Their
releases must be accounted for together. The method requires at least
four initial pair rows, corresponding to at least eight observations.

The Huber tuning parameters are:

| Parameter | Meaning and constraint |
|:---|:---|
| r\_{\min},r\_{\max} | Public integer bounds on the dyadic scale indices, with r\_{\min}\<r\_{\max}. |
| \rho | Fraction of scree epsilon assigned to private scale estimation, 0\<\rho\<1. |
| \mu_0 | Public initialization of the mean iterate. |
| \eta_0 | Positive gradient-descent step size. |
| T | Number of noisy gradient steps. |
| M | Number of blocks for private scale estimation. |

The iteration and block counts may be chosen as functions of public
sample size; the corresponding automatic choices used by the package are
described in
[`huber_control()`](https://yejinjo0220.github.io/dppca/reference/huber_control.md)

## 5. Nonnegative scree and PVE

Random noise can produce negative or nonmonotone raw scree estimates.
The reported estimates are truncated at zero and may additionally be
adjusted by nonincreasing isotonic regression. Writing the resulting
nonnegative sequence as
\widetilde{\lambda}\_1,\ldots,\widetilde{\lambda}\_k, the reported PVE
is

\widetilde{\mathrm{PVE}}\_\ell= \begin{cases}
\displaystyle\frac{\widetilde{\lambda}\_\ell}
{\sum\_{j=1}^k\widetilde{\lambda}\_j},&
\sum\_{j=1}^k\widetilde{\lambda}\_j\>0,\\\[6pt\] 0,&\text{otherwise}.
\end{cases}

## 6. DP scree/PVE plots in `dppca`

If the sum is positive, the k reported PVE values sum to one. When k\<p,
this is the **share of the estimated variation among the retained k
components**, not necessarily a fraction of the total variance in the
original p-dimensional data. If all reported scree values are zero,
every PVE is set to zero. Truncation, isotonic adjustment, and
normalization are post-processing of private estimates and do not
require additional privacy budget.

The resulting estimates
\widetilde{\lambda}\_1,\ldots,\widetilde{\lambda}\_k and
\widetilde{\mathrm{PVE}}\_1,\ldots,\widetilde{\mathrm{PVE}}\_k are
displayed in the **DP scree plot** and **DP PVE plot**, respectively.
For implementation details, plotting options, and examples, see the [DP
Scree Plot
Tutorial](https://yejinjo0220.github.io/dppca/articles/tutorials-scree-plots.md).

## References

Ramsay, K. and Spicker, D. (2025). *Improved subsample-and-aggregate via
the private modified winsorized mean*.
<https://arxiv.org/abs/2501.14095>

Durfee, D. (2023). *Unbounded differentially private quantile and
maximum estimation*. Advances in Neural Information Processing Systems,
36, 77691–77712. <https://doi.org/10.52202/075280-3394>

Yu, M., Ren, Z., and Zhou, W.-X. (2024). *Gaussian differentially
private robust mean estimation and inference*. Bernoulli, 30(4),
3059–3088. <https://doi.org/10.3150/23-BEJ1706>
