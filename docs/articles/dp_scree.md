# DP scree in dppca

This article describes the scree estimation used in the `dppca` package.
The main idea is to rewrite each scree value as a mean estimation
problem.

The package considers three main methods as follows.

1.  **Clipped mean using Gaussian noise**
2.  **Huber GDP mean based on noisy gradient descent**
3.  **Private modified winsorized mean (PMWM)**

## Overview

Let
``` math

X \in \mathbb{R}^{n \times p}
```

be the preprocessed data matrix. This means that the data have been
centered, and possibly standardized.

Let
``` math

V_k = [v_1,\ldots,v_k] \in \mathbb{R}^{p \times k}
```

be the matrix of PC directions. Depending on the analysis, $`V_k`$ may
be non-private or differentially private.

The goal is to estimate the scree vector
$`\lambda = (\lambda_1,\ldots,\lambda_k)`$ in a differentially private
way.

In the sample PCA case, the $`\ell`$-th scree value is

``` math

\hat\lambda_\ell
=
v_\ell^\top \hat\Sigma v_\ell
\quad \text{where} \quad
\hat\Sigma = \frac{1}{n-1}X^\top X.
```

In `dppca`, we construct private estimates
$`\widetilde\lambda_1,\ldots,\widetilde\lambda_k`$ and use them to build
a differentially private scree plot.

### Scree Estimation in dppca

For the $`\ell`$-th principal component direction $`v_\ell`$, define the
score vector $`z_\ell = X v_\ell`$.

The sample variance of the $`\ell`$-th score vector is

``` math

\hat\lambda_\ell
=
\frac{1}{n-1}
\sum_{i=1}^n
z_{i\ell}^2
=
\frac{n}{n-1}
\left(
\frac{1}{n}\sum_{i=1}^n w_{i\ell}
\right)
\quad \text{where} \quad 
w_{i\ell} = z_{i\ell}^2.
```

Therefore, for each component $`\ell`$, scree estimation can be viewed
as a mean estimation problem for $`w_{1\ell},\ldots,w_{n\ell}`$.

In `dppca`, we want to construct a differentially private estimate of
$`\widehat{\lambda}_\ell`$, denoted by $`\widetilde{\lambda}_\ell`$. So,
we first privately estimate the mean of $`W_{1\ell},\ldots,W_{n\ell}`$,
and then rescale it by $`n/(n-1)`$.

``` math

\widetilde{\lambda}_\ell
=
\frac{n}{n-1}
\cdot
\operatorname{DP-Mean}\left(w_{1\ell},\ldots,w_{n\ell}\right).
```

## Method 1: Clipped scree estimation

The simplest way to privately estimate the mean of
$`w_{1\ell},\ldots,w_{n\ell}`$ is to first clip the values to a bounded
interval and then add Gaussian noise.

### Clipping

Choose a clipping threshold
``` math

C > 0.
```

Define the clipped observations by
$`w_{i\ell}^{\mathrm{clip}} = \min(w_{i\ell}, C)`$.

Then
``` math

0 \leq w_{i\ell}^{\mathrm{clip}} \leq C.
```

The clipped empirical mean is

``` math

\hat\mu_\ell^{\mathrm{clip}}
=
\frac{1}{n}
\sum_{i=1}^n
w_{i\ell}^{\mathrm{clip}}
=
\frac{1}{n}
\sum_{i=1}^n
\min(w_{i\ell}, C).
```

The corresponding non-private clipped scree estimate is
``` math

\hat\lambda_\ell^{\mathrm{clip}} = \frac{n}{n-1} \hat\mu_\ell^{\mathrm{clip}}.
```

### Sensitivity

Because each clipped observation lies in $`[0,C]`$, changing one
observation can change the mean by at most $`\frac{C}{n}`$.

After multiplying by $`n/(n-1)`$, the sensitivity of the clipped scree
estimate is

``` math

\Delta_\ell
\leq
\frac{n}{n-1}
\cdot
\frac{C}{n}
=
\frac{C}{n-1}.
```

### DP clipped scree estimate

For privacy parameters $`(\epsilon_\ell,\delta_\ell)`$, the Gaussian
mechanism releases

``` math

\widetilde\lambda_\ell^{\mathrm{clip}}
=
\hat\lambda_\ell^{\mathrm{clip}}
+
\xi_\ell,
```

where

``` math

\xi_\ell \sim N(0,\sigma_\ell^2)
```
and
``` math

\sigma_\ell
=
\frac{\Delta_\ell\sqrt{2\log(1.25/\delta_\ell)}}{\epsilon_\ell}
=
\frac{(C/(n-1))\sqrt{2\log(1.25/\delta_\ell)}}{\epsilon_\ell}.
```

### Parameter

Clipped scree estimation is simple and easy to implement, but it depends
on the choice of the clipping threshold $`C`$.

- If $`C`$ is too small, the scree values may be underestimated.
- If $`C`$ is too large, the required privacy noise may increase.

### Example usage

In the `dppca`, we set clipped-mean parameter by

``` math

\texttt{scree\_clipped\_control(C\_clip = C)}.
```

``` r
library(dppca)

# x: numeric data matrix or data frame
# epsilon, delta: privacy parameters

out <- dp_scree_plot(
  x,
  k = 5,
  epsilon = 1,
  delta = 1e-6,
  method = "clipped"
)

out
```

## Method 2: Huber scree estimation

The Huber method is based on the differentially private robust mean
estimator of [Yu, Ren, and Zhou (2024)](#ref-Yu2024). For each component
$`\ell`$, it estimates the mean of $`w_{1\ell},\ldots,w_{n\ell}`$ using
a robust loss function instead of the ordinary sample mean.

### Huber loss

For a robustification parameter $`\tau > 0`$, the Huber loss is

``` math

\rho_\tau(r)
=
\begin{cases}
\frac{1}{2}r^2, & |r| \leq \tau,\\
\tau |r| - \frac{1}{2}\tau^2, & |r| > \tau.
\end{cases}
```

For small residuals, the Huber loss behaves like squared error loss. For
large residuals, it grows linearly, which reduces the influence of
extreme observations.

The derivative of the Huber loss is the score function

``` math

\psi_\tau(r)
=
\rho_\tau'(r)
=
\begin{cases}
-\tau, & r < -\tau,\\
r, & |r| \leq \tau,\\
\tau, & r > \tau.
\end{cases}
```

Thus, $`\psi_\tau(r)`$ clips the residual $`r`$ to the interval
$`[-\tau,\tau]`$.

### Huber mean estimator

For the $`\ell`$-th component, the Huber mean estimator is defined as

``` math

\hat\mu_{\tau,\ell}
\in
\arg\min_{\mu \in \mathbb{R}}
\frac{1}{n}
\sum_{i=1}^n
\rho_\tau(w_{i\ell}-\mu).
```

The first-order condition is

``` math

\frac{1}{n}
\sum_{i=1}^n
\psi_\tau(w_{i\ell}-\mu)
=
0.
```

This equation shows that large residuals are clipped through the score
function. As a result, the estimator is less sensitive to outliers than
the ordinary mean.

### Noisy gradient descent

The Huber objective is optimized using noisy gradient descent.

Given the current estimate $`\mu_\ell^{(t)}`$ at iteration $`t`$, we
define the corresponding residuals as

``` math

r_{i\ell}^{(t)}
=
w_{i\ell}
-
\mu_\ell^{(t)},
\qquad i=1,\ldots,n.
```

The clipped score values are $`\psi_\tau(r_{i\ell}^{(t)})`$,

and the average gradient type quantity is

``` math

g_\ell^{(t)}
=
\frac{1}{n}
\sum_{i=1}^n
\psi_\tau(r_{i\ell}^{(t)}).
```

A non-private update would take the form

``` math

\mu_\ell^{(t+1)}
=
\mu_\ell^{(t)}
+
\eta_t g_\ell^{(t)}
\quad \text{where} \quad
\eta_t > 0 \quad \text{is the step size.}
```

To ensure privacy, Gaussian noise is added at each iteration.

``` math

\mu_\ell^{(t+1)} = \mu_\ell^{(t)} + \eta_t g_\ell^{(t)} + \xi_\ell^{(t)}
\quad \text{where} \quad
\xi_\ell^{(t)} \sim N(0,\sigma_{t,\ell}^2).
```

### Sensitivity of one gradient step

Since

``` math

\psi_\tau(r) \in [-\tau,\tau],
```

changing one observation can change the average score by at most
$`\frac{2\tau}{n}`$.

After multiplying by the step size $`\eta_t`$, the one-step sensitivity
is

``` math

\Delta_{\mathrm{step}} = \frac{2\eta_t \tau}{n}.
```

This sensitivity determines the scale of the Gaussian noise added at
each gradient step.

### DP Huber scree estimate

After $`T`$ noisy gradient descent steps, let the final private Huber
mean estimate be

``` math

\mu_\ell^{(T)}.
```

The final Huber scree estimate is

``` math

\widetilde\lambda_\ell^{\mathrm{Huber}}
=
\frac{n}{n-1}
\mu_\ell^{(T)}.
```

The Huber method is useful when the squared scores $`w_{i\ell}`$ may be
heavy-tailed or affected by outliers.

### Robustification parameter

The robustification parameter $`\tau`$ controls the trade-off between
bias, robustness, and privacy.

A small $`\tau`$ gives stronger clipping and more robustness, but may
increase the bias of the estimator. A large $`\tau`$ behaves more like
the ordinary mean, but may be less robust and may require more noise.

In the Huber DP mean, $`\tau`$ is chosen using a scale quantity related
to the second moment. A typical theoretical form is

``` math

\tau
\asymp
\sqrt{m_2}
\left(
\frac{\epsilon n}
{\sqrt{(d+\log n)\log n}}
\right)^{1/2},
```

where

``` math

m_2 = \mathbb{E}\|X-\mu\|_2^2 \quad \text{is a second-moment scale.}
```

Because $`m_2`$ is usually unknown and can be sensitive to outliers, it
may need to be estimated robustly and privately. we use a [private
robust estimator for
$`m_2`$](https://yejinjo0220.github.io/dppca/articles/algorithms.html#alg-m2)
described in Algorithm 2.

### Other control parameters

The Huber scree method also uses additional control parameters for noisy
gradient descent and for the private scale-proxy step used to estimate
$`m_2`$. The default values are chosen based on the recommendations in
the [Yu, Ren, and Zhou (2024)](#ref-Yu2024)

- `mu0`: Initial value for noisy gradient descent. (default:
  $`\mu_\ell^{(0)} = 0`$)

- `eta0`: Step size for noisy gradient descent. (default:
  $`\eta_0 = 1`$)

- `T`: Number of noisy gradient descent iterations. (default:
  $`T = \lfloor \log n \rfloor`$)

- `M`: Number of blocks used in the private estimator for $`m_2`$.
  (default: $`M = \lfloor \sqrt{n}/2 \rfloor`$)

- `k_min_m2` and `k_max_m2`: Lower and upper dyadic bin indices used in
  the private histogram step for estimating $`m_2`$. The histogram
  searches over scale levels $`2^k`$ for

  ``` math

   k_{\min} \le k \le k_{\max}.
   
  ```

- `m2_frac`: Fraction of the Huber scree privacy budget used to
  privately estimate $`m_2`$.

  If $`(\epsilon_{\mathrm{scree}},\delta_{\mathrm{scree}})`$ is the
  budget for Huber scree estimation, then
  $`(\epsilon_{m_2},\delta_{m_2}) = \text{m2_frac}\cdot
  (\epsilon_{\text{scree}},\delta_{\text{scree}})`$, while the remaining
  budget
  $`(\epsilon_{\text{gd}},\delta_{\text{gd}}) = (1-\text{m2_frac})\cdot (\epsilon_{\text{scree}},\delta_{\text{scree}})`$
  is used for Huber noisy gradient descent.

### Example usage

In `dppca`, the Huber scree estimator is controlled by

``` math

\texttt{scree\_huber\_control(
mu0,\ eta0,\ T,\ M,\ k\_min\_m2,\ k\_max\_m2,\ m2\_frac)}.
```

``` r
library(dppca)

out <- dp_scree_plot(
  x,
  k = 5,
  epsilon = 1,
  delta = 1e-6,
  method = "huber",
  control = scree_huber_control(
    mu0 = 0,
    eta0 = 1,
    T = NULL,
    M = NULL,
    k_min_m2 = -40,
    k_max_m2 = 40,
    m2_frac = 1/4
  )
)
```

## Method 3: Private modified winsorized mean scree estimation

The PMWM method is based on the private modified winsorized mean of
[Ramsay and Spicker (2025)](#ref-Ramsay2025). For each component
$`\ell`$, it estimates the mean of $`w_{1\ell},\ldots,w_{n\ell}`$ by
privately estimating tail quantiles, winsorizing the data, and then
adding noise to the winsorized mean.

### Non-private modified winsorized mean

The PMWM method builds on the non-private modified winsorized mean of
[Lugosi and Mendelson (2021)](#ref-Lugosi2021).

The non-private modified winsorized mean starts with a clipping
proportion $`0 < p < \frac{1}{2}`$, which determines the lower and upper
tail quantiles used for winsorization.

Let $`\hat\xi_p \quad\text{and}\quad \hat\xi_{1-p}`$ be empirical lower
and upper quantiles.

Define the function

``` math

\phi_{a,b}(x)
=
\begin{cases}
a, & x < a,\\
x, & a \leq x \leq b,\\
b, & x > b.
\end{cases}
```

The non-private modified winsorized mean is

``` math

\hat\mu_p
=
\frac{1}{n}
\sum_{i=1}^n
\phi_{\hat\xi_p,\hat\xi_{1-p}}(X_i).
```

The idea is to first estimate the lower and upper tail cutoffs, and then
winsorize the data by replacing extreme values with the corresponding
cutoffs.

### Sample splitting notation

In the theoretical description, the data may be split into two subsets:

- a quantile estimation subset with size $`n_q`$,
- a mean estimation subset with size $`n_m`$.

If the full sample size is $`n`$, a simple split is
``` math

n_q = n_m = \frac{n}{2}.
```

In practice, the implementation may use all available observations at
each step instead of splitting the sample.

### Private quantile estimation

PMWM uses private quantile estimates instead of non-private empirical
quantiles. For component $`\ell`$, let the clipping proportion be
``` math

p_\ell = \zeta_\ell.
```

A theoretical choice has the form

``` math

\zeta_\ell
=
16\eta
+
\frac{112}{3}
\frac{
\log\left(32(\beta(u-\ell)/(\beta-1)\vee 1)/\delta\right)
}
{n_q},
```

where

- $`\eta`$: Contamination level,
- $`\beta > 1`$: Log-grid parameter in the private quantile estimator,
- $`\ell, u`$: Lower and upper search bounds used in private quantile
  estimation.
- $`\delta`$: Confidence parameter controlling the high-probability
  accuracy statement of the PMWM estimator.
- $`n_q`$: Number of observations used for quantile estimation.

For practical implementation, the paper suggests using the clipping
proportion

``` math

p
=
\frac{C}{n_q} \vee \eta,
```

where $`C`$ is a user-chosen trimming constant. This is also the choice
used in the `dppca` package.

### Winsorized mean

Let $`L_\ell`$ and $`U_\ell`$ be private estimates of the lower and
upper quantiles.

``` math

L_\ell \approx Q_{\zeta_\ell}(w_{1\ell},\ldots,w_{n\ell}),
\quad \text{and} \quad
U_\ell \approx Q_{1-\zeta_\ell}(w_{1\ell},\ldots,w_{n\ell}).
```

Define the winsorized observations by

``` math

w_{i\ell}^{\mathrm{win}}
=
\min\left\{
\max(w_{i\ell},L_\ell),
U_\ell
\right\}.
```

Equivalently,

- if $`w_{i\ell} < L_\ell`$, replace it by $`L_\ell`$;
- if $`L_\ell \leq w_{i\ell} \leq U_\ell`$, keep it unchanged;
- if $`w_{i\ell} > U_\ell`$, replace it by $`U_\ell`$.

Thus,
``` math

L_\ell
\leq
w_{i\ell}^{\mathrm{win}}
\leq
U_\ell.
```

Using the mean estimation subset $`I_m`$, define
$`\hat\mu_\ell^{\mathrm{win}} = \frac{1}{n_m} \sum_{i\in I_m} w_{i\ell}^{\mathrm{win}}`$.

The corresponding non-private winsorized scree estimate is

``` math

\hat\lambda_\ell^{\mathrm{PMWM,np}}
=
\frac{n}{n-1}
\hat\mu_\ell^{\mathrm{win}}.
```

### Sensitivity

Because all winsorized observations lie in $`[L_\ell,U_\ell]`$, the
sensitivity of the winsorized mean is

``` math

\Delta_\ell^{(\mu)}
=
\frac{U_\ell-L_\ell}{n_m}.
```

After multiplying by $`n/(n-1)`$, the scree sensitivity is
$`\Delta_\ell^{(\lambda)} = \frac{n}{n-1} \cdot \frac{U_\ell-L_\ell}{n_m}`$.

### DP PMWM scree estimate

The final PMWM scree estimate can be written as

``` math

\widetilde\lambda_\ell^{\mathrm{PMWM}}
=
\hat\lambda_\ell^{\mathrm{PMWM,np}} + Z_\ell 
\quad \text{where} \quad
Z_\ell \sim N(0,\sigma_\ell^2).
```

For privacy parameters $`(\epsilon_{M,\ell},\delta_{M,\ell})`$, the
noise scale is

``` math

\sigma_\ell
=
\frac{
\Delta_\ell^{(\lambda)}
\sqrt{2\log(1.25/\delta_{M,\ell})}
}
{\epsilon_{M,\ell}}
=
\frac{
\left(
\frac{n}{n-1}
\cdot
\frac{U_\ell-L_\ell}{n_m}
\right)
\sqrt{2\log(1.25/\delta_{M,\ell})}
}
{\epsilon_{M,\ell}}.
```

### Privacy budget splitting

PMWM uses privacy budget for both private quantile estimation and
private mean estimation.

For component $`\ell`$, let the total component-level budget be
$`(\epsilon_\ell,\delta_\ell)`$.

This can be split as
``` math

\epsilon_\ell
=
\epsilon_{Q,\ell}
+
\epsilon_{M,\ell},
\qquad
\delta_\ell
=
\delta_{Q,\ell}
+
\delta_{M,\ell}.
```

The quantile budget itself is used to estimate two quantiles, so it can
be split again.

``` math

\epsilon_{Q,\ell}
=
\epsilon_{q1,\ell}
+
\epsilon_{q2,\ell},
\qquad
\delta_{Q,\ell}
=
\delta_{q1,\ell}
+
\delta_{q2,\ell}.
```

### Parameters

The PMWM scree estimator uses additional parameters for private quantile
estimation and winsorization.

- `beta`: Log-binning base used in the private quantile estimator.

  It determines the spacing of the geometric search grid and must
  satisfy $`\beta > 1`$.

- `a`, `b`: Lower and upper search bounds supplied to the private
  quantile. The private lower and upper clipping cutoffs are searched
  within this range.

- `trim_const`, `eta`: Parameters used to set the practical clipping
  proportion
  ``` math

  p
  =
  \min\left\{
  \max\left(\frac{\texttt{trim\_const}}{n_q}, \eta\right),
  0.49
  \right\}.
  ```
  Here, `trim_const / n_q` controls the baseline clipping level, while
  `eta` gives a lower bound reflecting the expected contamination level.

- `split_mode`: Logical value indicating whether the sample is split
  into two parts.

  If `TRUE`, one part is used for private quantile estimation and the
  other part is used for the winsorized mean step. If `FALSE`, all
  observations are used in both steps.

- `max_extra_bins`: Maximum number of additional log-grid bins searched
  beyond the largest occupied bin in the private quantile.

### Example usage

In `dppca`, the PMWM-specific parameters can be specified through

``` math

\text{scree_pmw_control(
beta = beta,a = a,b = b, trim_const = C, eta = eta split_mode = split)}
```

``` r
library(dppca)

out <- dp_scree_plot(
  x,
  k = 5,
  epsilon = 1,
  delta = 1e-6,
  method = "pmw",
  control = scree_pmw_control(
    beta = 1.01,
    a = 0,
    b = 100,
    trim_const = 10,
    eta = 0.01,
    split_mode = TRUE,
    max_extra_bins = 1000
  )
)

out
```

## Post-processing

Because of the added privacy noise, the raw DP scree estimates

``` math

(\widetilde\lambda_1,\ldots,\widetilde\lambda_k)
```

may not have the usual scree shape. They may not be decreasing, and some
values may be negative.

In ordinary PCA, scree values satisfy

``` math

\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_k \ge 0.
```

Therefore, `dppca` can apply post-processing to make the DP scree
estimates nonnegative and decreasing.

``` math

\widetilde\lambda_1^{\mathrm{mono}}
\ge
\widetilde\lambda_2^{\mathrm{mono}}
\ge
\cdots
\ge
\widetilde\lambda_k^{\mathrm{mono}}
\ge
0.
```

If monotone post-processing is used, the PVE can be computed from the
post-processed scree values.

``` math

\widetilde{\operatorname{PVE}}_\ell^{\mathrm{mono}}
=
\frac{\widetilde\lambda_\ell^{\mathrm{mono}}}
{\sum_{j=1}^k \widetilde\lambda_j^{\mathrm{mono}}}.
```

This step only modifies the already released DP estimates, so it does
not use any additional privacy budget.

## References

Myeonghun Yu. Zhao Ren. Wen-Xin Zhou. “Gaussian differentially private
robust mean estimation and inference”. Bernoulli 30 (4) 3059 - 3088,
November 2024. <https://doi.org/10.3150/23-BEJ1706>

Kelly Ramsay and Dylan Spicker. (2025). “Improved
subsample-and-aggregate via the private modified winsorized mean”. arXiv
preprint. <https://arxiv.org/abs/2501.14095>

Gábor Lugosi. Shahar Mendelson. “Robust multivariate mean estimation:
The optimality of trimmed mean.” Ann. Statist. 49 (1) 393 - 410,
February 2021. <https://doi.org/10.1214/20-AOS1961>
