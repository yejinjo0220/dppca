# Differentially Private Scree Estimation

## Overview

A scree plot is one of the most common summaries of Principal Component
Analysis (PCA). It visualizes the amount of variance explained by each
principal component.

In ordinary PCA, the scree values are the eigenvalues of the sample
covariance matrix. In a privacy-preserving analysis, however, directly
releasing these values may reveal information about individual
observations. The goal of differentially private scree estimation is to
release noisy but useful estimates of PCA scree values under a formal
differential privacy guarantee.

This article describes the scree estimation framework used in the
`dppca` package. The main idea is to rewrite each scree value as a mean
estimation problem and then privatize that mean using one of several
differentially private mean estimators.

The package considers three main methods:

1.  **Clipped mean using Gaussian noise**.
2.  **Huber GDP mean based on noisy gradient descent**.
3.  **Private modified winsorized mean (PMWM)**.

## Goal

Let
``` math

X_{\mathrm{proc}} \in \mathbb{R}^{n \times p}
```
be the preprocessed data matrix. In most PCA workflows, this means that
the data have been centered, and possibly standardized.

Let
``` math

V_k = [v_1,\ldots,v_k] \in \mathbb{R}^{p \times k}
```
be the matrix of principal component directions. Depending on the
analysis, $`V_k`$ may be non-private or differentially private.

The score matrix is
``` math

Y = X_{\mathrm{proc}} V_k.
```

The goal is to estimate the scree vector
``` math

\lambda = (\lambda_1,\ldots,\lambda_k)
```
in a differentially private way.

In the sample PCA case, the $`\ell`$-th scree value is
``` math

\hat\lambda_\ell
=
v_\ell^\top \hat\Sigma v_\ell,
```
where
``` math

\hat\Sigma
=
\frac{1}{n-1}X_{\mathrm{proc}}^\top X_{\mathrm{proc}}.
```

In `dppca`, we construct private estimates
``` math

\widetilde\lambda_1,\ldots,\widetilde\lambda_k
```
and use them to build a differentially private scree plot.

## From PCA scores to scree values

For the $`\ell`$-th principal component direction $`v_\ell`$, define the
score vector
``` math

y_\ell = X_{\mathrm{proc}} v_\ell.
```

Let $`y_{i\ell}`$ be the $`i`$-th score on component $`\ell`$, and
define
``` math

\bar y_\ell
=
\frac{1}{n}
\sum_{i=1}^n y_{i\ell}.
```

The sample variance of the $`\ell`$-th score vector is
``` math

\hat\lambda_\ell
=
\frac{1}{n-1}
\sum_{i=1}^n
(y_{i\ell}-\bar y_\ell)^2.
```

Define
``` math

w_{i\ell}
=
(y_{i\ell}-\bar y_\ell)^2.
```

Then
``` math

\hat\lambda_\ell
=
\frac{1}{n-1}
\sum_{i=1}^n w_{i\ell}
=
\frac{n}{n-1}
\left(
\frac{1}{n}\sum_{i=1}^n w_{i\ell}
\right).
```

Therefore, scree estimation can be viewed as a mean estimation problem:
``` math

\text{estimate}
\qquad
\frac{1}{n}\sum_{i=1}^n w_{i\ell}
\qquad
\text{privately}.
```

After obtaining a private mean estimate for the transformed values
``` math

w_{1\ell},\ldots,w_{n\ell},
```
we multiply by $`n/(n-1)`$ to obtain a private estimate of the scree
value.

## Common workflow

The differentially private scree estimation procedures share the
following workflow.

1.  Preprocess the data to obtain $`X_{\mathrm{proc}}`$.
2.  Compute principal component directions $`V_k`$, either non-privately
    or privately.
3.  Compute the score matrix
    ``` math

    Y = X_{\mathrm{proc}}V_k.
    ```
4.  For each component $`\ell = 1,\ldots,k`$, define
    ``` math

    w_{i\ell} = (y_{i\ell}-\bar y_\ell)^2.
    ```
5.  Estimate the mean of
    ``` math

    w_{1\ell},\ldots,w_{n\ell}
    ```
    using a differentially private mean estimator.
6.  Convert the private mean estimate into a private scree estimate by
    multiplying by $`n/(n-1)`$.
7.  Optionally apply post-processing so that the final scree vector is
    nonnegative and monotone decreasing.

This common structure makes it possible to compare multiple private mean
estimators within the same PCA scree estimation framework.

## Method 1: Clipped scree estimation

The simplest way to privately estimate the mean of
$`w_{1\ell},\ldots,w_{n\ell}`$ is to first clip the values to a bounded
interval and then add Gaussian noise.

### Clipping

Choose a clipping threshold
``` math

C_\ell > 0.
```

Define the clipped observations by
``` math

w_{i\ell}^{\mathrm{clip}}
=
\min(w_{i\ell}, C_\ell).
```

Then
``` math

0 \leq w_{i\ell}^{\mathrm{clip}} \leq C_\ell.
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
\min(w_{i\ell}, C_\ell).
```

The corresponding non-private clipped scree estimate is
``` math

\hat\lambda_\ell^{\mathrm{clip}}
=
\frac{n}{n-1}
\hat\mu_\ell^{\mathrm{clip}}.
```

### Sensitivity

Because each clipped observation lies in $`[0,C_\ell]`$, changing one
observation can change the mean by at most
``` math

\frac{C_\ell}{n}.
```

After multiplying by $`n/(n-1)`$, the sensitivity of the clipped scree
estimate is
``` math

\Delta_\ell
\leq
\frac{n}{n-1}
\cdot
\frac{C_\ell}{n}
=
\frac{C_\ell}{n-1}.
```

### Gaussian mechanism

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
\frac{(C_\ell/(n-1))\sqrt{2\log(1.25/\delta_\ell)}}{\epsilon_\ell}.
```

This produces an $`(\epsilon_\ell,\delta_\ell)`$-differentially private
estimate for the $`\ell`$-th scree value.

### Interpretation

Clipped scree estimation is simple and easy to implement. However, the
choice of $`C_\ell`$ is important.

- If $`C_\ell`$ is too small, the estimator can be biased downward.
- If $`C_\ell`$ is too large, the sensitivity and therefore the privacy
  noise become large.

Thus, clipped estimation is useful as a baseline but can be sensitive to
tuning and to heavy-tailed score distributions.

## Method 2: Huber scree estimation

The Huber method estimates the mean of $`w_{1\ell},\ldots,w_{n\ell}`$
using a robust loss function instead of the ordinary sample mean.

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

Thus, $`\psi_\tau(r)`$ clips the residual contribution to the interval
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

Let $`\mu_\ell^{(t)}`$ be the current estimate at iteration $`t`$.
Define the residuals
``` math

r_{i\ell}^{(t)}
=
w_{i\ell}
-
\mu_\ell^{(t)}.
```

The clipped score values are
``` math

\psi_\tau(r_{i\ell}^{(t)}),
```
and the average gradient-type quantity is
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
\eta_t g_\ell^{(t)},
```
where $`\eta_t > 0`$ is the step size.

To ensure privacy, Gaussian noise is added at each iteration:
``` math

\mu_\ell^{(t+1)}
=
\mu_\ell^{(t)}
+
\eta_t g_\ell^{(t)}
+
\xi_\ell^{(t)},
```
where
``` math

\xi_\ell^{(t)}
\sim
N(0,\sigma_{t,\ell}^2).
```

### Sensitivity of one gradient step

Since
``` math

\psi_\tau(r) \in [-\tau,\tau],
```
changing one observation can change the average score by at most
``` math

\frac{2\tau}{n}.
```

After multiplying by the step size $`\eta_t`$, the one-step sensitivity
is
``` math

\Delta_{\mathrm{step}}
=
\frac{2\eta_t \tau}{n}.
```

This sensitivity is used to calibrate the Gaussian noise in each
gradient step.

### Choice of the robustification parameter

The robustification parameter $`\tau`$ controls the trade-off between
bias, robustness, and privacy.

A small $`\tau`$ gives stronger clipping and more robustness, but may
introduce bias. A large $`\tau`$ behaves more like the ordinary mean,
but may be less robust and may require more noise.

In the Huber GDP mean framework, $`\tau`$ is chosen using a scale
quantity related to the second moment. A typical theoretical form is
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

m_2 = \mathbb{E}\|X-\mu\|_2^2
```
is a second-moment scale.

Because $`m_2`$ is usually unknown and can be sensitive to outliers, it
may need to be estimated robustly and privately. A private histogram
learner and a private robust estimator for $`m_2`$ can be used for this
purpose.

### Final DP Huber scree estimate

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

## Method 3: Private modified winsorized mean scree estimation

The private modified winsorized mean (PMWM) method estimates the mean of
$`w_{1\ell},\ldots,w_{n\ell}`$ by privately estimating tail quantiles,
winsorizing the data, and then releasing a noisy winsorized mean.

This method is designed to improve robustness against extreme values
while maintaining differential privacy.

### Non-private modified winsorized mean

The non-private modified winsorized mean starts with a clipping
proportion
``` math

0 < p < \frac{1}{2}.
```

Let
``` math

\hat\xi_p
\quad\text{and}\quad
\hat\xi_{1-p}
```
be empirical lower and upper quantiles.

Define the projection function
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

The idea is to estimate the tails and then project extreme observations
back into the estimated central interval.

### Sample splitting notation

In the theoretical description, the data may be split into two subsets:

- a quantile estimation subset with size $`n_q`$,
- a mean estimation subset with size $`n_m`$.

If the full sample size is $`n`$, a simple split is
``` math

n_q = n_m = \frac{n}{2}.
```

In practice, variants may use all available observations in each step,
depending on the implementation.

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

- $`\eta`$ is an upper bound on the contamination level,
- $`\beta > 1`$ is the log-grid parameter in the private quantile
  estimator,
- $`\ell`$ and $`u`$ are lower and upper anchors,
- $`\delta`$ is a probability parameter,
- $`n_q`$ is the quantile estimation sample size.

In practical implementations, a simpler clipping proportion such as
``` math

\zeta_\ell
=
\frac{C}{n_q}\vee \eta
```
can be used, where $`C`$ is a user-chosen clipping count.

### Winsorization

Let $`L_\ell`$ and $`U_\ell`$ be private estimates of the lower and
upper quantiles:
``` math

L_\ell
\approx
Q_{\zeta_\ell}(w_{1\ell},\ldots,w_{n\ell}),
```
and
``` math

U_\ell
\approx
Q_{1-\zeta_\ell}(w_{1\ell},\ldots,w_{n\ell}).
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

### Winsorized mean

Using the mean estimation subset $`I_m`$, define
``` math

\hat\mu_\ell^{\mathrm{win}}
=
\frac{1}{n_m}
\sum_{i\in I_m}
w_{i\ell}^{\mathrm{win}}.
```

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
``` math

\Delta_\ell^{(\lambda)}
=
\frac{n}{n-1}
\cdot
\frac{U_\ell-L_\ell}{n_m}.
```

### Gaussian mechanism for the winsorized mean

The final PMWM scree estimate can be written as
``` math

\widetilde\lambda_\ell^{\mathrm{PMWM}}
=
\max\left\{
\hat\lambda_\ell^{\mathrm{PMWM,np}}
+
Z_\ell,
0
\right\},
```
where
``` math

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

The truncation at zero is a post-processing step and therefore does not
weaken the privacy guarantee.

### Privacy budget splitting

PMWM uses privacy budget for both private quantile estimation and
private mean estimation.

For component $`\ell`$, let the total component-level budget be
``` math

(\epsilon_\ell,\delta_\ell).
```

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
be split again:
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

By sequential composition, the full procedure satisfies the sum of the
privacy budgets used across the quantile and mean estimation steps.

## Post-processing

The raw differentially private scree estimates
``` math

(\widetilde\lambda_1,\ldots,\widetilde\lambda_k)
```
may fail to satisfy the usual scree constraints because of the added
privacy noise.

In ordinary PCA, scree values satisfy
``` math

\lambda_1
\geq
\lambda_2
\geq
\cdots
\geq
\lambda_k
\geq
0.
```

Therefore, `dppca` may apply monotone post-processing to obtain
``` math

\widetilde\lambda_1^{\mathrm{mono}}
\geq
\widetilde\lambda_2^{\mathrm{mono}}
\geq
\cdots
\geq
\widetilde\lambda_k^{\mathrm{mono}}
\geq
0.
```

This post-processing can be viewed as projecting the noisy scree vector
onto the cone of nonnegative decreasing sequences. Since post-processing
is applied after the private release, it does not consume additional
privacy budget.

## Explained variance ratio

After obtaining private scree values, one may also report the private
explained variance ratio:
``` math

\widetilde{\operatorname{EVR}}_\ell
=
\frac{\widetilde\lambda_\ell}
{\sum_{j=1}^k \widetilde\lambda_j}.
```

If monotone post-processing is used, the EVR can be computed from the
post-processed scree values:
``` math

\widetilde{\operatorname{EVR}}_\ell^{\mathrm{mono}}
=
\frac{\widetilde\lambda_\ell^{\mathrm{mono}}}
{\sum_{j=1}^k \widetilde\lambda_j^{\mathrm{mono}}}.
```

This is also a post-processing operation.

## Summary of methods

| Method | Main idea | Strength | Main tuning parameters |
|----|----|----|----|
| Clipped mean | Clip $`w_{i\ell}`$ and add Gaussian noise | Simple baseline | $`C_\ell`$, $`\epsilon_\ell`$, $`\delta_\ell`$ |
| Huber GDP mean | Optimize Huber loss with noisy gradient descent | Robust to outliers | $`\tau`$, step size, number of iterations |
| PMWM | Private quantiles, winsorization, and noisy mean | Robust to tails and contamination | $`\zeta_\ell`$, $`\beta`$, $`\ell,u`$, budget split |

## Example usage

The exact function arguments may depend on the installed version of
`dppca`. A typical workflow is as follows.

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

A Huber-based version may look like:

``` r
out_huber <- dp_scree_plot(
  x,
  k = 5,
  epsilon = 1,
  delta = 1e-6,
  method = "huber",
  tau = NULL
)
```

A PMWM-based version may look like:

``` r
out_pmwm <- dp_scree_plot(
  x,
  k = 5,
  epsilon = 1,
  delta = 1e-6,
  method = "pmwm",
  beta = 1.001,
  eta = 0
)
```

These examples are intentionally not evaluated in this vignette so that
the page can be built even when users have not loaded a specific
dataset.

## Practical recommendations

The following points are useful when applying private scree estimation.

1.  **Start with a moderate privacy budget.** Very small $`\epsilon`$
    can make the scree plot unstable because the added noise is large.
2.  **Use clipping as a baseline.** The clipped mean method is simple
    and useful for checking whether the pipeline works.
3.  **Use robust methods for heavy-tailed data.** If the squared scores
    have large outliers, Huber or PMWM may be preferable.
4.  **Check monotone post-processing.** Post-processing often improves
    the interpretability of noisy scree vectors.
5.  **Interpret small components carefully.** Later components usually
    have smaller scree values and are more easily dominated by privacy
    noise.

## References

Dwork, C. and Roth, A. (2014). *The Algorithmic Foundations of
Differential Privacy*. Foundations and Trends in Theoretical Computer
Science.

Lugosi, G. and Mendelson, S. (2021). Robust multivariate mean
estimation: The optimality of trimmed mean. *The Annals of Statistics*,
49, 393–410.

Ramsay, K. and Spicker, D. (2025). *Improved subsample-and-aggregate via
the private modified winsorized mean*.

Yu, M., Ren, Z., and Zhou, W.-X. (2024). Gaussian differentially private
robust mean estimation and inference. *Bernoulli*, 30, 3059–3088.

Kim, M. and Jung, S. (2025). *Robust and Differentially Private
Principal Component Analysis*.
