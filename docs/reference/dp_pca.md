# Estimate shared PCA quantities once

Summarizes an existing fit without new estimation or random sampling.

## Usage

``` r
dp_pca(
  X,
  privacy,
  scree,
  score = list(),
  axes = c(1, 2),
  cpp.option = TRUE,
  non_private = FALSE,
  V_dp = NULL,
  group = NULL
)

# S3 method for class 'dp_pca'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'dp_pca'
summary(object, ...)

# S3 method for class 'summary.dp_pca'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- X:

  A numeric matrix or data frame with observations in rows and variables
  in columns. Non-numeric data-frame columns are excluded with a
  warning; numeric values must be finite. At least two observations and
  two numeric features are required. If directions must be estimated,
  the pairwise-separation condition of
  [`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
  also applies.

- privacy:

  A
  [`privacy_control()`](https://yejinjo0220.github.io/dppca/reference/privacy_control.md)
  object allocating budget to `loading`, `scree`, and `score` (equal
  thirds by default).

- scree:

  Required named list with `k = NULL` (all components),
  `method = "clipped"`, `control` (required), and `mono = TRUE`. `k`,
  when supplied, must be between `2` and the number of numeric features.
  Use
  [`clipped_control()`](https://yejinjo0220.github.io/dppca/reference/clipped_control.md),
  [`pmwm_control()`](https://yejinjo0220.github.io/dppca/reference/pmwm_control.md),
  or
  [`huber_control()`](https://yejinjo0220.github.io/dppca/reference/huber_control.md)
  for tuning.

- score:

  Named list specifying `method = "sparse"` or `"add"` and
  `bins = c(20, 20)` by default. Exactly one histogram method is used.

- axes:

  Two distinct PC indices between `1` and the fitted `k` (`c(1, 2)` by
  default).

- cpp.option:

  Whether to use the compiled spherical-Kendall implementation (`TRUE`
  by default) when estimating private directions; ignored if `V_dp` is
  supplied.

- non_private:

  Whether to compute and store ordinary PCA reference results (`FALSE`
  by default).

- V_dp:

  Optional `p` by `p` private PC direction matrix (`NULL` by default),
  where `p` is the number of numeric PCA features. Columns must be
  orthonormal and rows must follow the feature order in `X`; matching
  row names are required when both sets of names exist.

- group:

  Optional group column name in `X` or vector of group labels (`NULL` by
  default). A selected group column is excluded from PCA features, even
  when numeric; group labels and membership are assumed public.

- x:

  A `summary.dp_pca` object.

- digits:

  Significant digits used when printing tables.

- ...:

  Additional arguments. The summary calculation takes no extra options;
  print methods pass formatting arguments to the printed tables.

- object:

  A `dp_pca()` result.

## Value

A `dp_pca` object containing private `directions`, `scree`, `pve`,
`hist`, `frame`, and optional `score_groups`. `nonprivate` contains the
optional ordinary PCA reference, or `NULL` when `non_private = FALSE`.
`metadata` contains data dimensions, privacy allocations, methods,
controls, grouping information and axes. The print method returns the
fitted object invisibly.

A `summary.dp_pca` object. `importance` contains DP eigenvalues, PVE and
cumulative PVE, with a column for each estimated component. `budget`
contains the allocated eps and delta for loading, scree, score frame and
score histogram, plus the total. Other fields describe data dimensions,
selected axes, methods and whether a non-private reference is stored.
Raw scores and loading matrices are not copied into the summary. The
print method returns its input invisibly.

## Details

A single private PC direction matrix is estimated or supplied and reused
for scree and score estimation. Scree values estimate variance along
those directions using disjoint paired differences; an odd final
observation is unused for scree, while the score histogram uses all
rows. PVE is normalized over the retained `k` components, not
necessarily all original variables.

The selected scree method uses the entire scree allocation. Of the score
epsilon budget, 35% goes to the private frame and 65% to the histogram;
all score delta goes to the histogram. With supplied `V_dp`, the loading
budget records its previous release rather than being reassigned.

Grouped data use one direction estimate and private frame. Their
disjoint histograms use parallel composition when group labels and
membership are public and fixed. The pooled histogram is post-processed
from the grouped histograms.

Input is not automatically centered or standardized. Data-dependent
preprocessing needs separate privacy accounting, and sample size,
schema, selected components, and tuning settings are assumed public. The
loading conditions in
[`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
also apply. The optional non-private reference uses sample-mean-centered
ordinary PCA scores and is not private.

PVE is normalized over the `k` retained components; a total of zero
produces zero PVE. When `k` is smaller than the number of variables, a
cumulative PVE of one does not describe all original variation.

The budget table reports allocated budgets, not cumulative privacy
spending. With supplied `V_dp`, the loading allocation records its prior
cost. The optional non-private comparison is excluded from the
importance table.

## Examples

``` r
data(gau, package = "dppca")
set.seed(123)
fit <- dp_pca(
  gau[1:200, ],
  privacy = privacy_control(eps = 5, delta = 1e-4),
  scree = list(k = 2, control = clipped_control(C_clip = 3)),
  score = list(method = "sparse", bins = c(10, 10)),
  non_private = FALSE
)
fit
#> Differentially Private PCA
#> Data: 200 observations x 20 variables
#> Components: 20 PC directions; 2 scree estimates
#> Score axes: PC1, PC2 
#> Methods: scree = clipped ; score = sparse 
#> Total privacy budget: eps = 5 ; delta = 1e-04 
#> Non-private reference: not stored 
summary(fit)
#> DP PCA Summary
#> Data: 200 observations x 20 variables
#> Components: 20 PC directions; 2 scree estimates
#> Score axes: PC1, PC2 
#> Methods: scree = clipped ; score = sparse 
#> 
#> DP component importance:
#>                   PC1    PC2
#> Eigenvalue     0.9521 0.8882
#> PVE            0.5174 0.4826
#> Cumulative PVE 0.5174 1.0000
#> 
#> Proportions are normalized over 2 estimated components.
#> 
#> Allocated privacy budget:
#>                    eps     delta
#> Loading         1.6667 3.333e-05
#> Scree           1.6667 3.333e-05
#> Score frame     0.5833 0.000e+00
#> Score histogram 1.0833 3.333e-05
#> Total           5.0000 1.000e-04
#> 
#> Non-private reference: not stored 
fit$directions
#>             PC1          PC2         PC3         PC4          PC5         PC6
#> V1   0.19413754  0.257797730  0.04284186  0.04243370 -0.180986738  0.01031595
#> V2   0.11370700 -0.316668014 -0.09719936 -0.04104061 -0.051496659 -0.39437698
#> V3   0.03944075 -0.008402996 -0.56728045 -0.18160434  0.263973413  0.10997240
#> V4  -0.10639338  0.187437548 -0.07491462 -0.03641717  0.127752997  0.11420546
#> V5  -0.10669240 -0.216073257 -0.08186225  0.31926062  0.136641992 -0.05545757
#> V6  -0.32776539  0.177324372 -0.22849591  0.37152230  0.377287804  0.33032518
#> V7  -0.22221606 -0.277439785 -0.07938759 -0.23823756  0.221347153 -0.42131868
#> V8   0.40209813 -0.120061102  0.04825535  0.23888581  0.027810689 -0.07064820
#> V9   0.17523776  0.167307252  0.21189517 -0.34755632  0.208195520 -0.04333997
#> V10 -0.10742646 -0.232271653 -0.15039194 -0.01458754 -0.407721944  0.28317202
#> V11 -0.34598251  0.361965776  0.06602439  0.07273690 -0.339745625 -0.08800841
#> V12  0.01168949 -0.045508426 -0.12241860  0.45072498 -0.208720218 -0.34877728
#> V13 -0.30590591 -0.015398251  0.03620102 -0.15858824  0.304988304 -0.23591729
#> V14  0.30581684  0.041770263  0.10862067 -0.26894362  0.003417028  0.21134391
#> V15  0.21157477  0.101964588 -0.26594436  0.09882708 -0.100652539 -0.10144672
#> V16  0.06024710  0.409644962 -0.50981428 -0.24922919 -0.160088582 -0.18729867
#> V17  0.22213938  0.209131165 -0.20229657  0.09581129  0.054359938 -0.26663251
#> V18 -0.09569474  0.175655433  0.20587776 -0.07235810  0.024401079 -0.15433075
#> V19 -0.22197050 -0.335419405 -0.22713018 -0.30324872 -0.367207461  0.14536093
#> V20 -0.31823010  0.214994429  0.16807229 -0.05691178 -0.179024144 -0.22658499
#>             PC7          PC8          PC9        PC10         PC11
#> V1   0.09632647  0.007459308  0.350744312  0.20240800  0.201319529
#> V2  -0.12377897  0.280058004 -0.040612595  0.19382301 -0.215128480
#> V3   0.47880682 -0.303460936 -0.003158392  0.23080015  0.072262162
#> V4  -0.11416819 -0.081030328 -0.605096536  0.23796653 -0.317889232
#> V5  -0.29629347 -0.430493133  0.294142472  0.06977654 -0.030503769
#> V6  -0.20321276  0.324189049  0.082609141 -0.07527589 -0.113296170
#> V7   0.05518077  0.170793053  0.024183483 -0.15224224 -0.306194295
#> V8   0.05514612  0.055491412  0.050065573 -0.06722717 -0.047614894
#> V9  -0.13567729 -0.160172535  0.019129946  0.09263512  0.036758138
#> V10 -0.32172285 -0.173415784  0.037922744 -0.15623127 -0.123438781
#> V11  0.16808379  0.127570589  0.289880775  0.39707401 -0.365736929
#> V12  0.06171612 -0.418127570 -0.175217132  0.05030340 -0.105754115
#> V13 -0.10063612 -0.178302301  0.392304244  0.12230237  0.093999094
#> V14 -0.18187376 -0.218666906  0.157812928  0.13769438 -0.579040921
#> V15  0.24438094  0.260985905  0.151170607 -0.21658568 -0.221938230
#> V16 -0.28177015 -0.099159609  0.077838612 -0.41699846 -0.002004277
#> V17 -0.44047980  0.182621152 -0.141340211  0.44038140  0.261584168
#> V18 -0.18467804  0.004382302 -0.003231587 -0.22651652  0.070646563
#> V19 -0.09939377  0.156126454 -0.008815188  0.25635145  0.245253860
#> V20  0.15373256 -0.195072878 -0.261575307 -0.12029344  0.092946029
#>              PC12         PC13         PC14        PC15        PC16
#> V1   0.5668866313 -0.133404971 -0.291934422 -0.14715535 -0.19143711
#> V2   0.1644770016  0.274694819  0.408814730 -0.03394769 -0.03553233
#> V3   0.0638320743  0.140988741  0.158401571  0.03276296  0.16858380
#> V4   0.1488432168 -0.396052426 -0.002116551 -0.28906063 -0.02440537
#> V5  -0.1027713119  0.010176414 -0.169597791 -0.45614443  0.38369387
#> V6   0.2450786055  0.171084732 -0.101245882  0.33574991  0.04636009
#> V7   0.1068348319 -0.008201181 -0.519647191 -0.02103049 -0.01197135
#> V8   0.1207780823 -0.135096158  0.085643501  0.30125016  0.51703790
#> V9  -0.2295190710  0.172079127 -0.363853834  0.23104333  0.05663727
#> V10 -0.1427755026 -0.166584923  0.041949359  0.21758481 -0.08933684
#> V11 -0.2438775116  0.176901537  0.054447621 -0.05645084  0.14579442
#> V12  0.1889509152  0.081269740 -0.148423686  0.28233633 -0.32575583
#> V13  0.0008207009 -0.489246685  0.345085124  0.26635183 -0.21695517
#> V14  0.2337933205  0.012856402  0.030903125  0.20054413  0.07060153
#> V15 -0.1986294157 -0.538298803 -0.103939033 -0.04113782  0.10299726
#> V16  0.0124837317  0.181357302  0.063132469 -0.07781205 -0.02763371
#> V17 -0.2070058692 -0.107501508 -0.083746160  0.13315009  0.08387004
#> V18  0.4005474829 -0.014582246  0.262801967 -0.22127368  0.25428692
#> V19  0.2377216740 -0.097660632 -0.188553016  0.08245285  0.23984546
#> V20  0.0801137219 -0.059729093 -0.023776348  0.32334395  0.43014143
#>             PC17          PC18        PC19        PC20
#> V1  -0.051573915 -0.3612887358 -0.05991851  0.16134489
#> V2   0.138171570 -0.2550964360 -0.37908268  0.19570525
#> V3  -0.060787146  0.0195494377  0.12624319  0.28881436
#> V4   0.234888522 -0.2205648977  0.04429999 -0.06978389
#> V5  -0.003649037 -0.0535827789 -0.21890492 -0.02447647
#> V6   0.076629919  0.0003913972 -0.15165340  0.06612594
#> V7  -0.202576186 -0.0832092080  0.31726989  0.06025109
#> V8   0.245725091 -0.2790186075  0.39537545 -0.22477158
#> V9   0.549308047 -0.0390164354 -0.16180677  0.27442248
#> V10 -0.039556598 -0.2793965657  0.21405006  0.50983260
#> V11  0.133064460 -0.0065398413  0.24584230 -0.01004403
#> V12  0.205579269  0.3122324201  0.02830726 -0.02583897
#> V13  0.102272795 -0.0507838755  0.01095784 -0.14769852
#> V14 -0.333010186  0.2588696433 -0.13103962 -0.14076626
#> V15  0.112520781  0.2557827008 -0.34939817  0.21114325
#> V16  0.094282927 -0.1400452836  0.02538594 -0.32195027
#> V17 -0.365833358  0.1322774231  0.13939097  0.12692062
#> V18  0.128056622  0.4383242287  0.25771979  0.42335664
#> V19  0.247419590  0.2822400236 -0.08729676 -0.25314761
#> V20 -0.307327055 -0.2160725816 -0.37440964  0.05447205
fit$nonprivate # NULL by default
#> NULL
fit$metadata$privacy
#> $total
#>   eps delta 
#> 5e+00 1e-04 
#> 
#> $budgets
#> $budgets$loading
#>          eps        delta 
#> 1.666667e+00 3.333333e-05 
#> 
#> $budgets$scree
#>          eps        delta 
#> 1.666667e+00 3.333333e-05 
#> 
#> $budgets$score
#>          eps        delta 
#> 1.666667e+00 3.333333e-05 
#> 
#> 
#> attr(,"class")
#> [1] "dppca_privacy_control"

set.seed(123)
plots <- plot(fit) # DP panels only
#> Warning: The private histogram contains no positive probability mass; no synthetic score points were generated.

plots$private$plots$biplot

```
