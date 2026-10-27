# Package index

## Integrated DP PCA

Estimate shared private PCA quantities and display five visualizations.

- [`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
  [`print(`*`<dp_pca>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
  [`summary(`*`<dp_pca>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
  [`print(`*`<summary.dp_pca>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
  : Estimate shared PCA quantities once
- [`plot(`*`<dp_pca>`*`)`](https://yejinjo0220.github.io/dppca/reference/plot.dp_pca.md)
  : Plot an integrated PCA fit

## DP Loading

Estimate and visualize private principal component directions.

- [`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
  [`print(`*`<dp_loading>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
  [`summary(`*`<dp_loading>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
  [`print(`*`<summary.dp_loading>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
  : Estimate differentially private PC directions
- [`plot(`*`<dp_loading>`*`)`](https://yejinjo0220.github.io/dppca/reference/plot.dp_loading.md)
  [`print(`*`<dp_loading_plot>`*`)`](https://yejinjo0220.github.io/dppca/reference/plot.dp_loading.md)
  [`summary(`*`<dp_loading_plot>`*`)`](https://yejinjo0220.github.io/dppca/reference/plot.dp_loading.md)
  [`plot(`*`<dp_loading_plot>`*`)`](https://yejinjo0220.github.io/dppca/reference/plot.dp_loading.md)
  : Plot stored private PC directions
- [`loading_plot_control()`](https://yejinjo0220.github.io/dppca/reference/loading_plot_control.md)
  : Control loading plot appearance

## DP Scree and PVE

Estimate private scree values and proportions of variance explained.

- [`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
  [`print(`*`<dp_scree>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
  [`summary(`*`<dp_scree>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
  [`print(`*`<summary.dp_scree>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
  : Differentially private scree values
- [`plot(`*`<dp_scree>`*`)`](https://yejinjo0220.github.io/dppca/reference/plot.dp_scree.md)
  : Plot stored private scree estimates
- [`scree_plot_control()`](https://yejinjo0220.github.io/dppca/reference/scree_plot_control.md)
  : Control options for scree and PVE plots
- [`clipped_control()`](https://yejinjo0220.github.io/dppca/reference/clipped_control.md)
  : Control options for clipped scree estimation
- [`pmwm_control()`](https://yejinjo0220.github.io/dppca/reference/pmwm_control.md)
  : Control options for private modified winsorized scree estimation
- [`huber_control()`](https://yejinjo0220.github.io/dppca/reference/huber_control.md)
  : Control options for Huber scree estimation

## DP Score Plots

Private histograms, sampled score plots, and grouped scores.

- [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
  [`summary(`*`<dp_score>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
  [`summary(`*`<dp_score_group>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
  [`print(`*`<summary.dp_score>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
  [`print(`*`<dp_score>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
  [`print(`*`<dp_score_group>`*`)`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
  : Differentially private score histograms
- [`plot(`*`<dp_score>`*`)`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score.md)
  [`plot(`*`<dp_score_plot>`*`)`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score.md)
  [`summary(`*`<dp_score_plot>`*`)`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score.md)
  [`print(`*`<dp_score_plot>`*`)`](https://yejinjo0220.github.io/dppca/reference/plot.dp_score.md)
  : Plot a stored private score histogram
- [`score_plot_control()`](https://yejinjo0220.github.io/dppca/reference/score_plot_control.md)
  : Control options for score plots
- [`sampling_control()`](https://yejinjo0220.github.io/dppca/reference/sampling_control.md)
  : Control options for synthetic score sampling

## Privacy Controls

Allocate privacy budgets for the estimation steps.

- [`privacy_control()`](https://yejinjo0220.github.io/dppca/reference/privacy_control.md)
  : Allocate privacy budgets for shared PCA estimation

## Shiny Application

Launch the optional interactive app (separate from the tutorials).

- [`dppca_app()`](https://yejinjo0220.github.io/dppca/reference/dppca_app.md)
  : Launch the dppca Shiny app

## Example Datasets

- [`gau`](https://yejinjo0220.github.io/dppca/reference/gau.md) : Five
  Gaussian Cluster Dataset
- [`gau_g`](https://yejinjo0220.github.io/dppca/reference/gau_g.md) :
  Five Gaussian Cluster Dataset with Group Labels
- [`vowel`](https://yejinjo0220.github.io/dppca/reference/vowel.md) :
  Vowel Recognition Dataset
