# Package index

## DP PCA

Core functions for computing principal component directions and
differentially private scree and score summaries.

- [`dp_pc_dir()`](https://yejinjo0220.github.io/dppca/reference/dp_pc_dir.md)
  : Estimate principal component directions
- [`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
  : Differentially private scree values
- [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
  : Differentially private score histograms

## DP PCA visualization

Plotting functions for differentially private PCA scree curves and score
histograms.

- [`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
  : Plot differentially private scree estimates
- [`dp_score_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot.md)
  : Plot differentially private score histograms

## Group score visualization

Functions for group-wise differentially private PCA score summaries and
visualizations.

- [`dp_score_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_group.md)
  : Group-wise differentially private score histograms
- [`dp_score_plot_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot_group.md)
  : Plot group-wise differentially private score histograms

## Method controls

Control objects for method-specific tuning parameters used in private
scree estimation.

- [`clipped_control()`](https://yejinjo0220.github.io/dppca/reference/clipped_control.md)
  : Control options for clipped scree estimation
- [`pmwm_control()`](https://yejinjo0220.github.io/dppca/reference/pmwm_control.md)
  : Control options for private modified winsorized scree estimation
- [`huber_control()`](https://yejinjo0220.github.io/dppca/reference/huber_control.md)
  : Control options for Huber scree estimation

## Shiny app

Interactive Shiny app for exploring differentially private PCA scree and
score plots.

- [`dppca_app()`](https://yejinjo0220.github.io/dppca/reference/dppca_app.md)
  : Launch the dppca Shiny app

## Example datasets

Datasets included for examples and demonstrations.

- [`adult`](https://yejinjo0220.github.io/dppca/reference/adult.md) :
  Adult numeric data
- [`gau`](https://yejinjo0220.github.io/dppca/reference/gau.md) : Five
  Gaussian clusters
- [`gau_g`](https://yejinjo0220.github.io/dppca/reference/gau_g.md) :
  Five Gaussian clusters with group labels
