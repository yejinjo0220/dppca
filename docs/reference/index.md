# Package index

## Shiny app

Interactive Shiny app for exploring differentially private PCA scree and
score plots.

- [`dppca_app()`](https://yejinjo0220.github.io/dppca/reference/dppca_app.md)
  : Launch the dppca Shiny app

## DP PCA Visualization

Functions for constructing differentially private PCA visualizations.

- [`dp_score_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot.md)
  : Plot differentially private score histograms
- [`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
  : Plot differentially private scree curves

## DP PCA

Core functions for computing differentially private PCA score and scree
quantities.

- [`dp_pc_dir()`](https://yejinjo0220.github.io/dppca/reference/dp_pc_dir.md)
  : Compute principal component directions
- [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
  : Differentially private score histograms
- [`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
  : Compute differentially private scree estimates

## Method controls

Helper functions for specifying method-specific tuning parameters.

- [`clipped_control()`](https://yejinjo0220.github.io/dppca/reference/clipped_control.md)
  : Control parameters for the clipped DP scree estimator
- [`pmwm_control()`](https://yejinjo0220.github.io/dppca/reference/pmwm_control.md)
  : Control parameters for the PMWM DP scree estimator
- [`huber_control()`](https://yejinjo0220.github.io/dppca/reference/huber_control.md)
  : Control parameters for the Huber DP scree estimator

## Group PCA Visualization

Functions for group-wise differentially private PCA score visualization.

- [`dp_score_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_group.md)
  : Group-wise DP score histograms
- [`dp_score_plot_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot_group.md)
  : Plot group-wise DP score histograms

## Example datasets

Datasets included for examples and demonstrations.

- [`adult`](https://yejinjo0220.github.io/dppca/reference/adult.md) :
  Adult data example
- [`eur_map`](https://yejinjo0220.github.io/dppca/reference/eur_map.md)
  : Europe map PCA data
- [`eur_map_g`](https://yejinjo0220.github.io/dppca/reference/eur_map_g.md)
  : Group labels for Europe map data
- [`gau`](https://yejinjo0220.github.io/dppca/reference/gau.md) :
  Gaussian example data
- [`gau_g`](https://yejinjo0220.github.io/dppca/reference/gau_g.md) :
  Group labels for Gaussian example data
