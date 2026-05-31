# Package index

## Direction

Functions for estimating private principal component directions.

- [`dp_pc_dir()`](https://yejinjo0220.github.io/dppca/reference/dp_pc_dir.md)
  : Estimate principal component directions

## Scree

Functions for estimating and visualizing differentially private scree
values. The control helpers specify method-specific tuning options for
the clipped, private modified winsorized mean, and Huber-type scree
estimators.

- [`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md)
  : Differentially private scree values
- [`dp_scree_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_scree_plot.md)
  : Plot differentially private scree estimates
- [`clipped_control()`](https://yejinjo0220.github.io/dppca/reference/clipped_control.md)
  : Control options for clipped scree estimation
- [`pmwm_control()`](https://yejinjo0220.github.io/dppca/reference/pmwm_control.md)
  : Control options for private modified winsorized scree estimation
- [`huber_control()`](https://yejinjo0220.github.io/dppca/reference/huber_control.md)
  : Control options for Huber scree estimation

## Score

Functions for constructing differentially private PCA score summaries,
score plots, and group-wise score visualizations.

- [`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
  : Differentially private score histograms
- [`dp_score_plot()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot.md)
  : Plot differentially private score histograms
- [`dp_score_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_group.md)
  : Group-wise differentially private score histograms
- [`dp_score_plot_group()`](https://yejinjo0220.github.io/dppca/reference/dp_score_plot_group.md)
  : Plot group-wise differentially private score histograms

## shiny app

shiny app for exploring differentially private PCA scree and score plots
through a graphical interface.

- [`dppca_app()`](https://yejinjo0220.github.io/dppca/reference/dppca_app.md)
  : Launch the dppca Shiny app

## Data

Datasets included for examples and demonstrations.

- [`adult`](https://yejinjo0220.github.io/dppca/reference/adult.md) :
  Adult numeric data
- [`gau`](https://yejinjo0220.github.io/dppca/reference/gau.md) : Five
  Gaussian clusters
- [`gau_g`](https://yejinjo0220.github.io/dppca/reference/gau_g.md) :
  Five Gaussian clusters with group labels
