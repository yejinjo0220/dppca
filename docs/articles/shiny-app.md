# Interactive Shiny App

The **`dppca` Shiny app** provides an interactive interface for
estimating and visualizing differentially private principal component
analysis (PCA) results. You can explore the example datasets, choose
privacy budgets and estimation methods, and compare available
visualizations without writing R code.

**[Launch the online Shiny
app](https://yejinjo0220-dppca.share.connect.posit.cloud/)**

The online application is hosted on **Posit Connect Cloud**. It runs
separately from this documentation website, so you do not need to
install R or `dppca` to try the example analyses.

## 1. Available analyses

Use the tabs at the top of the results panel to switch between four
analyses.

| Tab | What you can explore |
|:---|:---|
| **PCA** | An integrated analysis with scree/PVE, score, loading, variable, and biplot views based on a shared private PCA fit. |
| **Loading** | Private PC directions displayed as vector, point, or heatmap plots. |
| **Scree** | Private scree values and proportions of variance explained (PVE), using the clipped, PMWM, or Huber estimators. |
| **Score** | Private two-dimensional score histograms and synthetic sampled points, using additive and/or sparse histogram methods. Grouped displays are available when group labels are selected. |

For more detail on the methods, see the [DP
Loading](https://yejinjo0220.github.io/dppca/articles/theory-dp-loading.md),
[DP
Scree](https://yejinjo0220.github.io/dppca/articles/theory-dp-scree.md),
and [DP
Score](https://yejinjo0220.github.io/dppca/articles/theory-dp-score.md)
articles.

## 2. Getting started

1.  **Open the [online
    app](https://yejinjo0220-dppca.share.connect.posit.cloud/).** Under
    **Data**, select **Example data** and choose *Gaussian*, *Gaussian
    with groups*, or *Vowel*.
2.  **Choose an analysis.** The **PCA** tab gives an overview of five
    complementary visualizations; the other tabs let you focus on one
    part of the analysis.
3.  **Set the privacy parameters.** Enter the desired epsilon and delta
    values. For integrated analyses, configure the privacy budget
    allocation using either total-budget proportions or individual
    component budgets.
4.  **Adjust the analysis settings.** Select the PC directions, scree or
    histogram method, and other parameters relevant to the active tab.
    Group labels can be used for PCA and Score analyses.
5.  **Run the analysis.** Click **Run PCA**, **Run Loading**, **Run
    Scree**, or **Run Score**, depending on the selected tab, to view
    the results.

After an analysis has run, you can adjust supported **plot settings** to
change the display without re-estimating the underlying private
quantities. For PCA, the **View** menu lets you inspect individual plots
or the integrated overview.

## 3. Data and comparison options

The app includes three example datasets: `gau` (synthetic Gaussian
data), `gau_g` (Gaussian data with group labels), and `vowel` (vowel
measurements). It also supports **CSV upload** through the **Data**
panel.

Use **Include non-private comparison** only when a non-private reference
is appropriate. Ordinary PCA results are **not** differentially private
and should not be released from sensitive data without authorization.

For the R interfaces behind the app, see the [Getting Started
tutorial](https://yejinjo0220.github.io/dppca/articles/tutorials-getting-started.md)
and [DP PCA Visualization
tutorial](https://yejinjo0220.github.io/dppca/articles/tutorials-dp-visualization.md).

## 4. Run the app locally

If you prefer to run the application on your own computer, install
`dppca` and `shiny`, then launch the bundled app:

``` r
install.packages("shiny")
# Install the current development version of dppca if needed:
# remotes::install_github("yejinjo0220/dppca")

library(dppca)
dppca_app()
```

You can also supply your own data using `dppca_app(X)` or
`dppca_app(X, group = "group")`; see
[`dppca_app()`](https://yejinjo0220.github.io/dppca/reference/dppca_app.md)
for details.

**Online application:**
[yejinjo0220-dppca.share.connect.posit.cloud](https://yejinjo0220-dppca.share.connect.posit.cloud/)
