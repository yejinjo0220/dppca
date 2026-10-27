# dppca: Differentially Private Principal Component Analysis Visualization

Tools for estimating principal component directions and constructing
differentially private PCA visualizations. The package provides private
principal component directions; clipped, modified-winsorized, and Huber
scree and proportion-of-variance-explained summaries; additive and
sparse two-dimensional score histograms; group-wise score
visualizations; synthetic score sampling from released private
histograms; and an interactive Shiny application.

## Details

[`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
stores a shared set of private estimates for summaries and scree, score,
loading, variable, and biplot displays. Its default settings disable
sample-mean centering and omit non-private reference results.
[`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md),
[`dp_scree()`](https://yejinjo0220.github.io/dppca/reference/dp_scree.md),
and
[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
also support print, summary and plot methods for their stored results.
For example, call
[`dp_loading()`](https://yejinjo0220.github.io/dppca/reference/dp_loading.md)
with a privacy budget, then plot the fit. The Shiny app also supports
diagnostic comparisons that include non-private data summaries.

Data frames may include numeric and non-numeric columns. PCA uses the
numeric feature columns and warns when non-numeric columns are excluded.
Set `group` to a column name in
[`dp_score()`](https://yejinjo0220.github.io/dppca/reference/dp_score.md)
or [`dp_pca()`](https://yejinjo0220.github.io/dppca/reference/dp_pca.md)
to obtain grouped score results; the selected group column is excluded
from the PCA features, even when it is numeric. Leave `group = NULL` for
an ungrouped analysis. Summaries identify the selected group variable
and its groups. Numeric features are used without automatic centering or
standardization. See the individual functions for preprocessing
assumptions and privacy-budget accounting.

## See also

Useful links:

- <https://github.com/yejinjo0220/dppca>

- <https://yejinjo0220.github.io/dppca/>

- Report bugs at <https://github.com/yejinjo0220/dppca/issues>

## Author

**Maintainer**: Yejin Jo <yejinjo0220@gmail.com>

Authors:

- Minwoo Kim <mwkim.stat@gmail.com>

- Sungkyu Jung <sungkyu@snu.ac.kr>
