#' Adult numeric data
#'
#' A numerical subset of the Adult dataset from the UCI Machine Learning
#' Repository. The original Adult dataset is based on data extracted from
#' the 1994 United States Census database.
#'
#' This package dataset retains five numerical variables:
#' `age`, `education_num`, `capital_gain`, `capital_loss`, and
#' `hours_per_week`. These selected numerical variables contain no missing
#' values in the original data file used here. The resulting dataset contains
#' 32,561 observations.
#'
#' Since the variables have substantially different units and scales, scaling
#' is recommended before applying PCA-based methods.
#'
#' @format A data frame with 32,561 rows and 5 columns:
#' \describe{
#'   \item{age}{Age of the individual.}
#'   \item{education_num}{Number of years of education.}
#'   \item{capital_gain}{Capital gain.}
#'   \item{capital_loss}{Capital loss.}
#'   \item{hours_per_week}{Number of working hours per week.}
#' }
#'
#' @source
#' UCI Machine Learning Repository, Adult dataset.
#' The Adult dataset is licensed under the Creative Commons Attribution 4.0
#' International (CC BY 4.0) license.
#'
#' @references
#' Becker, B. and Kohavi, R. (1996).
#' Adult dataset. UCI Machine Learning Repository.
#'
"adult"


#' Five Gaussian clusters
#'
#' A simulated 20-dimensional Gaussian cluster dataset used as an example for
#' principal component analysis and differentially private PCA visualization.
#'
#' The dataset contains 5,000 observations in 20 dimensions. It consists of
#' five groups, with 1,000 observations in each group. Each group is generated
#' from a multivariate normal distribution. The group mean vectors and
#' covariance matrices are chosen differently across groups so that the data
#' contain both separated and partially overlapping cluster structures.
#'
#' @format A data frame with 5,000 rows and 20 columns:
#' \describe{
#'   \item{V1}{Simulated numerical variable.}
#'   \item{V2}{Simulated numerical variable.}
#'   \item{V3}{Simulated numerical variable.}
#'   \item{V4}{Simulated numerical variable.}
#'   \item{V5}{Simulated numerical variable.}
#'   \item{V6}{Simulated numerical variable.}
#'   \item{V7}{Simulated numerical variable.}
#'   \item{V8}{Simulated numerical variable.}
#'   \item{V9}{Simulated numerical variable.}
#'   \item{V10}{Simulated numerical variable.}
#'   \item{V11}{Simulated numerical variable.}
#'   \item{V12}{Simulated numerical variable.}
#'   \item{V13}{Simulated numerical variable.}
#'   \item{V14}{Simulated numerical variable.}
#'   \item{V15}{Simulated numerical variable.}
#'   \item{V16}{Simulated numerical variable.}
#'   \item{V17}{Simulated numerical variable.}
#'   \item{V18}{Simulated numerical variable.}
#'   \item{V19}{Simulated numerical variable.}
#'   \item{V20}{Simulated numerical variable.}
#' }
#'
#' @source
#' Simulated by the authors of the package. The data-generating code is
#' available in `data-raw/make_datasets.R` in the package source repository.
#'
"gau"


#' Five Gaussian clusters with group labels
#'
#' A grouped version of [`gau`] with an additional group label column.
#'
#' The dataset contains 5,000 observations in 20 dimensions and one group label
#' column. There are five groups, and each group contains 1,000 observations.
#'
#' @format A data frame with 5,000 rows and 21 columns:
#' \describe{
#'   \item{V1}{Simulated numerical variable.}
#'   \item{V2}{Simulated numerical variable.}
#'   \item{V3}{Simulated numerical variable.}
#'   \item{V4}{Simulated numerical variable.}
#'   \item{V5}{Simulated numerical variable.}
#'   \item{V6}{Simulated numerical variable.}
#'   \item{V7}{Simulated numerical variable.}
#'   \item{V8}{Simulated numerical variable.}
#'   \item{V9}{Simulated numerical variable.}
#'   \item{V10}{Simulated numerical variable.}
#'   \item{V11}{Simulated numerical variable.}
#'   \item{V12}{Simulated numerical variable.}
#'   \item{V13}{Simulated numerical variable.}
#'   \item{V14}{Simulated numerical variable.}
#'   \item{V15}{Simulated numerical variable.}
#'   \item{V16}{Simulated numerical variable.}
#'   \item{V17}{Simulated numerical variable.}
#'   \item{V18}{Simulated numerical variable.}
#'   \item{V19}{Simulated numerical variable.}
#'   \item{V20}{Simulated numerical variable.}
#'   \item{group}{Group label. One of `group1`, `group2`, `group3`,
#'   `group4`, or `group5`.}
#' }
#'
#' @source
#' Simulated by the authors of the package. The data-generating code is
#' available in `data-raw/make_datasets.R` in the package source repository.
#'
"gau_g"
