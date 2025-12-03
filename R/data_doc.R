#' European map genotype dataset
#'
#' @format A data frame with 1387 rows and 21 variables:
#' \describe{
#'   \item{X.1, ..., X.20}{Numeric PCA scores for 20 principal components.}
#'   \item{color}{Group labels.}
#' }
#'
#' @usage data(eur_map)
#' @examples
#' data(eur_map)
#' str(eur_map)
"eur_map"


#' Gaussian mixture synthetic dataset
#'
#' @format A data frame with:
#' \describe{
#'   \item{X}{A numeric matrix (5000 × 20).}
#'   \item{color}{Group labels ("red", "orange", "green", "blue", "purple").}
#' }
#'
#' @usage data(gaussian_groups)
#' @examples
#' data(gaussian_groups)
#' table(gaussian_groups$color)
"gaussian_groups"


#' Numeric subset of Adult dataset
#'
#' @format A data frame with numeric variables:
#' \describe{
#'   \item{age}{age}
#'   \item{education_num}{years of education}
#'   \item{capital_gain}{capital gains}
#'   \item{capital_loss}{capital losses}
#'   \item{hours_per_week}{hours worked per week}
#' }
#'
#' @usage data(adult)
#' @examples
#' data(adult)
#' summary(adult)
"adult"
