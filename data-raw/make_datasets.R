# data-raw/make_datasets.R

# This script creates package datasets used in dppca.
# The processed datasets are saved in data/ as .rda files.
#
# To update package datasets, run this file from the package root.

# 1. Five Gaussian clusters ----------------------------------------------

set.seed(123)

n_groups    <- 5
n_per_group <- 1000
d           <- 20

means <- list(
  rep(0, d),
  c(rep(5, 5), rep(0, d - 5)),
  c(rep(-3, 5), rep(0, d - 5)),
  c(rep(0, 5), rep(-1, 5), rep(0, d - 10)),
  c(rep(0, 5), rep(3, 5), rep(0, d - 10))
)

Sigma_list <- list(
  diag(d),
  diag(d),
  5 * diag(d),
  5 * diag(d),
  8 * diag(d)
)

X_mix <- do.call(rbind, lapply(seq_len(n_groups), function(i) {
  MASS::mvrnorm(
    n = n_per_group,
    mu = means[[i]],
    Sigma = Sigma_list[[i]]
  )
}))

gau <- as.data.frame(X_mix)
names(gau) <- paste0("V", seq_len(d))

group_gau <- rep(paste0("group", seq_len(n_groups)), each = n_per_group)
gau_g <- data.frame(gau, group = group_gau)


# 2. Adult numeric data ---------------------------------------------------

colnames_adult <- c(
  "age", "workclass", "fnlwgt", "education", "education_num",
  "marital_status", "occupation", "relationship", "race", "sex",
  "capital_gain", "capital_loss", "hours_per_week",
  "native_country", "income"
)

adult_data <- read.table(
  "data-raw/adult.data",
  sep = ",",
  strip.white = TRUE,
  header = FALSE,
  col.names = colnames_adult,
  na.strings = "?",
  stringsAsFactors = FALSE
)

numeric_vars <- c(
  "age",
  "education_num",
  "capital_gain",
  "capital_loss",
  "hours_per_week"
)

adult <- adult_data[, numeric_vars]

# The selected five numerical variables do not contain missing values
# in the original Adult data file used here.
stopifnot(!anyNA(adult))


# 3. Save package data ----------------------------------------------------

usethis::use_data(
  gau, gau_g, adult,
  overwrite = TRUE,
  compress = "xz"
)


