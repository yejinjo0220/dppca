# 1. eur_map -------------------------------------------------------------

load("data-raw/popres.RData")
X_eur_map <- read.table("data-raw/eurmap_data.txt")

color <- gsub("[0-9]+", "", popres$x$color)

eur_map <- X_eur_map
eur_map_g <- data.frame(X_eur_map, color = color)

# 2. gaussian_groups -----------------------------------------------------

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

X_mix <- do.call(rbind, lapply(1:5, function(i) {
  MASS::mvrnorm(n = n_per_group, mu = means[[i]], Sigma = Sigma_list[[i]])
}))

colors <- c("red", "orange", "green", "blue", "purple")
color_gau <- rep(colors, each = n_per_group)

gau <- X_mix
gau_g <- data.frame(X_mix, color = color_gau)

# 3. adult (numeric subset) ----------------------------------------------

colnames_adult <- c(
  "age", "workclass", "fnlwgt", "education", "education_num",
  "marital_status", "occupation", "relationship", "race", "sex",
  "capital_gain", "capital_loss", "hours_per_week",
  "native_country", "income"
)

adult_data <- read.table(
  "data-raw/adult.data",
  sep           = ",",
  strip.white   = TRUE,
  header        = FALSE,
  col.names     = colnames_adult,
  na.strings    = "?",
  stringsAsFactors = FALSE
)

numeric_vars <- c("age", "education_num",
                  "capital_gain", "capital_loss",
                  "hours_per_week")

adult <- adult_data[, numeric_vars]

usethis::use_data(eur_map, eur_map_g,
                  gau, gau_g,
                  adult,
                  overwrite = TRUE)


