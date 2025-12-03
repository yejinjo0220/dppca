# 1. eur_map -------------------------------------------------------------

load("data-raw/popres.RData")
X_eur_map <- read.table("data-raw/eurmap_data.txt")

color <- gsub("[0-9]+", "", popres$x$color)
eur_map <- data.frame(X = X_eur_map, color = color)

# 2. gaussian_groups -----------------------------------------------------

set.seed(123)

n_groups    <- 5
n_per_group <- 1000
d           <- 20
n_total     <- n_groups * n_per_group

means <- list(
  rep(0, d),                                        # Group 1
  c(rep(5, 5),  rep(0, d - 5)),                     # Group 2
  c(rep(-3, 5), rep(0, d - 5)),                     # Group 3
  c(rep(0, 5),  rep(-1, 5), rep(0, d - 10)),        # Group 4
  c(rep(0, 5),  rep(3, 5),  rep(0, d - 10))         # Group 5
)

Sigma_group1 <- diag(d)
Sigma_group2 <- diag(d)
Sigma_group3 <- 5 * diag(d)
Sigma_group4 <- 5 * diag(d)
Sigma_group5 <- 8 * diag(d)

X_mix <- do.call(rbind, lapply(1:5, function(i) {
  Sigma <- switch(i,
                  Sigma_group1,
                  Sigma_group2,
                  Sigma_group3,
                  Sigma_group4,
                  Sigma_group5
  )
  MASS::mvrnorm(n = n_per_group, mu = means[[i]], Sigma = Sigma)
}))

colors <- c("red", "orange", "green", "blue", "purple")
color  <- rep(colors, each = n_per_group)

gaussian_groups <- data.frame(X = X_mix, color = color)

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
