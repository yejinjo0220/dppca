## code to prepare `vowel` dataset goes here

usethis::use_data(vowel, overwrite = TRUE)



# ==================================================
# Vowel Recognition Dataset
# Source: UCI Machine Learning Repository
# File: vowel-context.data
# ==================================================

# 1. Read original data
vowel <- read.table(
  "data-raw/vowel-context.data",
  header = FALSE
)

# 2. Assign variable names
names(vowel) <- c(
  "Split",
  "Speaker",
  "Sex",
  paste0("F", 0:9),
  "Class"
)

# 3. Convert categorical variables
vowel$Split <- factor(
  vowel$Split,
  levels = c(0, 1),
  labels = c("Train", "Test")
)

vowel$Speaker <- factor(vowel$Speaker)

vowel$Sex <- factor(
  vowel$Sex,
  levels = c(0, 1),
  labels = c("Male", "Female")
)

vowel$Class <- factor(
  vowel$Class,
  levels = 0:10
)

# 4. Verify the data
stopifnot(
  nrow(vowel) == 990,
  ncol(vowel) == 14,
  !anyNA(vowel)
)

str(vowel)
dim(vowel)
summary(vowel)

# 5. Save as package data
usethis::use_data(
  vowel,
  overwrite = TRUE,
  compress = "xz"
)


