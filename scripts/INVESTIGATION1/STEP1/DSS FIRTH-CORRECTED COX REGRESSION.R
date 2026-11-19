# ==================================================
# STAGE 1A — DSS FIRTH-CORRECTED COX REGRESSION
# ==================================================

# Install once if needed
if (!requireNamespace("coxphf", quietly = TRUE)) {
  install.packages("coxphf")
}

library(survival)
library(coxphf)

# Load cleaned data and original validation splits
dat <- readRDS(
  "outputs/colorectal_primary_expression_clinical_merged.rds"
)

split_records <- readRDS(
  "outputs/Stage1A/split_records.rds"
)

# Prepare clinical variables
stage_group <- function(x) {
  x <- toupper(trimws(as.character(x)))
  grp <- rep(NA_character_, length(x))
  
  grp[grepl("^(STAGE )?IV[A-C]?$", x)] <- "IV"
  grp[grepl("^(STAGE )?III[A-C]?$", x)] <- "III"
  grp[grepl("^(STAGE )?II[A-C]?$", x)] <- "II"
  grp[grepl("^(STAGE )?I[A-C]?$", x)] <- "I"
  
  factor(grp, levels = c("I", "II", "III", "IV"))
}

dat$stage <- stage_group(dat$ajcc_pathologic_tumor_stage)
dat$sex <- factor(dat$gender)
dat$age <- as.numeric(dat$age_at_initial_pathologic_diagnosis)

# Select DSS-complete patients
use <- complete.cases(
  dat[, c("age", "sex", "stage", "DSS", "DSS.time")]
)

rows <- which(use)
d <- dat[use, ]

d$time <- as.numeric(as.character(d$DSS.time))
d$status <- as.numeric(as.character(d$DSS))

stopifnot(nrow(d) == 338)
stopifnot(sum(d$status) == 37)

# Fit full-cohort Firth model
d$age_z <- as.numeric(scale(d$age))

fit_full <- coxphf(
  Surv(time, status) ~ age_z + sex + stage,
  data = d
)

print(summary(fit_full))

# ==================================================
# REPEAT EXISTING 100 VALIDATION SPLITS
# ==================================================

results <- data.frame()

for (i in 1:100) {
  
  split <- split_records[[paste0("DSS_", i)]]
  
  # Original split indices refer to the 358-row dataset
  train_idx <- match(split$train, rows)
  test_idx <- match(split$test, rows)
  
  stopifnot(
    !anyNA(train_idx),
    !anyNA(test_idx)
  )
  
  train <- d[train_idx, ]
  hold <- d[test_idx, ]
  
  # Standardise age using training data only
  mu <- mean(train$age)
  s <- sd(train$age)
  
  train$age_z <- (train$age - mu) / s
  hold$age_z <- (hold$age - mu) / s
  
  # Fit Firth-corrected Cox
  fit <- tryCatch(
    coxphf(
      Surv(time, status) ~ age_z + sex + stage,
      data = train,
      pl = FALSE
    ),
    error = function(e) NULL
  )
  
  c_index <- NA_real_
  
  if (!is.null(fit)) {
    
    # Construct the held-out design matrix
    X <- model.matrix(
      ~ age_z + sex + stage,
      data = hold
    )[, -1, drop = FALSE]
    
    beta <- coef(fit)
    
    if (all(is.finite(beta)) &&
        all(names(beta) %in% colnames(X))) {
      
      # Calculate linear risk score
      risk <- as.vector(
        X[, names(beta), drop = FALSE] %*% beta
      )
      
      c_index <- tryCatch(
        concordance(
          Surv(time, status) ~ risk,
          data = hold,
          reverse = TRUE
        )$concordance,
        error = function(e) NA_real_
      )
    }
  }
  
  results <- rbind(
    results,
    data.frame(
      split = i,
      c_index = c_index,
      valid = is.finite(c_index),
      train_events = sum(train$status),
      test_events = sum(hold$status)
    )
  )
}

# ==================================================
# RESULTS
# ==================================================

cat("Valid splits:", sum(results$valid), "/100\n")

cat(
  "Mean held-out DSS C-index:",
  mean(results$c_index, na.rm = TRUE),
  "\n"
)

cat(
  "SD of held-out DSS C-index:",
  sd(results$c_index, na.rm = TRUE),
  "\n"
)

write.csv(
  results,
  "outputs/Stage1A/DSS_firth_cindex.csv",
  row.names = FALSE
)


# does it fix it 
cat("Valid splits:", sum(results$valid), "/100\n")

cat("Mean DSS C-index:",
    mean(results$c_index, na.rm = TRUE), "\n")

cat("Failed splits:", sum(!results$valid), "\n")