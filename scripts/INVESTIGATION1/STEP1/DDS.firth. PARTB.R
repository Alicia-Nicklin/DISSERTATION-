
# ============================================================
# STAGE 1B — DSS FIRTH-CORRECTED COX REGRESSION
# ============================================================
#
# Patients: 338 with complete DSS data
# DSS events: 37
# Splits: Original 100 DSS splits from Stage 1A
#
# Model 1: age + sex + stage
# Model 2: age + sex + stage + CES1
#
# Expected clinical C-index: 0.7831906
# ============================================================

library(survival)
library(coxphf)


# ============================================================
# 1. LOAD DATA AND ORIGINAL SPLITS
# ============================================================

df <- readRDS(
  "outputs/colorectal_primary_expression_clinical_merged.rds"
)

saved_splits <- readRDS(
  "outputs/Stage1A/split_records.rds"
)

stopifnot(
  all(paste0("DSS_", 1:100) %in% names(saved_splits))
)

cat("Original DSS splits loaded.\n")


# ============================================================
# 2. PREPARE CLINICAL VARIABLES AND CES1
# ============================================================

ces1_gene <- "ENSG00000198848.12"

stopifnot(ces1_gene %in% names(df))

df$age <- as.numeric(
  as.character(df$age_at_initial_pathologic_diagnosis)
)

df$sex <- factor(df$gender)

df$ces1 <- as.numeric(
  as.character(df[[ces1_gene]])
)


# Collapse stage into I, II, III, IV

stage_raw <- toupper(trimws(
  as.character(df$ajcc_pathologic_tumor_stage)
))

stage_raw <- gsub(
  "^STAGE[[:space:]]*",
  "",
  stage_raw
)

stage_group <- rep(NA_character_, length(stage_raw))

stage_group[grepl("^I[A-C]?$", stage_raw)] <- "I"
stage_group[grepl("^II[A-C]?$", stage_raw)] <- "II"
stage_group[grepl("^III[A-C]?$", stage_raw)] <- "III"
stage_group[grepl("^IV[A-C]?$", stage_raw)] <- "IV"

keep_stage <- !is.na(stage_group)

df <- df[keep_stage, , drop = FALSE]

df$stage <- factor(
  stage_group[keep_stage],
  levels = c("I", "II", "III", "IV")
)

stopifnot(nrow(df) == 358)


# ============================================================
# 3. PREPARE DSS DATA — COMPLETE CASES ONLY
# ============================================================

dss_time <- as.numeric(as.character(df$DSS.time))
dss_status <- as.numeric(as.character(df$DSS))

# These are row positions in the original 358-patient
# clinical cohort, before DSS filtering.

dss_keep <- complete.cases(
  dss_time,
  dss_status,
  df$age,
  df$sex,
  df$stage,
  df$ces1
)

original_rows <- which(dss_keep)

dat <- data.frame(
  time = dss_time[dss_keep],
  status = dss_status[dss_keep],
  age = df$age[dss_keep],
  sex = df$sex[dss_keep],
  stage = df$stage[dss_keep],
  ces1 = df$ces1[dss_keep]
)

# Verify the original Stage 1A DSS cohort

stopifnot(
  nrow(dat) == 338,
  sum(dat$status) == 37,
  all(dat$status %in% c(0, 1)),
  all(dat$time >= 0),
  all(complete.cases(dat))
)

cat("\n====================================\n")
cat("DSS DATA CHECK\n")
cat("====================================\n")

cat("Patients:", nrow(dat), "\n")
cat("DSS events:", sum(dat$status), "\n")
cat("Missing CES1:", sum(is.na(dat$ces1)), "\n")

print(table(dat$stage))


# ============================================================
# 4. FIRTH RISK PREDICTION FUNCTION
# ============================================================

firth_risk <- function(fit, newdata, rhs_formula) {
  
  X <- model.matrix(
    rhs_formula,
    data = newdata
  )
  
  # Remove intercept
  X <- X[, colnames(X) != "(Intercept)",
         drop = FALSE]
  
  beta <- fit$coefficients
  
  if (
    any(!is.finite(beta)) ||
    !setequal(colnames(X), names(beta))
  ) {
    stop("Predictor matrix and coefficients do not match.")
  }
  
  # Match coefficient ordering
  X <- X[, names(beta), drop = FALSE]
  
  as.numeric(X %*% beta)
}


# ============================================================
# 5. SET UP RESULTS
# ============================================================

n_splits <- 100

results <- data.frame(
  split = 1:n_splits,
  clinical_c = NA_real_,
  ces1_c = NA_real_,
  difference = NA_real_,
  valid = FALSE
)

clinical_rhs <- ~ age_z + sex + stage
ces1_rhs <- ~ age_z + sex + stage + ces1_z


# ============================================================
# 6. RUN ORIGINAL 100 DSS SPLITS
# ============================================================

for (i in 1:n_splits) {
  
  # ----------------------------------------------------------
  # LOAD SAVED SPLIT
  # ----------------------------------------------------------
  
  split <- saved_splits[[paste0("DSS_", i)]]
  
  # The saved indices refer to the 358-patient cohort.
  # Convert them into positions in the 338-patient
  # DSS-complete dataset using match().
  
  train_idx <- match(split$train, original_rows)
  test_idx <- match(split$test, original_rows)
  
  # Check mapping and complete coverage of 338 patients
  
  stopifnot(
    !anyNA(train_idx),
    !anyNA(test_idx),
    !anyDuplicated(train_idx),
    !anyDuplicated(test_idx),
    all(train_idx %in% seq_len(nrow(dat))),
    all(test_idx %in% seq_len(nrow(dat))),
    length(intersect(train_idx, test_idx)) == 0,
    setequal(
      c(train_idx, test_idx),
      seq_len(nrow(dat))
    )
  )
  
  train <- dat[train_idx, , drop = FALSE]
  test <- dat[test_idx, , drop = FALSE]
  
  
  # ==========================================================
  # 7. STANDARDISE AGE AND CES1
  # ==========================================================
  
  age_mean <- mean(train$age)
  age_sd <- sd(train$age)
  
  ces1_mean <- mean(train$ces1)
  ces1_sd <- sd(train$ces1)
  
  if (
    any(!is.finite(c(age_sd, ces1_sd))) ||
    age_sd == 0 ||
    ces1_sd == 0
  ) {
    stop(paste("Invalid scaling on DSS split", i))
  }
  
  train$age_z <- (train$age - age_mean) / age_sd
  test$age_z <- (test$age - age_mean) / age_sd
  
  train$ces1_z <- (train$ces1 - ces1_mean) / ces1_sd
  test$ces1_z <- (test$ces1 - ces1_mean) / ces1_sd
  
  
  # ==========================================================
  # 8. FIT BOTH FIRTH-CORRECTED COX MODELS
  # ==========================================================
  
  scores <- tryCatch({
    
    # Clinical-only model
    
    clinical_fit <- coxphf(
      formula = Surv(time, status) ~
        age_z + sex + stage,
      data = train,
      firth = TRUE,
      pl = FALSE
    )
    
    # Clinical + CES1 model
    
    ces1_fit <- coxphf(
      formula = Surv(time, status) ~
        age_z + sex + stage + ces1_z,
      data = train,
      firth = TRUE,
      pl = FALSE
    )
    
    
    # ========================================================
    # 9. PREDICT HELD-OUT RISK
    # ========================================================
    
    clinical_risk <- firth_risk(
      clinical_fit,
      test,
      clinical_rhs
    )
    
    ces1_risk <- firth_risk(
      ces1_fit,
      test,
      ces1_rhs
    )
    
    
    # ========================================================
    # 10. HARRELL'S C-INDEX
    # ========================================================
    
    clinical_c <- concordance(
      Surv(test$time, test$status) ~ clinical_risk,
      reverse = TRUE,
      timewt = "n"
    )$concordance
    
    ces1_c <- concordance(
      Surv(test$time, test$status) ~ ces1_risk,
      reverse = TRUE,
      timewt = "n"
    )$concordance
    
    c(clinical_c, ces1_c)
    
  }, error = function(e) {
    
    warning(
      paste("DSS split", i, e$message)
    )
    
    c(NA_real_, NA_real_)
    
  })
  
  
  # ==========================================================
  # 11. STORE PAIRED RESULTS
  # ==========================================================
  
  if (all(is.finite(scores))) {
    
    results$clinical_c[i] <- scores[1]
    results$ces1_c[i] <- scores[2]
    
    results$difference[i] <- scores[2] - scores[1]
    
    results$valid[i] <- TRUE
  }
  
}


# ============================================================
# 12. PRINT RESULTS
# ============================================================

cat("\n====================================\n")
cat("STAGE 1B — DSS RESULTS\n")
cat("====================================\n")

cat(
  "Valid splits:",
  sum(results$valid),
  "/100\n"
)

cat(
  "Clinical mean C-index:",
  mean(results$clinical_c, na.rm = TRUE),
  "\n"
)

cat(
  "CES1 mean C-index:",
  mean(results$ces1_c, na.rm = TRUE),
  "\n"
)

cat(
  "Mean paired improvement:",
  mean(results$difference, na.rm = TRUE),
  "\n"
)


# ============================================================
# 13. VERIFY STAGE 1A DSS BASELINE
# ============================================================

expected_baseline <- 0.7831906

observed_baseline <- mean(
  results$clinical_c,
  na.rm = TRUE
)

# Every split must be valid

if (sum(results$valid) != 100) {
  
  stop(
    "Not all 100 DSS splits were valid. ",
    "Do not interpret CES1 results."
  )
  
}

# Must reproduce Stage 1A baseline

if (
  !is.finite(observed_baseline) ||
  abs(observed_baseline - expected_baseline) > 1e-5
) {
  
  stop(
    "DSS baseline mismatch! Expected 0.7831906, got ",
    observed_baseline,
    ". Do not use CES1 results."
  )
  
}

cat("\nDSS BASELINE VALIDATED!\n")


# ============================================================
# 14. SAVE RESULTS ONLY AFTER VALIDATION
# ============================================================

write.csv(
  results,
  "outputs/Stage1B/DSS_CES1_results.csv",
  row.names = FALSE
)

saveRDS(
  results,
  "outputs/Stage1B/DSS_CES1_results.rds"
)

summary_table <- data.frame(
  endpoint = "DSS",
  clinical_mean = mean(results$clinical_c),
  ces1_mean = mean(results$ces1_c),
  mean_improvement = mean(results$difference),
  valid_splits = sum(results$valid)
)

write.csv(
  summary_table,
  "outputs/Stage1B/DSS_CES1_summary.csv",
  row.names = FALSE
)


# ============================================================
# 15. FINAL OUTPUT
# ============================================================

cat("\n====================================\n")
cat("STAGE 1B — DSS COMPLETE\n")
cat("====================================\n")

print(summary_table, row.names = FALSE)
