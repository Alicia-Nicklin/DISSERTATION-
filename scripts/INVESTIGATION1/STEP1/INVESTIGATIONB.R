
# ============================================================
# STAGE 1B — CLINICAL MODEL + CES1
# ============================================================
#
# Endpoints: OS and PFI
# Model: Ordinary Cox proportional hazards
# Patients: 358
# Splits: SAME 100 splits as Stage 1A
#
# Compare:
# Model 1: age + sex + stage
# Model 2: age + sex + stage + CES1
#
# Evaluation: Harrell's C-index
# ============================================================


library(survival)


# ============================================================
# 1. LOAD DATA AND ORIGINAL SPLITS
# ============================================================

df <- readRDS(
  "outputs/colorectal_primary_expression_clinical_merged.rds"
)

saved_splits <- readRDS(
  "outputs/Stage1A/split_records.rds"
)

dir.create(
  "outputs/Stage1B",
  recursive = TRUE,
  showWarnings = FALSE
)


# Check original split names

expected_splits <- c(
  paste0("OS_", 1:100),
  paste0("PFI_", 1:100)
)

stopifnot(
  all(expected_splits %in% names(saved_splits))
)

cat("Original Stage 1A splits loaded.\n")


# ============================================================
# 2. PREPARE CLINICAL VARIABLES AND CES1
# ============================================================

ces1_gene <- "ENSG00000198848.12"

stopifnot(ces1_gene %in% names(df))


# Age

df$age <- as.numeric(
  as.character(
    df$age_at_initial_pathologic_diagnosis
  )
)


# Sex

df$sex <- factor(df$gender)


# CES1 expression

df$ces1 <- as.numeric(
  as.character(df[[ces1_gene]])
)


# ------------------------------------------------------------
# Collapse pathological stages into I, II, III and IV
# ------------------------------------------------------------

stage_raw <- toupper(trimws(
  as.character(df$ajcc_pathologic_tumor_stage)
))

stage_raw <- gsub(
  "^STAGE[[:space:]]*",
  "",
  stage_raw
)

stage_group <- rep(
  NA_character_,
  length(stage_raw)
)

stage_group[
  grepl("^I[A-C]?$", stage_raw)
] <- "I"

stage_group[
  grepl("^II[A-C]?$", stage_raw)
] <- "II"

stage_group[
  grepl("^III[A-C]?$", stage_raw)
] <- "III"

stage_group[
  grepl("^IV[A-C]?$", stage_raw)
] <- "IV"


# Remove patients without a valid stage
# Preserve original row ordering

keep <- !is.na(stage_group)

df <- df[keep, , drop = FALSE]

df$stage <- factor(
  stage_group[keep],
  levels = c("I", "II", "III", "IV")
)

# Stage I is the reference

df$stage <- relevel(
  df$stage,
  ref = "I"
)


# ============================================================
# 3. CHECK DATA
# ============================================================

cat("\nPatients:", nrow(df), "\n")

cat("\nStage distribution:\n")
print(table(df$stage))

cat("\nMissing CES1 values:",
    sum(is.na(df$ces1)), "\n")


# Must match Stage 1A cohort

stopifnot(nrow(df) == 358)

stopifnot(
  all(complete.cases(
    df[, c("age", "sex", "stage", "ces1")]
  ))
)


# ============================================================
# 4. DEFINE ENDPOINTS
# ============================================================

endpoints <- list(
  OS  = c("OS.time", "OS"),
  PFI = c("PFI.time", "PFI")
)

n_splits <- 100

all_results <- list()


# ============================================================
# 5. LOOP THROUGH OS AND PFI
# ============================================================

for (endpoint in names(endpoints)) {
  
  cat("\n====================================\n")
  cat("Running:", endpoint, "\n")
  cat("====================================\n")
  
  
  time_col <- endpoints[[endpoint]][1]
  event_col <- endpoints[[endpoint]][2]
  
  
  dat <- data.frame(
    time = as.numeric(as.character(df[[time_col]])),
    status = as.numeric(as.character(df[[event_col]])),
    age = df$age,
    sex = df$sex,
    stage = df$stage,
    ces1 = df$ces1
  )
  
  
  # No further patient removal
  
  stopifnot(
    nrow(dat) == 358,
    all(complete.cases(dat)),
    all(dat$status %in% c(0, 1)),
    all(dat$time >= 0)
  )
  
  
  cat("Patients:", nrow(dat), "\n")
  cat("Events:", sum(dat$status), "\n")
  
  
  # ----------------------------------------------------------
  # 6. SET UP RESULTS
  # ----------------------------------------------------------
  
  results <- data.frame(
    split = 1:n_splits,
    clinical_c = NA_real_,
    ces1_c = NA_real_,
    difference = NA_real_,
    valid = FALSE
  )
  
  
  # ==========================================================
  # 7. RUN THE ORIGINAL 100 SPLITS
  # ==========================================================
  
  for (i in 1:n_splits) {
    
    
    # --------------------------------------------------------
    # LOAD ORIGINAL STAGE 1A SPLIT
    # --------------------------------------------------------
    
    split <- saved_splits[[
      paste0(endpoint, "_", i)
    ]]
    
    train_idx <- split$train
    test_idx <- split$test
    
    
    # Verify split indices
    
    stopifnot(
      is.numeric(train_idx),
      is.numeric(test_idx),
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
    
    
    train <- dat[
      train_idx, ,
      drop = FALSE
    ]
    
    test <- dat[
      test_idx, ,
      drop = FALSE
    ]
    
    
    # ========================================================
    # 8. STANDARDISE AGE AND CES1
    # ========================================================
    #
    # Calculate mean and SD from TRAINING data only.
    # Apply same values to held-out test patients.
    # ========================================================
    
    age_mean <- mean(train$age)
    age_sd <- sd(train$age)
    
    ces1_mean <- mean(train$ces1)
    ces1_sd <- sd(train$ces1)
    
    
    if (
      any(!is.finite(c(age_sd, ces1_sd))) ||
      age_sd == 0 ||
      ces1_sd == 0
    ) {
      warning(
        paste("Invalid scaling on", endpoint, "split", i)
      )
      next
    }
    
    
    # Age
    
    train$age_z <- (
      train$age - age_mean
    ) / age_sd
    
    test$age_z <- (
      test$age - age_mean
    ) / age_sd
    
    
    # CES1
    
    train$ces1_z <- (
      train$ces1 - ces1_mean
    ) / ces1_sd
    
    test$ces1_z <- (
      test$ces1 - ces1_mean
    ) / ces1_sd
    
    
    # ========================================================
    # 9. FIT THE TWO COX MODELS
    # ========================================================
    
    scores <- tryCatch({
      
      # MODEL 1: Clinical baseline
      
      clinical_fit <- coxph(
        Surv(time, status) ~
          age_z + sex + stage,
        data = train,
        ties = "efron"
      )
      
      
      # MODEL 2: Clinical + CES1
      
      ces1_fit <- coxph(
        Surv(time, status) ~
          age_z + sex + stage + ces1_z,
        data = train,
        ties = "efron"
      )
      
      
      # Check coefficient estimates
      
      if (
        any(!is.finite(coef(clinical_fit))) ||
        any(!is.finite(coef(ces1_fit)))
      ) {
        stop("Non-finite Cox coefficients")
      }
      
      
      # ======================================================
      # 10. PREDICT HELD-OUT RISK
      # ======================================================
      
      clinical_risk <- predict(
        clinical_fit,
        newdata = test,
        type = "lp"
      )
      
      ces1_risk <- predict(
        ces1_fit,
        newdata = test,
        type = "lp"
      )
      
      
      # ======================================================
      # 11. HARRELL'S C-INDEX
      # ======================================================
      #
      # Higher Cox linear predictor = greater risk.
      # reverse = TRUE corrects concordance direction.
      # ======================================================
      
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
        paste(
          endpoint,
          "split",
          i,
          e$message
        )
      )
      
      c(NA_real_, NA_real_)
      
    })
    
    
    # ========================================================
    # 12. STORE PAIRED RESULTS
    # ========================================================
    
    if (all(is.finite(scores))) {
      
      results$clinical_c[i] <- scores[1]
      
      results$ces1_c[i] <- scores[2]
      
      results$difference[i] <-
        scores[2] - scores[1]
      
      results$valid[i] <- TRUE
    }
    
  }
  
  
  # ==========================================================
  # 13. SAVE ENDPOINT RESULTS
  # ==========================================================
  
  all_results[[endpoint]] <- results
  
  
  write.csv(
    results,
    paste0(
      "outputs/Stage1B/",
      endpoint,
      "_CES1_results.csv"
    ),
    row.names = FALSE
  )
  
  
  cat("\nRESULTS FOR", endpoint, "\n")
  
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
  
}


# ============================================================
# 14. CREATE FINAL SUMMARY
# ============================================================

summary_table <- data.frame(
  
  endpoint = names(all_results),
  
  clinical_mean = sapply(
    all_results,
    function(x) mean(
      x$clinical_c,
      na.rm = TRUE
    )
  ),
  
  ces1_mean = sapply(
    all_results,
    function(x) mean(
      x$ces1_c,
      na.rm = TRUE
    )
  ),
  
  mean_improvement = sapply(
    all_results,
    function(x) mean(
      x$difference,
      na.rm = TRUE
    )
  ),
  
  valid_splits = sapply(
    all_results,
    function(x) sum(x$valid)
  )
  
)


# ============================================================
# 15. SAVE ALL RESULTS
# ============================================================

saveRDS(
  all_results,
  "outputs/Stage1B/CES1_results.rds"
)

write.csv(
  summary_table,
  "outputs/Stage1B/CES1_summary.csv",
  row.names = FALSE
)


# ============================================================
# 16. CHECK AGAINST STAGE 1A BASELINE
# ============================================================

baseline <- c(
  OS = 0.696,
  PFI = 0.694
)

for (ep in names(baseline)) {
  
  observed <- summary_table$clinical_mean[
    summary_table$endpoint == ep
  ]
  
  if (
    length(observed) == 1 &&
    is.finite(observed) &&
    abs(observed - baseline[ep]) > 0.001
  ) {
    
    warning(
      ep,
      ": Clinical baseline differs from Stage 1A. ",
      "Check cohort ordering and preprocessing."
    )
    
  }
  
}


# ============================================================
# 17. FINAL OUTPUT
# ============================================================

cat("\n====================================\n")
cat("STAGE 1B COMPLETE — OS AND PFI\n")
cat("====================================\n")

print(
  summary_table,
  row.names = FALSE
)


