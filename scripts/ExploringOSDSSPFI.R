# ============================================================
# CES1 SURVIVAL ANALYSIS — INITIAL CHECKS
# ============================================================

# 1. Confirm CES1 is present
ces1_col <- grep(
  "^ENSG00000198848(\\.|$)",
  names(crc_data),
  value = TRUE
)

print(ces1_col)

# 2. Find available survival columns
survival_cols <- grep(
  "OS|DSS|PFI|survival|time|event",
  names(crc_data),
  value = TRUE,
  ignore.case = TRUE
)

cat("\nSurvival-related columns:\n")
print(survival_cols)

# 3. Check CES1 expression distribution
if (length(ces1_col) == 1) {
  ces1 <- as.numeric(crc_data[[ces1_col]])
  
  cat("\nCES1 expression summary:\n")
  print(summary(ces1))
  
  cat("\nCES1 quartile thresholds:\n")
  print(quantile(
    ces1,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    na.rm = TRUE
  ))
  
  cat("\nMissing CES1 values:\n")
  print(sum(is.na(ces1)))
}

# 4. Check dimensions
cat("\nDataset dimensions:\n")
print(dim(crc_data))


# ============================================================
# 03 — CES1 KAPLAN–MEIER SURVIVAL ANALYSIS
# Colorectal cancer (COAD + READ)
# ============================================================

library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

# ------------------------------------------------------------
# 1. PREPARE DATA
# ------------------------------------------------------------

ces1_col <- "ENSG00000198848.12"

km_data <- data.frame(
  CES1 = as.numeric(crc_data[[ces1_col]]),
  OS = as.numeric(as.character(crc_data$OS)),
  OS.time = as.numeric(as.character(crc_data$OS.time)),
  DSS = as.numeric(as.character(crc_data$DSS)),
  DSS.time = as.numeric(as.character(crc_data$DSS.time)),
  PFI = as.numeric(as.character(crc_data$PFI)),
  PFI.time = as.numeric(as.character(crc_data$PFI.time))
)

# ------------------------------------------------------------
# 2. CREATE EXPRESSION GROUPS
# ------------------------------------------------------------

median_ces1 <- median(km_data$CES1, na.rm = TRUE)

km_data$median_group <- factor(
  ifelse(km_data$CES1 <= median_ces1, "Low", "High"),
  levels = c("Low", "High")
)

# Rank-based quartiles give approximately equal group sizes
km_data$quartile <- factor(
  paste0("Q", dplyr::ntile(km_data$CES1, 4)),
  levels = paste0("Q", 1:4)
)

km_data$extreme_group <- factor(
  ifelse(
    km_data$quartile == "Q1", "Q1 Low",
    ifelse(km_data$quartile == "Q4", "Q4 High", NA)
  ),
  levels = c("Q1 Low", "Q4 High")
)

cat("\nMedian groups:\n")
print(table(km_data$median_group))

cat("\nQuartiles:\n")
print(table(km_data$quartile))

# ------------------------------------------------------------
# 3. OUTPUT FOLDER
# ------------------------------------------------------------

out_dir <- file.path(
  "outputs", "CES1_Kaplan_Meier"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 4. FUNCTION FOR KM ANALYSIS
# ------------------------------------------------------------

run_km <- function(endpoint, group_col, analysis_name) {
  
  time_col <- paste0(endpoint, ".time")
  
  df <- km_data[, c(time_col, endpoint, group_col)]
  names(df) <- c("time", "event", "group")
  
  df <- df[complete.cases(df), ]
  df <- df[
    df$time >= 0 &
      df$event %in% c(0, 1),
  ]
  
  df$group <- droplevels(df$group)
  
  if (nlevels(df$group) < 2) {
    warning(paste("Insufficient groups:", endpoint, analysis_name))
    return(NULL)
  }
  
  # TCGA CDR survival times are in days
  df$time_years <- df$time / 365.25
  
  fit <- survfit(
    Surv(time_years, event) ~ group,
    data = df
  )
  
  logrank <- survdiff(
    Surv(time_years, event) ~ group,
    data = df
  )
  
  p_value <- pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1,
    lower.tail = FALSE
  )
  
  plot <- ggsurvplot(
    fit,
    data = df,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = FALSE,
    xlab = "Time (years)",
    ylab = "Survival probability",
    title = paste(
      "CES1 —", endpoint, "—", analysis_name
    ),
    legend.title = "CES1 expression",
    ggtheme = theme_bw()
  )
  
  pdf(
    file.path(
      out_dir,
      paste0(endpoint, "_", analysis_name, ".pdf")
    ),
    width = 9,
    height = 7
  )
  
  print(plot)
  dev.off()
  
  result <- data.frame(
    endpoint = endpoint,
    analysis = analysis_name,
    n = nrow(df),
    events = sum(df$event),
    logrank_p = p_value
  )
  
  return(result)
}

# ------------------------------------------------------------
# 5. RUN ALL NINE ANALYSES
# ------------------------------------------------------------

endpoints <- c("OS", "DSS", "PFI")

results <- list()

for (ep in endpoints) {
  
  results[[paste0(ep, "_median")]] <- run_km(
    ep, "median_group", "Median_Split"
  )
  
  results[[paste0(ep, "_quartiles")]] <- run_km(
    ep, "quartile", "Four_Quartiles"
  )
  
  results[[paste0(ep, "_extremes")]] <- run_km(
    ep, "extreme_group", "Extreme_Quartiles"
  )
}

# ------------------------------------------------------------
# 6. SAVE SUMMARY
# ------------------------------------------------------------

summary_table <- bind_rows(results)

print(summary_table)

write.csv(
  summary_table,
  file.path(out_dir, "CES1_KM_summary.csv"),
  row.names = FALSE
)

cat("\nKaplan–Meier analysis complete!\n")


# ============================================================
# CES1 + MGLL — PROOF OF CONCEPT
# Step 1: Check gene availability
# ============================================================

# MGLL Ensembl gene ID
mgll_id <- "ENSG00000074416"

mgll_col <- grep(
  paste0("^", mgll_id, "(\\.|$)"),
  names(crc_data),
  value = TRUE
)

cat("MGLL matching column:\n")
print(mgll_col)

if (length(mgll_col) == 1) {
  
  cat("\nMGLL expression summary:\n")
  print(summary(crc_data[[mgll_col]]))
  
  cat("\nCorrelation with CES1:\n")
  print(cor(
    as.numeric(crc_data[["ENSG00000198848.12"]]),
    as.numeric(crc_data[[mgll_col]]),
    use = "complete.obs"
  ))
  
}


