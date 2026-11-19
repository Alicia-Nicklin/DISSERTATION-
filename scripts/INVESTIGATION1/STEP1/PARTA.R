library(survival)

dat <- readRDS("outputs/colorectal_primary_expression_clinical_merged.rds")

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
print(table(dat$ajcc_pathologic_tumor_stage, dat$stage, useNA = "ifany"))

endpoints <- list(
  OS  = c("OS",  "OS.time"),
  PFI = c("PFI", "PFI.time"),
  DSS = c("DSS", "DSS.time")
)

set.seed(1)
n_splits <- 100
scores <- data.frame()
split_records <- list()

for (nm in names(endpoints)) {
  status <- endpoints[[nm]][1]
  time <- endpoints[[nm]][2]
  use <- complete.cases(dat[, c("age", "sex", "stage", time, status)])
  rows <- which(use)
  d <- dat[use, ]
  d$time <- as.numeric(as.character(d[[time]]))
  d$status <- as.numeric(as.character(d[[status]]))
  stopifnot(all(is.finite(d$time)), all(d$time >= 0), all(d$status %in% c(0, 1)))
  
  fit <- coxph(Surv(time, status) ~ age_z + sex + stage,
               data = transform(d, age_z = as.numeric(scale(age))))
  print(nm)
  print(summary(fit))
  print(cox.zph(fit))
  
  for (i in seq_len(n_splits)) {
    test <- sample.int(nrow(d), size = floor(0.2 * nrow(d)))
    train <- d[-test, ]
    hold <- d[test, ]
    mu <- mean(train$age)
    s <- sd(train$age)
    train$age_z <- (train$age - mu) / s
    hold$age_z <- (hold$age - mu) / s
    
    warns <- character()
    m <- tryCatch(
      withCallingHandlers(
        coxph(
          Surv(time, status) ~ age_z + sex + stage,
          data = train,
          control = coxph.control(iter.max = 50)
        ),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )
    
    bad_warn <- any(grepl("infinite|did not converge|Ran out of iterations",
                          warns, ignore.case = TRUE))
    ok <- !is.null(m) && all(is.finite(coef(m))) &&
      sum(hold$status) > 0 && !bad_warn
    
    c_index <- NA_real_
    if (ok) {
      risk <- predict(m, newdata = hold, type = "lp")
      c_index <- tryCatch(
        concordance(Surv(time, status) ~ risk,
                    data = hold, reverse = TRUE)$concordance,
        error = function(e) NA_real_
      )
    }
    
    scores <- rbind(scores, data.frame(
      endpoint = nm, split = i, c_index = c_index,
      train_events = sum(train$status),
      test_n = nrow(hold), test_events = sum(hold$status),
      valid = is.finite(c_index),
      warning = paste(warns, collapse = " | ")
    ))
    split_records[[paste(nm, i, sep = "_")]] <- list(
      train = rows[-test], test = rows[test]
    )
  }
}

dir.create("outputs/Stage1A", recursive = TRUE, showWarnings = FALSE)
write.csv(scores, "outputs/Stage1A/clinical_baseline_cindex.csv", row.names = FALSE)
saveRDS(split_records, "outputs/Stage1A/split_records.rds")
print(tapply(scores$c_index, scores$endpoint, mean, na.rm = TRUE))
print(tapply(scores$valid, scores$endpoint, sum))

#this is 


table(
  dat$ajcc_pathologic_tumor_stage,
  dat$stage,
  useNA = "ifany"
)

sum(is.na(dat$stage))

unique(
  dat$ajcc_pathologic_tumor_stage[is.na(dat$stage)]
)