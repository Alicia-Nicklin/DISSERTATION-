# ============================================================
# 02_merge_expression_clinical.R
#
# Purpose:
#   - Load clean primary tumour sample map
#   - Match those samples to TPM expression data
#   - Keep only matching expression columns
#   - Convert expression sample IDs to patient IDs
#   - Load COAD/READ TCGA-CDR clinical/survival data
#   - Match by patient ID
#   - Save merged analysis dataset
# ============================================================


# ------------------------------------------------------------
# 1. LOAD PACKAGES
# ------------------------------------------------------------

library(readr)
library(dplyr)
library(stringr)
library(readxl)
library(data.table)


# ------------------------------------------------------------
# 2. FILE PATHS
# ------------------------------------------------------------

data_dir <- "../../Data"

sample_map_path <- "outputs/colorectal_primary_sample_map.csv"

expression_path <- file.path(
  data_dir,
  "tcga_RSEM_gene_tpm.gz"
)

cdr_path <- file.path(
  data_dir,
  "TCGA-CDR-SupplementalTableS1.xlsx"
)


stopifnot(file.exists(sample_map_path))
stopifnot(file.exists(expression_path))
stopifnot(file.exists(cdr_path))


# ------------------------------------------------------------
# 3. LOAD CLEAN PRIMARY TUMOUR SAMPLE MAP
# ------------------------------------------------------------

primary_sample_map <- read_csv(
  sample_map_path,
  show_col_types = FALSE
)

cat("\n============================================\n")
cat("PRIMARY SAMPLE MAP\n")
cat("============================================\n")

cat("Rows in sample map: ",
    nrow(primary_sample_map), "\n")

cat("Unique patient IDs: ",
    n_distinct(primary_sample_map$patient_id), "\n")

cat("Unique sample IDs: ",
    n_distinct(primary_sample_map$sample), "\n")


# ------------------------------------------------------------
# 4. READ EXPRESSION FILE HEADER ONLY
# ------------------------------------------------------------

expression_header <- names(
  fread(
    expression_path,
    nrows = 0
  )
)

# First column is the gene ID column
expression_sample_ids <- expression_header[-1]


cat("\n============================================\n")
cat("EXPRESSION HEADER\n")
cat("============================================\n")

cat("Total expression sample columns: ",
    length(expression_sample_ids), "\n")


# ------------------------------------------------------------
# 5. MATCH PRIMARY TUMOUR SAMPLES TO TPM
# ------------------------------------------------------------

primary_samples_in_tpm <- primary_sample_map %>%
  filter(sample %in% expression_sample_ids)

missing_from_tpm <- setdiff(
  primary_sample_map$sample,
  expression_sample_ids
)


cat("\n============================================\n")
cat("PRIMARY TUMOUR -> TPM MATCH\n")
cat("============================================\n")

cat("Primary tumour samples expected: ",
    nrow(primary_sample_map), "\n")

cat("Primary tumour samples found in TPM: ",
    nrow(primary_samples_in_tpm), "\n")

cat("Primary tumour samples missing from TPM: ",
    length(missing_from_tpm), "\n")


if (length(missing_from_tpm) > 0) {
  
  cat("\nMissing sample IDs:\n")
  print(missing_from_tpm)
  
}


# ------------------------------------------------------------
# 6. LOAD ONLY MATCHED TPM COLUMNS
# ------------------------------------------------------------

columns_to_load <- c(
  "sample",
  primary_samples_in_tpm$sample
)

cat("\nLoading matched TPM expression columns...\n")

expression_filtered <- fread(
  expression_path,
  select = columns_to_load
)

cat("\nFiltered expression dimensions:\n")
print(dim(expression_filtered))


# Rename gene identifier column

setnames(
  expression_filtered,
  old = "sample",
  new = "gene_id"
)


# ------------------------------------------------------------
# 7. CHECK EXPRESSION SAMPLE ORDER
# ------------------------------------------------------------

expression_loaded_samples <- names(expression_filtered)[-1]

sample_map_expression <- primary_samples_in_tpm %>%
  slice(
    match(
      expression_loaded_samples,
      sample
    )
  )


stopifnot(
  identical(
    sample_map_expression$sample,
    expression_loaded_samples
  )
)


cat("\nPASS: expression columns match sample map order.\n")


# ------------------------------------------------------------
# 8. CONVERT EXPRESSION SAMPLE IDs TO PATIENT IDs
# ------------------------------------------------------------

expression_patient_ids <- sample_map_expression$patient_id


cat("\n============================================\n")
cat("EXPRESSION PATIENT CHECK\n")
cat("============================================\n")

cat("Expression samples: ",
    length(expression_loaded_samples), "\n")

cat("Expression patients: ",
    length(expression_patient_ids), "\n")

cat("Unique expression patients: ",
    n_distinct(expression_patient_ids), "\n")


# Because script 01 found no duplicate primary tumour patients,
# this should also be zero here.

duplicate_expression_patients <- expression_patient_ids[
  duplicated(expression_patient_ids)
]

cat("Duplicate expression patient IDs: ",
    length(duplicate_expression_patients), "\n")


if (length(duplicate_expression_patients) > 0) {
  
  stop(
    "Unexpected duplicate patient IDs found in expression data."
  )
  
}


# ------------------------------------------------------------
# 9. RENAME EXPRESSION COLUMNS FROM SAMPLE ID -> PATIENT ID
# ------------------------------------------------------------

setnames(
  expression_filtered,
  old = expression_loaded_samples,
  new = expression_patient_ids
)


# ------------------------------------------------------------
# 10. LOAD TCGA-CDR CLINICAL DATA
# ------------------------------------------------------------

cdr <- read_excel(
  cdr_path,
  sheet = "TCGA-CDR"
)


# Remove blank first Excel column if present

if ("...1" %in% names(cdr)) {
  
  cdr <- cdr %>%
    select(-`...1`)
  
}


cat("\n============================================\n")
cat("TCGA-CDR\n")
cat("============================================\n")

cat("CDR dimensions:\n")
print(dim(cdr))


# ------------------------------------------------------------
# 11. KEEP COAD + READ PATIENTS IN CDR
# ------------------------------------------------------------

cdr_colorectal <- cdr %>%
  filter(
    type %in% c("COAD", "READ")
  )


cat("\nCOAD/READ rows in CDR: ",
    nrow(cdr_colorectal), "\n")

cat("Unique COAD/READ patients in CDR: ",
    n_distinct(cdr_colorectal$bcr_patient_barcode), "\n")


# ------------------------------------------------------------
# 12. MATCH EXPRESSION PATIENTS TO CDR
# ------------------------------------------------------------

cdr_patient_ids <- unique(
  cdr_colorectal$bcr_patient_barcode
)

matched_patient_ids <- intersect(
  expression_patient_ids,
  cdr_patient_ids
)

expression_missing_from_cdr <- setdiff(
  expression_patient_ids,
  cdr_patient_ids
)

cdr_missing_from_expression <- setdiff(
  cdr_patient_ids,
  expression_patient_ids
)


cat("\n============================================\n")
cat("EXPRESSION -> CDR MATCH\n")
cat("============================================\n")

cat("Expression patients: ",
    length(expression_patient_ids), "\n")

cat("CDR colorectal patients: ",
    length(cdr_patient_ids), "\n")

cat("Patients present in BOTH: ",
    length(matched_patient_ids), "\n")

cat("Expression patients absent from CDR: ",
    length(expression_missing_from_cdr), "\n")

cat("CDR patients absent from expression cohort: ",
    length(cdr_missing_from_expression), "\n")


if (length(expression_missing_from_cdr) > 0) {
  
  cat("\nExpression patients missing from CDR:\n")
  print(expression_missing_from_cdr)
  
}

if (length(cdr_missing_from_expression) > 0) {
  
  cat("\nCDR patients missing from expression cohort:\n")
  print(cdr_missing_from_expression)
  
}


# ------------------------------------------------------------
# 13. RESTRICT EXPRESSION TO MATCHED PATIENTS
# ------------------------------------------------------------

expression_matched <- expression_filtered[
  ,
  c("gene_id", matched_patient_ids),
  with = FALSE
]


# ------------------------------------------------------------
# 14. RESTRICT CLINICAL DATA TO MATCHED PATIENTS
# ------------------------------------------------------------

clinical_matched <- cdr_colorectal %>%
  filter(
    bcr_patient_barcode %in% matched_patient_ids
  ) %>%
  distinct(
    bcr_patient_barcode,
    .keep_all = TRUE
  )


# Put clinical rows in same patient order as expression columns

clinical_matched <- clinical_matched %>%
  slice(
    match(
      matched_patient_ids,
      bcr_patient_barcode
    )
  )


stopifnot(
  identical(
    clinical_matched$bcr_patient_barcode,
    names(expression_matched)[-1]
  )
)


cat("\nPASS: clinical rows and expression columns are aligned.\n")


# ------------------------------------------------------------
# 15. TRANSPOSE EXPRESSION: PATIENTS x GENES
# ------------------------------------------------------------

expression_matrix <- as.matrix(
  expression_matched[, -1, with = FALSE]
)

rownames(expression_matrix) <- expression_matched$gene_id


expression_patient_df <- as.data.frame(
  t(expression_matrix)
)

expression_patient_df$patient_id <- rownames(
  expression_patient_df
)

expression_patient_df <- expression_patient_df %>%
  select(
    patient_id,
    everything()
  )


# ------------------------------------------------------------
# 16. MERGE CLINICAL + EXPRESSION
# ------------------------------------------------------------

clinical_for_merge <- clinical_matched %>%
  rename(
    patient_id = bcr_patient_barcode
  )


analysis_data <- clinical_for_merge %>%
  inner_join(
    expression_patient_df,
    by = "patient_id"
  )


cat("\n============================================\n")
cat("FINAL MERGED DATASET\n")
cat("============================================\n")

cat("Patients: ",
    nrow(analysis_data), "\n")

cat("Columns: ",
    ncol(analysis_data), "\n")


# ------------------------------------------------------------
# 17. SAVE OUTPUTS
# ------------------------------------------------------------

saveRDS(
  analysis_data,
  "outputs/colorectal_primary_expression_clinical_merged.rds"
)

saveRDS(
  clinical_matched,
  "outputs/colorectal_primary_clinical_matched.rds"
)

write_csv(
  tibble(
    patient_id = matched_patient_ids
  ),
  "outputs/colorectal_matched_patient_ids.csv"
)


cat("\n============================================\n")
cat("MERGE COMPLETE\n")
cat("============================================\n")

cat(
  "Saved merged dataset to:\n",
  "outputs/colorectal_primary_expression_clinical_merged.rds\n"
)


# AHHHHH CHECK WHY WE LOST SO MANY 
# ============================================================
# DIAGNOSTIC: WHY ARE 248 PRIMARY TUMOURS MISSING FROM TPM?
# ============================================================

# Create patient IDs from ALL TPM sample columns
tpm_patient_ids <- str_sub(expression_sample_ids, 1, 12)

# ------------------------------------------------------------
# 1. Do the missing samples at least have the SAME PATIENT
#    somewhere in the TPM dataset?
# ------------------------------------------------------------

missing_sample_patients <- str_sub(
  missing_from_tpm,
  1,
  12
)

patients_found_with_different_sample <- intersect(
  missing_sample_patients,
  tpm_patient_ids
)

cat("\n============================================\n")
cat("MISSING TPM DIAGNOSTIC\n")
cat("============================================\n")

cat(
  "Missing exact primary-tumour sample IDs: ",
  length(missing_from_tpm),
  "\n"
)

cat(
  "Of these, patients appearing elsewhere in TPM: ",
  length(unique(patients_found_with_different_sample)),
  "\n"
)


# ------------------------------------------------------------
# 2. Show examples where patient exists but sample ID differs
# ------------------------------------------------------------

if (length(patients_found_with_different_sample) > 0) {
  
  examples <- head(
    unique(patients_found_with_different_sample),
    20
  )
  
  cat("\nExamples of alternative TPM samples for these patients:\n")
  
  for (id in examples) {
    
    cat("\n", id, ":\n", sep = "")
    
    print(
      expression_sample_ids[
        str_sub(expression_sample_ids, 1, 12) == id
      ]
    )
    
  }
  
}


# ------------------------------------------------------------
# 3. Count missing patients by cancer type
# ------------------------------------------------------------

missing_details <- primary_sample_map %>%
  filter(sample %in% missing_from_tpm) %>%
  count(`_primary_disease`, sort = TRUE)

cat("\nMissing samples by disease:\n")
print(missing_details)

#another check why but just found out we chilling- even tho we lost 200 patientrs lol 
# ============================================================
# CHECK WHETHER MISSING PRIMARY TUMOURS EXIST IN OTHER
# EXPRESSION FILES
# ============================================================

# ------------------------------------------------------------
# 1. READ HEADERS ONLY
# ------------------------------------------------------------

tpm_ids <- names(
  fread(
    "../../Data/tcga_RSEM_gene_tpm.gz",
    nrows = 0
  )
)[-1]

hugo_ids <- names(
  fread(
    "../../Data/tcga_RSEM_Hugo_norm_count.gz",
    nrows = 0
  )
)[-1]

expected_count_ids <- names(
  fread(
    "../../Data/tcga_gene_expected_count.gz",
    nrows = 0
  )
)[-1]


# ------------------------------------------------------------
# 2. CHECK THE 628 PRIMARY TUMOUR SAMPLES AGAINST EACH FILE
# ------------------------------------------------------------

expression_file_check <- tibble(
  
  dataset = c(
    "TPM",
    "Hugo normalised count",
    "Expected count"
  ),
  
  total_expression_samples = c(
    length(tpm_ids),
    length(hugo_ids),
    length(expected_count_ids)
  ),
  
  colorectal_primary_found = c(
    sum(primary_sample_map$sample %in% tpm_ids),
    sum(primary_sample_map$sample %in% hugo_ids),
    sum(primary_sample_map$sample %in% expected_count_ids)
  )
  
) %>%
  
  mutate(
    
    colorectal_primary_missing =
      nrow(primary_sample_map) -
      colorectal_primary_found,
    
    percent_found =
      round(
        100 * colorectal_primary_found /
          nrow(primary_sample_map),
        1
      )
  )


cat("\n============================================\n")
cat("EXPRESSION FILE COMPARISON\n")
cat("============================================\n")

print(expression_file_check)


# ------------------------------------------------------------
# 3. CHECK WHETHER TPM-MISSING SAMPLES ARE PRESENT ELSEWHERE
# ------------------------------------------------------------

missing_tpm_found_hugo <- intersect(
  missing_from_tpm,
  hugo_ids
)

missing_tpm_found_expected <- intersect(
  missing_from_tpm,
  expected_count_ids
)


cat("\nTPM-missing primary tumours found in Hugo file: ",
    length(missing_tpm_found_hugo), "\n")

cat("TPM-missing primary tumours found in expected-count file: ",
    length(missing_tpm_found_expected), "\n")


if (length(missing_tpm_found_hugo) > 0) {
  
  cat("\nExamples found in Hugo but not TPM:\n")
  
  print(
    head(
      missing_tpm_found_hugo,
      20
    )
  )
  
}


if (length(missing_tpm_found_expected) > 0) {
  
  cat("\nExamples found in expected counts but not TPM:\n")
  
  print(
    head(
      missing_tpm_found_expected,
      20
    )
  )
  
}

#finalcheck for OS DSS and PFI 
# ============================================================
# 03_audit_final_cohort.R
#
# Purpose:
#   - Create cohort filtering counts table
#   - Audit clinical variable availability / missingness
#   - Audit OS, DSS and PFI missingness and event counts
#
# Uses outputs created by scripts 01 and 02
# ============================================================


# ------------------------------------------------------------
# 1. LOAD PACKAGES
# ------------------------------------------------------------

library(readr)
library(dplyr)


# ------------------------------------------------------------
# 2. LOAD SAVED OUTPUTS
# ------------------------------------------------------------

primary_sample_map <- read_csv(
  "outputs/colorectal_primary_sample_map.csv",
  show_col_types = FALSE
)

matched_patient_ids <- read_csv(
  "outputs/colorectal_matched_patient_ids.csv",
  show_col_types = FALSE
)

clinical_matched <- readRDS(
  "outputs/colorectal_primary_clinical_matched.rds"
)

analysis_data <- readRDS(
  "outputs/colorectal_primary_expression_clinical_merged.rds"
)


# ------------------------------------------------------------
# 3. BUILD COHORT FILTERING COUNTS TABLE
# ------------------------------------------------------------

# We already established:
# 732 total COAD/READ phenotype samples
# 628 unique COAD/READ patients
# 628 primary tumour samples / patients
# 380 primary tumour samples found in TPM
# 377 patients matched to CDR

counts_table <- tibble(
  
  step = c(
    "All COAD/READ phenotype samples",
    "Unique COAD/READ phenotype patients",
    "Primary tumour samples",
    "Unique primary tumour patients",
    "Primary tumour patients with TPM expression",
    "Patients with TPM expression + TCGA-CDR"
  ),
  
  n = c(
    732,
    628,
    628,
    628,
    380,
    nrow(analysis_data)
  )
  
)


cat("\n============================================\n")
cat("COHORT FILTERING COUNTS\n")
cat("============================================\n")

print(counts_table)


# Save it

write_csv(
  counts_table,
  "outputs/cohort_filter_counts.csv"
)


# ------------------------------------------------------------
# 4. CONFIRM FINAL COHORT SIZE
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FINAL COHORT CHECK\n")
cat("============================================\n")

cat("Final merged patients: ",
    nrow(analysis_data), "\n")

cat("Clinical matched patients: ",
    nrow(clinical_matched), "\n")

cat("Saved matched patient IDs: ",
    nrow(matched_patient_ids), "\n")


# ------------------------------------------------------------
# 5. CLINICAL VARIABLE AUDIT
# ------------------------------------------------------------

clinical_variables <- c(
  
  "age_at_initial_pathologic_diagnosis",
  "gender",
  "race",
  "ajcc_pathologic_tumor_stage",
  "clinical_stage",
  "histological_type",
  "histological_grade",
  "vital_status",
  "tumor_status"
  
)


clinical_variable_audit <- tibble(
  variable = clinical_variables
) %>%
  
  rowwise() %>%
  
  mutate(
    
    available =
      variable %in% names(clinical_matched),
    
    n_total =
      nrow(clinical_matched),
    
    n_non_missing =
      if (available) {
        
        sum(
          !is.na(
            clinical_matched[[variable]]
          )
        )
        
      } else {
        
        NA_integer_
        
      },
    
    n_missing =
      if (available) {
        
        sum(
          is.na(
            clinical_matched[[variable]]
          )
        )
        
      } else {
        
        NA_integer_
        
      },
    
    percent_missing =
      if (available) {
        
        round(
          100 *
            mean(
              is.na(
                clinical_matched[[variable]]
              )
            ),
          1
        )
        
      } else {
        
        NA_real_
        
      }
    
  ) %>%
  
  ungroup()


cat("\n============================================\n")
cat("CLINICAL VARIABLE AUDIT\n")
cat("============================================\n")

print(
  clinical_variable_audit,
  n = Inf
)


write_csv(
  clinical_variable_audit,
  "outputs/clinical_variable_audit.csv"
)


# ------------------------------------------------------------
# 6. SHOW VALUES / CODING FOR KEY CLINICAL VARIABLES
# ------------------------------------------------------------

cat("\n============================================\n")
cat("CLINICAL VARIABLE VALUES\n")
cat("============================================\n")


cat("\nGENDER:\n")
print(
  clinical_matched %>%
    count(gender, sort = TRUE)
)


cat("\nPATHOLOGIC STAGE:\n")
print(
  clinical_matched %>%
    count(
      ajcc_pathologic_tumor_stage,
      sort = TRUE
    )
)


cat("\nCLINICAL STAGE:\n")
print(
  clinical_matched %>%
    count(
      clinical_stage,
      sort = TRUE
    )
)


cat("\nHISTOLOGICAL GRADE:\n")
print(
  clinical_matched %>%
    count(
      histological_grade,
      sort = TRUE
    )
)


cat("\nHISTOLOGICAL TYPE:\n")
print(
  clinical_matched %>%
    count(
      histological_type,
      sort = TRUE
    )
)


# ------------------------------------------------------------
# 7. AGE SUMMARY
# ------------------------------------------------------------

cat("\n============================================\n")
cat("AGE SUMMARY\n")
cat("============================================\n")

print(
  summary(
    clinical_matched$
      age_at_initial_pathologic_diagnosis
  )
)


# ------------------------------------------------------------
# 8. SURVIVAL ENDPOINT AUDIT FUNCTION
# ------------------------------------------------------------

audit_survival_endpoint <- function(data,
                                    endpoint,
                                    time_variable) {
  
  status <- data[[endpoint]]
  time   <- data[[time_variable]]
  
  tibble(
    
    endpoint = endpoint,
    
    status_variable = endpoint,
    
    time_variable = time_variable,
    
    total_patients =
      nrow(data),
    
    status_available =
      sum(
        !is.na(status)
      ),
    
    status_missing =
      sum(
        is.na(status)
      ),
    
    time_available =
      sum(
        !is.na(time)
      ),
    
    time_missing =
      sum(
        is.na(time)
      ),
    
    events =
      sum(
        status == 1,
        na.rm = TRUE
      ),
    
    censored =
      sum(
        status == 0,
        na.rm = TRUE
      ),
    
    complete_cases =
      sum(
        !is.na(status) &
          !is.na(time)
      ),
    
    percent_complete =
      round(
        100 *
          mean(
            !is.na(status) &
              !is.na(time)
          ),
        1
      )
    
  )
  
}


# ------------------------------------------------------------
# 9. AUDIT OS, DSS AND PFI
# ------------------------------------------------------------

survival_endpoint_audit <- bind_rows(
  
  audit_survival_endpoint(
    clinical_matched,
    "OS",
    "OS.time"
  ),
  
  audit_survival_endpoint(
    clinical_matched,
    "DSS",
    "DSS.time"
  ),
  
  audit_survival_endpoint(
    clinical_matched,
    "PFI",
    "PFI.time"
  )
  
)


cat("\n============================================\n")
cat("SURVIVAL ENDPOINT AUDIT\n")
cat("============================================\n")

print(
  survival_endpoint_audit,
  n = Inf
)


write_csv(
  survival_endpoint_audit,
  "outputs/survival_endpoint_audit.csv"
)


# ------------------------------------------------------------
# 10. SHOW RAW STATUS CODING
# ------------------------------------------------------------

cat("\n============================================\n")
cat("SURVIVAL STATUS CODING\n")
cat("============================================\n")


cat("\nOS:\n")
print(
  clinical_matched %>%
    count(OS, sort = TRUE)
)


cat("\nDSS:\n")
print(
  clinical_matched %>%
    count(DSS, sort = TRUE)
)


cat("\nPFI:\n")
print(
  clinical_matched %>%
    count(PFI, sort = TRUE)
)


# ------------------------------------------------------------
# 11. ENDPOINT TIME SUMMARIES
# ------------------------------------------------------------

cat("\n============================================\n")
cat("SURVIVAL TIME SUMMARIES\n")
cat("============================================\n")


cat("\nOS.time:\n")
print(
  summary(
    clinical_matched$OS.time
  )
)


cat("\nDSS.time:\n")
print(
  summary(
    clinical_matched$DSS.time
  )
)


cat("\nPFI.time:\n")
print(
  summary(
    clinical_matched$PFI.time
  )
)


# ------------------------------------------------------------
# 12. FINAL CHECK
# ------------------------------------------------------------

cat("\n============================================\n")
cat("AUDIT COMPLETE\n")
cat("============================================\n")

cat(
  "Final cohort size: ",
  nrow(analysis_data),
  "\n"
)

cat(
  "Saved:\n",
  "- outputs/cohort_filter_counts.csv\n",
  "- outputs/clinical_variable_audit.csv\n",
  "- outputs/survival_endpoint_audit.csv\n"
)

