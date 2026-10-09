library(readr)
library(dplyr)
library(stringr)
library(readxl)

# ============================================================
# COAD PATIENT COUNT AUDIT
# Aim: investigate why some sources give ~245 COAD cases
# while TCGA-CDR gives ~459
# ============================================================


# ------------------------------------------------------------
# 1. DATA FOLDER
# ------------------------------------------------------------

data_dir <- "C:/Users/alici/OneDrive/Documents/University/year 4/Technical Project/Data"

cat("\nFILES IN DATA FOLDER:\n")
print(list.files(data_dir))


# ------------------------------------------------------------
# 2. LOAD TCGA PHENOTYPE DATA
# ------------------------------------------------------------

phenotype <- read_tsv(
  file.path(data_dir, "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"),
  show_col_types = FALSE
)

cat("\nPHENOTYPE DIMENSIONS:\n")
print(dim(phenotype))

cat("\nPHENOTYPE COLUMN NAMES:\n")
print(names(phenotype))


# ------------------------------------------------------------
# 3. LOOK AT DISEASE LABELS
# ------------------------------------------------------------

cat("\nPRIMARY DISEASE VALUES CONTAINING 'COLON':\n")

print(
  sort(
    unique(
      phenotype$`_primary_disease`[
        str_detect(
          phenotype$`_primary_disease`,
          regex("colon", ignore_case = TRUE)
        )
      ]
    )
  )
)


# ------------------------------------------------------------
# 4. SELECT COAD SAMPLES
# ------------------------------------------------------------

coad_all <- phenotype %>%
  filter(`_primary_disease` == "colon adenocarcinoma") %>%
  mutate(
    patient_id = str_sub(sample, 1, 12)
  )

cat("\nALL COAD SAMPLES IN PHENOTYPE DATA:\n")
print(nrow(coad_all))

cat("\nUNIQUE COAD PATIENTS IN PHENOTYPE DATA:\n")
print(n_distinct(coad_all$patient_id))


# ------------------------------------------------------------
# 5. LOOK AT SAMPLE TYPES WITHIN COAD
# ------------------------------------------------------------

cat("\nCOAD SAMPLE TYPES:\n")

print(
  coad_all %>%
    count(sample_type, sort = TRUE)
)


# ------------------------------------------------------------
# 6. KEEP PRIMARY TUMOURS ONLY
# ------------------------------------------------------------

coad_primary <- coad_all %>%
  filter(sample_type == "Primary Tumor")

cat("\nCOAD PRIMARY TUMOUR SAMPLES:\n")
print(nrow(coad_primary))

cat("\nUNIQUE COAD PATIENTS WITH PRIMARY TUMOUR SAMPLES:\n")
print(n_distinct(coad_primary$patient_id))


# ------------------------------------------------------------
# 7. CHECK FOR DUPLICATE PRIMARY TUMOUR SAMPLES
# ------------------------------------------------------------

duplicate_primary <- coad_primary %>%
  count(patient_id, name = "n_primary_samples") %>%
  filter(n_primary_samples > 1) %>%
  arrange(desc(n_primary_samples))

cat("\nNUMBER OF PATIENTS WITH >1 PRIMARY TUMOUR SAMPLE:\n")
print(nrow(duplicate_primary))

cat("\nDUPLICATE PRIMARY TUMOUR PATIENTS:\n")
print(duplicate_primary)


# ------------------------------------------------------------
# 8. MATCH PRIMARY TUMOURS TO TCGA-CDR
# ------------------------------------------------------------

phenotype_ids <- unique(coad_primary$patient_id)
cdr_ids <- unique(cdr_coad$bcr_patient_barcode)

cat("\n--------------------------------------------\n")
cat("COAD MATCHING SUMMARY\n")
cat("--------------------------------------------\n")

cat("\nUnique COAD patients in phenotype data:\n")
print(n_distinct(coad_all$patient_id))

cat("\nUnique COAD patients with primary tumour data:\n")
print(length(phenotype_ids))

cat("\nUnique COAD patients in TCGA-CDR:\n")
print(length(cdr_ids))

cat("\nPatients present in BOTH:\n")
print(length(intersect(phenotype_ids, cdr_ids)))

cat("\nCDR patients without a primary tumour phenotype sample:\n")
print(length(setdiff(cdr_ids, phenotype_ids)))

cat("\nPrimary tumour phenotype patients absent from CDR:\n")
print(length(setdiff(phenotype_ids, cdr_ids)))


# ------------------------------------------------------------
# 5. LOOK AT SAMPLE TYPES WITHIN COAD
# ------------------------------------------------------------

cat("\nCOAD SAMPLE TYPES:\n")

print(
  coad_all %>%
    count(sample_type, sort = TRUE)
)


# ------------------------------------------------------------
# 6. KEEP PRIMARY TUMOURS ONLY
# ------------------------------------------------------------

coad_primary <- coad_all %>%
  filter(sample_type == "Primary Tumor")

cat("\nCOAD PRIMARY TUMOUR SAMPLES:\n")
print(nrow(coad_primary))

cat("\nUNIQUE COAD PATIENTS WITH PRIMARY TUMOUR SAMPLES:\n")
print(n_distinct(coad_primary$patient_id))


# ------------------------------------------------------------
# 7. CHECK FOR DUPLICATE PRIMARY TUMOUR SAMPLES
# ------------------------------------------------------------

duplicate_primary <- coad_primary %>%
  count(patient_id, name = "n_primary_samples") %>%
  filter(n_primary_samples > 1) %>%
  arrange(desc(n_primary_samples))

cat("\nNUMBER OF PATIENTS WITH >1 PRIMARY TUMOUR SAMPLE:\n")
print(nrow(duplicate_primary))

cat("\nDUPLICATE PRIMARY TUMOUR PATIENTS:\n")
print(duplicate_primary)


# ------------------------------------------------------------
# 8. LOAD TCGA CLINICAL DATA RESOURCE (CDR)
# ------------------------------------------------------------

cdr_path <- file.path(
  data_dir,
  "TCGA-CDR-SupplementalTableS1.xlsx"
)

cat("\nCDR EXCEL SHEETS:\n")
print(excel_sheets(cdr_path))

cdr <- read_excel(
  cdr_path,
  sheet = 1
)

cat("\nCDR DIMENSIONS:\n")
print(dim(cdr))

cat("\nCDR COLUMN NAMES:\n")
print(names(cdr))


# ------------------------------------------------------------
# 9. COUNT COAD PATIENTS IN CDR
# ------------------------------------------------------------

if ("type" %in% names(cdr)) {
  
  cdr_coad <- cdr %>%
    filter(type == "COAD")
  
  cat("\nCOAD ROWS IN CDR:\n")
  print(nrow(cdr_coad))
  
  if ("bcr_patient_barcode" %in% names(cdr_coad)) {
    
    cat("\nUNIQUE COAD PATIENTS IN CDR:\n")
    print(n_distinct(cdr_coad$bcr_patient_barcode))
    
  }
  
} else {
  
  cat("\nWARNING: No column called 'type' found in the CDR.\n")
  cat("Look at the column names printed above.\n")
  
}


# ------------------------------------------------------------
# 10. MATCH XENA PRIMARY-TUMOUR PATIENTS TO CDR
# ------------------------------------------------------------

if (
  exists("cdr_coad") &&
  "bcr_patient_barcode" %in% names(cdr_coad)
) {
  
  phenotype_ids <- unique(coad_primary$patient_id)
  
  cdr_ids <- unique(cdr_coad$bcr_patient_barcode)
  
  cat("\n--------------------------------------------\n")
  cat("COAD MATCHING SUMMARY\n")
  cat("--------------------------------------------\n")
  
  cat("\nUnique COAD patients with primary tumour data:\n")
  print(length(phenotype_ids))
  
  cat("\nUnique COAD patients in TCGA-CDR:\n")
  print(length(cdr_ids))
  
  cat("\nPatients present in BOTH datasets:\n")
  print(length(intersect(phenotype_ids, cdr_ids)))
  
  cat("\nCDR COAD patients NOT found among Xena primary tumours:\n")
  print(length(setdiff(cdr_ids, phenotype_ids)))
  
  cat("\nXena primary-tumour patients NOT found in CDR:\n")
  print(length(setdiff(phenotype_ids, cdr_ids)))
  
  
  # ----------------------------------------------------------
  # 11. SHOW THE UNMATCHED PATIENT IDS
  # ----------------------------------------------------------
  
  missing_from_xena <- setdiff(
    cdr_ids,
    phenotype_ids
  )
  
  missing_from_cdr <- setdiff(
    phenotype_ids,
    cdr_ids
  )
  
  cat("\nCDR patients missing from Xena primary-tumour cohort:\n")
  print(missing_from_xena)
  
  cat("\nXena primary-tumour patients missing from CDR:\n")
  print(missing_from_cdr)
  
}


# ------------------------------------------------------------
# 12. OPTIONAL: CHECK SURVIVAL ENDPOINT AVAILABILITY
# ------------------------------------------------------------

if (exists("cdr_coad")) {
  
  cat("\n--------------------------------------------\n")
  cat("SURVIVAL ENDPOINT AVAILABILITY\n")
  cat("--------------------------------------------\n")
  
  if ("OS" %in% names(cdr_coad)) {
    
    cat("\nPatients with OS status available:\n")
    print(sum(!is.na(cdr_coad$OS)))
    
  }
  
  if ("OS.time" %in% names(cdr_coad)) {
    
    cat("\nPatients with OS time available:\n")
    print(sum(!is.na(cdr_coad$OS.time)))
    
  }
  
  if ("PFI" %in% names(cdr_coad)) {
    
    cat("\nPatients with PFI status available:\n")
    print(sum(!is.na(cdr_coad$PFI)))
    
  }
  
  if ("PFI.time" %in% names(cdr_coad)) {
    
    cat("\nPatients with PFI time available:\n")
    print(sum(!is.na(cdr_coad$PFI.time)))
    
  }
  
  if ("DSS" %in% names(cdr_coad)) {
    
    cat("\nPatients with DSS status available:\n")
    print(sum(!is.na(cdr_coad$DSS)))
    
  }
  
  if ("DSS.time" %in% names(cdr_coad)) {
    
    cat("\nPatients with DSS time available:\n")
    print(sum(!is.na(cdr_coad$DSS.time)))
    
  }
  
}


# ------------------------------------------------------------
# 13. FINAL QUICK SUMMARY
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FINAL COAD AUDIT SUMMARY\n")
cat("============================================\n")

cat("\nAll COAD phenotype samples: ",
    nrow(coad_all),
    "\n")

cat("Unique COAD phenotype patients: ",
    n_distinct(coad_all$patient_id),
    "\n")

cat("COAD primary tumour samples: ",
    nrow(coad_primary),
    "\n")

cat("Unique COAD primary tumour patients: ",
    n_distinct(coad_primary$patient_id),
    "\n")

cat("Patients with duplicate primary tumour samples: ",
    nrow(duplicate_primary),
    "\n")

if (exists("cdr_coad")) {
  
  cat("COAD rows in TCGA-CDR: ",
      nrow(cdr_coad),
      "\n")
  
}

if (
  exists("cdr_coad") &&
  "bcr_patient_barcode" %in% names(cdr_coad)
) {
  
  cat("Matched primary-tumour patients + CDR: ",
      length(intersect(phenotype_ids, cdr_ids)),
      "\n")
  
}

cat("============================================\n")