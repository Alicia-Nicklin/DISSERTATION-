# ============================================================
# 00_data_audit.R
# TCGA Pan-Cancer dissertation
# Purpose: inspect raw data before cleaning
# ============================================================

library(data.table)
library(dplyr)
library(readxl)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

DATA_DIR <- "../../Data"

files <- list(
  tpm = file.path(DATA_DIR, "tcga_RSEM_gene_tpm.gz"),
  expected = file.path(DATA_DIR, "tcga_gene_expected_count.gz"),
  norm = file.path(DATA_DIR, "tcga_RSEM_Hugo_norm_count.gz"),
  phenotype = file.path(DATA_DIR, "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"),
  survival = file.path(DATA_DIR, "Survival_SupplementalTable_S1_20171025_xena_sp"),
  clinical = file.path(DATA_DIR, "TCGA-CDR-SupplementalTableS1.xlsx")
)

# Check all files exist
data.frame(
  dataset = names(files),
  exists = file.exists(unlist(files))
)

# ------------------------------------------------------------
# 2. Quick file sizes
# ------------------------------------------------------------

file_info <- file.info(unlist(files))

data.frame(
  dataset = names(files),
  size_MB = round(file_info$size / 1024^2, 2)
)

# ------------------------------------------------------------
# 3. Phenotype data audit
# ------------------------------------------------------------

library(data.table)
library(dplyr)

phenotype <- fread(files$phenotype)

# Basic structure
dim(phenotype)
names(phenotype)
head(phenotype)

# Missing values
phenotype_missing <- data.frame(
  variable = names(phenotype),
  n_missing = sapply(phenotype, function(x) sum(is.na(x) | x == "")),
  percent_missing = round(
    sapply(phenotype, function(x) mean(is.na(x) | x == "") * 100),
    2
  )
)

phenotype_missing

# Sample type distribution
table(phenotype$sample_type, useNA = "ifany")

# Primary tumour count
sum(phenotype$sample_type == "Primary Tumor", na.rm = TRUE)

# Cancer type distribution
sort(table(phenotype$`_primary_disease`), decreasing = TRUE)

# Duplicate sample IDs
sum(duplicated(phenotype$sample))

# Create patient-level TCGA ID
phenotype <- phenotype %>%
  mutate(patient_id = substr(sample, 1, 12))

# Count patients with multiple samples
patient_duplicates <- phenotype %>%
  count(patient_id, sort = TRUE) %>%
  filter(n > 1)

nrow(patient_duplicates)
head(patient_duplicates, 20)

# ------------------------------------------------------------
# 4. Duplicate patients among primary tumour samples only
# ------------------------------------------------------------

primary_phenotype <- phenotype %>%
  filter(sample_type == "Primary Tumor")

primary_patient_duplicates <- primary_phenotype %>%
  count(patient_id, sort = TRUE) %>%
  filter(n > 1)

nrow(primary_patient_duplicates)

head(primary_patient_duplicates, 20)

# ------------------------------------------------------------
# 5. COAD phenotype audit
# ------------------------------------------------------------

coad_phenotype <- phenotype %>%
  filter(`_primary_disease` == "colon adenocarcinoma")

dim(coad_phenotype)

table(coad_phenotype$sample_type, useNA = "ifany")

coad_primary <- coad_phenotype %>%
  filter(sample_type == "Primary Tumor")

nrow(coad_primary)

coad_primary_duplicates <- coad_primary %>%
  count(patient_id, sort = TRUE) %>%
  filter(n > 1)

nrow(coad_primary_duplicates)

head(coad_primary_duplicates, 20)

# ------------------------------------------------------------
# 6. Clinical data audit
# ------------------------------------------------------------

library(readxl)

excel_sheets(files$clinical)

clinical <- read_excel(
  files$clinical,
  sheet = "TCGA-CDR"
)

dim(clinical)

names(clinical)

head(clinical)