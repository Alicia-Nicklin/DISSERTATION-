# Inspect the first few lines without loading the full files

readLines(gzfile("../../Data/tcga_RSEM_gene_tpm.gz"), n = 5)

readLines(gzfile("../../Data/tcga_RSEM_Hugo_norm_count.gz"), n = 5)

readLines(gzfile("../../Data/tcga_gene_expected_count.gz"), n = 5)

library(data.table)

tpm_preview <- fread("../../Data/tcga_RSEM_gene_tpm.gz", nrows = 5)
hugo_preview <- fread("../../Data/tcga_RSEM_Hugo_norm_count.gz", nrows = 5)
count_preview <- fread("../../Data/tcga_gene_expected_count.gz", nrows = 5)

# TPM
dim(tpm_preview)
names(tpm_preview)[1:10]
tpm_preview[, 1:6]

# Hugo normalised counts
dim(hugo_preview)
names(hugo_preview)[1:10]
hugo_preview[, 1:6]

# Expected counts
dim(count_preview)
names(count_preview)[1:10]
count_preview[, 1:6]

#this is basically the filtering, phentotyope, COAD/READ section, checks for duplicates, save clean sample 

# ============================================================
# 01_filter_primary_tumours.R
#
# Purpose:
#   - Load TCGA phenotype data
#   - Keep COAD + READ
#   - Keep primary tumours only
#   - Remove normal / metastatic / recurrent samples
#   - Create patient IDs
#   - Check for duplicate primary tumour samples
#   - Save clean sample-patient map for later merging
# ============================================================


# ------------------------------------------------------------
# 1. LOAD PACKAGES
# ------------------------------------------------------------

library(readr)
library(dplyr)
library(stringr)


# ------------------------------------------------------------
# 2. FILE PATHS
# ------------------------------------------------------------

data_dir <- "../../Data"

phenotype_path <- file.path(
  data_dir,
  "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"
)

stopifnot(file.exists(phenotype_path))


# ------------------------------------------------------------
# 3. LOAD PHENOTYPE DATA
# ------------------------------------------------------------

phenotype <- read_tsv(
  phenotype_path,
  show_col_types = FALSE
)

cat("\nPhenotype dimensions:\n")
print(dim(phenotype))

cat("\nPhenotype columns:\n")
print(names(phenotype))


# ------------------------------------------------------------
# 4. KEEP COAD + READ SAMPLES
# ------------------------------------------------------------

colorectal_all <- phenotype %>%
  filter(
    `_primary_disease` %in% c(
      "colon adenocarcinoma",
      "rectum adenocarcinoma"
    )
  ) %>%
  mutate(
    patient_id = str_sub(sample, 1, 12)
  )


cat("\n============================================\n")
cat("COLORECTAL COHORT BEFORE FILTERING\n")
cat("============================================\n")

cat("All COAD/READ samples: ",
    nrow(colorectal_all), "\n")

cat("Unique COAD/READ patients: ",
    n_distinct(colorectal_all$patient_id), "\n")


# ------------------------------------------------------------
# 5. INSPECT SAMPLE TYPES
# ------------------------------------------------------------

sample_type_counts <- colorectal_all %>%
  count(sample_type, sort = TRUE)

cat("\nSample types before filtering:\n")
print(sample_type_counts)


# ------------------------------------------------------------
# 6. KEEP PRIMARY TUMOURS ONLY
# ------------------------------------------------------------

colorectal_primary <- colorectal_all %>%
  filter(sample_type == "Primary Tumor")


cat("\n============================================\n")
cat("PRIMARY TUMOUR COHORT\n")
cat("============================================\n")

cat("Primary tumour samples: ",
    nrow(colorectal_primary), "\n")

cat("Unique primary tumour patients: ",
    n_distinct(colorectal_primary$patient_id), "\n")


# ------------------------------------------------------------
# 7. VERIFY NON-PRIMARY SAMPLES ARE REMOVED
# ------------------------------------------------------------

remaining_types <- colorectal_primary %>%
  count(sample_type)

cat("\nSample types remaining:\n")
print(remaining_types)

if (!all(colorectal_primary$sample_type == "Primary Tumor")) {
  stop("ERROR: non-primary samples are still present.")
}

cat("\nPASS: only Primary Tumor samples remain.\n")


# ------------------------------------------------------------
# 8. CHECK FOR DUPLICATE PRIMARY TUMOUR PATIENTS
# ------------------------------------------------------------

duplicate_primary <- colorectal_primary %>%
  count(
    patient_id,
    name = "n_primary_samples"
  ) %>%
  filter(n_primary_samples > 1) %>%
  arrange(desc(n_primary_samples))


cat("\n============================================\n")
cat("DUPLICATE PATIENT CHECK\n")
cat("============================================\n")

cat("Patients with >1 primary tumour sample: ",
    nrow(duplicate_primary), "\n")

print(duplicate_primary)


# ------------------------------------------------------------
# 9. CREATE CLEAN SAMPLE-PATIENT MAP
# ------------------------------------------------------------

primary_sample_map <- colorectal_primary %>%
  select(
    sample,
    patient_id,
    sample_type,
    `_primary_disease`
  )


# ------------------------------------------------------------
# 10. SAVE OUTPUTS
# ------------------------------------------------------------

dir.create(
  "outputs",
  showWarnings = FALSE
)

write_csv(
  sample_type_counts,
  "outputs/colorectal_sample_type_counts.csv"
)

write_csv(
  duplicate_primary,
  "outputs/duplicate_primary_tumour_patients.csv"
)

write_csv(
  primary_sample_map,
  "outputs/colorectal_primary_sample_map.csv"
)


cat("\n============================================\n")
cat("PRIMARY TUMOUR FILTERING COMPLETE\n")
cat("============================================\n")

cat(
  "Saved primary tumour sample map to:\n",
  "outputs/colorectal_primary_sample_map.csv\n"
)