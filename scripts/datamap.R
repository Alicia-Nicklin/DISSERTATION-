# =====================================================================
# 00_data_audit.R
# TCGA Pan-Cancer Dissertation
#
# PURPOSE:
# Build a simple, evidence-based data map for the raw TCGA files:
#
#   file | source | rows | columns | ID key | units | purpose | issues
#
# The script also performs supporting QC checks so that the data map is
# based on what is actually present in the files rather than assumptions.
#
# IMPORTANT:
# DRY_RUN = TRUE  -> quick test using first 2,000 genes.
#                    Creates a PREVIEW only.
#
# DRY_RUN = FALSE -> full scan of expression matrices.
#                    Creates the FINAL data map.
#
# Raw data are NEVER modified.
# =====================================================================


# =====================================================================
# 0. PACKAGES
# =====================================================================

required_packages <- c(
  "readr",
  "readxl",
  "dplyr",
  "tidyr",
  "stringr",
  "tibble",
  "matrixStats"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(matrixStats)
})


# =====================================================================
# 1. CONFIG
# =====================================================================

# Raw data folder
DATA_DIR <- "../../Data"

# Output folder inside GitHub/R project
OUT_DIR <- "outputs/data_map"

# ---------------------------------------------------------
# START WITH TRUE
# Change to FALSE only after the dry run works.
# ---------------------------------------------------------

Sys.setenv(VROOM_CONNECTION_SIZE = 1048576)

DRY_RUN <- TRUE

# Number of genes to inspect during dry run
DRY_N <- 2000

# Rows processed at once during full scan
CHUNK_SIZE <- 5000


dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# =====================================================================
# 2. FILE PATHS
# =====================================================================

FILES <- c(
  
  tpm =
    "tcga_RSEM_gene_tpm.gz",
  
  expected =
    "tcga_gene_expected_count.gz",
  
  norm =
    "tcga_RSEM_Hugo_norm_count.gz",
  
  phenotype =
    "TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
  
  cdr =
    "TCGA-CDR-SupplementalTableS1.xlsx",
  
  survival =
    "Survival_SupplementalTable_S1_20171025_xena_sp"
)


P <- setNames(
  file.path(DATA_DIR, FILES),
  names(FILES)
)


# Check every file exists

file_check <- tibble(
  dataset = names(P),
  file = basename(P),
  exists = file.exists(P)
)

print(file_check)


if (!all(file_check$exists)) {
  
  stop(
    "One or more files are missing. Check DATA_DIR and filenames."
  )
}


# =====================================================================
# 3. HELPER FUNCTIONS
# =====================================================================

section <- function(title) {
  
  cat(
    "\n",
    paste(rep("=", 70), collapse = ""),
    "\n",
    title,
    "\n",
    paste(rep("=", 70), collapse = ""),
    "\n",
    sep = ""
  )
}


# ---------------------------------------------------------------------
# Count rows in a .gz expression file without loading it into memory.
#
# This means even during DRY_RUN we still know the REAL number of rows.
# ---------------------------------------------------------------------

count_gz_rows <- function(path) {
  
  con <- gzfile(path, open = "rt")
  
  on.exit(close(con))
  
  # remove header
  readLines(con, n = 1)
  
  total <- 0
  
  repeat {
    
    x <- readLines(
      con,
      n = 100000
    )
    
    if (length(x) == 0) {
      break
    }
    
    total <- total + length(x)
  }
  
  total
}


# ---------------------------------------------------------------------
# Missing-value tokens found in TCGA-style clinical data
# ---------------------------------------------------------------------

NA_TOKENS <- c(
  "",
  "NA",
  "N/A",
  "#N/A",
  "n/a",
  "NULL",
  "null"
)


is_missing_token <- function(x) {
  
  is.na(x) |
    x %in% NA_TOKENS |
    (!is.na(x) & str_detect(x, "^\\[.*\\]$"))
}


# =====================================================================
# 4. FILE INVENTORY
# =====================================================================

section("1. FILE INVENTORY")


file_inventory <- tibble(
  
  dataset = names(P),
  
  file = basename(P),
  
  size_MB = round(
    file.size(P) / 1024^2,
    2
  ),
  
  md5 = unname(
    tools::md5sum(P)
  )
)


print(
  file_inventory,
  n = Inf
)


write_csv(
  file_inventory,
  file.path(
    OUT_DIR,
    "file_inventory.csv"
  )
)


# =====================================================================
# 5. EXPRESSION MATRIX SCANNER
# =====================================================================

# This function:
#
# - gets exact rows
# - gets exact columns
# - checks gene IDs
# - checks sample IDs
# - counts NAs
# - counts literal zeros
# - counts negative values
# - checks constant genes
# - records observed range
#
# In DRY_RUN only first 2,000 genes undergo QC.
# Exact dimensions are still calculated from the full file.


scan_expression <- function(
    path,
    dataset_name,
    target_gene = NULL,
    target_prefix = FALSE) {
  
  
  cat(
    "\nReading:",
    basename(path),
    "\n"
  )
  
  
  # ---------------------------------------------------------
  # Header
  # ---------------------------------------------------------
  
  header <- names(
    read_tsv(
      path,
      n_max = 0,
      col_types = cols(.default = col_character()),
      progress = FALSE
    )
  )
  
  
  n_samples <- length(header) - 1
  
  
  # ---------------------------------------------------------
  # Exact number of gene rows
  # ---------------------------------------------------------
  
  cat(
    "Counting exact rows...\n"
  )
  
  
  exact_rows <- count_gz_rows(path)
  
  
  cat(
    "Exact dimensions:",
    exact_rows,
    "rows x",
    length(header),
    "columns\n"
  )
  
  
  # first column = gene ID
  col_spec <- paste0(
    "c",
    strrep(
      "d",
      n_samples
    )
  )
  
  
  # =========================================================
  # DRY RUN
  # =========================================================
  
  if (DRY_RUN) {
    
    
    df <- read_tsv(
      path,
      n_max = DRY_N,
      col_types = col_spec,
      progress = FALSE
    )
    
    
    ids <- df[[1]]
    
    m <- as.matrix(
      df[, -1]
    )
    
    storage.mode(m) <- "double"
    
    
    gene_qc <- tibble(
      
      id = ids,
      
      n_na =
        rowSums(
          is.na(m)
        ),
      
      n_exact_zero =
        rowSums(
          m == 0,
          na.rm = TRUE
        ),
      
      min =
        matrixStats::rowMins(
          m,
          na.rm = TRUE
        ),
      
      max =
        matrixStats::rowMaxs(
          m,
          na.rm = TRUE
        ),
      
      mean =
        rowMeans(
          m,
          na.rm = TRUE
        ),
      
      sd =
        matrixStats::rowSds(
          m,
          na.rm = TRUE
        )
    )
    
    
    sample_qc <- tibble(
      
      sample =
        colnames(m),
      
      n_na =
        colSums(
          is.na(m)
        ),
      
      n_exact_zero =
        colSums(
          m == 0,
          na.rm = TRUE
        ),
      
      n_negative =
        colSums(
          m < 0,
          na.rm = TRUE
        )
    )
    
    
    target_values <- NULL
    
    
    if (!is.null(target_gene)) {
      
      hit <- if (target_prefix) {
        
        startsWith(
          ids,
          target_gene
        )
        
      } else {
        
        ids == target_gene
      }
      
      
      if (any(hit)) {
        
        target_values <- setNames(
          
          m[
            which(hit)[1],
          ],
          
          colnames(m)
        )
      }
    }
    
    
    return(
      
      list(
        
        dataset = dataset_name,
        
        header = header,
        
        exact_rows = exact_rows,
        
        exact_columns = length(header),
        
        n_samples = n_samples,
        
        scanned_rows = nrow(df),
        
        gene_qc = gene_qc,
        
        sample_qc = sample_qc,
        
        global_min = min(
          m,
          na.rm = TRUE
        ),
        
        global_max = max(
          m,
          na.rm = TRUE
        ),
        
        target_values = target_values,
        
        full_scan = FALSE
      )
    )
  }
  
  
  # =========================================================
  # FULL SCAN
  # =========================================================
  
  
  sample_na <-
    numeric(n_samples)
  
  
  sample_zero <-
    numeric(n_samples)
  
  
  sample_negative <-
    numeric(n_samples)
  
  
  # Used only as transformation diagnostics.
  # They are NOT automatically interpreted as the true units.
  
  implied_sum_offset_0001 <-
    numeric(n_samples)
  
  
  implied_sum_offset_1 <-
    numeric(n_samples)
  
  
  gene_blocks <- list()
  
  
  global_min <- Inf
  
  global_max <- -Inf
  
  
  target_values <- NULL
  
  
  callback_function <- function(df, pos) {
    
    
    ids <- df[[1]]
    
    
    m <- as.matrix(
      df[, -1]
    )
    
    
    storage.mode(m) <- "double"
    
    
    # -----------------------------------------------------
    # Global range
    # -----------------------------------------------------
    
    global_min <<-
      min(
        global_min,
        min(
          m,
          na.rm = TRUE
        )
      )
    
    
    global_max <<-
      max(
        global_max,
        max(
          m,
          na.rm = TRUE
        )
      )
    
    
    # -----------------------------------------------------
    # Sample-level QC
    # -----------------------------------------------------
    
    sample_na <<-
      sample_na +
      colSums(
        is.na(m)
      )
    
    
    sample_zero <<-
      sample_zero +
      colSums(
        m == 0,
        na.rm = TRUE
      )
    
    
    sample_negative <<-
      sample_negative +
      colSums(
        m < 0,
        na.rm = TRUE
      )
    
    
    # -----------------------------------------------------
    # Back-transform diagnostics
    # -----------------------------------------------------
    
    implied_sum_offset_0001 <<-
      implied_sum_offset_0001 +
      colSums(
        2^m - 0.001,
        na.rm = TRUE
      )
    
    
    implied_sum_offset_1 <<-
      implied_sum_offset_1 +
      colSums(
        2^m - 1,
        na.rm = TRUE
      )
    
    
    # -----------------------------------------------------
    # Gene-level QC
    # -----------------------------------------------------
    
    gene_blocks[[length(gene_blocks) + 1]] <<- tibble(
      
      id = ids,
      
      n_na = rowSums(
        is.na(m)
      ),
      
      n_exact_zero = rowSums(
        m == 0,
        na.rm = TRUE
      ),
      
      min = matrixStats::rowMins(
        m,
        na.rm = TRUE
      ),
      
      max = matrixStats::rowMaxs(
        m,
        na.rm = TRUE
      ),
      
      mean = rowMeans(
        m,
        na.rm = TRUE
      ),
      
      sd = matrixStats::rowSds(
        m,
        na.rm = TRUE
      )
    )
    
    
    # -----------------------------------------------------
    # Capture CES1 row if requested
    # -----------------------------------------------------
    
    if (!is.null(target_gene)) {
      
      
      hit <- if (target_prefix) {
        
        startsWith(
          ids,
          target_gene
        )
        
      } else {
        
        ids == target_gene
      }
      
      
      if (
        any(hit) &&
        is.null(target_values)
      ) {
        
        
        target_values <<-
          setNames(
            
            m[
              which(hit)[1],
            ],
            
            colnames(m)
          )
      }
    }
  }
  
  
  read_tsv_chunked(
    
    path,
    
    callback =
      SideEffectChunkCallback$new(
        callback_function
      ),
    
    chunk_size =
      CHUNK_SIZE,
    
    col_types =
      col_spec,
    
    progress =
      TRUE
  )
  
  
  gene_qc <-
    bind_rows(
      gene_blocks
    )
  
  
  sample_qc <-
    tibble(
      
      sample =
        header[-1],
      
      n_na =
        sample_na,
      
      n_exact_zero =
        sample_zero,
      
      n_negative =
        sample_negative,
      
      implied_sum_offset_0001 =
        implied_sum_offset_0001,
      
      implied_sum_offset_1 =
        implied_sum_offset_1
    )
  
  
  list(
    
    dataset =
      dataset_name,
    
    header =
      header,
    
    exact_rows =
      exact_rows,
    
    exact_columns =
      length(header),
    
    n_samples =
      n_samples,
    
    scanned_rows =
      nrow(gene_qc),
    
    gene_qc =
      gene_qc,
    
    sample_qc =
      sample_qc,
    
    global_min =
      global_min,
    
    global_max =
      global_max,
    
    target_values =
      target_values,
    
    full_scan =
      TRUE
  )
}


# =====================================================================
# 6. SCAN THE THREE EXPRESSION MATRICES
# =====================================================================

section(
  "2. EXPRESSION MATRICES"
)


# CES1 Ensembl ID
#
# This is used only as a QC cross-check.
# We will verify the ID before relying on it scientifically.

CES1_ENSEMBL <-
  "ENSG00000198848"


SC <- list(
  
  
  tpm =
    scan_expression(
      
      P["tpm"],
      
      dataset_name =
        "TPM",
      
      target_gene =
        CES1_ENSEMBL,
      
      target_prefix =
        TRUE
    ),
  
  
  expected =
    scan_expression(
      
      P["expected"],
      
      dataset_name =
        "Expected count",
      
      target_gene =
        CES1_ENSEMBL,
      
      target_prefix =
        TRUE
    ),
  
  
  norm =
    scan_expression(
      
      P["norm"],
      
      dataset_name =
        "Normalised count",
      
      target_gene =
        "CES1",
      
      target_prefix =
        FALSE
    )
)


# =====================================================================
# 7. EXPRESSION QC SUMMARY
# =====================================================================

section(
  "2A. EXPRESSION QC SUMMARY"
)


for (nm in names(SC)) {
  
  
  s <- SC[[nm]]
  
  
  cat(
    "\n---",
    nm,
    "---\n"
  )
  
  
  cat(
    "Exact rows:",
    s$exact_rows,
    "\n"
  )
  
  
  cat(
    "Exact columns:",
    s$exact_columns,
    "\n"
  )
  
  
  cat(
    "Sample columns:",
    s$n_samples,
    "\n"
  )
  
  
  cat(
    "Rows scanned for QC:",
    s$scanned_rows,
    "\n"
  )
  
  
  cat(
    "Observed range:",
    s$global_min,
    "to",
    s$global_max,
    "\n"
  )
  
  
  cat(
    "NA cells in scanned data:",
    sum(
      s$sample_qc$n_na
    ),
    "\n"
  )
  
  
  cat(
    "Exact zero cells in scanned data:",
    sum(
      s$sample_qc$n_exact_zero
    ),
    "\n"
  )
  
  
  cat(
    "Negative cells in scanned data:",
    sum(
      s$sample_qc$n_negative
    ),
    "\n"
  )
  
  
  cat(
    "Duplicate sample IDs:",
    sum(
      duplicated(
        s$header[-1]
      )
    ),
    "\n"
  )
  
  
  cat(
    "Duplicate gene IDs in scanned data:",
    sum(
      duplicated(
        s$gene_qc$id
      )
    ),
    "\n"
  )
  
  
  cat(
    "Constant genes in scanned data:",
    sum(
      s$gene_qc$sd == 0,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  cat(
    "First three gene IDs:",
    paste(
      head(
        s$gene_qc$id,
        3
      ),
      collapse = " | "
    ),
    "\n"
  )
}


# =====================================================================
# 8. EXPRESSION MATRIX AGREEMENT
# =====================================================================

section(
  "2B. DO THE EXPRESSION MATRICES AGREE?"
)


same_sample_set <- function(a, b) {
  
  setequal(
    a$header[-1],
    b$header[-1]
  )
}


same_sample_order <- function(a, b) {
  
  identical(
    a$header[-1],
    b$header[-1]
  )
}


cat(
  "TPM vs expected counts - same sample set:",
  same_sample_set(
    SC$tpm,
    SC$expected
  ),
  "\n"
)


cat(
  "TPM vs expected counts - same sample order:",
  same_sample_order(
    SC$tpm,
    SC$expected
  ),
  "\n"
)


cat(
  "TPM vs normalised counts - same sample set:",
  same_sample_set(
    SC$tpm,
    SC$norm
  ),
  "\n"
)


cat(
  "TPM vs normalised counts - same sample order:",
  same_sample_order(
    SC$tpm,
    SC$norm
  ),
  "\n"
)


# Ensembl version suffix check

cat(
  "TPM gene IDs with .version suffix in scanned data:",
  sum(
    str_detect(
      SC$tpm$gene_qc$id,
      "\\.\\d+$"
    )
  ),
  "\n"
)


# =====================================================================
# 9. CES1 CROSS-CHECK
# =====================================================================

section(
  "2C. CES1 CROSS-FILE CHECK"
)


if (
  !DRY_RUN &&
  !is.null(SC$tpm$target_values) &&
  !is.null(SC$expected$target_values) &&
  !is.null(SC$norm$target_values)
) {
  
  
  common_samples <-
    Reduce(
      intersect,
      list(
        names(
          SC$tpm$target_values
        ),
        names(
          SC$expected$target_values
        ),
        names(
          SC$norm$target_values
        )
      )
    )
  
  
  cor_tpm_norm <-
    cor(
      SC$tpm$target_values[
        common_samples
      ],
      SC$norm$target_values[
        common_samples
      ],
      use = "complete.obs"
    )
  
  
  cor_tpm_expected <-
    cor(
      SC$tpm$target_values[
        common_samples
      ],
      SC$expected$target_values[
        common_samples
      ],
      use = "complete.obs"
    )
  
  
  cat(
    "Shared samples:",
    length(
      common_samples
    ),
    "\n"
  )
  
  
  cat(
    "CES1 correlation TPM vs normalised:",
    round(
      cor_tpm_norm,
      4
    ),
    "\n"
  )
  
  
  cat(
    "CES1 correlation TPM vs expected count:",
    round(
      cor_tpm_expected,
      4
    ),
    "\n"
  )
  
  
} else {
  
  
  cat(
    "CES1 cross-file correlation is only performed after the FULL scan.\n"
  )
}


# =====================================================================
# 10. UNIT / TRANSFORMATION DIAGNOSTICS
# =====================================================================

section(
  "2D. UNIT / TRANSFORMATION DIAGNOSTICS"
)


infer_tpm_units <- function(s) {
  
  
  if (DRY_RUN) {
    
    return(
      paste(
        "CHECK after full scan.",
        "Dry-run gene subset cannot be used to test a TPM total."
      )
    )
  }
  
  
  total_0001 <-
    median(
      s$sample_qc$implied_sum_offset_0001,
      na.rm = TRUE
    )
  
  
  total_1 <-
    median(
      s$sample_qc$implied_sum_offset_1,
      na.rm = TRUE
    )
  
  
  near_million <- function(x) {
    
    abs(
      x - 1e6
    ) / 1e6 < 0.02
  }
  
  
  if (
    near_million(
      total_0001
    )
  ) {
    
    return(
      paste0(
        "Consistent with log2(TPM + 0.001). ",
        "Median back-transformed sample total = ",
        signif(
          total_0001,
          4
        )
      )
    )
  }
  
  
  if (
    near_million(
      total_1
    )
  ) {
    
    return(
      paste0(
        "Consistent with log2(TPM + 1). ",
        "Median back-transformed sample total = ",
        signif(
          total_1,
          4
        )
      )
    )
  }
  
  
  paste0(
    "TPM-labelled matrix, but neither tested back-transformation ",
    "gave a median total near 1e6. Verify against Xena metadata."
  )
}


TPM_UNIT_RESULT <-
  infer_tpm_units(
    SC$tpm
  )


cat(
  "TPM:",
  TPM_UNIT_RESULT,
  "\n"
)


cat(
  "Expected count:",
  "RSEM expected count matrix according to source/file label;",
  "stored transformation must be confirmed from Xena metadata.",
  "\n"
)


cat(
  "Normalised count:",
  "RSEM normalised count matrix according to source/file label;",
  "stored transformation must be confirmed from Xena metadata.",
  "\n"
)


# =====================================================================
# 11. PHENOTYPE FILE
# =====================================================================

section(
  "3. PHENOTYPE FILE"
)


phenotype <-
  read_tsv(
    
    P["phenotype"],
    
    col_types =
      cols(
        .default =
          col_character()
      ),
    
    na =
      c(
        "",
        "NA"
      ),
    
    progress =
      FALSE
  )


cat(
  "Dimensions:",
  nrow(
    phenotype
  ),
  "rows x",
  ncol(
    phenotype
  ),
  "columns\n"
)


cat(
  "Columns:",
  paste(
    names(
      phenotype
    ),
    collapse = " | "
  ),
  "\n"
)


# Add barcode information

phenotype <-
  phenotype %>%
  
  mutate(
    
    patient_id =
      substr(
        sample,
        1,
        12
      ),
    
    sample_code =
      substr(
        sample,
        14,
        15
      )
  )


cat(
  "Duplicate sample IDs:",
  sum(
    duplicated(
      phenotype$sample
    )
  ),
  "\n"
)


cat(
  "Missing sample_type:",
  sum(
    is.na(
      phenotype$sample_type
    )
  ),
  "\n"
)


cat(
  "\nSample types:\n"
)


print(
  table(
    phenotype$sample_type,
    useNA = "ifany"
  )
)


cat(
  "\nPrimary diseases:\n"
)


print(
  sort(
    table(
      phenotype$`_primary_disease`
    ),
    decreasing = TRUE
  )
)


# ---------------------------------------------------------
# Primary tumour duplicate check
# ---------------------------------------------------------

primary_phenotype <-
  phenotype %>%
  
  filter(
    sample_type ==
      "Primary Tumor"
  )


primary_duplicates <-
  primary_phenotype %>%
  
  count(
    patient_id
  ) %>%
  
  filter(
    n > 1
  )


cat(
  "\nPrimary Tumor samples:",
  nrow(
    primary_phenotype
  ),
  "\n"
)


cat(
  "Patients with >1 Primary Tumor sample:",
  nrow(
    primary_duplicates
  ),
  "\n"
)


# ---------------------------------------------------------
# COAD information
# ---------------------------------------------------------

coad_pheno <-
  phenotype %>%
  
  filter(
    `_primary_disease` ==
      "colon adenocarcinoma"
  )


cat(
  "\nCOAD total phenotype samples:",
  nrow(
    coad_pheno
  ),
  "\n"
)


cat(
  "COAD Primary Tumor samples:",
  sum(
    coad_pheno$sample_type ==
      "Primary Tumor",
    na.rm = TRUE
  ),
  "\n"
)


# =====================================================================
# 12. TCGA-CDR CLINICAL FILE
# =====================================================================

section(
  "4. TCGA-CDR CLINICAL FILE"
)


cdr_sheets <-
  excel_sheets(
    P["cdr"]
  )


cat(
  "Sheets:",
  paste(
    cdr_sheets,
    collapse = " | "
  ),
  "\n"
)


cdr <-
  suppressWarnings(
    suppressMessages(
      
      read_excel(
        
        P["cdr"],
        
        sheet =
          "TCGA-CDR",
        
        col_types =
          "text"
      )
    )
  )


cat(
  "Dimensions:",
  nrow(
    cdr
  ),
  "rows x",
  ncol(
    cdr
  ),
  "columns\n"
)


cat(
  "Duplicate patient IDs:",
  sum(
    duplicated(
      cdr$bcr_patient_barcode
    )
  ),
  "\n"
)


# ---------------------------------------------------------
# Missingness and placeholder audit
# ---------------------------------------------------------

cdr_missingness <-
  bind_rows(
    
    lapply(
      names(
        cdr
      ),
      
      function(column_name) {
        
        
        x <-
          cdr[
            [column_name]
          ]
        
        
        placeholder <-
          !is.na(x) &
          is_missing_token(x)
        
        
        tibble(
          
          variable =
            column_name,
          
          true_NA =
            sum(
              is.na(x)
            ),
          
          placeholders =
            sum(
              placeholder,
              na.rm = TRUE
            ),
          
          percent_missing =
            round(
              
              100 *
                (
                  sum(
                    is.na(x)
                  ) +
                    sum(
                      placeholder,
                      na.rm = TRUE
                    )
                ) /
                length(x),
              
              2
            )
        )
      }
    )
  )


print(
  cdr_missingness %>%
    
    arrange(
      desc(
        percent_missing
      )
    ),
  
  n = Inf
)


write_csv(
  
  cdr_missingness,
  
  file.path(
    OUT_DIR,
    "cdr_missingness.csv"
  )
)


# =====================================================================
# 13. SURVIVAL ENDPOINT CHECKS
# =====================================================================

section(
  "4A. CDR SURVIVAL ENDPOINTS"
)


endpoint_columns <-
  c(
    "OS",
    "OS.time",
    "DSS",
    "DSS.time",
    "DFI",
    "DFI.time",
    "PFI",
    "PFI.time"
  )


endpoint_columns <-
  intersect(
    endpoint_columns,
    names(
      cdr
    )
  )


cdr_numeric <-
  cdr %>%
  
  mutate(
    
    across(
      
      all_of(
        endpoint_columns
      ),
      
      ~ suppressWarnings(
        
        as.numeric(
          
          ifelse(
            is_missing_token(.x),
            NA,
            .x
          )
        )
      )
    )
  )


for (
  endpoint in
  c(
    "OS",
    "DSS",
    "DFI",
    "PFI"
  )
) {
  
  
  time_col <-
    paste0(
      endpoint,
      ".time"
    )
  
  
  if (
    endpoint %in%
    names(
      cdr_numeric
    ) &&
    time_col %in%
    names(
      cdr_numeric
    )
  ) {
    
    
    events <-
      cdr_numeric[
        [endpoint]
      ]
    
    
    times <-
      cdr_numeric[
        [time_col]
      ]
    
    
    cat(
      "\n",
      endpoint,
      "\n",
      sep = ""
    )
    
    
    cat(
      "Known event values:",
      sum(
        !is.na(
          events
        )
      ),
      "\n"
    )
    
    
    cat(
      "Known times:",
      sum(
        !is.na(
          times
        )
      ),
      "\n"
    )
    
    
    cat(
      "Events not coded 0/1:",
      sum(
        !events %in%
          c(
            0,
            1
          ) &
          !is.na(
            events
          )
      ),
      "\n"
    )
    
    
    cat(
      "Time <= 0:",
      sum(
        times <= 0,
        na.rm = TRUE
      ),
      "\n"
    )
    
    
    cat(
      "Maximum time:",
      max(
        times,
        na.rm = TRUE
      ),
      "\n"
    )
  }
}


# =====================================================================
# 14. XENA SURVIVAL FILE
# =====================================================================

section(
  "5. XENA SURVIVAL FILE"
)


survival_xena <-
  tryCatch(
    
    
    read_tsv(
      
      P["survival"],
      
      col_types =
        cols(
          .default =
            col_character()
        ),
      
      progress =
        FALSE
    ),
    
    
    error = function(e) {
      
      
      cat(
        "Could not parse survival file as TSV.\n"
      )
      
      
      cat(
        "Error:",
        conditionMessage(e),
        "\n"
      )
      
      
      NULL
    }
  )


if (
  !is.null(
    survival_xena
  )
) {
  
  
  cat(
    "Dimensions:",
    nrow(
      survival_xena
    ),
    "rows x",
    ncol(
      survival_xena
    ),
    "columns\n"
  )
  
  
  cat(
    "Columns:",
    paste(
      names(
        survival_xena
      ),
      collapse = " | "
    ),
    "\n"
  )
  
  
  cat(
    "\nFirst rows:\n"
  )
  
  
  print(
    head(
      survival_xena
    )
  )
}


# =====================================================================
# 15. EXPRESSION <-> PHENOTYPE MATCH
# =====================================================================

section(
  "6. CROSS-FILE SAMPLE MATCH"
)


expression_samples <-
  tibble(
    
    sample =
      SC$tpm$header[-1]
  ) %>%
  
  mutate(
    
    patient_id =
      substr(
        sample,
        1,
        12
      ),
    
    sample_code =
      substr(
        sample,
        14,
        15
      )
  )


cat(
  "Expression samples:",
  nrow(
    expression_samples
  ),
  "\n"
)


cat(
  "Expression sample IDs found in phenotype:",
  sum(
    expression_samples$sample %in%
      phenotype$sample
  ),
  "of",
  nrow(
    expression_samples
  ),
  "\n"
)


cat(
  "\nExpression sample barcode codes:\n"
)


print(
  table(
    expression_samples$sample_code
  )
)


# =====================================================================
# 16. SUPPORTING QC FILES
# =====================================================================

section(
  "7. SAVE SUPPORTING QC"
)


for (
  nm in names(
    SC
  )
) {
  
  
  write_csv(
    
    SC[[nm]]$gene_qc,
    
    file.path(
      OUT_DIR,
      paste0(
        "gene_qc_",
        nm,
        ifelse(
          DRY_RUN,
          "_DRY_RUN",
          ""
        ),
        ".csv"
      )
    )
  )
  
  
  write_csv(
    
    SC[[nm]]$sample_qc,
    
    file.path(
      OUT_DIR,
      paste0(
        "sample_qc_",
        nm,
        ifelse(
          DRY_RUN,
          "_DRY_RUN",
          ""
        ),
        ".csv"
      )
    )
  )
}


# =====================================================================
# 17. BUILD SIMPLE DATA MAP
# =====================================================================

section(
  "8. SIMPLE DATA MAP"
)


# ---------------------------------------------------------
# Expression issues helper
# ---------------------------------------------------------

expression_issues <- function(s) {
  
  
  if (DRY_RUN) {
    
    
    paste0(
      
      "Large expression matrix. QC currently based on first ",
      
      s$scanned_rows,
      
      " genes only; run full scan before finalising issues."
    )
    
    
  } else {
    
    
    paste(
      
      paste0(
        sum(
          s$sample_qc$n_na
        ),
        " NA cells"
      ),
      
      paste0(
        sum(
          duplicated(
            s$gene_qc$id
          )
        ),
        " duplicated gene IDs"
      ),
      
      paste0(
        sum(
          s$gene_qc$sd == 0,
          na.rm = TRUE
        ),
        " constant genes"
      ),
      
      sep = "; "
    )
  }
}


# ---------------------------------------------------------
# Units
# ---------------------------------------------------------

tpm_units <-
  if (
    DRY_RUN
  ) {
    
    "RSEM gene-level TPM; stored transformation CHECK after full scan"
    
  } else {
    
    TPM_UNIT_RESULT
  }


expected_units <-
  paste(
    "RSEM expected counts;",
    "stored transformation to be confirmed from Xena metadata"
  )


norm_units <-
  paste(
    "RSEM normalised gene counts;",
    "HUGO gene symbols;",
    "stored transformation to be confirmed from Xena metadata"
  )


# ---------------------------------------------------------
# Survival-file ID
# ---------------------------------------------------------

if (
  !is.null(
    survival_xena
  )
) {
  
  
  candidate_survival_ids <-
    intersect(
      
      c(
        "sample",
        "patient",
        "bcr_patient_barcode",
        "sampleID",
        "patient_id"
      ),
      
      names(
        survival_xena
      )
    )
  
  
  survival_id <-
    if (
      length(
        candidate_survival_ids
      ) > 0
    ) {
      
      candidate_survival_ids[1]
      
    } else {
      
      names(
        survival_xena
      )[1]
    }
  
  
  survival_rows <-
    nrow(
      survival_xena
    )
  
  
  survival_cols <-
    ncol(
      survival_xena
    )
  
  
} else {
  
  
  survival_id <-
    "CHECK"
  
  
  survival_rows <-
    NA_integer_
  
  
  survival_cols <-
    NA_integer_
}


# ---------------------------------------------------------
# Phenotype issues
# ---------------------------------------------------------

phenotype_issue_text <-
  paste0(
    
    sum(
      is.na(
        phenotype$sample_type
      )
    ),
    
    " missing sample_type values; ",
    
    "contains primary, normal, metastatic and recurrent samples; ",
    
    nrow(
      primary_duplicates
    ),
    
    " patients have >1 labelled Primary Tumor sample."
  )


# ---------------------------------------------------------
# CDR issues
# ---------------------------------------------------------

cdr_issue_text <-
  paste0(
    
    "Clinical variables have variable missingness; ",
    
    sum(
      cdr_missingness$placeholders > 0
    ),
    
    " columns contain text-style missing-value placeholders; ",
    
    "Excel sheet also contains an index-style first column."
  )


# ---------------------------------------------------------
# Build map
# ---------------------------------------------------------

DATA_MAP <-
  bind_rows(
    
    
    tibble(
      
      file =
        FILES[["tpm"]],
      
      source =
        "UCSC Xena; Toil RNA-seq Recompute; TCGA Pan-Cancer (PANCAN)",
      
      rows =
        SC$tpm$exact_rows,
      
      columns =
        SC$tpm$exact_columns,
      
      `ID key` =
        paste0(
          "Gene ID = ",
          SC$tpm$gene_qc$id[1],
          " (Ensembl-style); sample IDs = TCGA sample barcodes"
        ),
      
      units =
        tpm_units,
      
      purpose =
        paste(
          "Primary gene-expression matrix for gene-level,",
          "pathway and survival/prognostic analyses."
        ),
      
      issues =
        paste(
          expression_issues(
            SC$tpm
          ),
          "Gene identifiers must be mapped to symbols for pathway work."
        )
    ),
    
    
    tibble(
      
      file =
        FILES[["expected"]],
      
      source =
        "UCSC Xena; Toil RNA-seq Recompute; TCGA Pan-Cancer (PANCAN)",
      
      rows =
        SC$expected$exact_rows,
      
      columns =
        SC$expected$exact_columns,
      
      `ID key` =
        paste0(
          "Gene ID = ",
          SC$expected$gene_qc$id[1],
          " (Ensembl-style); sample IDs = TCGA sample barcodes"
        ),
      
      units =
        expected_units,
      
      purpose =
        paste(
          "Alternative expression representation for methods",
          "requiring count-like expression data."
        ),
      
      issues =
        paste(
          expression_issues(
            SC$expected
          ),
          "Transformation and suitability depend on modelling method."
        )
    ),
    
    
    tibble(
      
      file =
        FILES[["norm"]],
      
      source =
        "UCSC Xena; Toil RNA-seq Recompute; TCGA Pan-Cancer (PANCAN)",
      
      rows =
        SC$norm$exact_rows,
      
      columns =
        SC$norm$exact_columns,
      
      `ID key` =
        paste0(
          "Gene ID = ",
          SC$norm$gene_qc$id[1],
          " (HUGO/symbol-style); sample IDs = TCGA sample barcodes"
        ),
      
      units =
        norm_units,
      
      purpose =
        paste(
          "Alternative normalised expression matrix;",
          "useful for gene-symbol-based cross-checking."
        ),
      
      issues =
        paste(
          expression_issues(
            SC$norm
          ),
          "Gene symbols may need checking for duplicates/aliases."
        )
    ),
    
    
    tibble(
      
      file =
        FILES[["phenotype"]],
      
      source =
        "UCSC Xena; TCGA Pan-Cancer phenotype/sample annotation",
      
      rows =
        nrow(
          phenotype
        ),
      
      columns =
        4,
      
      `ID key` =
        paste(
          "sample = TCGA sample barcode;",
          "patient ID = first 12 characters;",
          "sample type encoded in barcode characters 14-15"
        ),
      
      units =
        "Categorical sample metadata",
      
      purpose =
        paste(
          "Identify cancer type and sample type;",
          "used to restrict analysis to primary tumours."
        ),
      
      issues =
        phenotype_issue_text
    ),
    
    
    tibble(
      
      file =
        FILES[["cdr"]],
      
      source =
        "Liu et al. (2018), TCGA Pan-Cancer Clinical Data Resource, Supplementary Table S1",
      
      rows =
        nrow(
          cdr
        ),
      
      columns =
        ncol(
          cdr
        ),
      
      `ID key` =
        "bcr_patient_barcode (12-character TCGA patient barcode)",
      
      units =
        paste(
          "Clinical variables plus survival endpoints;",
          "OS/DSS/DFI/PFI events coded 0/1;",
          "time variables expected in days and checked separately."
        ),
      
      purpose =
        paste(
          "Clinical covariates and survival/progression outcomes",
          "for patient-level prognostic modelling."
        ),
      
      issues =
        cdr_issue_text
    ),
    
    
    tibble(
      
      file =
        FILES[["survival"]],
      
      source =
        "UCSC Xena version of TCGA Pan-Cancer survival supplementary data",
      
      rows =
        survival_rows,
      
      columns =
        survival_cols,
      
      `ID key` =
        survival_id,
      
      units =
        "Survival/clinical variables; exact fields and time units to verify from file/source",
      
      purpose =
        paste(
          "Convenient Xena-formatted survival resource;",
          "used as a cross-check against TCGA-CDR before deciding whether it is needed."
        ),
      
      issues =
        paste(
          "May duplicate information already available in TCGA-CDR;",
          "retain only if it adds useful variables or simplifies matching."
        )
    )
  )


print(
  DATA_MAP,
  n = Inf,
  width = Inf
)


# =====================================================================
# 18. SAVE DATA MAP
# =====================================================================

if (
  DRY_RUN
) {
  
  
  preview_path <-
    file.path(
      OUT_DIR,
      "data_map_DRY_RUN_PREVIEW.csv"
    )
  
  
  write_csv(
    DATA_MAP,
    preview_path
  )
  
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "DRY RUN COMPLETE\n"
  )
  
  cat(
    "Preview saved to:\n",
    preview_path,
    "\n"
  )
  
  cat(
    "\nDO NOT use this as the final data map yet.\n"
  )
  
  cat(
    "If there are no errors, change:\n\n"
  )
  
  cat(
    "DRY_RUN <- TRUE\n"
  )
  
  cat(
    "to:\n"
  )
  
  cat(
    "DRY_RUN <- FALSE\n"
  )
  
  cat(
    "\nand run the whole script again.\n"
  )
  
  cat(
    "============================================================\n"
  )
  
  
} else {
  
  
  final_path <-
    file.path(
      OUT_DIR,
      "data_map.csv"
    )
  
  
  write_csv(
    DATA_MAP,
    final_path
  )
  
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "FULL DATA AUDIT COMPLETE\n"
  )
  
  cat(
    "Final data map saved to:\n",
    final_path,
    "\n"
  )
  
  cat(
    "============================================================\n"
  )
}