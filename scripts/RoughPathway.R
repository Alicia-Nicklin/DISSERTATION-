# ============================================================
# PART 1: PATHWAYS FOR CES1, CPT1A, ACSL3 AND MGLL
# ============================================================

library(msigdbr)
library(dplyr)

genes <- c("CES1", "CPT1A", "ACSL3", "MGLL")

# Download full collections
c2 <- msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "C2"
)

hallmark <- msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)

# Keep Reactome and KEGG from C2
canonical <- c2 %>%
  filter(gs_subcollection %in% c(
    "CP:REACTOME",
    "CP:KEGG_LEGACY",
    "CP:KEGG_MEDICUS"
  ))

all_sets <- bind_rows(canonical, hallmark)

# Identify gene-pathway memberships
gene_pathways <- all_sets %>%
  filter(gene_symbol %in% genes) %>%
  distinct(gs_collection, gs_subcollection,
           gs_name, gene_symbol)

# Count how many of our four genes are in each pathway
pathway_comparison <- gene_pathways %>%
  group_by(gs_collection, gs_subcollection, gs_name) %>%
  summarise(
    n_target_genes = n_distinct(gene_symbol),
    genes_present = paste(
      sort(unique(gene_symbol)),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(n_target_genes), gs_name)

# Display pathways containing at least two target genes
shared_pathways <- pathway_comparison %>%
  filter(n_target_genes >= 2)

print(shared_pathways, n = Inf)

# Save complete results
dir.create("outputs", showWarnings = FALSE)

write.csv(
  gene_pathways,
  "outputs/four_genes_all_pathways.csv",
  row.names = FALSE
)

write.csv(
  pathway_comparison,
  "outputs/four_genes_pathway_comparison.csv",
  row.names = FALSE
)

write.csv(
  shared_pathways,
  "outputs/four_genes_shared_pathways.csv",
  row.names = FALSE
)

cat("\nNumber of pathways per gene:\n")
print(gene_pathways %>%
        count(gene_symbol, name = "n_pathways"))


# next step 


# ============================================================
# PART 2: WHICH GENES SHARE PATHWAYS WITH CES1?
# ============================================================

ces1_sets <- gene_pathways %>%
  filter(gene_symbol == "CES1") %>%
  select(gs_collection, gs_subcollection, gs_name)

ces1_overlaps <- gene_pathways %>%
  filter(gene_symbol != "CES1") %>%
  semi_join(
    ces1_sets,
    by = c("gs_collection",
           "gs_subcollection",
           "gs_name")
  ) %>%
  arrange(gene_symbol, gs_name)

print(ces1_overlaps, n = Inf)

write.csv(
  ces1_overlaps,
  "outputs/CES1_shared_pathways_with_other_genes.csv",
  row.names = FALSE
)

cat("\nShared pathway counts with CES1:\n")

print(
  data.frame(gene_symbol = c("CPT1A", "ACSL3", "MGLL")) %>%
    left_join(
      ces1_overlaps %>%
        count(gene_symbol, name = "shared_pathways"),
      by = "gene_symbol"
    ) %>%
    mutate(shared_pathways = coalesce(shared_pathways, 0L))
)


# ============================================================
# CES1 EXPRESSION EXPLORATION
# Week 2 - Colorectal cohort (COAD/READ)
# ============================================================

# Check working directory
getwd()

# Find existing cleaned datasets
files <- list.files(
  path = ".",
  pattern = "\\.(csv|rds|tsv)(\\.gz)?$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

print(files)


# ============================================================
# PART 2: INSPECT CLEANED COLORECTAL DATA
# ============================================================

# Load the previously cleaned and merged dataset
crc_data <- readRDS(
  "outputs/colorectal_primary_expression_clinical_merged.rds"
)

# What type of object is it?
class(crc_data)

# What are its dimensions?
dim(crc_data)

# Look at its structure (without printing everything)
str(crc_data, max.level = 1, list.len = 15)

# Check its variable names
if (is.data.frame(crc_data)) {
  print(names(crc_data))
} else if (is.list(crc_data)) {
  print(names(crc_data))
} else if (is.matrix(crc_data)) {
  print(head(rownames(crc_data), 10))
  print(head(colnames(crc_data), 10))
}


# ============================================================
# PART 3: FIND CES1 IN CLEANED COLORECTAL DATA
# ============================================================

# CES1 Ensembl gene ID
ces1_id <- "ENSG00000198848"

# Find matching column (ignoring version suffix)
ces1_column <- grep(
  paste0("^", ces1_id, "(\\.|$)"),
  names(crc_data),
  value = TRUE
)

cat("\nCES1 matching columns:\n")
print(ces1_column)

# Check number of patients and columns
cat("\nDataset dimensions:\n")
print(dim(crc_data))

# Check whether the cancer-type variable is present
cat("\nCancer type columns:\n")
print(
  grep("type|disease|cohort",
       names(crc_data),
       value = TRUE,
       ignore.case = TRUE)
)


# ============================================================
# PART 4: EXPLORE CES1 EXPRESSION
# ============================================================

library(ggplot2)

# Extract CES1 expression
ces1_expression <- data.frame(
  CES1 = as.numeric(crc_data[[ces1_column]]),
  type = as.character(crc_data$type)
)

# Check cancer-type labels
cat("\nCancer types:\n")
print(table(ces1_expression$type, useNA = "ifany"))

# Check CES1 values
cat("\nCES1 expression summary:\n")
print(summary(ces1_expression$CES1))

cat("\nMissing CES1 values:",
    sum(is.na(ces1_expression$CES1)), "\n")

cat("\nStandard deviation:",
    sd(ces1_expression$CES1, na.rm = TRUE), "\n")

# Plot 1: Distribution of CES1 expression
p1 <- ggplot(ces1_expression, aes(x = CES1)) +
  geom_histogram(
    bins = 30,
    fill = "#4F8C83",
    colour = "white"
  ) +
  labs(
    title = "CES1 expression in colorectal cancer",
    x = "CES1 expression (scale to confirm)",
    y = "Number of patients"
  ) +
  theme_minimal()

print(p1)

# Plot 2: Compare expression by available type labels
p2 <- ggplot(
  ces1_expression,
  aes(x = type, y = CES1)
) +
  geom_boxplot(fill = "#89B6AD") +
  labs(
    title = "CES1 expression by cancer type",
    x = "Cancer type",
    y = "CES1 expression (scale to confirm)"
  ) +
  theme_minimal()

print(p2)

# Save plots
dir.create("figures", showWarnings = FALSE)

ggsave(
  "figures/CES1_expression_distribution.png",
  p1, width = 7, height = 5, dpi = 300
)

ggsave(
  "figures/CES1_expression_by_type.png",
  p2, width = 7, height = 5, dpi = 300
)



