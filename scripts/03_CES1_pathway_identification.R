# ============================================================
# CES1 PATHWAY IDENTIFICATION
# ============================================================

# Install packages (only needed once)
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}

library(msigdbr)

# 1. Download human gene sets from MSigDB
pathways <- msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "C2"
)

# 2. Keep Reactome and KEGG pathways
pathways <- pathways[
  pathways$gs_subcollection %in% c(
    "CP:REACTOME",
    "CP:KEGG_LEGACY",
    "CP:KEGG_MEDICUS"
  ),
]

# 3. Find which pathways contain CES1
ces1 <- pathways[pathways$gene_symbol == "CES1", ]

# 4. Display just the pathway names
results <- unique(
  ces1[, c("gs_subcollection", "gs_name")]
)

print(results, row.names = FALSE)

cat("\nTotal pathways containing CES1:",
    nrow(results), "\n")

# 5. Save results
dir.create("outputs", showWarnings = FALSE)

write.csv(
  results,
  "outputs/CES1_pathways_initial.csv",
  row.names = FALSE
)

#PART2

# ============================================================
# COUNT GENES IN EACH CES1 PATHWAY
# ============================================================

# Keep all genes belonging to CES1-containing pathways
ces1_genes <- pathways[
  pathways$gs_name %in% results$gs_name,
]

# Count unique genes in each pathway
pathway_sizes <- aggregate(
  gene_symbol ~ gs_name,
  data = ces1_genes,
  FUN = function(x) length(unique(x))
)

# Rename column
names(pathway_sizes)[2] <- "n_genes"

# Sort from smallest to largest
pathway_sizes <- pathway_sizes[
  order(pathway_sizes$n_genes),
]

# Display results
print(pathway_sizes, row.names = FALSE)

# Save
write.csv(
  pathway_sizes,
  "outputs/CES1_pathway_sizes.csv",
  row.names = FALSE
)

#MAKEACSV


# ============================================================
# PART 3: EXTRACT GENES IN CES1 PATHWAYS
# ============================================================

# Create a table of each pathway and its member genes
ces1_gene_table <- unique(
  ces1_genes[, c(
    "gs_subcollection",
    "gs_name",
    "gene_symbol"
  )]
)

# Sort by pathway and gene
ces1_gene_table <- ces1_gene_table[
  order(
    ces1_gene_table$gs_name,
    ces1_gene_table$gene_symbol
  ),
]

# Show the genes in the smallest pathway
smallest <- pathway_sizes$gs_name[1]

cat("\nGenes in:", smallest, "\n")

print(
  ces1_gene_table$gene_symbol[
    ces1_gene_table$gs_name == smallest
  ]
)

# Count unique genes across ALL pathways
cat("\nTotal unique genes:",
    length(unique(ces1_gene_table$gene_symbol)),
    "\n")

# Save complete gene membership table
write.csv(
  ces1_gene_table,
  "outputs/CES1_pathway_gene_membership.csv",
  row.names = FALSE
)





