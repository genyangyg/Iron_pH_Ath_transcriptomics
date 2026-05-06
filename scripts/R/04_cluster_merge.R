# ================================
# 04_cluster_merge.R
# ================================

  library(readxl)
  library(openxlsx)
  library(dplyr)


# -----------------------------
# 1. Paths
# -----------------------------
input_dir <- "results/03_functional_annotation"
cluster_dir <- "data/clustering"
output_dir <- "results/04_expression_clusters_mfuzz"

input_file <- file.path(
  input_dir,
  "universal_DEG_table_with_annotation.xlsx"
)

cluster_file <- file.path(
  cluster_dir,
  "clustering_mfuzz.csv"
)

output_file <- file.path(
  output_dir,
  "universal_DEG_table_with_annotation_cluster.xlsx"
)

# -----------------------------
# 2. Read input tables
# -----------------------------
deg_df <- read_excel(input_file)

cluster_df <- read.csv(cluster_file, stringsAsFactors = FALSE)

# Basic checks
if (!"Gene_ID" %in% colnames(deg_df)) {
  stop("Input DEG table must contain a 'Gene_ID' column.")
}
if (!all(c("gene", "cluster", "membership") %in% colnames(cluster_df))) {
  stop("Cluster file must contain columns: gene, cluster, membership.")
}

# -----------------------------
# 3. Prepare cluster annotation
# -----------------------------
cluster_annot <- cluster_df %>%
  dplyr::select(gene, cluster, membership) %>%
  distinct(gene, .keep_all = TRUE)

# -----------------------------
# 4. Merge cluster information
# -----------------------------
merged_df <- deg_df %>%
  left_join(cluster_annot, by = c("Gene_ID" = "gene"))

# Rename columns
colnames(merged_df)[colnames(merged_df) == "cluster"] <- "Cluster_Mfuzz"
colnames(merged_df)[colnames(merged_df) == "membership"] <- "Membership_Mfuzz"

# -----------------------------
# 5. Reorder columns
# -----------------------------
cols <- colnames(merged_df)
cols <- setdiff(cols, c("Cluster_Mfuzz", "Membership_Mfuzz"))

if ("Gene_Description" %in% cols) {
  desc_pos <- which(cols == "Gene_Description")
  cols_new <- c(
    cols[1:desc_pos],
    "Cluster_Mfuzz",
    "Membership_Mfuzz",
    cols[(desc_pos + 1):length(cols)]
  )
} else if ("Gene_Symbol" %in% cols) {
  sym_pos <- which(cols == "Gene_Symbol")
  cols_new <- c(
    cols[1:sym_pos],
    "Cluster_Mfuzz",
    "Membership_Mfuzz",
    cols[(sym_pos + 1):length(cols)]
  )
} else {
  gene_pos <- which(cols == "Gene_ID")
  cols_new <- c(
    cols[1:gene_pos],
    "Cluster_Mfuzz",
    "Membership_Mfuzz",
    cols[(gene_pos + 1):length(cols)]
  )
}

cols_new <- unique(cols_new)
merged_df <- merged_df[, cols_new]

# -----------------------------
# 6. Export
# -----------------------------
write.xlsx(merged_df, output_file, overwrite = TRUE)

cat("Merged file saved to:", output_file, "\n")