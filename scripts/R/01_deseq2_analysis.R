# ================================
# 01_DESeq2_core_and_classification.R
# ================================

  library(DESeq2)
  library(dplyr)
  library(writexl)


# ---- Parameters and file paths ----
counts_file <- "data/raw/countData_root_Ath_clean.csv"
coldata_file <- "data/raw/coldata_root_Ath_clean.csv"

output_dir <- "results/01_DESeq2_model"
output_file <- file.path(
  output_dir,
  "deseq2_deg_classification.xlsx"
)

# Intermediate files for downstream scripts
dds_rds_file      <- file.path(output_dir, "dds_deseq2.rds")
res_list_rds_file  <- file.path(output_dir, "deg_results.rds")
gene_sets_rds_file <- file.path(output_dir, "gene_sets_classified.rds")
norm_counts_rds_file <- file.path(output_dir, "norm_counts.rds")

padj_cut <- 0.05
lfc_cut  <- log2(1.5)  # ~0.585

# ---- 1. Read data ----
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
coldata <- read.csv(coldata_file, stringsAsFactors = FALSE)

# Basic safety checks
if (!"Sample" %in% colnames(coldata)) {
  stop("coldata must contain a 'Sample' column.")
}
if (!all(coldata$Sample %in% colnames(counts))) {
  missing_samples <- setdiff(coldata$Sample, colnames(counts))
  stop("These samples are in coldata but missing from counts: ",
       paste(missing_samples, collapse = ", "))
}
if (!all(colnames(counts) %in% coldata$Sample)) {
  extra_samples <- setdiff(colnames(counts), coldata$Sample)
  stop("These samples are in counts but missing from coldata: ",
       paste(extra_samples, collapse = ", "))
}

# Set sample rownames in coldata for DESeq2 consistency
rownames(coldata) <- coldata$Sample

# Factor settings
coldata$Fe <- factor(coldata$Fe, levels = c("WFe", "LFe"))
coldata$pH <- factor(coldata$pH, levels = c("5.5", "7.5"))
coldata$batch <- factor(coldata$batch)
coldata$Group <- factor(
  paste0(coldata$Fe, "_pH", coldata$pH),
  levels = c("WFe_5.5", "LFe_5.5", "WFe_7.5", "LFe_7.5")
)

# Reorder counts to match sample order in coldata
counts <- counts[, coldata$Sample, drop = FALSE]
if (!identical(colnames(counts), coldata$Sample)) {
  stop("Counts column order does not exactly match coldata$Sample after reordering.")
}

# ---- 2. Build DESeq2 dataset with interaction design ----
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ batch + Fe * pH
)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run DESeq once
dds <- DESeq(dds)

# ---- 3. Extract contrasts ----
res_Fe_5.5 <- results(dds, name = "Fe_LFe_vs_WFe")
res_Fe_7.5 <- results(dds, list(c("Fe_LFe_vs_WFe", "FeLFe.pH7.5")))

res_pH_WFe <- results(dds, name = "pH_7.5_vs_5.5")
res_pH_LFe <- results(dds, list(c("pH_7.5_vs_5.5", "FeLFe.pH7.5")))

res_combined <- results(dds, contrast = list(c("Fe_LFe_vs_WFe", "pH_7.5_vs_5.5", "FeLFe.pH7.5")))
res_interaction <- results(dds, name = "FeLFe.pH7.5")

res_list <- list(
  Fe_5.5 = res_Fe_5.5,
  Fe_7.5 = res_Fe_7.5,
  pH_WFe = res_pH_WFe,
  pH_LFe = res_pH_LFe,
  combined = res_combined,
  interaction = res_interaction
)

# ---- 4. Helper function to get significant up/down gene sets ----
get_sig_genes <- function(res) {
  up <- rownames(res)[!is.na(res$padj) & res$padj < padj_cut & res$log2FoldChange >= lfc_cut]
  down <- rownames(res)[!is.na(res$padj) & res$padj < padj_cut & res$log2FoldChange <= -lfc_cut]
  list(up = up, down = down)
}

Fe_5.5_genes <- get_sig_genes(res_Fe_5.5)
Fe_7.5_genes <- get_sig_genes(res_Fe_7.5)
pH_WFe_genes <- get_sig_genes(res_pH_WFe)
pH_LFe_genes <- get_sig_genes(res_pH_LFe)
combined_genes <- get_sig_genes(res_combined)

int_sig_genes <- rownames(res_interaction)[
  !is.na(res_interaction$padj) & res_interaction$padj < padj_cut
]

# ---- 5. Classify genes ----
# Fe_core_up/down: significant and same direction at both pH
Fe_core_up <- intersect(Fe_5.5_genes$up, Fe_7.5_genes$up)
Fe_core_down <- intersect(Fe_5.5_genes$down, Fe_7.5_genes$down)
Fe_core <- union(Fe_core_up, Fe_core_down)

# Fe_flip: significant in both pH but opposite directions
common_Fe_sig <- intersect(
  union(Fe_5.5_genes$up, Fe_5.5_genes$down),
  union(Fe_7.5_genes$up, Fe_7.5_genes$down)
)

Fe_flip <- common_Fe_sig[
  (res_Fe_5.5[common_Fe_sig, "log2FoldChange"] * res_Fe_7.5[common_Fe_sig, "log2FoldChange"]) < 0
]

# Fe_only_up/down: significant in Fe but not in pH union and not core/flip
pH_union_up <- union(pH_WFe_genes$up, pH_LFe_genes$up)
pH_union_down <- union(pH_WFe_genes$down, pH_LFe_genes$down)

Fe_union_up <- union(Fe_5.5_genes$up, Fe_7.5_genes$up)
Fe_union_down <- union(Fe_5.5_genes$down, Fe_7.5_genes$down)

Fe_only_up <- setdiff(Fe_union_up, union(pH_union_up, union(Fe_core_up, Fe_flip)))
Fe_only_down <- setdiff(Fe_union_down, union(pH_union_down, union(Fe_core_down, Fe_flip)))
Fe_only <- union(Fe_only_up, Fe_only_down)

# pH_core_up/down: significant and same direction at both Fe
pH_core_up <- intersect(pH_WFe_genes$up, pH_LFe_genes$up)
pH_core_down <- intersect(pH_WFe_genes$down, pH_LFe_genes$down)
pH_core <- union(pH_core_up, pH_core_down)

# pH_flip: significant in both Fe but opposite directions
common_pH_sig <- intersect(
  union(pH_WFe_genes$up, pH_WFe_genes$down),
  union(pH_LFe_genes$up, pH_LFe_genes$down)
)

pH_flip <- common_pH_sig[
  (res_pH_WFe[common_pH_sig, "log2FoldChange"] * res_pH_LFe[common_pH_sig, "log2FoldChange"]) < 0
]

# pH_only_up/down: significant in pH but not in Fe union and not core/flip
pH_only_up <- setdiff(pH_union_up, union(Fe_union_up, union(pH_core_up, pH_flip)))
pH_only_down <- setdiff(pH_union_down, union(Fe_union_down, union(pH_core_down, pH_flip)))
pH_only <- union(pH_only_up, pH_only_down)

# Crossed single effect
crossed_single_effect_up <- intersect(
  setdiff(Fe_union_up, union(Fe_core_up, Fe_flip)),
  setdiff(pH_union_up, union(pH_core_up, pH_flip))
)
crossed_single_effect_down <- intersect(
  setdiff(Fe_union_down, union(Fe_core_down, Fe_flip)),
  setdiff(pH_union_down, union(pH_core_down, pH_flip))
)
crossed_single_effect <- union(crossed_single_effect_up, crossed_single_effect_down)

single_factor_up <- union(
  union(Fe_5.5_genes$up, Fe_7.5_genes$up),
  union(pH_WFe_genes$up, pH_LFe_genes$up)
)
single_factor_down <- union(
  union(Fe_5.5_genes$down, Fe_7.5_genes$down),
  union(pH_WFe_genes$down, pH_LFe_genes$down)
)

combined_up <- combined_genes$up
combined_down <- combined_genes$down

combined_shared_up <- setdiff(intersect(combined_up, single_factor_up), crossed_single_effect)
combined_shared_down <- setdiff(intersect(combined_down, single_factor_down), crossed_single_effect)

combined_specific_up <- setdiff(combined_up, union(combined_shared_up, crossed_single_effect))
combined_specific_down <- setdiff(combined_down, union(combined_shared_down, crossed_single_effect))

# Interaction_only: significant in interaction contrast only (exclude all other groups)
all_other_genes <- Reduce(union, list(
  Fe_core_up, Fe_core_down, Fe_only_up, Fe_only_down, Fe_flip,
  pH_core_up, pH_core_down, pH_only_up, pH_only_down, pH_flip,
  crossed_single_effect,
  combined_shared_up, combined_shared_down,
  combined_specific_up, combined_specific_down
))

interaction_only <- setdiff(int_sig_genes, all_other_genes)

gene_sets <- list(
  Fe_core_up = Fe_core_up,
  Fe_core_down = Fe_core_down,
  Fe_only_up = Fe_only_up,
  Fe_only_down = Fe_only_down,
  Fe_flip = Fe_flip,
  pH_core_up = pH_core_up,
  pH_core_down = pH_core_down,
  pH_only_up = pH_only_up,
  pH_only_down = pH_only_down,
  pH_flip = pH_flip,
  crossed_single_effect = crossed_single_effect,
  combined_shared_up = combined_shared_up,
  combined_shared_down = combined_shared_down,
  combined_specific_up = combined_specific_up,
  combined_specific_down = combined_specific_down,
  interaction_only = interaction_only
)

# ---- 6. Save results for downstream scripts ----
norm_counts <- counts(dds, normalized = TRUE)

saveRDS(dds, dds_rds_file)
saveRDS(res_list, res_list_rds_file)
saveRDS(gene_sets, gene_sets_rds_file)
saveRDS(norm_counts, norm_counts_rds_file)

# ---- 7. Save DEG classification table ----
sheets <- list(
  Fe_core_up = data.frame(gene = sort(Fe_core_up)),
  Fe_core_down = data.frame(gene = sort(Fe_core_down)),
  Fe_only_up = data.frame(gene = sort(Fe_only_up)),
  Fe_only_down = data.frame(gene = sort(Fe_only_down)),
  Fe_flip = data.frame(gene = sort(Fe_flip)),
  pH_core_up = data.frame(gene = sort(pH_core_up)),
  pH_core_down = data.frame(gene = sort(pH_core_down)),
  pH_only_up = data.frame(gene = sort(pH_only_up)),
  pH_only_down = data.frame(gene = sort(pH_only_down)),
  pH_flip = data.frame(gene = sort(pH_flip)),
  crossed_single_effect = data.frame(gene = sort(crossed_single_effect)),
  Combined_shared_up = data.frame(gene = sort(combined_shared_up)),
  Combined_shared_down = data.frame(gene = sort(combined_shared_down)),
  Combined_specific_up = data.frame(gene = sort(combined_specific_up)),
  Combined_specific_down = data.frame(gene = sort(combined_specific_down)),
  Interaction_only = data.frame(gene = sort(interaction_only))
)

write_xlsx(sheets, path = output_file)

# ---- 8. Print summary ----
cat("Summary of gene categories:\n")
cat("Fe_core_up:", length(Fe_core_up), "\n")
cat("Fe_core_down:", length(Fe_core_down), "\n")
cat("Fe_only_up:", length(Fe_only_up), "\n")
cat("Fe_only_down:", length(Fe_only_down), "\n")
cat("Fe_flip:", length(Fe_flip), "\n")
cat("pH_core_up:", length(pH_core_up), "\n")
cat("pH_core_down:", length(pH_core_down), "\n")
cat("pH_only_up:", length(pH_only_up), "\n")
cat("pH_only_down:", length(pH_only_down), "\n")
cat("pH_flip:", length(pH_flip), "\n")
cat("crossed_single_effect:", length(crossed_single_effect), "\n")
cat("Combined_shared_up:", length(combined_shared_up), "\n")
cat("Combined_shared_down:", length(combined_shared_down), "\n")
cat("Combined_specific_up:", length(combined_specific_up), "\n")
cat("Combined_specific_down:", length(combined_specific_down), "\n")
cat("Interaction_only:", length(interaction_only), "\n")
