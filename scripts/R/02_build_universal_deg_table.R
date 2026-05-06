# ================================
# 02_build_universal_deg_table.R
# ================================

  library(dplyr)
  library(tidyr)
  library(readxl)
  library(writexl)
  library(DESeq2)

# ---- Paths ----
input_dir <- "results/01_DESeq2_model"
output_dir <- "results/02_DEG_classification"


dds_rds_file         <- file.path(input_dir, "dds_deseq2.rds")
res_list_rds_file    <- file.path(input_dir, "deg_results.rds")
gene_sets_rds_file   <- file.path(input_dir, "gene_sets_classified.rds")
norm_counts_rds_file <- file.path(input_dir, "norm_counts.rds")

annot_file <- "data/annotation/internal/Ath_TAIR10_gene_annotation_processed.xlsx"

output_file <- file.path(
  output_dir,
  "universal_DEG_table.xlsx"
)

# ---- Load objects from script 1 ----
dds <- readRDS(dds_rds_file)
res_list <- readRDS(res_list_rds_file)
gene_sets <- readRDS(gene_sets_rds_file)
norm_counts <- readRDS(norm_counts_rds_file)

# ---- Load gene annotation ----
annot_df <- read_excel(annot_file)
colnames(annot_df)[1:3] <- c("Gene_ID", "Gene_Symbol", "Gene_Description")
annot_df <- annot_df %>% distinct(Gene_ID, .keep_all = TRUE)

# ---- Helper: annotation table builder ----
annotate_joint <- function(gene_vec, res_list, annot_df, norm_counts, conditions, cond_order) {
  df <- data.frame(Gene_ID = gene_vec, stringsAsFactors = FALSE)
  
  # Basic annotation
  df <- left_join(df, annot_df, by = "Gene_ID")
  
  # DESeq2 statistics
  res_label_map <- list(
    Fe_5.5 = c(
      "Log2FC (LFe vs. WFe) at pH 5.5",
      "P-value (LFe vs. WFe) at pH 5.5"
    ),
    Fe_7.5 = c(
      "Log2FC (LFe vs. WFe) at pH 7.5",
      "P-value (LFe vs. WFe) at pH 7.5"
    ),
    pH_WFe = c(
      "Log2FC (pH 7.5 vs. pH 5.5) at WFe",
      "P-value (pH 7.5 vs. pH 5.5) at WFe"
    ),
    pH_LFe = c(
      "Log2FC (pH 7.5 vs. pH 5.5) at LFe",
      "P-value (pH 7.5 vs. pH 5.5) at LFe"
    ),
    combined = c(
      "Log2FC (LFe_pH 7.5 vs. WFe_pH 5.5)",
      "P-value (LFe_pH 7.5 vs. WFe_pH 5.5)"
    ),
    interaction = c(
      "Log2FC (Interaction: Fe × pH)",
      "P-value (Interaction: Fe × pH)"
    )
  )
  
  for (nm in names(res_list)) {
    r <- as.data.frame(res_list[[nm]])
    idx <- match(gene_vec, rownames(r))
    
    df[[res_label_map[[nm]][1]]] <- r$log2FoldChange[idx]
    df[[res_label_map[[nm]][2]]] <- r$padj[idx]
  }
  
  # Normalized expression summary
  existing <- intersect(rownames(norm_counts), gene_vec)
  norm_sub <- norm_counts[existing, , drop = FALSE]
  
  missing <- setdiff(gene_vec, rownames(norm_sub))
  if (length(missing) > 0) {
    na_mat <- matrix(
      NA_real_,
      nrow = length(missing),
      ncol = ncol(norm_counts),
      dimnames = list(missing, colnames(norm_counts))
    )
    norm_sub <- rbind(norm_sub, na_mat)
  }
  
  norm_sub <- norm_sub[gene_vec, , drop = FALSE]
  
  norm_long <- as.data.frame(t(norm_sub))
  norm_long$Sample <- rownames(norm_long)
  norm_long$Condition <- conditions[norm_long$Sample]
  
  norm_long_gather <- norm_long %>%
    pivot_longer(
      cols = -c(Sample, Condition),
      names_to = "Gene_ID",
      values_to = "Expression"
    )
  
  summary_expr <- norm_long_gather %>%
    group_by(Gene_ID, Condition) %>%
    summarise(
      Mean_Normalized_Count = mean(Expression, na.rm = TRUE),
      SD = sd(Expression, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Condition,
      values_from = c(Mean_Normalized_Count, SD),
      names_sep = "_"
    )
  
  summary_expr <- summary_expr %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))
  
  df <- left_join(df, summary_expr, by = "Gene_ID")
  
  # Column order
  annot_cols <- c("Gene_ID", intersect(c("Gene_Symbol", "Gene_Description"), colnames(df)))
  stat_cols <- unlist(lapply(names(res_list), function(nm) {
    c(
      res_label_map[[nm]][1],
      res_label_map[[nm]][2]
    )
  }))
  stat_cols <- stat_cols[stat_cols %in% colnames(df)]
  
  expr_cols <- unlist(lapply(cond_order, function(co) {
    c(
      paste0("Mean_Normalized_Count_", co),
      paste0("SD_", co)
    )
  }))
  expr_cols <- expr_cols[expr_cols %in% colnames(df)]
  
  df <- df[, c(annot_cols, stat_cols, expr_cols)]
  df
}

# ---- Sample condition labels ----
# Must match the order of columns in norm_counts
conditions <- c(
  "WFe_pH5.5", "WFe_pH5.5", "WFe_pH5.5",
  "LFe_pH5.5", "LFe_pH5.5", "LFe_pH5.5",
  "LFe_pH7.5", "LFe_pH7.5", "LFe_pH7.5",
  "WFe_pH7.5", "WFe_pH7.5", "WFe_pH7.5"
)
names(conditions) <- colnames(norm_counts)
cond_order <- c("WFe_pH5.5", "LFe_pH5.5", "LFe_pH7.5", "WFe_pH7.5")

# ---- Build the universal table ----
all_genes <- rownames(norm_counts)
joint_df <- annotate_joint(
  gene_vec = all_genes,
  res_list = res_list,
  annot_df = annot_df,
  norm_counts = norm_counts,
  conditions = conditions,
  cond_order = cond_order
)

# ---- Add regulation assignment columns ----
joint_df[["Fe-regulated"]] <- NA_character_
joint_df[["pH-regulated"]] <- NA_character_
joint_df[["Combined/Interaction"]] <- NA_character_

joint_df[["Fe-regulated"]][joint_df$Gene_ID %in% gene_sets$Fe_core_up]   <- "Fe_core_up"
joint_df[["Fe-regulated"]][joint_df$Gene_ID %in% gene_sets$Fe_core_down] <- "Fe_core_down"
joint_df[["Fe-regulated"]][joint_df$Gene_ID %in% gene_sets$Fe_only_up]    <- "Fe_only_up"
joint_df[["Fe-regulated"]][joint_df$Gene_ID %in% gene_sets$Fe_only_down]  <- "Fe_only_down"
joint_df[["Fe-regulated"]][joint_df$Gene_ID %in% gene_sets$Fe_flip]       <- "Fe_flip"

joint_df[["pH-regulated"]][joint_df$Gene_ID %in% gene_sets$pH_core_up]   <- "pH_core_up"
joint_df[["pH-regulated"]][joint_df$Gene_ID %in% gene_sets$pH_core_down] <- "pH_core_down"
joint_df[["pH-regulated"]][joint_df$Gene_ID %in% gene_sets$pH_only_up]    <- "pH_only_up"
joint_df[["pH-regulated"]][joint_df$Gene_ID %in% gene_sets$pH_only_down]  <- "pH_only_down"
joint_df[["pH-regulated"]][joint_df$Gene_ID %in% gene_sets$pH_flip]       <- "pH_flip"

joint_df[["Combined/Interaction"]][joint_df$Gene_ID %in% gene_sets$crossed_single_effect]    <- "crossed_single_effect"
joint_df[["Combined/Interaction"]][joint_df$Gene_ID %in% gene_sets$combined_shared_up]       <- "Combined_shared_up"
joint_df[["Combined/Interaction"]][joint_df$Gene_ID %in% gene_sets$combined_shared_down]     <- "Combined_shared_down"
joint_df[["Combined/Interaction"]][joint_df$Gene_ID %in% gene_sets$combined_specific_up]     <- "Combined_specific_up"
joint_df[["Combined/Interaction"]][joint_df$Gene_ID %in% gene_sets$combined_specific_down]   <- "Combined_specific_down"
joint_df[["Combined/Interaction"]][joint_df$Gene_ID %in% gene_sets$interaction_only]         <- "Interaction_only"

# ---- Move the assignment columns close to the gene info ----
joint_df <- joint_df %>%
  relocate(`Fe-regulated`, `pH-regulated`, `Combined/Interaction`, .after = "Gene_Description")


# ---- Write output ----
write_xlsx(list(All_DEGs = joint_df), path = output_file)

message("Joint table written to: ", output_file)