# ================================
# Complete Subgroup classification and export (A, A+, A+exp, Aexp, B, Bexp)
# ================================

library(dplyr)
library(tidyr)
library(readxl)
library(writexl)

# ---- Input ----
input_dir <- "results/04_expression_clusters_mfuzz"
output_dir <- "results/05_regulatory_subgroup_classification"

input_file <- file.path(
  input_dir, "universal_DEG_table_with_annotation_cluster.xlsx")
output_file <- file.path(
  output_dir, "universal_DEG_table_with_annotation_cluster_regulatory_subgroup.xlsx")

joint_df <- read_excel(input_file)

# ---- thresholds (reassign to be safe) ----
padj_cut <- 0.05
lfc_cut  <- log2(1.5)
padj_generous <- 0.2
lfc_generous <- log2(2.5)

# ---- Directional Kim / Ferrome helpers ----

# Kim et al.
Kim_up <- function(df_row) {
  !is.na(df_row$`Fold Changes (-Fe/+Fe) [Kim et al., 2019]`) &
    df_row$`Fold Changes (-Fe/+Fe) [Kim et al., 2019]` > 1.5
}

Kim_down <- function(df_row) {
  !is.na(df_row$`Fold Changes (-Fe/+Fe) [Kim et al., 2019]`) &
    df_row$`Fold Changes (-Fe/+Fe) [Kim et al., 2019]` < 0.67
}

# Ferrome
Ferrome_up <- function(df_row) {
  !is.na(df_row$`Direction R Ferrome [McInturf et al., 2022]`) &
    df_row$`Direction R Ferrome [McInturf et al., 2022]` == 1
}

Ferrome_down <- function(df_row) {
  !is.na(df_row$`Direction R Ferrome [McInturf et al., 2022]`) &
    df_row$`Direction R Ferrome [McInturf et al., 2022]` == -1
}

# ---- Safety: ensure columns we use exist ----
needed_cols <- c(
  "Gene_ID",
  "Log2FC (LFe vs. WFe) at pH 5.5", "P-value (LFe vs. WFe) at pH 5.5",
  "Log2FC (LFe vs. WFe) at pH 7.5", "P-value (LFe vs. WFe) at pH 7.5",
  "Log2FC (pH 7.5 vs. pH 5.5) at WFe", "P-value (pH 7.5 vs. pH 5.5) at WFe",
  "Log2FC (pH 7.5 vs. pH 5.5) at LFe", "P-value (pH 7.5 vs. pH 5.5) at LFe",
  "Fold Changes (-Fe/+Fe) [Kim et al., 2019]", "Direction R Ferrome [McInturf et al., 2022]",
  "Combined/Interaction"
)

missing_cols <- setdiff(needed_cols, colnames(joint_df))
if(length(missing_cols) > 0){
  stop("Missing expected columns in joint_df: ", paste(missing_cols, collapse = ", "))
}

# -----------------------------
# 1) Subgroup A: genes responding to low Fe at both pH (up & down)
# -----------------------------
A_up <- joint_df %>%
  filter(
    (`Log2FC (LFe vs. WFe) at pH 5.5` > lfc_cut &
       `P-value (LFe vs. WFe) at pH 5.5` < padj_cut) &
      ((`Log2FC (LFe vs. WFe) at pH 7.5` > lfc_cut &
          `P-value (LFe vs. WFe) at pH 7.5` < padj_cut) |
         (`Log2FC (LFe vs. WFe) at pH 7.5` > lfc_generous &
            `P-value (LFe vs. WFe) at pH 7.5` < padj_generous))
  ) %>% pull(Gene_ID)

A_down <- joint_df %>%
  filter(
    (`Log2FC (LFe vs. WFe) at pH 5.5` < -lfc_cut &
       `P-value (LFe vs. WFe) at pH 5.5` < padj_cut) &
      (
        (`Log2FC (LFe vs. WFe) at pH 7.5` < -lfc_cut &
           `P-value (LFe vs. WFe) at pH 7.5` < padj_cut) |
          (`Log2FC (LFe vs. WFe) at pH 7.5` < -lfc_generous &
             `P-value (LFe vs. WFe) at pH 7.5` < padj_generous)
      )
  ) %>% pull(Gene_ID)

Subgroup_A <- union(A_up, A_down)

# ---- A+ ----
Aplus_up <- joint_df %>%
  filter(
    Gene_ID %in% A_up,
    abs(`Log2FC (pH 7.5 vs. pH 5.5) at WFe`) < lfc_cut
  ) %>% pull(Gene_ID)

Aplus_down <- joint_df %>%
  filter(
    Gene_ID %in% A_down,
    abs(`Log2FC (pH 7.5 vs. pH 5.5) at WFe`) < lfc_cut
  ) %>% pull(Gene_ID)

Subgroup_Aplus <- union(Aplus_up, Aplus_down)

# -----------------------------
# 2) Subgroup Aexp
# -----------------------------
A_exp_up <- joint_df %>%
  filter(
    (Kim_up(.) | Ferrome_up(.)) &
      !(`P-value (LFe vs. WFe) at pH 5.5` < padj_cut &
          `Log2FC (LFe vs. WFe) at pH 5.5` > lfc_cut) &
      (
        (`Log2FC (LFe vs. WFe) at pH 7.5` > lfc_cut &
           `P-value (LFe vs. WFe) at pH 7.5` < padj_cut) |
          (`Log2FC (LFe vs. WFe) at pH 7.5` > lfc_generous &
             `P-value (LFe vs. WFe) at pH 7.5` < padj_generous)
      )
  ) %>% pull(Gene_ID)

A_exp_down <- joint_df %>%
  filter(
    (Kim_down(.) | Ferrome_down(.)) &
      !(`P-value (LFe vs. WFe) at pH 5.5` < padj_cut &
          `Log2FC (LFe vs. WFe) at pH 5.5` < -lfc_cut) &
      (
        (`Log2FC (LFe vs. WFe) at pH 7.5` < -lfc_cut &
           `P-value (LFe vs. WFe) at pH 7.5` < padj_cut) |
          (`Log2FC (LFe vs. WFe) at pH 7.5` < -lfc_generous &
             `P-value (LFe vs. WFe) at pH 7.5` < padj_generous)
      )
  ) %>% pull(Gene_ID)

Subgroup_Aexp <- union(A_exp_up, A_exp_down)

# -----------------------------
# 3) Subgroup A+exp
# -----------------------------
Aplus_exp_up <- joint_df %>%
  filter(
    (Kim_up(.) | Ferrome_up(.)) &
      !(`P-value (LFe vs. WFe) at pH 5.5` < padj_cut &
          `Log2FC (LFe vs. WFe) at pH 5.5` > lfc_cut) &
      (
        (
          (`Log2FC (LFe vs. WFe) at pH 7.5` > lfc_cut &
             `P-value (LFe vs. WFe) at pH 7.5` < padj_cut) |
            (`Log2FC (LFe vs. WFe) at pH 7.5` > lfc_generous &
               `P-value (LFe vs. WFe) at pH 7.5` < padj_generous)
        ) &
          abs(`Log2FC (pH 7.5 vs. pH 5.5) at WFe`) < lfc_cut
      )
  ) %>% pull(Gene_ID)

Aplus_exp_down <- joint_df %>%
  filter(
    (Kim_down(.) | Ferrome_down(.)) &
      (
        !(`P-value (LFe vs. WFe) at pH 5.5` < padj_cut &
            `Log2FC (LFe vs. WFe) at pH 5.5` < -lfc_cut) &
          (
            (`Log2FC (LFe vs. WFe) at pH 7.5` < -lfc_cut &
               `P-value (LFe vs. WFe) at pH 7.5` < padj_cut) |
              (`Log2FC (LFe vs. WFe) at pH 7.5` < -lfc_generous &
                 `P-value (LFe vs. WFe) at pH 7.5` < padj_generous)
          ) &
          abs(`Log2FC (pH 7.5 vs. pH 5.5) at WFe`) < lfc_cut
      )
  ) %>% pull(Gene_ID)

Subgroup_Aplus_exp <- union(Aplus_exp_up, Aplus_exp_down)

# -----------------------------
# 4) Subgroup B and Bexp
# -----------------------------
B_up <- joint_df %>%
  filter(
    (`Log2FC (LFe vs. WFe) at pH 5.5` > lfc_cut &
       `P-value (LFe vs. WFe) at pH 5.5` < padj_cut) &
      !(
        (`Log2FC (LFe vs. WFe) at pH 7.5` > lfc_cut &
           `P-value (LFe vs. WFe) at pH 7.5` < padj_cut) |
          (`Log2FC (LFe vs. WFe) at pH 7.5` > lfc_generous &
             `P-value (LFe vs. WFe) at pH 7.5` < padj_generous)
      ) &
      (`Log2FC (pH 7.5 vs. pH 5.5) at WFe` > lfc_cut &
         `P-value (pH 7.5 vs. pH 5.5) at WFe` < padj_cut)
  ) %>% pull(Gene_ID)

B_down <- joint_df %>%
  filter(
    (`Log2FC (LFe vs. WFe) at pH 5.5` < -lfc_cut &
       `P-value (LFe vs. WFe) at pH 5.5` < padj_cut) &
      !(
        (`Log2FC (LFe vs. WFe) at pH 7.5` < -lfc_cut &
           `P-value (LFe vs. WFe) at pH 7.5` < padj_cut) |
          (`Log2FC (LFe vs. WFe) at pH 7.5` < -lfc_generous &
             `P-value (LFe vs. WFe) at pH 7.5` < padj_generous)
      ) &
      (`Log2FC (pH 7.5 vs. pH 5.5) at WFe` < -lfc_cut &
         `P-value (pH 7.5 vs. pH 5.5) at WFe` < padj_cut)
  ) %>% pull(Gene_ID)

Subgroup_B <- union(B_up, B_down)

Bexp_up <- joint_df %>%
  filter(
    (Kim_up(.) | Ferrome_up(.)) &
      (
        !(`Log2FC (LFe vs. WFe) at pH 5.5` > lfc_cut &
            `P-value (LFe vs. WFe) at pH 5.5` < padj_cut) &
          !(
            (`Log2FC (LFe vs. WFe) at pH 7.5` > lfc_cut &
               `P-value (LFe vs. WFe) at pH 7.5` < padj_cut) |
              (`Log2FC (LFe vs. WFe) at pH 7.5` > lfc_generous &
                 `P-value (LFe vs. WFe) at pH 7.5` < padj_generous)
          )
      ) &
      (`Log2FC (pH 7.5 vs. pH 5.5) at WFe` > lfc_cut &
         `P-value (pH 7.5 vs. pH 5.5) at WFe` < padj_cut)
  ) %>% pull(Gene_ID)

Bexp_down <- joint_df %>%
  filter(
    (Kim_down(.) | Ferrome_down(.)) &
      (
        !(`Log2FC (LFe vs. WFe) at pH 5.5` < -lfc_cut &
            `P-value (LFe vs. WFe) at pH 5.5` < padj_cut) &
          !(
            (`Log2FC (LFe vs. WFe) at pH 7.5` < -lfc_cut &
               `P-value (LFe vs. WFe) at pH 7.5` < padj_cut) |
              (`Log2FC (LFe vs. WFe) at pH 7.5` < -lfc_generous &
                 `P-value (LFe vs. WFe) at pH 7.5` < padj_generous)
          )
      ) &
      (`Log2FC (pH 7.5 vs. pH 5.5) at WFe` < -lfc_cut &
         `P-value (pH 7.5 vs. pH 5.5) at WFe` < padj_cut)
  ) %>% pull(Gene_ID)

Subgroup_Bexp <- union(Bexp_up, Bexp_down)

# -----------------------------
# 5) Safety
# -----------------------------
if(!exists("Subgroup_A")) Subgroup_A <- character(0)
if(!exists("Subgroup_Aplus")) Subgroup_Aplus <- character(0)
if(!exists("Subgroup_Aexp")) Subgroup_Aexp <- character(0)
if(!exists("Subgroup_Aplus_exp")) Subgroup_Aplus_exp <- character(0)
if(!exists("Subgroup_B")) Subgroup_B <- character(0)
if(!exists("Subgroup_Bexp")) Subgroup_Bexp <- character(0)

# -----------------------------
# 6) Assemble Subgroup column
# -----------------------------
gene_ids <- joint_df$Gene_ID

subgroup_tag_vec <- sapply(gene_ids, function(gid) {
  tags <- c()
  if (gid %in% Subgroup_Aplus) tags <- c(tags, "A+")
  else if (gid %in% Subgroup_A) tags <- c(tags, "A")
  if (gid %in% Subgroup_Aplus_exp) tags <- c(tags, "A+exp")
  else if (gid %in% Subgroup_Aexp) tags <- c(tags, "Aexp")
  if (gid %in% Subgroup_B) tags <- c(tags, "B")
  if (gid %in% Subgroup_Bexp) tags <- c(tags, "Bexp")
  if (length(tags) == 0) return(NA_character_)
  paste(tags, collapse = "_")
}, USE.NAMES = FALSE)

joint_df$Subgroup <- subgroup_tag_vec

# -----------------------------
# 7) Reorder columns
# -----------------------------
cols <- colnames(joint_df)
pos <- which(cols == "Combined/Interaction")
new_order <- c(cols[1:pos], "Subgroup", cols[(pos+1):length(cols)])
new_order <- unique(new_order)
joint_df <- joint_df[, new_order]

# -----------------------------
# 8) QC
# -----------------------------
cat("Counts per subgroup:\n")
cat("A:", length(Subgroup_A),
    "A+:", length(Subgroup_Aplus),
    "Aexp:", length(Subgroup_Aexp),
    "A+exp:", length(Subgroup_Aplus_exp), "\n")
cat("B:", length(Subgroup_B),
    "Bexp:", length(Subgroup_Bexp), "\n")

# -----------------------------
# 9) Export
# -----------------------------
write_xlsx(joint_df, output_file)
cat("Wrote output to:", output_file, "\n")