# ================================
# 03_External_Annotation_Enrichment.R
# ================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(writexl)
})

# ---- Paths ----
output_dir <- "E:/桌面/Iron Project_Gen/!Transcriptomic experiment/!Data Analysis/Ath_featurecount_Ath_root/!For TRR paper/Gene lists"

input_file <- file.path(
  output_dir,
  "02_Ath_root_universal_DEG_table_with_annotation_and_expression_20260505.xlsx"
)

output_file <- file.path(
  output_dir,
  "03_Ath_root_universal_DEG_table_with_annotation_and_expression_external_annotations_20260505.xlsx"
)

# ---- Helper: load gene IDs from .txt (first column) ----
load_gene_ids <- function(path) {
  df <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  ids <- toupper(trimws(df[[1]]))
  unique(ids)
}

# ---- Load main DEG table ----
deg_df <- read_excel(input_file)

if (!"Gene_ID" %in% colnames(deg_df)) {
  stop("Input table must contain a 'Gene_ID' column.")
}

deg_df <- deg_df %>% distinct(Gene_ID, .keep_all = TRUE)

# ================================
# 1. Kim et al., 2019
# ================================
kim_file <- file.path(output_dir, "Kim et al_2019.xlsx")
kim_df <- read_excel(kim_file)

deg_df <- deg_df %>%
  left_join(
    kim_df %>% select(Gene_ID, `Fold Changes (-Fe/+Fe)`, `Raw p value`),
    by = "Gene_ID"
  ) %>%
  rename(
    `Fold Changes (-Fe/+Fe) [Kim et al., 2019]` = `Fold Changes (-Fe/+Fe)`,
    `P-value [Kim et al., 2019]` = `Raw p value`
  )

# ================================
# 2. Internal Metal Homeostasis
# ================================
metal_file <- file.path(output_dir, "!metal_homeostasis_2024_03_18_v17_deleted_redundancy.xlsx")
metal_df <- read_excel(metal_file, sheet = "Metal_homeostasis_generous_upd") %>%
  mutate(AGI_Number = toupper(trimws(AGI_Number)))

deg_df <- deg_df %>%
  left_join(
    metal_df %>% select(AGI_Number, `Short Name`, Function, DOI),
    by = c("Gene_ID" = "AGI_Number")
  ) %>%
  rename(
    `Metal Homeostasis (short name)` = `Short Name`,
    `Metal Homeostasis (function)` = Function,
    `Metal Homeostasis Citation(s)` = DOI
  )

# ================================
# 3. Genes affecting ionome
# ================================
ionome_file <- file.path(output_dir, "!genes_affecting_ionome.xlsx")
ionome_df <- read_excel(ionome_file, sheet = "A.thaliana")

deg_df <- deg_df %>%
  left_join(
    ionome_df %>% select(GeneID, Elements, `Citation(s) - DOI only`),
    by = c("Gene_ID" = "GeneID")
  ) %>%
  rename(
    `Genes Affecting Ionome - Elements [Whitt et al., 2020]` = Elements,
    `Genes Affecting Ionome - Citation(s) [Whitt et al., 2020]` = `Citation(s) - DOI only`
  )

# ================================
# 4. Fe Metalloprotein
# ================================
fe_files <- list(
  "Fe cation" = "2023_07_26_Fe_cation_binding_genes_ZhangFPLS.txt",
  "heme/cytochrome" = "2023_07_26_Heme_containing_proteins_ZhangFPLS.txt",
  "FeS cluster use" = "2023_07_26_Fe-S_cluster_containing_proteins_Zhang.txt",
  "FeS cluster biogenesis" = "2023_07_26_FE-S_assembly_proteins_Zhang.txt"
)

deg_df$`Fe Metalloprotein` <- NA_character_

for (type in names(fe_files)) {
  genes <- load_gene_ids(file.path(output_dir, fe_files[[type]]))
  
  deg_df$`Fe Metalloprotein` <- ifelse(
    deg_df$Gene_ID %in% genes,
    ifelse(
      is.na(deg_df$`Fe Metalloprotein`),
      type,
      paste(deg_df$`Fe Metalloprotein`, type, sep = "; ")
    ),
    deg_df$`Fe Metalloprotein`
  )
}

# ================================
# 5. Metal-containing proteins
# ================================
metal_contain_files <- list(
  "Mn-Containing" = "2023_07_26_Mn_containing_proteins_ZhangFPLS.txt",
  "Zn-Containing" = "2023_07_26_Zn_containing_proteins_ZhangFPLS.txt",
  "Cu-Containing" = "2023_07_26_Cu_containing_proteins_ZhangFPLS.txt"
)

for (col_name in names(metal_contain_files)) {
  genes <- load_gene_ids(file.path(output_dir, metal_contain_files[[col_name]]))
  deg_df[[col_name]] <- ifelse(deg_df$Gene_ID %in% genes, "Yes", "")
}

# ================================
# 6. Casparian Strip / Suberin
# ================================
casparian_file <- file.path(output_dir, "Casparian_strip_and_suberin_Baohai2019.xlsx")
casparian_df <- read_excel(casparian_file)

deg_df <- deg_df %>%
  left_join(
    casparian_df %>% select(Gene_ID, `Casparian Strip/Suberin`),
    by = "Gene_ID"
  )

# ================================
# 7. Meiosis
# ================================
meiosis_file <- file.path(output_dir, "2023_07_26_Meiosis_2019.txt")
meiosis_genes <- load_gene_ids(meiosis_file)

deg_df$Meiosis <- ifelse(deg_df$Gene_ID %in% meiosis_genes, "Yes", "")

# ================================
# 8. DNA Repair
# ================================
dna_file <- file.path(output_dir, "2023_08_03_DDR_and_recombination_Fullset_AP2023_edited.txt")
dna_genes <- load_gene_ids(dna_file)

deg_df$`DNA Repair` <- ifelse(deg_df$Gene_ID %in% dna_genes, "Yes", "")

# ================================
# 9. Root Hair / Trichoblast
# ================================
root_file <- file.path(output_dir, "Root hair and trichoblast _GO_TAIR_Gen_20250824.xlsx")
root_df <- read_excel(root_file, sheet = "Root hair and trichoblast") %>%
  select(Gene_ID, GO_Term_Description)

root_summary <- root_df %>%
  group_by(Gene_ID) %>%
  summarise(
    `Root Hair/Trichoblast` = paste(unique(GO_Term_Description), collapse = "; "),
    .groups = "drop"
  )

deg_df <- deg_df %>%
  left_join(root_summary, by = "Gene_ID")

# ================================
# 10. Root Ferrome
# ================================
ferrome_file <- file.path(
  output_dir,
  "!leaf_root_ferrome_McInturf_Mondoza-Cozatl_JExpBot_2021_suppl_supplementary_table_s1.xlsx"
)
ferrome_df <- read_excel(ferrome_file, sheet = "Root")

colnames(ferrome_df)[colnames(ferrome_df) %in% c("AGI R Ferrome", "SN R Ferrome", "Direction R Ferrome")] <-
  c("Gene_ID_Ferrome", "SN R Ferrome", "Direction R Ferrome")

deg_df <- deg_df %>%
  left_join(
    ferrome_df %>% select(Gene_ID_Ferrome, `SN R Ferrome`, `Direction R Ferrome`),
    by = c("Gene_ID" = "Gene_ID_Ferrome")
  )%>%
  rename(
    `SN R Ferrome [McInturf et al., 2022]` = `SN R Ferrome`,
    `Direction R Ferrome [McInturf et al., 2022]` = `Direction R Ferrome`
  )

# ================================
# 11. FIT / PYE target
# ================================
FIT_PYE_file <- file.path(output_dir, "Schmidt_Buckhout_FIT_PYE_target.xlsx")
FIT_PYE_df <- read_excel(FIT_PYE_file)

deg_df <- deg_df %>%
  left_join(
    FIT_PYE_df %>% select(ATG, `Ferrome [Schmidt and Buckhout, 2011]`),
    by = c("Gene_ID" = "ATG")
  )

# ================================
# 12. Final cleanup
# ================================
deg_df <- deg_df %>%
  distinct(Gene_ID, .keep_all = TRUE)

# ================================
# 13. Save enriched table
# ================================
write_xlsx(deg_df, output_file)

message("External annotation table written to: ", output_file)