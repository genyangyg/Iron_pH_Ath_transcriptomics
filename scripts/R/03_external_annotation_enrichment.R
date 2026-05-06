# ================================
# 03_external_annotation_enrichment.R
# ================================

  library(dplyr)
  library(readxl)
  library(writexl)


# ---- Paths ----
input_dir <- "results/02_DEG_classification"
internal_dir <- "data/annotation/internal"
external_dir <- "data/annotation/external"
output_dir <- "results/03_functional_annotation"

input_file <- file.path(
  input_dir,
  "universal_DEG_table.xlsx"
)

output_file <- file.path(
  output_dir,
  "universal_DEG_table_with_annotation.xlsx"
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
kim_file <- file.path(external_dir, "Kim_et_al_2019_Fe_response.xlsx")
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
metal_file <- file.path(internal_dir, "metal_homeostasis_gene.xlsx")
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
ionome_file <- file.path(external_dir, "Whitt_et_al_2020_genes_affecting_ionome.xlsx")
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
  "Fe cation" = "Zhang_2018_Fe_cation_binding_genes.txt",
  "heme/cytochrome" = "Zhang_2018_Heme_containing_proteins.txt",
  "FeS cluster use" = "Zhang_2018_Fe-S_cluster_containing_proteins.txt",
  "FeS cluster biogenesis" = "Zhang_2018_FE-S_assembly_proteins.txt"
)

deg_df$`Fe Metalloprotein` <- NA_character_

for (type in names(fe_files)) {
  genes <- load_gene_ids(file.path(internal_dir, fe_files[[type]]))
  
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
  "Mn-Containing" = "Zhang_2018_Mn_containing_proteins.txt",
  "Zn-Containing" = "Zhang_2018_Zn_containing_proteins.txt",
  "Cu-Containing" = "Zhang_2018_Cu_containing_proteins.txt"
)

for (col_name in names(metal_contain_files)) {
  genes <- load_gene_ids(file.path(internal_dir, metal_contain_files[[col_name]]))
  deg_df[[col_name]] <- ifelse(deg_df$Gene_ID %in% genes, "Yes", "")
}

# ================================
# 6. Casparian Strip / Suberin
# ================================
casparian_file <- file.path(internal_dir, "Casparian_strip_and_suberin.xlsx")
casparian_df <- read_excel(casparian_file)

deg_df <- deg_df %>%
  left_join(
    casparian_df %>% select(Gene_ID, `Casparian Strip/Suberin`),
    by = "Gene_ID"
  )

# ================================
# 7. Meiosis
# ================================
meiosis_file <- file.path(internal_dir, "Meiosis.txt")
meiosis_genes <- load_gene_ids(meiosis_file)

deg_df$Meiosis <- ifelse(deg_df$Gene_ID %in% meiosis_genes, "Yes", "")

# ================================
# 8. DNA Repair
# ================================
dna_file <- file.path(internal_dir, "DDR_and_recombination.txt")
dna_genes <- load_gene_ids(dna_file)

deg_df$`DNA Repair` <- ifelse(deg_df$Gene_ID %in% dna_genes, "Yes", "")

# ================================
# 9. Root Hair / Trichoblast
# ================================
root_file <- file.path(internal_dir, "Root_hair_and_trichoblast_GO_TAIR_processed.xlsx")
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
  external_dir,
  "McInturf_et_al_2022_leaf_root_ferrome.xlsx"
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
FIT_PYE_file <- file.path(external_dir, "Schmidt_and_Buckhout_2011_FIT_PYE_target.xlsx")
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