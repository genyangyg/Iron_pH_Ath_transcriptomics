# ================================
# GO enrichment per cluster (BP/MF/CC)
# ================================

library(clusterProfiler)
library(org.At.tair.db)
library(readr)
library(dplyr)
library(openxlsx)

# ---- Input ----
input_file <- "data/clustering/clustering_mfuzz.csv"

data <- read_csv(input_file, show_col_types = FALSE)

clusters <- sort(unique(data$cluster))
clusters <- clusters[!is.na(clusters)]

# ---- Output ----
out_dir <- "results/06_GO_enrichment"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- GO function ----
run_GO <- function(gene_list, ont){
  
  ego <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  
  df <- as.data.frame(ego)
  df$Ontology <- ont
  return(df)
}

run_GO_cluster <- function(cl){
  
  genes <- data$gene[data$cluster == cl]
  genes <- unique(na.omit(genes))
  
  if (length(genes) == 0) return(NULL)
  
  bp <- run_GO(genes, "BP")
  mf <- run_GO(genes, "MF")
  cc <- run_GO(genes, "CC")
  
  bind_rows(bp, mf, cc)
}

# ---- workbook ----
wb <- createWorkbook()

for (cl in clusters){
  cat("GO Cluster", cl, "\n")
  
  res <- run_GO_cluster(cl)
  
  if (!is.null(res) && nrow(res) > 0){
    addWorksheet(wb, paste0("Cluster", cl))
    writeData(wb, paste0("Cluster", cl), res)
  }
}

saveWorkbook(
  wb,
  file.path(out_dir, "GO_clusters.xlsx"),
  overwrite = TRUE
)

message("GO done")