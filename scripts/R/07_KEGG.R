# ================================
# KEGG enrichment per cluster
# ================================

library(clusterProfiler)
library(org.At.tair.db)
library(AnnotationDbi)
library(readr)
library(dplyr)
library(openxlsx)

# ---- Input ----
input_file <- "data/clustering/clustering_mfuzz.csv"
data <- read_csv(input_file, show_col_types = FALSE)

clusters <- sort(unique(data$cluster))
clusters <- clusters[!is.na(clusters)]

# ---- Output ----
out_dir <- "results/07_KEGG_enrichment"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- helper: TAIR → Entrez ----
map_to_entrez <- function(genes){
  gene_entrez <- bitr(
    genes,
    fromType = "TAIR",
    toType   = "ENTREZID",
    OrgDb    = org.At.tair.db
  )
  
  gene_entrez <- gene_entrez %>%
    filter(!is.na(ENTREZID))
  
  unique(gene_entrez$ENTREZID)
}

run_KEGG_cluster <- function(cl){
  
  genes <- data$gene[data$cluster == cl]
  genes <- unique(na.omit(genes))
  
  if (length(genes) == 0) return(NULL)
  
  entrez <- map_to_entrez(genes)
  
  if (length(entrez) == 0) return(NULL)
  
  kk <- enrichKEGG(
    gene          = entrez,
    organism      = "ath",
    keyType       = "ncbi-geneid",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 1
  )
  
  if (is.null(kk) || nrow(as.data.frame(kk)) == 0) return(NULL)
  
  as.data.frame(kk)
}

# ---- workbook ----
wb <- createWorkbook()

for (cl in clusters){
  cat("KEGG Cluster", cl, "\n")
  
  res <- run_KEGG_cluster(cl)
  
  if (!is.null(res) && nrow(res) > 0){
    sheet <- paste0("Cluster", cl)
    addWorksheet(wb, sheet)
    writeData(wb, sheet, res)
  }
}

saveWorkbook(
  wb,
  file.path(out_dir, "KEGG_clusters.xlsx"),
  overwrite = TRUE
)

message("KEGG done")