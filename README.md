# Iron_pH_Ath_transcriptomics
This repository contains the RNA-seq data analysis pipeline for Arabidopsis thaliana root samples under different iron (Fe) and pH conditions.


---

## Read mapping

Clean paired-end reads (150 bp) were mapped to the Arabidopsis thaliana TAIR10.1 reference genome using STAR v2.7.10b.

Gene-level read counts were generated using featureCounts v2.0.0 with the following parameters:
- feature type: exon (-t exon)
- strand-specific: yes (-s 2)
- paired-end mode: -p

---

## Differential expression analysis

Filtering:
- genes with counts ≥ 10 in at least 3 samples were retained

DE analysis:
- performed in R v4.5.1 using DESeq2 v1.48.2
- design: ~ batch + Fe * pH

Significance thresholds:
- |log2FC| > log2(1.5) ≈ 0.585
- adjusted p-value < 0.05

---

## Gene clustering

Significantly differentially expressed genes were clustered using Mfuzz (soft clustering).

Comparisons excluded:
- LFe pH 5.5 vs WFe pH 7.5

---

## Functional enrichment

For each cluster:
- GO enrichment analysis (adjusted p-value < 0.05)
- KEGG pathway enrichment analysis (adjusted p-value < 0.1)

Performed using clusterProfiler v4.16.0 with Benjamini–Hochberg correction.


---

## Outputs

- DESeq2 results (normalized counts, DEGs, statistics)
- Mfuzz clustering results
- GO enrichment tables
- KEGG enrichment tables
