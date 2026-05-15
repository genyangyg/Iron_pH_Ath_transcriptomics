# External Annotation Data

## Overview

This directory contains gene annotation datasets obtained from published studies and public resources.  
These datasets are integrated into the analysis pipeline to functionally annotate differentially expressed genes (DEGs).

All files are used in:

scripts/03_external_annotation_enrichment.R

---

## Data Organization

Each dataset corresponds to a published study or curated external resource.  
File names are standardized for consistency and reproducibility.

---

## Datasets and Sources

### 1. Iron Deficiency Response (Kim et al., 2019)
- File: Kim_et_al_2019_Fe_response.xlsx  
- Description: Gene expression changes under iron deficiency (-Fe vs. +Fe)  
- Key columns used:
  - Fold Changes (-Fe/+Fe)
  - Raw p value  
- Source:
  Kim et al., 2019  

---

### 2. Genes Affecting Ionome (Whitt et al., 2020)
- File: Whitt_et_al_2020_genes_affecting_ionome.xlsx  
- Description: Genes influencing elemental composition (ionome) in Arabidopsis  
- Key columns used:
  - Elements
  - Citation(s) (DOI only)  
- Source:
  Whitt et al., 2020  

---

### 3. Root Ferrome (Mclnturf et al., 2022)
- File: Mclnturf_et_al_2022_leaf_root_ferrome.xlsx  
- Description: Root ferrome gene set and expression direction  
- Key columns used:
  - SN R Ferrome
  - Direction R Ferrome  
- Source:
  McInturf et al., 2022  

---

### 4. FIT and PYE Target Genes (Schmidt and Buckhout, 2011)
- File: Schmidt_and_Buckhout_2011_FIT_PYE_target.xlsx  
- Description: Known transcriptional targets of FIT and PYE  
- Key columns used:
  - Ferrome classification  
- Source:
  Schmidt and Buckhout, 2011  

---

### 5. Cell-type specific Fe deficiency response (Dinneny et al., 2008)
- File: Dinneny_et_al_2008_cell_type_Fe_response.xls
- Description: Cell-type specific differential gene expression under iron deficiency in Arabidopsis root
-Key column used:
  - Cell type differentially expressed in
- Source:
Dinneny et al., 2008

---

## Data Processing Notes

- Gene identifiers were standardized to uppercase AGI format (e.g., AT1G01010).
- Column names were simplified or renamed for consistency across datasets.
- In some cases, only a subset of columns was retained for integration.
- No biological values were altered.

---

## Usage in Pipeline

These datasets are integrated using left joins based on Gene_ID:

- Matching key: Gene_ID (AGI format)
- Integration method: dplyr::left_join()

Annotated features include:
- Iron deficiency response metrics
- Ionome-related gene annotations
- Ferrome classification
- Known transcriptional targets

---

## Reproducibility

To ensure reproducibility:

- File names must remain unchanged
- File structure should not be modified
- Scripts assume exact column names as used in the pipeline

If any file is modified, corresponding script updates are required.

---

## Licensing and Data Use

- These datasets originate from published studies.
- Copyright and data ownership remain with the original authors and publishers.
- Files are included here solely for:
  - Reproducibility
  - Academic research use

If you intend to reuse these datasets independently, please consult the original publications.

---

## Citation

Please cite the original publications when using these datasets:

- DINNENY, J. R., LONG, T. A., WANG, J. Y., JUNG, J. W., MACE, D., POINTER, S., BARRON, C., BRADY, S. M., SCHIEFELBEIN, J. & BENFEY, P. N. 2008. Cell identity mediates the response of Arabidopsis roots to abiotic stress. Science, 320(5878), 942-945.
- KIM, S. A., LACROIX, I. S., GERBER, S. A. & GUERINOT, M. L. 2019. The iron deficiency response in Arabidopsis thaliana requires the phosphorylated transcription factor URI. Proc Natl Acad Sci USA, 116, 24933-24942.
- WHITT, L., RICACHENEVSKY, F. K., ZIEGLER, G. Z., CLEMENS, S., WALKER, E., MAATHUIS, F. J. M., KEAR, P. & BAXTER, I. 2020. A curated list of genes that affect the plant ionome. Plant Direct, 4, e00272.
- MCINTURF, S. A., KHAN, M. A., LI, J., MARJAULT, H.-B., FICHMAN, Y., KUNZ, H.-H., GOKUL, A., CASTRO-GUERRERO, N. A., HÖHNER, R., GOGGIN, F. L., KEYSTER, M., NECHUSTHAI, R., MITTLER, R. & MENDOZA-CÓZATL, D. G. 2022. Cadmium interference with iron sensing reveals transcriptional programs sensitive and insensitive to reactive oxygen species. Journal of Experimental Botany, 73, 324-338.
- SCHMIDT, W. & BUCKHOUT, T. J. 2011. A hitchhiker's guide to the Arabidopsis ferrome. Plant Physiol Biochem, 49, 462-70.

---

## Disclaimer

This repository redistributes processed versions of publicly available datasets.  
While care has been taken to preserve data integrity, users are encouraged to verify against the original sources.
