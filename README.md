# stromal_subclasses
code for data processing and visualisation associated with the Wu et al. manuscript "Single-cell analysis reveals diverse stromal subsets associated with immune evasion in triple-negative breast cancer"

### 01 cellranger count processing
job submission script for single-cell RNA-Seq processing using cellranger v2.1.1. 

### 02 seurat v2 processing of individual datasets
job submission script and R script for processing individual seurat objects

### 03 seurat v2 data integration
job submission script and R script for integrating seurat objects

### 04 cluster annotation and reclustering 
job submission script and R script for cluster annotations and reclustering of epithelial cells, stromal-immune cells and individual stromal subpopulations (CAFs and VDSCs). This R script also includes the generation of stromal gene signatures and export of gene expression matrices for downstream SCENIC (step 06).

### 05 AUCell gene signature scoring
job submission script and R script for scoring stromal-immune cells with cell type signatures (XCell database) using the AUCell method

### 06 pySCENIC transcription factor enrichment of reclustered CAFs and VDSCs
job submission script for running pySCENIC (python based command line version) using the CAFs and VDSC raw count expression matrix as input. This also includes R script for filtering top TF candidates for clustering and visualisation

### 07 stromal cell signalling predictions
R script visualisation scripts for filtering stromal cell-cell signalling predictions. For processing ligand-receptor analysis, please visit (enter git link for Ruis github link)

### Other analytical tools used
#### Immune evasion using TIDE
For computing T-cell dysfunction and exclusion analysis, please see the tumour immune dysfunction and exclusion (TIDE) method

- paper: https://www.nature.com/articles/s41591-018-0136-1
- link to online tool: http://tide.dfci.harvard.edu/


