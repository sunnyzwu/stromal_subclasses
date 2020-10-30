# stromal_subclasses
code for data processing and visualisation associated with the Wu et al. (2020) study "Stromal cell diversity associated with immune evasion in human triple‐negative breast cancer"
- https://www.embopress.org/doi/full/10.15252/embj.2019104063

### Data availability
Processed data, including count matrices (raw and normalised) and cell metadata, can be found at the Broad Single Cell Portal at the following link:
- https://singlecell.broadinstitute.org/single_cell/study/SCP1106/stromal-cell-diversity-associated-with-immune-evasion-in-human-triple-negative-breast-cancer

### Code Summary
#### 01 cellranger count processing
job submission script for single-cell RNA-Seq processing using cellranger v2.1.1. 

#### 02 seurat v2 processing of individual datasets
job submission script and R script for processing individual seurat objects

#### 03 seurat v2 data integration
job submission script and R script for integrating seurat objects

#### 04 cluster annotation and reclustering 
job submission script and R script for cluster annotations and reclustering of epithelial cells, stromal-immune cells and individual stromal subpopulations (CAFs and PVL cells). This R script also includes the generation of stromal gene signatures and export of gene expression matrices for downstream SCENIC (step 06).

#### 05 gene signature analysis
##### AUCell
job submission script and R scripts for scoring of stromal-immune cells with cell type signatures (XCell database), T-cell exhaustion signatures (Blackburn et al. 2008) and the pancreatic ductal adenocarcinoma CAF signatures (David Tuveson's lab) using the AUCell method

##### clusterProfiler
job submission script and R scripts for gene ontology enrichment of the gene signatures for each stromal subcluster. Top 250 DEGs are used. 

#### 06 pySCENIC transcription factor enrichment of reclustered CAFs and PVL cells
job submission script and R scripts for exporting, filtering gene expression matrices, and running pySCENIC (python based command line version) with the CAFs and PVL raw count expression matrix as input. This also includes R script for filtering top TF candidates for clustering and visualisation

#### 07 stromal cell signalling predictions
R script visualisation scripts for filtering stromal cell-cell signalling predictions. 


### Other analytical tools used
#### Immune evasion using TIDE
For computing T-cell dysfunction and exclusion analysis, please see the tumour immune dysfunction and exclusion (TIDE) method

- paper: https://www.nature.com/articles/s41591-018-0136-1
- link to online tool: http://tide.dfci.harvard.edu/

#### Cell signalling using NATMI
For processing ligand-receptor analysis, please see the NATMI documentation at:
- https://github.com/asrhou/NATMI

### Contact
Please email s.wu@garvan.org.au or a.swarbrick@garvan.org.au for any additional questions about the analytical methods used in this paper. All other relevant data are available from the authors upon request.
