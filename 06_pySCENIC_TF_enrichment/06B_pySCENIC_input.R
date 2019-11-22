# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 06B EXPORTING AND FILTERING GENE MATRICES FOR PYSCENIC INPUT
# R 3.5.0
#
### inputs (in order); 
#     - seurat object
# 
### REQUIRES R v.3.5.0
### QSUB ARGUMENTS
# qsub \
# -cwd \
# -pe smp 8 \
# -l mem_requested=10G \
# -P TumourProgression \
# -b y \
# -j y \
# -V \
# -N ${JOBNAME}\
# "${R} CMD BATCH \
# --no-save '--args \
# ${SEURATOBJECTPATH}' \
# ${SCRIPT}"
#
#
#
# SCENIC NORMALISED GENE EXPRESSION MATRIX + GENE FILTERING 
#
## from Lambrechts et al. 2018
## Phenotype molding of stromal cells in the lung tumor microenvironment
#
#
# SCENIC analysis. 
# The SCENIC analysis was run as described24 on the 52,698 
# cells that passed the filtering, using the 20-thousand motifs database for 
# RcisTarget and GRNboost (SCENIC version 0.1.5, which corresponds to RcisTarget 
# 0.99.0 and AUCell 0.99.5; with RcisTarget.hg19.motifDatabases.20k). The input 
# matrix was the normalized expression matrix, output from Seurat, from which 9,919 
# genes passed the filtering (sum of expression >3 × 0.005 × 52,698 and detected in 
# at least 0.5% of the cells).
#
#
# Analysis of differential pathway or regulon activities. 
# To assess differential activities 
# of pathways (GSVA) or regulons (SCENIC) between sets of cells
# (for example, derived from tumor or normal samples, or belonging to different subclusters), 
# we contrasted the activity scores for each cell using a generalized
# linear model. To avoid inflating signals because of interindividual differences
# (for example, in the relative frequencies of cells from different patients), 
# we always included the patient of origin as a categorical variable. Results of these linear models
# were visualized using bar plots or heatmaps. For the latter, pathways or regulons that did not show
# significant changes (Benjamini–Hochberg-corrected P value >0.05) in any of the sets of cells contrasted
# in one analysis were not visualized.
# 01 LOAD PACKAGES -----------------------------------------------------------

library(Seurat)

# 02 COMMAND LINE ARGUMENTS  ------------------------------------------------------------

temp_start_time <- date()
temp_args <-
  commandArgs(trailingOnly = T)

# input path to seurat object
temp_seurat_object_path <- 
  temp_args[1]

# 03 SET UP AND FUNCTIONS ---------------------------------------------------------------

dir.create("data")

# 04 LOAD DATA ---------------------------------------------------------------

seurat_10X <- 
  readRDS(temp_seurat_object_path)

# 05 EXPORT AND FILTER NORMALIZED EXPRESSION MATRIX --------------------------------------

  temp_normalised.data.frame <- as.matrix(x = seurat_10X@data)
  temp_cell_number <- length(seurat_10X@ident)
  temp_threshold <- 3*0.005*temp_cell_number

  temp_genes_to_filer <- data.frame(gene = row.names(temp_normalised.data.frame),
                                    threshold_value = rowSums(temp_normalised.data.frame))

  temp_genes_to_filer <- 
    subset(temp_genes_to_filer,
           threshold_value > temp_threshold)
  
  temp_normalised.data.frame_filtered <- temp_normalised.data.frame[row.names(temp_normalised.data.frame) %in% temp_genes_to_filer$gene ,]

  # export matrix
  write.csv(
    temp_normalised.data.frame_filtered,
    file = "data/Stromal_normalized_expression_matrix_filtered.csv",
    quote = F,
    row.names = T
  )
  
  # export metadata
  write.csv(
          data.frame(
            barcode = row.names(seurat_10X@meta.data),
            Cluster_ID = seurat_10X@meta.data$stromal_celltype),
    file = "data/Stromal_metadata.csv", 
    row.names = F,
    quote = F
  )
