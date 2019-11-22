# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 05B AUCELL SCORING SCRIPT HPC
# R 3.5.0
#
### inputs (in order); 
#     - seurat object
#     - gene sets
# 
### REQUIRES R v.3.5.0
### QSUB ARGUMENTS
# qsub \
# -cwd \
# -pe smp 16 \
# -l mem_requested=10G \
# -P TumourProgression \
# -b y \
# -j y \
# -V \
# -N ${JOBNAME}\
# "${R} CMD BATCH \
# --no-save '--args \
# ${SEURATOBJECTPATH} \
# ${GENESETS}' \
# ${SCRIPT}"
#
#
#         Input gene set should be in the format of a csv file
#         e.g.
#         gene_set_1,gene_set_2,gene_set_3
#         GENE1,GENE1,GENE1
#         GENE2,GENE2,GENE2
#         GENE3,GENE23,GENE3
#         ,GENE4,GENE4
#         ,GENE5,
#
#
# LOAD PACKAGES -----------------------------------------------------------

library(Seurat)
library(readr)
library(GSEABase)
library(AUCell)
library(reshape2)
library(NMF)


# COMMAND LINE ARGUMENTS  ------------------------------------------------------------

temp_start_time <- date()
temp_args <-
  commandArgs(trailingOnly = T)

# input path to seurat object
temp_seurat_object_path <- 
  temp_args[1]

# input gene sets
temp_genesets <- 
  temp_args[2]

# SET UP AND FUNCTIONS ---------------------------------------------------------------

dir.create("Output")
dir.create("Output/Figures")
dir.create("Output/Rdata")

# INPUT MATRIX ------------------------------------------------------------

seurat_10X <- 
  readRDS(temp_seurat_object_path)

#If using subset data or certain cell barcodes, set ident
temp_ident <- 
  NULL

temp_exprMatrix <-
  as.matrix(x = seurat_10X@raw.data[, WhichCells(object = seurat_10X,
                                                 ident = temp_ident,
                                                 cells.use = NULL)])

temp_first_geneset <- 
  (ncol(seurat_10X@meta.data) + 1)

# INPUT GENE SETS ---------------------------------------------------------

temp_xcell_genesets <- 
  read.csv(temp_genesets)

temp_GeneSetCollection_list <- NULL
for(i in c(1:ncol(temp_xcell_genesets))) {
  
  temp_set_name <- colnames(temp_xcell_genesets[i])
  temp_set_name <- as.character(temp_set_name)
  
  temp_set <- na.omit(temp_xcell_genesets[i])
  temp_set <- as.vector(temp_set[,1])
  temp_set <- temp_set[! temp_set %in% ""]
  temp_set <- unique(temp_set)
  temp_set <- GeneSet(temp_set,
                      setName = temp_set_name)
  temp_GeneSetCollection_list <- append(temp_GeneSetCollection_list,
                                        temp_set)
  
}

temp_gene_set_collection <- 
  GeneSetCollection(temp_GeneSetCollection_list)

rm(list = ls(pattern = "temp_set_name"))


# RUN AUCELL --------------------------------------------------------------

# build rankings
temp_cells_rankings <- 
  AUCell_buildRankings(temp_exprMatrix, 
                       nCores = 1, 
                       plotStats = F)

# subset gene sets
temp_subsetgeneSets <- 
  subsetGeneSets(temp_gene_set_collection, 
                 rownames(temp_exprMatrix)) 

# calculate area under the curve
temp_cells_AUC <- 
  AUCell_calcAUC(geneSets = temp_subsetgeneSets, 
                 rankings = temp_cells_rankings, 
                 aucMaxRank = ceiling(0.05 * nrow(temp_cells_rankings)))


#transpose matrix for seurat metadata assignment
temp_cells_AUC_matrix <- 
  t(as.data.frame(getAUC(temp_cells_AUC)))

# ADD TO SEURAT OBJECT  -----------------------------------------------------

temp_cells_AUC_matrix_sorted <- 
  temp_cells_AUC_matrix[rownames(seurat_10X@meta.data),,drop=FALSE]

temp_cells_AUC_matrix_sorted <- 
  as.data.frame.matrix(temp_cells_AUC_matrix_sorted)

saveRDS(temp_cells_AUC_matrix_sorted,
        "Output/Rdata/Rdata_cells_AUC_matrix_sorted.Rdata")

seurat_10X <- AddMetaData(seurat_10X, 
                          metadata = temp_cells_AUC_matrix_sorted)

saveRDS(seurat_10X,
        "Output/Rdata/Rdata_seurat_object_AUCell.Rdata")

# PLOTTING GENESETS -----------------------------------

temp_last_geneset <-
  (ncol(seurat_10X@meta.data))

temp_gene_set_names <-
  colnames(seurat_10X@meta.data[temp_first_geneset:temp_last_geneset])

# png function
temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 8, 
      height = 8, 
      res = 300, 
      units = 'in'
    )
  }


for(i in c(1:length(temp_gene_set_names))) {
  
  temp_gene_set_name <- 
    (temp_gene_set_names[i])
    
  temp_png_function(paste0("Output/Figures/0",
                           i,
                           "_featureplot_",
                           temp_gene_set_name,
                           ".png"))
  FeaturePlot(
    seurat_10X,
    features.plot = temp_gene_set_name,
    pt.size = 1,
    cols.use = c("light blue", "red3"), 
    no.axes = T, 
    min.cutoff = "q50",
    reduction.use = "umap_16_Stromal_Immune"
  )
  dev.off()
  
  temp_VlnPlot <- VlnPlot(
    seurat_10X,
    features.plot = temp_gene_set_name,
    group.by = "celltype",
    y.log = T,
    point.size.use = 0, 
    x.lab.rot=T
  )
  temp_png_function(paste0("Output/Figures/0",
                           i,
                           "_vlnplot_",
                           temp_gene_set_name,
                           ".png"))
  print(temp_VlnPlot)
  dev.off()
}

