# AUCell signature analysis
# inputs (in order); 
#     - working directory
#     - seurat_object
#     - gene_signatures
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
# outputs; 
#     - seurat_object with AUCell appended as Rdata
#     - plots of all AUcell signature scores overlayed to default TSNE generated in the seurat_object
#
# REQUIRES R v3.5.0
# 
# QSUB ONE LINER
# qsub -cwd -pe smp 16 -l h_vmem=100G -P TumourProgression -b y -j y -V -N AUCell "/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.5.0/bin/R CMD BATCH --no-save '--args /share/ScratchGeneral/sunwu/AUCell/CID4463_forDan/test_script_from_seurat_object/ ../CID4463.seurat_after_clustering.RData ../../gene_sets/CTP_BREAST_CANCER_signatures.gmx.csv' ./AUCell_signature_script.R" 
# 
# QSUB ARGUMENTS
# qsub 
# -cwd 
# -pe smp 16 
# -l h_vmem=100G 
# -P TumourProgression 
# -b y 
# -j y 
# -V 
# -N AUCell
# "/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.5.0/bin/R CMD BATCH 
# --no-save 
# '--args 
# /working/directory/path
# /path/to/seurat_object.Rdata 
# /path/to/gene/sets'
# /path/to/this/script.R" 
#
#
# LOAD PACKAGES -----------------------------------------------------------

if(!"Seurat" %in% installed.packages()){
  library(BiocInstaller)
  biocLite("Seurat")
}
if(!"readr" %in% installed.packages()){
  library(BiocInstaller)
  biocLite("readr")
}
if(!"GSEABase" %in% installed.packages()){
  library(BiocInstaller)
  biocLite("GSEABase")
}
if(!"AUCell" %in% installed.packages()){
  library(BiocInstaller)
  biocLite("AUCell")
}

if(!"reshape2" %in% installed.packages()){
  library(BiocInstaller)
  biocLite("reshape2")
}
if(!"NMF" %in% installed.packages()){
  library(BiocInstaller)
  biocLite("NMF")
}


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

# temp_wd 
temp_wd <- 
  temp_args[1]

# input path to seurat object
temp_seurat_object_path <- 
  temp_args[2]

# input gene sets
temp_genesets <- 
  temp_args[3]

# dimensional reduction to use for plotting
# UMAP_# e.g. UMAP_CC20
temp_dr_use <- 
  temp_args[4]
if(is.na(temp_dr_use)) {
  temp_dr_use <- 
    "tsne"
}

# threshold for top gene percentage to take in to account
temp_aucMaxRank <- 
  0.05

# FUNCTIONS ---------------------------------------------------------------

setwd(temp_wd)

# png function
temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 1189, 
      height = 519, 
      pointsize = 24
    )
  }

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


# if from multiple species alignment
rownames(temp_exprMatrix) <- paste(
  gsub("GRCh38_", "", rownames(temp_exprMatrix)))

temp_metadata_colnum <-
  ncol(seurat_10X@meta.data)

# INPUT GENE SETS ---------------------------------------------------------

temp_xcell_genesets <- 
  read_csv(temp_genesets)

temp_GeneSetCollection_list <- NULL

for(i in c(1:ncol(temp_xcell_genesets))) {
  n <- paste0("temp_set_name",
              i)
  assign(n,
         colnames(temp_xcell_genesets[i]))
  
  temp_set_name <- get(paste0("temp_set_name", 
                              i))
  
  temp_set <- na.omit(temp_xcell_genesets[i])
  
  colnames(temp_set) <- "gene_set"
  
  temp_set <- GeneSet(temp_set$gene_set,
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
                 aucMaxRank = ceiling(temp_aucMaxRank * nrow(temp_cells_rankings)))

save(temp_cells_AUC, 
     file="01_cells_AUC_aucMaxRank.RData")

#transpose matrix for seurat metadata assignment
temp_cells_AUC_matrix <- 
  t(as.data.frame(getAUC(temp_cells_AUC)))

saveRDS(temp_cells_AUC_matrix, 
        "02_cells_AUC_matrix_transposed.Rdata")


# ADD TO SEURAT OBJECT  -----------------------------------------------------

temp_cells_AUC_matrix_sorted <- 
  temp_cells_AUC_matrix[rownames(seurat_10X@meta.data),,drop=FALSE]

temp_cells_AUC_matrix_sorted <- 
  as.data.frame.matrix(temp_cells_AUC_matrix_sorted)

saveRDS(temp_cells_AUC_matrix_sorted,
        "03_cells_AUC_matrix_sorted.Rdata")

seurat_10X <- AddMetaData(seurat_10X, 
                          metadata = temp_cells_AUC_matrix_sorted)

saveRDS(seurat_10X,
        "04_seurat_object_AUCell.Rdata")

# PLOTTING GENESETS -----------------------------------

temp_first_geneset <- 
  (temp_metadata_colnum + 1)

temp_last_geneset <-
  (ncol(seurat_10X@meta.data))

temp_gene_set_names <-
  colnames(seurat_10X@meta.data[temp_first_geneset:temp_last_geneset])

dir.create("FeaturePlots_q50cutoff_AUCell_default_AUC_threshold")

for(i in c(1:length(temp_gene_set_names))) {
  
  temp_gene_set_name <- 
    (temp_gene_set_names[i])
    
  temp_png_function(paste0("FeaturePlots_q50cutoff_AUCell_default_AUC_threshold/FeaturePlot_AUCell_AUC_",
                           i,
                           "_",
                           temp_gene_set_name,
                           ".png"))
  FeaturePlot(
    seurat_10X,
    features.plot = temp_gene_set_name,
    pt.size = 2,
    cols.use = c("light blue", "red3"), 
    no.axes = T, 
    min.cutoff = "q50",
    reduction.use = temp_dr_use
  )
  dev.off()
}


# SAVE OBJECT -------------------------------------------------------------

saveRDS(seurat_10X,
        "seurat_object_processed_AUCell.Rdata")


