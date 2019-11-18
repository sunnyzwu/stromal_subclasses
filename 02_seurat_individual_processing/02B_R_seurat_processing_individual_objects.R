# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 02B SEURAT V2 PROCESSING SCRIPT HPC
# R 3.5.0
# 
#
### PARAMS FROM JOB SUBMISSION (in order); 
#     - path to working directory
#     - sample ID
#     - species (human or mouse)
#     - path to cellranger output 
# 
### REQUIRES R v.3.5.0
### QSUB ARGUMENTS
#     qsub 
#     -cwd 
#     -pe smp 32 
#     -l h_vmem=200G 
#     -P TumourProgression 
#     -b y 
#     -j y 
#     -V 
#     -N s_sampleID
#     "R CMD BATCH 
#     --no-save 
#     '--args 
#     sampleID
#     human
#     /path_to_raw_matrix/GRCh38/'
#     /path/to/this/script.R" 
#
# 01 LOAD PACKAGES -----------------------------------------------------------

library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(cowplot)
library(tidyr)
library(eply)
library(readr)
library(ggplot2)
library(DropletUtils)

# 02 COMMAND LINE ARGUMENTS  ------------------------------------------------------------

temp_start_time <- date()
temp_args <-
  commandArgs(trailingOnly = T)
# project name
temp_project_name <- 
  temp_args[1]
# species type (human or mouse)
temp_species_type <- 
  temp_args[2]
# path to raw gene-barcode matrix files from cellranger
temp_raw_matrix <- 
  temp_args[3]

# 03 SET UP AND FUNCTIONS ------------------------------------------------------------------

#Sub-directory Outputs
dir.create("Output")
dir.create("Output/Figures")
dir.create("Output/Rdata")

# PNG function
temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 10, 
      height = 10, 
      res = 300, 
      units = 'in'
    )
  }

# 04 READ10X MATRIX ----------------------------------------------------

temp_seurat_raw <-
  Read10X(temp_raw_matrix)

# add unique barcode names for subsequent integration
colnames(x = temp_seurat_raw) <- 
  paste(temp_project_name, 
        colnames(x = temp_seurat_raw), 
        sep = '_')

seurat_10X <-  CreateSeuratObject(
  temp_seurat_raw,
  min.cells = 1,
  min.genes = 0,
  project = temp_project_name
)

# add mito content

temp_mito.genes <- grep(
  pattern = (if (temp_species_type == "human")
    "MT-"
    else
      "mt-"),
  x = rownames(x = seurat_10X@data),
  value = TRUE
)

temp_percent.mito <-
  Matrix::colSums(seurat_10X@raw.data[temp_mito.genes, ]) / Matrix::colSums(seurat_10X@raw.data)

seurat_10X <- AddMetaData(object = seurat_10X,
                          metadata = temp_percent.mito,
                          col.name = "percent.mito")


# 05 EMPTY DROPS TO ID REAL CELL BARCODES -----------------------------------------------------------

  set.seed(100)
  temp_emptyDrops_out <- emptyDrops(temp_seurat_raw,
                                    lower = 250)
  
  temp_is.cell <- 
    temp_emptyDrops_out$FDR <= 0.01
  
  sum(temp_is.cell, 
      na.rm=TRUE)
  
  temp_emptyDrops_out_subset <- 
    subset(temp_emptyDrops_out, 
           !is.na(FDR))
  
  temp_cell_ids <-
    rownames(temp_emptyDrops_out_subset[temp_emptyDrops_out_subset$FDR <= 0.01,])
  
  # emptydrops filtering
  seurat_10X <- 
    SubsetData(seurat_10X,
               cells.use = temp_cell_ids)
  
  # additional conservative gene filter
  seurat_10X <- FilterCells(
    object = seurat_10X,
    subset.names = c("nGene", 
                     "percent.mito"),
    low.thresholds = c(
      200,
      0
    ),
    high.thresholds = c(
      Inf,
      0.1
    )
  )

# 06 NORMALISATION AND VARIABLE GENES ---------------------------------------------------

#log normalisation
seurat_10X <-
  NormalizeData(
    object = seurat_10X,
    normalization.method = "LogNormalize",
    scale.factor = 10000, 
    display.progress = F
  )


seurat_10X <- FindVariableGenes(
  object = seurat_10X,
  mean.function = ExpMean,
  dispersion.function = LogVMR,
  x.low.cutoff = 0.05,
  x.high.cutoff = 8,
  y.cutoff = 0.5,
  y.high.cutoff = Inf, 
  do.plot = F
)

seurat_10X <- ScaleData(
  object = seurat_10X,
  vars.to.regress = c("nUMI", "nGene"),
  display.progress = F
)

# 07 COMPUTE PCA, TSNE and UMAP ------------------------------------------------------------
  
  # PCA
  seurat_10X <- RunPCA(
    object = seurat_10X,
    pc.genes = seurat_10X@var.genes,
    do.print = F,
    pcs.compute = 20,
    do.fast =  T
  )

  # TSNE
  seurat_10X <-
    RunTSNE(object = seurat_10X,
            dims.use =  1:20,
            do.fast = TRUE)
  # UMAP
  seurat_10X <- 
    RunUMAP(seurat_10X, 
            reduction.use = "pca", 
            dims.use = 1:20
    )
  
  # CLUSTERING
  seurat_10X <-
    FindClusters(
      object = seurat_10X,
      reduction.type = "pca",
      dims.use = 1:20,
      resolution = 0.8,
      print.output = F,
      save.SNN = T
    )
  
  # PLOT
  for(reduction in c("tsne","umap")){  
  
  temp_png_function(paste0("Output/Figures/01_",reduction,"_Plot.png"))
  DimPlot(
    object = seurat_10X,
    do.label = T,
    label.size = 4,
    pt.size = 1,
    group.by = "ident",
    reduction.use = reduction
  )
  dev.off()
}


# 08 CELL TYPE MARKERS -----------------------------------------------------------------

  temp_markers_to_use <- 
    c("ACTB",	# housekeeper
      "EPCAM",	"KRT8",	"KRT5", # epithelial
      "COL1A1",	"PDGFRA",	"PDGFRB",	"ACTA2",	"PECAM1",	# stromal
      "PTPRC",	"CD68",	"MS4A1","JCHAIN", "CD3D",	"CD8A", "FOXP3" # immune
      )
    for(reduction in c("tsne","umap")){  
    temp_png_function(paste0("Output/Figures/02_FeaturePlot_",reduction,".png"))
    FeaturePlot(
      seurat_10X,
      features.plot = temp_markers_to_use,
      pt.size = 1,
      cols.use = c("light blue", "red3"), 
      reduction.use = reduction
    )
    dev.off()
    }
  
# 09 SAVE PROCESSED OBJECT ---------------------------------------------------

saveRDS(seurat_10X,
        "Output/Rdata/seurat_object_processed.RData")
