# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 03B CELL ANNOTATION AND SUBSET REPROCESSING SCRIPT SCRIPT HPC
# R 3.5.0
#
### inputs (in order); 
#     - project ID
#     - species (human or mouse)
#     - working directory
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
#     -N sCCA_sampleID
#     "R CMD BATCH 
#     --no-save 
#     '--args 
#     projectID
#     human'
#     /path/to/this/script.R" 
#
# 01 COMMAND LINE ARGUMENTS  ------------------------------------------------------------
# arguments from command line
temp_start_time <- date()
temp_args <-
  commandArgs(trailingOnly = T)

# project name
temp_project_name <- 
  temp_args[1]
# species type - "human" or "mouse" (must be lower case)
temp_species_type <- 
  temp_args[2]
# path to seurat objects
temp_objects_path <- 
  temp_args[3]

# sample ids
temp_sample_names <- 
  c("CID44041_P1", "CID44971_P2", "CID44991_P3", "CID4513_P4", "CID4515_P5")

# 02 LOAD PACKAGES -----------------------------------------------------------

library(Seurat)
library(plyr)
library(dplyr)
# library(Matrix)
# library(cowplot)
# library(tidyr)
# library(eply)

# 03 SET UP AND FUNCTIONS ------------------------------------------------------------------

# sub-directory outputs
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

# 04 LOAD DATA ------------------------------------------------------------

seurat_10X <- readRDS("Output/Rdata/seurat_CCA_aligned_processed.Rdata")


# 05 ANNOTATE DATASET -----------------------------------------------------

current.cluster.ids <- c(0:31)
new.cluster.ids <- c("CD8_T_Cells",
                     "Myeloid_1",
                     "Plasma_Cells",
                     "T_Regs",
                     "CD4_T_Cells",
                     "Myeloid_2",
                     "Epithelial_Basal",
                     "CD8_T_Cells",
                     "B_Cells",
                     "Fibroblasts_2",
                     "Epithelial_Basal",
                     "Myeloid_3",
                     "Fibroblasts_2",
                     "Epithelial_Basal_Cycling",
                     "Epithelial_Basal",
                     "Endothelial",
                     "CD8_T_Cells",
                     "Epithelial_Basal",
                     "CD8_T_Cells",
                     "Epithelial_Basal",
                     "Unassigned_1",
                     "Fibroblast_1",
                     "T_Cells_Cycling",
                     "Unassigned_2",
                     "Unassigned_3",
                     "Epithelial_Luminal_Mature",
                     "Myoepithelial",
                     "Plasma_Cells",
                     "Epithelial_Basal",
                     "Plasma_Cells",
                     "Myeloid_4")

seurat_10X@ident <- plyr::mapvalues(seurat_10X@ident, 
                                    from = current.cluster.ids, 
                                    to = new.cluster.ids)

temp_png_function(paste0("temp/TEMP_NEW_IDENTS.png"))
DimPlot(
  object = seurat_10X,
  do.label = T,
  label.size = 4,
  pt.size = 0.5,
  group.by = "ident",
  reduction.use = paste0("UMAP_CC",
                         15), 
  no.axes = T,
  plot.title = paste0("UMAP_CC",
                      15)
)
dev.off()

# add to metadata 
temp_celltype_df <- data.frame(row.names=row.names(seurat_10X@meta.data),
                               celltype=seurat_10X@ident)

seurat_10X <- AddMetaData(seurat_10X,
                          temp_celltype_df)



# 06 FILTER DATASET (Unassigned removal) ----------------------------------------------------------

temp_idents <- unique(seurat_10X@meta.data$celltype)

temp_idents <- temp_idents[!temp_idents %in% c("Unassigned_1", "Unassigned_2", "Unassigned_3")]

seurat_10X_subsetted <- SubsetData(seurat_10X,
                                   ident.use = temp_idents)

temp_png_function(paste0("temp/TEMP_NEW_NO_UNASSIGNED.png"))
DimPlot(
  object = seurat_10X_subsetted,
  do.label = T,
  label.size = 4,
  pt.size = 0.5,
  group.by = "ident",
  reduction.use = paste0("UMAP_CC",
                         15), 
  no.axes = T,
  plot.title = paste0("UMAP_CC",
                      15)
)
dev.off()

# FILTER DATASET (Manual doublet and clustering artifact removal removal) ----------------------------------------------------------

# Filtering EPCAM from all stromal and immune lineages
temp_cluster_ids <- unique(seurat_10X_subsetted@meta.data$celltype)
temp_original_barcodes <- row.names(seurat_10X_subsetted@meta.data)

temp_stromal_immune <- temp_cluster_ids[!grepl("Epithelial",temp_cluster_ids)]
temp_stromal_immune <- temp_stromal_immune[!grepl("Myoepithelial",temp_stromal_immune)]

temp_subset_to_filter <- SubsetData(seurat_10X_subsetted,
                                    ident.use = temp_stromal_immune)

temp_original_barcodes_subset <- row.names(temp_subset_to_filter@meta.data)
# number of cells original
print(length(temp_original_barcodes_subset))

for (i in c("EPCAM")) {
  temp_list_IDs <- WhichCells(
    temp_subset_to_filter,
    subset.name = i,
    accept.low = -Inf,
    accept.high = 0.1
  )
  
  temp_subset_to_filter <-
    SubsetData(temp_subset_to_filter,
               ident.use = NULL,
               cells.use = temp_list_IDs)
}

# number of cells remaining
temp_barcodes_filtered <- row.names(temp_subset_to_filter@meta.data)
print(length(temp_barcodes_filtered))
# number of cells filtered
temp_barcodes_to_filter <- temp_original_barcodes_subset[!temp_original_barcodes_subset %in% temp_barcodes_filtered]
print(length(temp_barcodes_to_filter))

# filter original dataset
temp_barcodes_to_keep <- temp_original_barcodes[!temp_original_barcodes %in% temp_barcodes_to_filter]
print(length(temp_barcodes_to_keep))

seurat_10X_subsetted_filtered <- SubsetData(seurat_10X_subsetted, 
                                            cells.use = temp_barcodes_to_keep)

# FILTER DATASET (Manual doublet and clustering artifact removal removal) ----------------------------------------------------------

# Filtering immune cells from epithelial and stromal lineages
temp_cluster_ids <- unique(seurat_10X_subsetted_filtered@meta.data$celltype)
temp_original_barcodes <- row.names(seurat_10X_subsetted_filtered@meta.data)

temp_non_immune <- as.vector(temp_cluster_ids[grep("Epithelial",temp_cluster_ids)])
temp_non_immune <- c(temp_non_immune, c("Myoepithelial",
                                        "Fibroblast_1",
                                        "Fibroblasts_2",
                                        "Endothelial"))

temp_subset_to_filter <- SubsetData(seurat_10X_subsetted_filtered,
                                    ident.use = temp_non_immune)

temp_original_barcodes_subset <- row.names(temp_subset_to_filter@meta.data)
# number of cells original
print(length(temp_original_barcodes_subset))

for (i in c("PTPRC", "CD3D", "CD3E", "CD3G", "CD19", "JCHAIN")) {
  temp_list_IDs <- WhichCells(
    temp_subset_to_filter,
    subset.name = i,
    accept.low = -Inf,
    accept.high = 0.1
  )
  
  temp_subset_to_filter <-
    SubsetData(temp_subset_to_filter,
               ident.use = NULL,
               cells.use = temp_list_IDs)
}

# number of cells remaining
temp_barcodes_filtered <- row.names(temp_subset_to_filter@meta.data)
print(length(temp_barcodes_filtered))
# number of cells filtered
temp_barcodes_to_filter <- temp_original_barcodes_subset[!temp_original_barcodes_subset %in% temp_barcodes_filtered]
print(length(temp_barcodes_to_filter))

# filter original dataset
temp_barcodes_to_keep <- temp_original_barcodes[!temp_original_barcodes %in% temp_barcodes_to_filter]
print(length(temp_barcodes_to_keep))

seurat_10X_subsetted_filtered <- SubsetData(seurat_10X_subsetted_filtered, 
                                            cells.use = temp_barcodes_to_keep)

