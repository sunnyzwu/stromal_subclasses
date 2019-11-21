# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 04B CELL ANNOTATION AND RECLUSTERING SCRIPT HPC
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
# temp_sample_names <- 
#   c("CID44041_P1", "CID44971_P2", "CID44991_P3", "CID4513_P4", "CID4515_P5")

# 02 LOAD PACKAGES -----------------------------------------------------------

library(Seurat)
library(plyr)
library(dplyr)

# 03 SET UP AND FUNCTIONS ------------------------------------------------------------------

# sub-directory outputs
dir.create("Output")
dir.create("Output/Figures")
dir.create("Output/Rdata")

# 04 LOAD DATA ------------------------------------------------------------

seurat_10X <- readRDS(paste0(temp_objects_path,
                             "Output/Rdata/seurat_CCA_aligned_processed.Rdata")
                      )

# load cell identities
temp_cell_ids <- 
  read.csv("cell_ids/01_cell_ids.csv", 
           row.names = "X")
# order cell barcodes
temp_cell_ids <- 
  temp_cell_ids[rownames(seurat_10X@meta.data),,drop=F]

# 05 LOAD MANUAL CELL IDs -----------------------------------------------------

  # check cell barcodes match perfectly
print(all.equal(rownames(seurat_10X@meta.data), 
                row.names(temp_cell_ids)))

seurat_10X <- AddMetaData(seurat_10X,
                          temp_cell_ids)

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

temp_png_function(paste0("Output/Figures/01_unfiltered_cell_ids.png"))
DimPlot(
  object = seurat_10X,
  do.label = T,
  label.size = 4,
  pt.size = 0.1,
  group.by = "celltype",
  reduction.use = "umap_20" 
)
dev.off()

# 06 ADDITIONAL FILTERING (UNASSIGNED - DOUBLETS) -------------------------------------------------------

  # filter small unassigned clusters
seurat_10X <- SetAllIdent(seurat_10X, 
                          "celltype")
temp_idents <- 
  unique(seurat_10X@meta.data$celltype)

temp_idents <- 
  temp_idents[!temp_idents %in% c("Unassigned_1", "Unassigned_2", "Unassigned_3")]

seurat_10X <- SubsetData(seurat_10X, 
                         ident.use = temp_idents)

  # Filtering EPCAM from all stromal and immune lineages
temp_cluster_ids <- 
  unique(seurat_10X@meta.data$celltype)
temp_original_barcodes <- 
  row.names(seurat_10X@meta.data)

temp_stromal_immune <- 
  temp_cluster_ids[!grepl("Epithelial",temp_cluster_ids)]
temp_stromal_immune <- 
  temp_stromal_immune[!grepl("Myoepithelial",temp_stromal_immune)]

temp_subset_to_filter <- SubsetData(seurat_10X,
                                    ident.use = temp_stromal_immune)

temp_original_barcodes_subset <- 
  row.names(temp_subset_to_filter@meta.data)

for (i in c("EPCAM")) {
  print(i)
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

temp_barcodes_filtered <- 
  row.names(temp_subset_to_filter@meta.data)

temp_barcodes_to_filter <- 
  temp_original_barcodes_subset[!temp_original_barcodes_subset %in% temp_barcodes_filtered]

temp_barcodes_to_keep <- 
  temp_original_barcodes[!temp_original_barcodes %in% temp_barcodes_to_filter]

seurat_10X_filtered <- SubsetData(seurat_10X,
                                  cells.use = temp_barcodes_to_keep)


  # Filtering immune cells from epithelial and stromal lineages
temp_cluster_ids <- 
  unique(seurat_10X_filtered@meta.data$celltype)
temp_original_barcodes <- 
  row.names(seurat_10X_filtered@meta.data)

temp_non_immune <- as.vector(temp_cluster_ids[grep("Epithelial",temp_cluster_ids)])
temp_non_immune <- c(temp_non_immune, c("Myoepithelial",
                                        "Fibroblasts_1",
                                        "Fibroblasts_2",
                                        "Endothelial"))

temp_subset_to_filter <- SubsetData(seurat_10X_filtered,
                                    ident.use = temp_non_immune)

temp_original_barcodes_subset <- row.names(temp_subset_to_filter@meta.data)

for (i in c("PTPRC", "CD3D", "CD3E", "CD3G", "CD19", "JCHAIN")) {
  print(i)
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

temp_barcodes_filtered <- 
  row.names(temp_subset_to_filter@meta.data)

temp_barcodes_to_filter <- 
  temp_original_barcodes_subset[!temp_original_barcodes_subset %in% temp_barcodes_filtered]

temp_barcodes_to_keep <- 
  temp_original_barcodes[!temp_original_barcodes %in% temp_barcodes_to_filter]

seurat_10X_filtered_2 <- SubsetData(seurat_10X_filtered,
                                  cells.use = temp_barcodes_to_keep)


  # Filtering epithelial and immune cells from fibroblast lineages
temp_cluster_ids <- unique(seurat_10X_filtered_2@meta.data$celltype)
temp_original_barcodes <- row.names(seurat_10X_filtered_2@meta.data)

temp_fibroblast <- c("Fibroblasts_1", "Fibroblasts_2")

temp_subset_to_filter <- SubsetData(seurat_10X_filtered_2,
                                    ident.use = temp_fibroblast)

temp_original_barcodes_subset <- row.names(temp_subset_to_filter@meta.data)

temp_genes_to_filter <- c("EPCAM",
                          "KRT5",
                          "KRT14",
                          "KRT8",
                          "KRT18",
                          "PECAM1",
                          "PTPRC",
                          "CD3D",
                          "CD3G",
                          "CD3E",
                          "CD19",
                          "JCHAIN")

for (i in temp_genes_to_filter) {
  print(i)
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

temp_barcodes_filtered <- 
  row.names(temp_subset_to_filter@meta.data)

temp_barcodes_to_filter <- 
  temp_original_barcodes_subset[!temp_original_barcodes_subset %in% temp_barcodes_filtered]

temp_barcodes_to_keep <- 
  temp_original_barcodes[!temp_original_barcodes %in% temp_barcodes_to_filter]

seurat_10X_filtered_3 <- SubsetData(seurat_10X_filtered_2,
                                    cells.use = temp_barcodes_to_keep)


  # plot filter object
temp_png_function(paste0("Output/Figures/02_filtered_cell_ids.png"))
DimPlot(
  object = seurat_10X_filtered_3,
  do.label = T,
  label.size = 4,
  pt.size = 0.1,
  group.by = "celltype",
  reduction.use = "umap_20" 
)
dev.off()

# 07 EPITHELIAL STROMAL SEPERATION AND PROCESSING -------------------------

# reprocess
temp_cluster_ids <- unique(seurat_10X_filtered_3@meta.data$celltype)
for(celltype in c("Epithelial", "Stromal_Immune")){
  print(celltype)
  if(celltype == "Epithelial"){
    temp_ids <- as.vector(temp_cluster_ids[grep("Epithelial",temp_cluster_ids)])
    temp_ids <- c(temp_ids, c("Myoepithelial"))
    temp_dims <- 15
  }
  if(celltype == "Stromal_Immune"){
    temp_ids <- temp_cluster_ids[!grepl("Epithelial",temp_cluster_ids)]
    temp_ids <- temp_ids[!grepl("Myoepithelial",temp_ids)]
    temp_dims <- 16
  }
  
  temp_subset <- SubsetData(seurat_10X_filtered_3,
                                   ident.use = temp_ids)
  
  temp_subset <- 
    RunUMAP(temp_subset, 
            reduction.use = "cca.aligned", 
            dims.use = 1:temp_dims, 
            reduction.name = paste0("umap_",temp_dims,"_",celltype)
    )
  
  n <- paste0("temp_subset_",celltype)
  assign(n, temp_subset)
  

  }

# plot
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

for(celltype in c("Epithelial", "Stromal_Immune")){
  temp_subset <- 
    get(paste0("temp_subset_",celltype))
  
  if(celltype == "Epithelial"){
    temp_dims <- 15
  }
  if(celltype == "Stromal_Immune"){
    temp_dims <- 16
  }
  temp_png_function(paste0("Output/Figures/03_",celltype,".png"))
  DimPlot(
    object = temp_subset,
    do.label = T,
    label.size = 4,
    pt.size = 0.1,
    group.by = "celltype",
    reduction.use = paste0("umap_",temp_dims,"_",celltype)
  )
  dev.off()
  }
  
# 08 SAVE REPROCESSED OBJECTS ---------------------------------------------

for(celltype in c("Epithelial", "Stromal_Immune")){
  temp_subset <- 
    get(paste0("temp_subset_",celltype))

  saveRDS(temp_subset,
          paste0("Output/Rdata/Rdata_",celltype,".png"))
  
  }

# 09 STROMAL RECLUSTERING: REALIGNMENT ----------------------------------------------

for(fibroblast in c("Fibroblasts_1","Fibroblasts_2")){
  print(fibroblast)
  
  temp_stromal_subset <- SubsetData(temp_subset_Stromal_Immune,
                               ident.use = fibroblast)

  temp_stromal_subset <- SetAllIdent(temp_stromal_subset, 
                                id = "orig.ident")
  
  temp_orig_ident <- unique(temp_stromal_subset@meta.data$orig.ident)
  
  temp_combined_genes <- c()
  for(i in temp_orig_ident){
    
    temp_stromal_subset_orig.ident <- SubsetData(temp_stromal_subset, 
                              ident.use = i)
    
    temp_stromal_subset_orig.ident <- FindVariableGenes(
      object = temp_stromal_subset_orig.ident,
      mean.function = ExpMean,
      dispersion.function = LogVMR,
      x.low.cutoff = 0.05,
      x.high.cutoff = 8,
      y.cutoff = 0.5,
      y.high.cutoff = Inf, 
      do.plot = F
    )
    
    temp_combined_genes <- c(temp_combined_genes,
                             head(
                               rownames(temp_stromal_subset_orig.ident@hvg.info),
                               1000
                             ))
    
    n <- paste0("temp_stromal_subset_orig.ident_",
                i)
    
    assign(n, 
           temp_stromal_subset_orig.ident)
    rm(temp_stromal_subset_orig.ident)
  }
  
  # make list of seurat objects for >2 sample processing
  temp_sample_list <-
    mget(ls(pattern = "temp_stromal_subset_orig.ident_*"))
  
  temp_combined_genes_unique <- names(which(table(temp_combined_genes) > 1))
  
  # conservative approach to remove all ribosomal & mitochondrial genes from integration to avoid artefacts
  temp_combined_genes_unique <- temp_combined_genes_unique[!grepl("RPL",temp_combined_genes_unique)]
  temp_combined_genes_unique <- temp_combined_genes_unique[!grepl("MT-",temp_combined_genes_unique)]
  temp_combined_genes_unique <- temp_combined_genes_unique[!grepl("RPS",temp_combined_genes_unique)]
  
  temp_stromal_subset_integrated <- RunMultiCCA(
    object.list = temp_sample_list,
    genes.use = temp_combined_genes_unique,
    niter = 25,
    num.ccs = 10,
    standardize = TRUE
  )
  
  temp_stromal_subset_integrated <- AlignSubspace(
    object = temp_stromal_subset_integrated,
    reduction.type = "cca",
    grouping.var = "orig.ident",
    verbose = T,
    dims.align = 1:10
  )

  if(fibroblast == "Fibroblasts_1"){
    temp_id <- "CAFs"
  }
  if(fibroblast == "Fibroblasts_2"){
    temp_id <- "VDSCs"
  }
  
  n <- paste0("temp_stromal_subset_integrated_",
              temp_id)
  assign(n, 
         temp_stromal_subset_integrated)
  
  rm(list = ls(pattern = "temp_stromal_subset_orig.ident_"))
  rm(temp_stromal_subset)
  rm(temp_sample_list)
  rm(temp_stromal_subset_integrated)
}

# 10 STROMAL RECLUSTERING: RERUN DIMENSIONAL REDUCTION ------------------------------------------

for(stromalsubset in c("CAFs", "VDSCs")){
  
  temp_stromal_subset_integrated <- get(paste0("temp_stromal_subset_integrated_",
                                               stromalsubset))
  if(stromalsubset == "CAFs"){
    temp_dims <- 20
  }
  if(stromalsubset == "VDSCs"){
    temp_dims <- 4
  }
  
  temp_stromal_subset_integrated <- RunTSNE(temp_stromal_subset_integrated,
                                            dims.use=1:temp_dims,
                                            reduction.use="cca.aligned")
  
  
  temp_stromal_subset_integrated <-
    FindClusters(
      object = temp_stromal_subset_integrated,
      reduction.type = "cca.aligned",
      dims.use = 1:temp_dims_use,
      resolution = c(0.2,0.3,0.4),
      print.output = F,
      save.SNN = T,
      force.recalc = T
    )
  
  n <- paste0("temp_stromal_subset_integrated_",
              stromalsubset)
  assign(n, 
         temp_stromal_subset_integrated)
}


# 11 STROMAL RECLUSTERING: VISUALISATION AND GENE EXPRESSION --------------

temp_GOIs <-
  c("PDGFRB","THY1","S100A4", "ITGB1",
    "PDGFRA","COL1A1", "PDPN", "FAP", 
    "CD34", "CXCL12",
    "ACTA2", "MCAM", "CAV1",
    "TAGLN", "MYH11", "MYLK",
    "CD36", "RGS5")

# TSNE

for(stromalsubset in c("CAFs", "VDSCs")){
  
  if(stromalsubset == "CAFs"){
    temp_cols_use <- c("#f4a582","#ca0020")
    temp_pt_size <- 4
  }
  if(stromalsubset == "VDSCs"){
    temp_cols_use <- c("#0571b0","#92c5de")
    temp_pt_size <- 5
  }
  
  temp_stromal_subset_integrated <- get(paste0("temp_stromal_subset_integrated_",
                                               stromalsubset))
  
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
  
  for(res in c(0.2,0.3,0.4)){
    temp_png_function(paste0("Output/Figures/04_",stromalsubset,"_",res,".png"))
    DimPlot(
      object = temp_stromal_subset_integrated,
      do.label = T,
      label.size = 4,
      pt.size = 0.1,
      group.by = paste0("res.",res),
      reduction.use = "tsne",
      cols = temp_cols_use
    )
    dev.off()
  }
  
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
  
  temp_featureplot <- FeaturePlot(
    temp_stromal_subset_integrated,
    features.plot = temp_GOIs,
    pt.size = 0.1,
    reduction.use = "tsne"
  )
  temp_png_function(paste0("Output/Figures/04_",stromalsubset,"_",res,".png"))
  print(temp_featureplot)
  dev.off()
  
}



# 12 ANNOTATION OF STROMAL CELLS ---------------------------------------------


# 13 GENERATION OF STROMAL GENE SIGNATURES --------------------------------



# 14 EXPORT MATRICES ---------------------------------------------------------


# 14 SAVE REPROCESSED OBJECTS ---------------------------------------------

for(stromalsubset in c("CAFs", "VDSCs")){
  
  temp_stromal_subset_integrated <- get(paste0("temp_stromal_subset_integrated_",
                                               stromalsubset))
  
  saveRDS(temp_stromal_subset_integrated,
          paste0("Output/Rdata/Rdata_",stromalsubset,".png"))
  
}

