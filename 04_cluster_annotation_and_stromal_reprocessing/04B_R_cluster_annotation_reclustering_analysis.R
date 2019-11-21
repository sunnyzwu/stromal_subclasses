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
dir.create("Output/Matrices_metadata")
dir.create("Output/Stromal_gene_signatures")


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
  
# 08 STROMAL RECLUSTERING: REALIGNMENT AND RECLUSTERING ----------------------------------------------

for(fibroblast in c("Fibroblasts_1","Fibroblasts_2")){
  print(fibroblast)
  
  temp_stromal_subset <- SubsetData(temp_subset_Stromal_Immune,
                               ident.use = fibroblast)
  # re-run CCA only the smaller FB2 cluster which showed no structure 
  if(fibroblast == "Fibroblasts_2"){

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
        num.ccs = 20,
        standardize = TRUE
      )
      
      temp_stromal_subset_integrated <- AlignSubspace(
        object = temp_stromal_subset_integrated,
        reduction.type = "cca",
        grouping.var = "orig.ident",
        verbose = T,
        dims.align = 1:20
      )

  }
  
  # new cell IDs
  if(fibroblast == "Fibroblasts_1"){
    temp_id <- "CAFs"
    temp_stromal_subset_integrated <- temp_stromal_subset
  }
  if(fibroblast == "Fibroblasts_2"){
    temp_id <- "VDSCs"
    rm(temp_sample_list)
    rm(list = ls(pattern = "temp_stromal_subset_orig.ident_"))
  }
  
  n <- paste0("temp_stromal_subset_integrated_",
              temp_id)
  assign(n, 
         temp_stromal_subset_integrated)
  
  rm(temp_stromal_subset)
  rm(temp_stromal_subset_integrated)
}

for(stromalsubset in c("CAFs", "VDSCs")){
  
  temp_stromal_subset_integrated <- get(paste0("temp_stromal_subset_integrated_",
                                               stromalsubset))
  if(stromalsubset == "CAFs"){
    temp_dims <- 20
    temp_dims_clustering <- 20
  }
  if(stromalsubset == "VDSCs"){
    temp_dims <- 4
    temp_dims_clustering <- 4
  }
  
  temp_stromal_subset_integrated <- RunTSNE(temp_stromal_subset_integrated,
                                            dims.use=1:temp_dims,
                                            reduction.use="cca.aligned")
  
  
  temp_stromal_subset_integrated <-
    FindClusters(
      object = temp_stromal_subset_integrated,
      reduction.type = "cca.aligned",
      dims.use = 1:temp_dims_clustering,
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


# 09 STROMAL RECLUSTERING: VISUALISATION AND GENE EXPRESSION --------------

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
    # temp_cols_use <- c("#f4a582","#ca0020")
    temp_pt_size <- 4
  }
  if(stromalsubset == "VDSCs"){
    # temp_cols_use <- c("#0571b0","#92c5de")
    temp_pt_size <- 5
  }
  
  temp_stromal_subset_integrated <- get(paste0("temp_stromal_subset_integrated_",
                                               stromalsubset))
  
  temp_png_function <-
    function(x) {
      png(
        file = (x), 
        width = 4, 
        height = 4, 
        res = 300, 
        units = 'in'
      )
    }
  
  for(res in c(0.2,0.3,0.4)){
    temp_png_function(paste0("Output/Figures/04_",stromalsubset,"_res.",res,".png"))
    DimPlot(
      object = temp_stromal_subset_integrated,
      do.label = T,
      label.size = 4,
      pt.size = temp_pt_size,
      group.by = paste0("res.",res),
      reduction.use = "tsne"
    )
    dev.off()
  }
  
  temp_png_function <-
    function(x) {
      png(
        file = (x), 
        width = 12, 
        height = 12, 
        res = 300, 
        units = 'in'
      )
    }
  
  temp_png_function(paste0("Output/Figures/05_featureplot_",stromalsubset,".png"))
  FeaturePlot(
    temp_stromal_subset_integrated,
    features.plot = temp_GOIs,
    pt.size = 2,
    reduction.use = "tsne", 
    do.return = T
  )
  dev.off()
  
}



# 10 ANNOTATION OF STROMAL/IMMUNE/EPITHELIAL CELLS AND COMBINE OBJECTS ---------------------------------------------

# annotate stromal cells
for(stromalsubset in c("CAFs", "VDSCs")){
  
  if(stromalsubset == "CAFs"){
    current.cluster.ids <- c(0,1)
    new.cluster.ids <- c("iCAFs", "myCAFs")
    res <-  0.2
  }
  if(stromalsubset == "VDSCs"){
    current.cluster.ids <- c(0,1)
    new.cluster.ids <- c("dVDSCs", "imVDSCs")
    res <-  0.2
  }
  
  temp_stromal_subset_integrated <- get(paste0("temp_stromal_subset_integrated_",
                                               stromalsubset))

  temp_stromal_subset_integrated <- SetAllIdent(temp_stromal_subset_integrated,
                                     id = paste0("res.",res))
  
  
  temp_stromal_subset_integrated@ident <- plyr::mapvalues(temp_stromal_subset_integrated@ident, 
                                               from = current.cluster.ids, 
                                               to = new.cluster.ids)
  
  temp_stromal_subset_integrated@meta.data$stromal_celltype <- 
    temp_stromal_subset_integrated@ident
  
  temp_png_function(paste0("Output/Figures/06_",stromalsubset,"_annotated.png"))
  DimPlot(
    object = temp_stromal_subset_integrated,
    do.label = T,
    label.size = 10,
    pt.size = temp_pt_size,
    group.by = "stromal_celltype",
    reduction.use = "tsne"
  )
  dev.off()
  
  n <- paste0("temp_stromal_subset_integrated_",
              stromalsubset)
  assign(n, 
         temp_stromal_subset_integrated)

  }

# combined objects
temp_merged_stromal <- MergeSeurat(temp_stromal_subset_integrated_CAFs,
                                   temp_stromal_subset_integrated_VDSCs)

temp_merged_stromal <- 
  SetAllIdent(temp_merged_stromal,
              id = "stromal_celltype")

# annotate stromal cells within complete stromal-immune object
temp_df <- data.frame(row.names = row.names(temp_subset_Stromal_Immune@meta.data),
                                  celltype = temp_subset_Stromal_Immune@meta.data$celltype)
temp_df <- 
  temp_df[!temp_df$celltype %in% c("Fibroblasts_1",
                                   "Fibroblasts_2"),,drop=F]

temp_df_stromal <- data.frame(row.names = row.names(temp_merged_stromal@meta.data),
                                  celltype = temp_merged_stromal@meta.data$stromal_celltype)

temp_df <- rbind(temp_df,
                 temp_df_stromal)

temp_df <- temp_df[row.names(temp_subset_Stromal_Immune@meta.data),,drop=F]

temp_subset_Stromal_Immune <- AddMetaData(temp_subset_Stromal_Immune, 
                                             metadata = temp_df)

# re-annotate stromal-immune cells
temp_subset_Stromal_Immune@meta.data$celltype <- factor(temp_subset_Stromal_Immune@meta.data$celltype,
                                                    levels=unique(temp_subset_Stromal_Immune@meta.data$celltype))

levels(temp_subset_Stromal_Immune@meta.data$celltype)[levels(temp_subset_Stromal_Immune@meta.data$celltype)=="T_Cells_Cycling"] <- "Cycling T-Cells"
levels(temp_subset_Stromal_Immune@meta.data$celltype)[levels(temp_subset_Stromal_Immune@meta.data$celltype)=="CD8_T_Cells"] <- "CD8+ T Cells"
levels(temp_subset_Stromal_Immune@meta.data$celltype)[levels(temp_subset_Stromal_Immune@meta.data$celltype)=="CD4_T_Cells"] <- "CD4+ T Cells"
levels(temp_subset_Stromal_Immune@meta.data$celltype)[levels(temp_subset_Stromal_Immune@meta.data$celltype)=="Plasma_Cells"] <- "Plasma Cells"
levels(temp_subset_Stromal_Immune@meta.data$celltype)[levels(temp_subset_Stromal_Immune@meta.data$celltype)=="T_Regs"] <- "T-Regs"
levels(temp_subset_Stromal_Immune@meta.data$celltype)[levels(temp_subset_Stromal_Immune@meta.data$celltype)=="B_Cells"] <- "B-Cells"
levels(temp_subset_Stromal_Immune@meta.data$celltype)[levels(temp_subset_Stromal_Immune@meta.data$celltype)=="Myeloid_4"] <- "Myeloid"
levels(temp_subset_Stromal_Immune@meta.data$celltype)[levels(temp_subset_Stromal_Immune@meta.data$celltype)=="Myeloid_3"] <- "Myeloid"
levels(temp_subset_Stromal_Immune@meta.data$celltype)[levels(temp_subset_Stromal_Immune@meta.data$celltype)=="Myeloid_2"] <- "Myeloid"
levels(temp_subset_Stromal_Immune@meta.data$celltype)[levels(temp_subset_Stromal_Immune@meta.data$celltype)=="Myeloid_1"] <- "Myeloid"

# re-annotate Epithelial cells
temp_subset_Epithelial@meta.data$celltype <- factor(temp_subset_Epithelial@meta.data$celltype,
                                             levels=unique(temp_subset_Epithelial@meta.data$celltype))
levels(temp_subset_Epithelial@meta.data$celltype)[levels(temp_subset_Epithelial@meta.data$celltype)=="Epithelial_Basal"] <- "Cancer Epithelial"
levels(temp_subset_Epithelial@meta.data$celltype)[levels(temp_subset_Epithelial@meta.data$celltype)=="Epithelial_Basal_Cycling"] <- "Cancer Proliferating"
levels(temp_subset_Epithelial@meta.data$celltype)[levels(temp_subset_Epithelial@meta.data$celltype)=="Epithelial_Luminal_Mature"] <- "Normal Epithelial Mature Luminal"
levels(temp_subset_Epithelial@meta.data$celltype)[levels(temp_subset_Epithelial@meta.data$celltype)=="Myoepithelial"] <- "Normal Epithelial Myoepithelial"

# merge stromal-immune and epithelial objects
temp_merged_all <- MergeSeurat(temp_subset_Stromal_Immune,
                               temp_subset_Epithelial)

# simpler patient annotations
temp_merged_all@meta.data$patientID <- NULL
temp_merged_all@meta.data$patientID[temp_merged_all@meta.data$orig.ident == "CID44041" ] <- "P1"
temp_merged_all@meta.data$patientID[temp_merged_all@meta.data$orig.ident == "CID44971" ] <- "P2"
temp_merged_all@meta.data$patientID[temp_merged_all@meta.data$orig.ident == "CID44991" ] <- "P3"
temp_merged_all@meta.data$patientID[temp_merged_all@meta.data$orig.ident == "CID4513" ] <- "P4"
temp_merged_all@meta.data$patientID[temp_merged_all@meta.data$orig.ident == "CID4515" ] <- "P5"

# 11 SAVE REPROCESSED OBJECTS ---------------------------------------------

# save individual objects
for(celltype in c("Epithelial", "Stromal_Immune")){
  temp_subset <- 
    get(paste0("temp_subset_",celltype))
  
  saveRDS(temp_subset,
          paste0("Output/Rdata/Rdata_",celltype,".png"))
  
}

for(stromalsubset in c("CAFs", "VDSCs")){
  
  temp_stromal_subset_integrated <- get(paste0("temp_stromal_subset_integrated_",
                                               stromalsubset))
  
  saveRDS(temp_stromal_subset_integrated,
          paste0("Output/Rdata/Rdata_",stromalsubset,".png"))
  
}


# save complete objects
saveRDS(temp_merged_all,
        paste0("Output/Rdata/Rdata_complete_dataset.png"))

saveRDS(temp_merged_stromal,
        paste0("Output/Rdata/Rdata_all_stromal_dataset.png"))


# 12 GENERATION OF STROMAL GENE SIGNATURES --------------------------------

temp_merged_stromal <- 
  SetAllIdent(temp_merged_stromal,
              id = "stromal_celltype")

temp_findallmarkers <- FindAllMarkers(temp_merged_stromal, 
                                      logfc.threshold = 0.1, 
                                      min.diff.pct = 0.1, 
                                      min.pct = 0.1,
                                      only.pos = T,
                                      test.use = 'MAST')

temp_findallmarkers <- dplyr::arrange(temp_findallmarkers,
                                      (cluster),
                                      desc(avg_logFC))
# number of DEGs per cluster
print(table(temp_findallmarkers$cluster))

write.csv(temp_findallmarkers,
          "Output/Stromal_gene_signatures/01_FindAllMarkers_stromal_celltype.csv", 
          quote = F, 
          row.names = T)

# 13 EXPORT CLUSTER AVERAGED LOG MATRICES ---------------------------------------------------------

# CLUSTER AVERAGED LOG NORMALISED EXPRESSION COMBINED STROMAL
print(unique(temp_merged_all@meta.data$celltype))
temp_merged_all <-
  SetAllIdent(temp_merged_all, 
              id = "celltype")

current.cluster.ids <- as.vector(unique(temp_merged_all@ident))
current.cluster.ids <- factor(current.cluster.ids,
                              levels = current.cluster.ids)
new.cluster.ids <- current.cluster.ids

levels(new.cluster.ids)[levels(new.cluster.ids)=="iCAFs"] <- "Stromal"
levels(new.cluster.ids)[levels(new.cluster.ids)=="myCAFs"] <- "Stromal"
levels(new.cluster.ids)[levels(new.cluster.ids)=="dVDSCs"] <- "Stromal"
levels(new.cluster.ids)[levels(new.cluster.ids)=="imVDSCs"] <- "Stromal"

temp_merged_all <-
  SetAllIdent(temp_merged_all, 
              id = "celltype")

temp_merged_all@ident <- plyr::mapvalues(temp_merged_all@ident, 
                                    from = as.vector(current.cluster.ids), 
                                    to = as.vector(new.cluster.ids)
                                    )

temp_merged_all@meta.data$celltype <-  temp_merged_all@ident

temp_cluster_averages <- 
  AverageExpression(temp_merged_all,
                    return.seurat = T,
                    show.progress = T)

temp_data_frame <- 
  as.matrix(temp_cluster_averages@data)

avg_log_df <- 
  as.data.frame(temp_data_frame)

write.csv(avg_log_df,
          "Output/Stromal_gene_signatures/ALL_CELLTYPES_merged_expression_values_log.csv",
          row.names = T)

# CLUSTER AVERAGED LOG NORMALISED EXPRESSION ONLY STROMAL
temp_merged_stromal <- 
  SetAllIdent(temp_merged_stromal,
              id = "stromal_celltype")

print(unique(temp_merged_stromal@meta.data$stromal_celltype))

current.cluster.ids <- as.vector(unique(temp_merged_stromal@ident))
current.cluster.ids <- factor(current.cluster.ids,
                              levels = current.cluster.ids)
new.cluster.ids <- current.cluster.ids

temp_merged_stromal@ident <- plyr::mapvalues(temp_merged_stromal@ident, 
                                         from = as.vector(current.cluster.ids), 
                                         to = as.vector(new.cluster.ids)
)

temp_merged_stromal@meta.data$celltype <-  temp_merged_stromal@ident

temp_cluster_averages <- 
  AverageExpression(temp_merged_stromal,
                    return.seurat = T,
                    show.progress = T)

temp_data_frame <- 
  as.matrix(temp_cluster_averages@data)

avg_log_df <- 
  as.data.frame(temp_data_frame)

write.csv(avg_log_df,
          "Output/Stromal_gene_signatures/ALL_STROMAL_merged_expression_values_log.csv",
          row.names = T)


