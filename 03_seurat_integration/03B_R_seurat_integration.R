# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 03B SEURAT V2 INTEGRATION SCRIPT HPC
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

# read Rdata for each input seurat object
for (samplename in temp_sample_names) {
  print(samplename)
  n <- paste0("temp_seurat_object_", samplename)
  assign(n,
         readRDS(paste0(temp_objects_path,
                        samplename,
                        "/Output/Rdata/seurat_object_processed.RData"))
         )
}

# additional conservative filters
for (samplename in temp_sample_names) {
  temp_seurat_object <- get(paste0("temp_seurat_object_", samplename))
  
  temp_seurat_object <- FilterCells(
    object = temp_seurat_object,
    subset.names = c("nGene", "percent.mito"),
    low.thresholds = c(200, 0.0),
    high.thresholds = c(Inf, 0.1))
  
  n <- paste0("temp_seurat_object_", samplename)
  assign(n,
         temp_seurat_object)
  
  rm(temp_seurat_object)
}


# 05 INTEGRATE DATA -----------------------------------------------------------------

  # make list of seurat objects 
  temp_sample_list <-
    mget(ls(pattern = "temp_seurat_object_*"))

  # variable genes for CCA
  ## take from pre-computed HVGs in each object
  temp_combined_genes <- c()
  for (i in c(1:length(temp_sample_list))) {
    temp_combined_genes <- c(temp_combined_genes,
                             head(
                               rownames(temp_sample_list[[i]]@hvg.info),
                               2000
                             ))
  }
  temp_combined_genes <- names(which(table(temp_combined_genes) > 1))
  ## filter for genes that are detected in every object
  for (i in c(1:length(temp_sample_list))) {
    temp_combined_genes <-
      temp_combined_genes[temp_combined_genes %in% rownames(temp_sample_list[[i]]@scale.data)]
  }

  # Run multiCCA
  temp_seurat_10X <- RunMultiCCA(
    object.list = temp_sample_list,
    temp_combined_genes,
    niter = 25,
    num.ccs = 20,
    standardize = TRUE
  )

  # CELL CYCLE REGRESSION
    # read cell cycle genes
    cc.genes <- c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1",
                  "UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7",
                  "POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2",
                  "USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8","HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2",
                  "TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2",
                  "CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","HJURP","CDCA3",
                  "HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23",
                  "HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA")
    s.genes <- cc.genes[1:43]
    g2m.genes <- cc.genes[44:98]
    
    temp_seurat_10X <- CellCycleScoring(temp_seurat_10X, 
                                        s.genes = s.genes, 
                                        g2m.genes = g2m.genes, 
                                        set.ident = T)
    
    temp_seurat_10X <- ScaleData(object = temp_seurat_10X, 
                                 vars.to.regress = c("S.Score", 
                                                     "G2M.Score"), 
                                 display.progress = F)
  
  
  # RUN ALIGNMENT
  temp_seurat_10X <- AlignSubspace(
    object = temp_seurat_10X,
    reduction.type = "cca",
    grouping.var = "orig.ident",
    verbose = T,
    dims.align = 1:20
  )

# 06 TSNE UMAP AND CLUSTERING ------------------------------------------------------

  temp_seurat_10X <-
  RunTSNE(
    object = temp_seurat_10X,
    reduction.use = "cca.aligned",
    dims.use = 1:20,
    do.fast = TRUE,
    reduction.name = paste0("tsne_",
                            20)
  )

  temp_seurat_10X <- 
    RunUMAP(temp_seurat_10X, 
            reduction.use = "cca.aligned", 
            dims.use = 1:20, 
            reduction.name = paste0("umap_",
                                    20)
    )
  
  temp_seurat_10X <-
    FindClusters(
      object = temp_seurat_10X,
      reduction.type = "cca.aligned",
      dims.use = 1:20,
      resolution = c(0.8,1.0,1.2),
      print.output = F,
      save.SNN = T
    )
  
# 07 DR PLOTS -------------------------------------------------------------------
  
  
  for(reduction in c("tsne","umap")){  
    
    for(res in c(c(0.8,1,1.2))){
      temp_png_function(paste0("Output/Figures/01_",reduction,"_Plot_res_",res,".png"))
      DimPlot(
        object = temp_seurat_10X,
        do.label = T,
        label.size = 4,
        pt.size = 1,
        group.by = paste0("res.",res),
        reduction.use = paste0(reduction,
                               "_",
                               20) 
      )
      dev.off()
    }
  }
  
  for(reduction in c("tsne","umap")){  
    
    temp_png_function(paste0("Output/Figures/02_",reduction,"_by_sampleID.png"))
    DimPlot(
      object = temp_seurat_10X,
      do.label = F,
      label.size = 4,
      pt.size = 1,
      group.by = "orig.ident",
      reduction.use = paste0(reduction,
                             "_",
                             20)  
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
    temp_png_function(paste0("Output/Figures/03_FeaturePlot_",reduction,".png"))
    FeaturePlot(
      temp_seurat_10X,
      features.plot = temp_markers_to_use,
      pt.size = 1,
      cols.use = c("light blue", "red3"), 
      reduction.use = paste0(reduction,
                             "_",
                             20)  
    )
    dev.off()
  }
  
  
# 09 SAVE CCA ALIGNED OBJECT --------------------------------------------------------------------

saveRDS(temp_seurat_10X,
        "Output/Rdata/seurat_CCA_aligned_processed.Rdata")
