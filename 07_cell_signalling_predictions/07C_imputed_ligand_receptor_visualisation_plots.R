# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 07C LIGAND RECEPTOR VISUALISATION OF IMPUTED GENE EXPRESSION PLOTS
# R 3.5.0
#
### inputs (in order); 
#     - seurat object
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
# --no-save \
# ${SCRIPT}"
#
# 01 LOAD PACKAGES -----------------------------------------------------------

library(Rmagic)
library(ggplot2)
library(Seurat)

library(readr)
library(ggplot2)
library(viridis)

library(dplyr)
library(stringr)


# 02 COMMAND LINE ARGUMENTS  ------------------------------------------------------------

temp_start_time <- date()
temp_args <-
  commandArgs(trailingOnly = T)

# input path to seurat object
temp_seurat_object_path <- 
  temp_args[1]

# 03 SET UP AND FUNCTIONS ---------------------------------------------------------------

dir.create("Output")
dir.create("Output/Figures")
dir.create("Output/Rdata")

# PNG function
temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 12, 
      height = 8, 
      res = 300, 
      units = 'in'
    )
  }

# 04 LOAD DATA ---------------------------------------------------------------

seurat_10X_original <-
  readRDS(temp_seurat_object_path)

# 05 RUN MAGIC ---------------------------------------------------------------

data <- as.matrix(t(seurat_10X_original@data))

MAGIC_data_allgenes <- magic(data, 
                             genes="all_genes")

seurat_10X_imputed <- seurat_10X_original

seurat_10X_imputed@data <- t(MAGIC_data_allgenes$result)

# 06 SUBSET OBJECTS FOR FINAL PLOT -------------------------------------------

seurat_10X_imputed <- SetAllIdent(seurat_10X_imputed, 
                                  "celltype")
temp_target_ids <- unique(seurat_10X_imputed@meta.data$celltype)

seurat_10X_imputed_CAFs <- SubsetData(seurat_10X_imputed,
                                      ident.use =  
                                        c("myCAFs","iCAFs", "dVDSCs","imVDSCs"))

seurat_10X_imputed_targets <- SubsetData(seurat_10X_imputed,
                                         ident.use = temp_target_ids[!temp_target_ids %in% c("myCAFs","iCAFs", "dVDSCs","imVDSCs")])

# 07 IMPUTED LIGAND RECEPTOR PLOTS: CANCER --------------------------------

# cancer interactions 
{temp_interactions_of_interest <-
    c("FGF7_FGFR2", "FGF10_FGFR1", "BMP4_BMPR1A", "BMP7_BMPR1B", "HGF_MET")
  
  # grep("HGF",rownames(seurat_10X@data),value=T)
  
  for(i in temp_interactions_of_interest) {
    
    temp_gene_table <- 
      str_split_fixed(i, "_", 2)  
    
    # ligand
    for(gene in c(temp_gene_table[1])) {
      
      tempVlnPLot <- VlnPlot(seurat_10X_imputed_CAFs,
                             features.plot = gene, 
                             group.by = "celltype",
                             y.log = F, 
                             point.size.use = 0,
                             x.lab.rot = F,
                             cols = rev(c("#ca0020","#f4a582","#92c5de","#0571b0")))
      
      tempVlnPLot <- tempVlnPLot + theme(
        plot.title = element_blank(),
        axis.text.x=element_text(size = 14,
                                 angle= 45,
                                 vjust= .8),
        plot.margin = unit(c(1.5, 0, 0, 0.25), "cm") # T, r, b, l
      ) + xlab(" ") + ylab(" ") + coord_flip() #+ geom_boxplot(outlier.size=0.1)
      
      n <- paste0("temp_vlnplot_",gene)
      assign(n, tempVlnPLot)
    }
    
    # receptor
    for(gene in c(temp_gene_table[2])) {
      
      tempVlnPLot <- VlnPlot(seurat_10X_imputed_targets,
                             features.plot = gene, 
                             group.by = "celltype",
                             y.log = F, 
                             point.size.use = 0,
                             x.lab.rot = F)
      
      tempVlnPLot <- tempVlnPLot + theme(
        plot.title = element_blank(),
        axis.text.x=element_text(size = 14,
                                 angle= 45,
                                 vjust= .8),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm") # T, r, b, l
      ) +
        xlab(" ") + ylab(" ") + coord_flip() #+ geom_boxplot(outlier.size=0.1)
      
      n <- paste0("temp_vlnplot_",gene)
      assign(n, tempVlnPLot)
    }
    
    temp_grid <- 
      plot_grid(plotlist=mget(paste0("temp_vlnplot_",c(temp_gene_table[1], temp_gene_table[2]))),
                labels = NULL,
                nrow = 1,
                ncol = 2)
    
    n <- paste0("temp_grid_", i)
    assign(n, temp_grid)
    
  }
  
  temp_grid_final <- 
    plot_grid(plotlist=mget(paste0("temp_grid_",temp_interactions_of_interest)),
              labels = paste0(gsub("_", "-", temp_interactions_of_interest)),
              label_size = 20,
              label_y = 1,
              label_x = .1,
              hjust = 0,
              nrow = 5,
              ncol = 1)
  
  # fucking labels wont align properly if you use hjust (default 0.5) as this is proportional to the physical width of each label. 
  # Just set label_y and label_x to the same coordinates for alignment.
  # https://github.com/wilkelab/cowplot/issues/32
  
  temp_png_function <-
    function(x) {
      png(
        file = (x), 
        width = 6, 
        height = 13, 
        res = 300, 
        units = 'in'
      )
    }
  
  temp_png_function(paste0("Output/Figures/03_plots_cancer.png"))
  print(temp_grid_final)
  dev.off()
  
  rm(list = ls(pattern = "temp_grid_"))
}

# 08 IMPUTED LIGAND RECEPTOR PLOTS: IMMUNE --------------------------------

# immune interactions
{temp_interactions_of_interest <-
    c("C5_C5AR1", "TGFB1_TGFBR1", "TGFB2_TGFBR1", "CCL8_CCR1", "IL6_IL6R", "CCL2_CCR1", # myeloid interactions
      "CXCL12_CXCR4", "CXCL13_CXCR5", "CXCL11_CXCR3","CD274_PDCD1", # T-cell interactions
      "CXCL9_CXCR3", "CCL21_CCR7"
      
    )
  
  library(stringr)
  for(i in temp_interactions_of_interest) {
    
    temp_gene_table <- 
      str_split_fixed(i, "_", 2)  
    
    # ligand
    for(gene in c(temp_gene_table[1])) {
      
      tempVlnPLot <- VlnPlot(seurat_10X_imputed_CAFs,
                             features.plot = gene, 
                             group.by = "celltype",
                             y.log = F, 
                             point.size.use = 0,
                             x.lab.rot = F,
                             cols = rev(c("#ca0020","#f4a582","#92c5de","#0571b0")))
      
      tempVlnPLot <- tempVlnPLot + theme(
        plot.title = element_blank(),
        axis.text.x=element_text(size = 14,
                                 angle= 45,
                                 vjust= .8),
        plot.margin = unit(c(1.5, 0, 0, 0.25), "cm") # T, r, b, l
      ) + xlab(" ") + ylab(" ") + coord_flip() #+ geom_boxplot(outlier.size=0.1)
      
      n <- paste0("temp_vlnplot_",gene)
      assign(n, tempVlnPLot)
    }
    
    # receptor
    for(gene in c(temp_gene_table[2])) {
      
      tempVlnPLot <- VlnPlot(seurat_10X_imputed_targets,
                             features.plot = gene, 
                             group.by = "celltype",
                             y.log = F, 
                             point.size.use = 0,
                             x.lab.rot = F)
      
      tempVlnPLot <- tempVlnPLot + theme(
        plot.title = element_blank(),
        axis.text.x=element_text(size = 14,
                                 angle= 45,
                                 vjust= .8),
        plot.margin = unit(c(0.5, 0, 0, 0), "cm") # T, r, b, l
      ) +
        xlab(" ") + ylab(" ") + coord_flip() #+ geom_boxplot(outlier.size=0.1)
      
      n <- paste0("temp_vlnplot_",gene)
      assign(n, tempVlnPLot)
    }
    
    temp_grid <- 
      plot_grid(plotlist=mget(paste0("temp_vlnplot_",c(temp_gene_table[1], temp_gene_table[2]))),
                labels = NULL,
                nrow = 1,
                ncol = 2)
    
    n <- paste0("temp_grid_", i)
    assign(n, temp_grid)
    
  }
  
  temp_grid_final <- 
    plot_grid(plotlist=mget(paste0("temp_grid_",temp_interactions_of_interest)),
              labels = paste0(gsub("_", "-", temp_interactions_of_interest)),
              label_size = 20,
              label_y = 1,
              label_x = .1,
              hjust = 0,
              nrow = 6,
              ncol = 2)
  
  temp_png_function <-
    function(x) {
      png(
        file = (x), 
        width = 12, 
        height = 13, 
        res = 300, 
        units = 'in'
      )
    }
  
  temp_png_function(paste0("Output/Figures/04_plots_immune.png"))
  print(temp_grid_final)
  dev.off()
}



# 09 SAVE OBJECT -------------------------------------------------------------

saveRDS(seurat_10X_imputed, 
        "Output/Rdata/Rdata_Stromal_Immune.Rdata_MAGIC.Rdata")
