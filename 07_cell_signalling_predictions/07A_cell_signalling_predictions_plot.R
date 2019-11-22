# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 07A CELL SIGNALLING PREDICTIONS AND VISUALISATION
# R 3.5.0
#
# w
# 01 LOAD PACKAGES -----------------------------------------------------------

library(circlize)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(plyr)
library(scales)


# 02 COMMAND LINE ARGUMENTS  ------------------------------------------------------------

# no command line arguments needed

# 03 SET UP AND FUNCTIONS ---------------------------------------------------------------

dir.create("Output")
dir.create("Output/Figures")

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

temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 10,
      height = 10,
      useDingbats=F
    )
  }


# 04 LOAD DATA ---------------------------------------------------------------

ad <-
  read.csv(
    paste0("data/cluster_interactions_signalType_paracrine_weightType_mean_specificityThreshold_0.1.csv"))

# ligands
temp_df_1 <- 
  read.csv("data/meanClusteredEM_paracrine_combined_ligand.csv")

# receptors
temp_df_2 <- 
  read.csv("data/meanClusteredEM_paracrine_combined_receptor.csv")

# new stromal IDs from raw inputs    
temp_CAF1A_short <- "myCAFs"
temp_CAF1B_short <- "iCAFs" 
temp_CAF2A_short <- "dVDSCs"
temp_CAF2B_short <- "imVDSCs" 

# 05 C2CCN CIRCOS PLOT --------------------------------------------

temp_cutoff <- 0.1

temp_backround_arrows <- "light grey"
grid.col = c(
  "B-cells" = temp_backround_arrows,
  "myCAFs" = "#ca0020",
  "iCAFs" = "#f4a582",
  "dVDSCs" = "#0571b0",
  "imVDSCs" = "#92c5de",
  "CD4+ T-cells" = temp_backround_arrows,
  "CD8+ T-cells" = temp_backround_arrows,
  "Endothelial" = temp_backround_arrows,
  "Cancer Epithelial" = temp_backround_arrows,
  "Proliferating Cancer" = temp_backround_arrows,
  "Epithelial_Luminal_Mature" = temp_backround_arrows,
  "Myeloid" = temp_backround_arrows,
  "Myoepithelial" = temp_backround_arrows,
  "Plasma_Cells" = temp_backround_arrows,
  "T-Regs" = temp_backround_arrows,
  "Proliferating T-cells" = temp_backround_arrows
)
order = c(
  'myCAFs',
  'iCAFs',
  'dVDSCs',
  'imVDSCs',
  'CD4+ T-cells',
  'CD8+ T-cells',
  'Proliferating T-cells',
  'T-Regs',
  'B-cells',
  'Myeloid',
  'Endothelial',
  'Cancer Epithelial',
  'Proliferating Cancer',
  'Epithelial_Luminal_Mature',
  'Myoepithelial',
  'Plasma_Cells'
)

  ad <- ad[! ad$sending.cluster.name %in% "Epithelial_Luminal_Mature",]
  ad <- ad[! ad$target.cluster.name %in% "Epithelial_Luminal_Mature",]
  ad <- ad[! ad$sending.cluster.name %in% "Myoepithelial",]
  ad <- ad[! ad$target.cluster.name %in% "Myoepithelial",]
  ad <- ad[! ad$sending.cluster.name %in% "Plasma_Cells",]
  ad <- ad[! ad$target.cluster.name %in% "Plasma_Cells",]
  
  # change IDs
  levels(ad$sending.cluster.name)[levels(ad$sending.cluster.name)=="B_Cells"] <- "B-cells"
  levels(ad$target.cluster.name)[levels(ad$target.cluster.name)=="B_Cells"] <- "B-cells"
  levels(ad$sending.cluster.name)[levels(ad$sending.cluster.name)=="CAF1A"] <- temp_CAF1A_short
  levels(ad$target.cluster.name)[levels(ad$target.cluster.name)=="CAF1A"] <- temp_CAF1A_short
  levels(ad$sending.cluster.name)[levels(ad$sending.cluster.name)=="CAF1B"] <- temp_CAF1B_short
  levels(ad$target.cluster.name)[levels(ad$target.cluster.name)=="CAF1B"] <- temp_CAF1B_short
  levels(ad$sending.cluster.name)[levels(ad$sending.cluster.name)=="CAF2A"] <- temp_CAF2A_short
  levels(ad$target.cluster.name)[levels(ad$target.cluster.name)=="CAF2A"] <- temp_CAF2A_short
  levels(ad$sending.cluster.name)[levels(ad$sending.cluster.name)=="CAF2B"] <- temp_CAF2B_short
  levels(ad$target.cluster.name)[levels(ad$target.cluster.name)=="CAF2B"] <- temp_CAF2B_short
  
  levels(ad$sending.cluster.name)[levels(ad$sending.cluster.name)=="T_Regs"] <- "T-Regs"
  levels(ad$target.cluster.name)[levels(ad$target.cluster.name)=="T_Regs"] <- "T-Regs"
  
  levels(ad$sending.cluster.name)[levels(ad$sending.cluster.name)=="CD8_T_Cells"] <- "CD8+ T-cells"
  levels(ad$target.cluster.name)[levels(ad$target.cluster.name)=="CD8_T_Cells"] <- "CD8+ T-cells"
  
  levels(ad$sending.cluster.name)[levels(ad$sending.cluster.name)=="CD4_T_Cells"] <- "CD4+ T-cells"
  levels(ad$target.cluster.name)[levels(ad$target.cluster.name)=="CD4_T_Cells"] <- "CD4+ T-cells"
  
  levels(ad$sending.cluster.name)[levels(ad$sending.cluster.name)=="T_Cells_Cycling"] <- "Proliferating T-cells"
  levels(ad$target.cluster.name)[levels(ad$target.cluster.name)=="T_Cells_Cycling"] <- "Proliferating T-cells"
  
  levels(ad$sending.cluster.name)[levels(ad$sending.cluster.name)=="Epithelial_Basal_Cycling"] <- "Proliferating Cancer"
  levels(ad$target.cluster.name)[levels(ad$target.cluster.name)=="Epithelial_Basal_Cycling"] <- "Proliferating Cancer"
  
  levels(ad$sending.cluster.name)[levels(ad$sending.cluster.name)=="Epithelial_Basal"] <- "Cancer Epithelial"
  levels(ad$target.cluster.name)[levels(ad$target.cluster.name)=="Epithelial_Basal"] <- "Cancer Epithelial"
  
  #draw specificity png
  for(type in c("edge.count", "total.weight", "total.specificity", "total.combined")) {
    temp_png_function(paste0("Output/Figures/01_",type,"_signalling_circos_plot_specificity_cutoff_",temp_cutoff, ".png"))
    circos.par(start.degree = 140)
    chordDiagram(
      ad[,colnames(ad) %in% c("sending.cluster.name", "target.cluster.name", type)],
      order = order,
      directional = 1,
      direction.type = c("diffHeight"),
      grid.col = grid.col,
      annotationTrack = c("grid") # remove axis labels
    )
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, 
                  CELL_META$ylim[1], 
                  CELL_META$sector.index, 
                  facing = "clockwise", 
                  niceFacing = TRUE, 
                  adj = c(-0.5, 0.5))
    }, 
    bg.border = NA) # here set bg.border to NA is important
    dev.off()
    circos.clear()
    }


# 06 LIGANDS RECEPTORS  ---------------------------------------------------------------------

for(i in c(1,2)) {
  
  temp_df <- get(paste0("temp_df_",i))
  temp_df_combined <- NULL
  if(i == 1){
    type <- "ligand"
  }
  if(i == 2){
    type <- "receptor"
  }
  for(col in c(2:length(colnames(temp_df)))) {
    
    temp_df_combined <- rbind(temp_df_combined,
                              (data.frame(celltype = colnames(temp_df)[col],
                                          value = sum(temp_df[,col]),
                                          type = type)))
  }
  
  print(i)
  print(temp_df_combined)
  n <- paste0("temp_df_",i, "colsums")
  assign(n, temp_df_combined)
  
  
}

temp_df_LR_combined <- rbind(temp_df_1colsums,
                             temp_df_2colsums)

levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="Epithelial_Basal"] <- "Cancer Epithelial"
levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="Epithelial_Basal_Cycling"] <- "Cancer Proliferating"
levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="Epithelial_Luminal_Mature"] <- "Normal Luminal Epithelial"
levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="Myoepithelial"] <- "Normal Myoepithelial"

levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="T_Cells_Cycling"] <- "CD8+ T Cells Cycling"
levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="CD8_T_Cells"] <- "CD8+ T Cells"
levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="CD4_T_Cells"] <- "CD4+ T Cells"
levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="Plasma_Cells"] <- "Plasma Cells"
levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="T_Regs"] <- "T-Regs"
levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="B_Cells"] <- "B-Cells"
levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="CAF1A"] <- temp_CAF1A_short
levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="CAF1B"] <- temp_CAF1B_short
levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="CAF2A"] <- temp_CAF2A_short
levels(temp_df_LR_combined$celltype)[levels(temp_df_LR_combined$celltype)=="CAF2B"] <- temp_CAF2B_short

temp_df_LR_combined$celltype <- factor(temp_df_LR_combined$celltype,
                                       levels = rev(c(temp_CAF1A_short,temp_CAF1B_short,temp_CAF2A_short, temp_CAF2B_short,
                                                      "Endothelial",
                                                      "CD8+ T Cells Cycling",
                                                      "CD8+ T Cells",
                                                      "CD4+ T Cells",
                                                      "T-Regs",
                                                      "B-Cells",
                                                      "Plasma Cells",
                                                      "Myeloid",
                                                      "Cancer Proliferating",
                                                      "Cancer Epithelial",
                                                      "Normal Luminal Epithelial",
                                                      "Normal Myoepithelial")))

temp_df_LR_combined$type <- factor(temp_df_LR_combined$type,
                                   levels = rev(c("ligand","receptor")))

temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 8, 
      height = 4, 
      res = 300, 
      units = 'in'
    )
  }


temp_ggplot <- ggplot(temp_df_LR_combined,
                      aes(x=celltype,
                          y=value,
                          group=type,
                          fill=type)) + 
  geom_bar(stat="identity",
           position=position_dodge(), 
           colour="black") + 
  scale_y_continuous(breaks = c(0,
                                10,
                                20,
                                30,
                                40,
                                50)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  ylab(label = "Product of Combined Sum") +
  xlab(label = "Cell Type") +
  coord_flip()


temp_png_function("Output/Figures/02_sum_productcombined_ligand.png")
print(temp_ggplot)
dev.off()





# 07 CELL SIGNALLING HEATMAPs ------------------------


# LOAD DATA
{temp_df_ALL_raw <- read.csv(paste0("data/meanClusteredEM_paracrine_signaling_interactions.csv"))
  
  names(temp_df_ALL_raw) <- c("sending_cluster", "target_cluster", "ligand", "receptor", "original_ligand",
                              "specified_ligand", "combined_ligand", "original_receptor", "specified_receptor",
                              "combined_receptor", "product_of_original", "product_of_specified", "product_of_combined")
  
  temp_df_ALL_raw$cluster_index <- paste0(temp_df_ALL_raw$sending_cluster,
                                          "__",
                                          temp_df_ALL_raw$target_cluster)
  
  temp_df_ALL_raw$interaction_index <- paste0(temp_df_ALL_raw$ligand,
                                              "__",
                                              temp_df_ALL_raw$receptor)
  
  temp_df_ALL_raw$index <- paste0(temp_df_ALL_raw$sending_cluster,
                                  "__",
                                  temp_df_ALL_raw$target_cluster,
                                  "__",
                                  temp_df_ALL_raw$ligand,
                                  "__",
                                  temp_df_ALL_raw$receptor)
  
  temp_df_ALL_raw_CAFs <- 
    temp_df_ALL_raw %>% 
    separate(index, c("sending_cluster", "target_cluster","ligand", "receptor"), 
             "__")
  
  temp_df_ALL_raw_CAFs$cluster_index <- paste0(temp_df_ALL_raw_CAFs$sending_cluster,
                                               "_",
                                               temp_df_ALL_raw_CAFs$target_cluster)
  
  temp_df_ALL_raw_CAFs$interaction_index <- paste0(temp_df_ALL_raw_CAFs$ligand,
                                                   "_",
                                                   temp_df_ALL_raw_CAFs$receptor)
  
  temp_df_ALL_raw_CAFs$index <- paste0(temp_df_ALL_raw_CAFs$sending_cluster,
                                       "_",
                                       temp_df_ALL_raw_CAFs$target_cluster,
                                       "_",
                                       temp_df_ALL_raw_CAFs$ligand,
                                       "_",
                                       temp_df_ALL_raw_CAFs$receptor)

  # rescale unfiltered dataset
  temp_df_ALL_raw_CAFs_rescaled <-
    transform(temp_df_ALL_raw_CAFs,
              rescale = ave(product_of_original,
                            interaction_index,
                            FUN = scale))
  
  # assigining NA solo interaction values a rescale figure
  temp_NA_df <- 
    temp_df_ALL_raw_CAFs_rescaled[temp_df_ALL_raw_CAFs_rescaled$rescale %in% NaN,]
  
  temp_NA_df_mirror <- NULL
  for(num in c(1:nrow(temp_NA_df))){
    temp_row <- 
      temp_NA_df[num,]
    temp_row[,1:9]<- 0 
    temp_row$cluster_index <- paste0("MIRROR")
    
    temp_NA_df_mirror <- rbind(temp_NA_df_mirror,temp_row)
  }
  temp_df_ALL_raw_CAFs_rescaled <- 
    rbind(temp_df_ALL_raw_CAFs_rescaled,temp_NA_df_mirror)
  
  # rerun rescale
  # rescale unfiltered dataset
  temp_df_ALL_raw_CAFs_rescaled <-
    transform(temp_df_ALL_raw_CAFs_rescaled,
              rescale = ave(product_of_original,
                            interaction_index,
                            FUN = scale))
  
  temp_df_ALL_raw_CAFs_rescaled <-
    temp_df_ALL_raw_CAFs_rescaled[!temp_df_ALL_raw_CAFs_rescaled$cluster_index %in% "MIRROR",]
  
  # export
  temp_df_ALL_raw_CAFs_rescaled_export <- 
    temp_df_ALL_raw_CAFs_rescaled[,colnames(temp_df_ALL_raw_CAFs_rescaled) 
                                  %in% c("interaction_index", 
                                         "rescale",
                                         "sending_cluster",
                                         "target_cluster",
                                         "product_of_original")]
  
}

# Filter CAFs and dataset
{temp_df_ALL_raw_CAFs_rescaled_CAFs <- NULL
  for(i in c("CAF1A", "CAF1B", "CAF2A", "CAF2B")) {
    
    temp_df <- temp_df_ALL_raw_CAFs_rescaled %>% filter(sending_cluster == i)
    
    temp_df_ALL_raw_CAFs_rescaled_CAFs <- rbind(temp_df_ALL_raw_CAFs_rescaled_CAFs,
                                                temp_df)
    rm(temp_df)
  }
  # # filter spreadsheet
  temp_df_ALL_raw_CAFs_rescaled_justCAFs_export_FILTERED <-
    temp_df_ALL_raw_CAFs_rescaled_CAFs[temp_df_ALL_raw_CAFs_rescaled_CAFs$product_of_original > 0,]
  
  # # filter spreadsheet
  temp_df_ALL_raw_CAFs_rescaled_justCAFs_export_FILTERED <-
    temp_df_ALL_raw_CAFs_rescaled_justCAFs_export_FILTERED[temp_df_ALL_raw_CAFs_rescaled_justCAFs_export_FILTERED$rescale > 0.1,]
  
  # remove all collagens
  temp_df_ALL_raw_CAFs_rescaled_justCAFs_export_FILTERED <-
    temp_df_ALL_raw_CAFs_rescaled_justCAFs_export_FILTERED[!grepl('^area$|COL', temp_df_ALL_raw_CAFs_rescaled_justCAFs_export_FILTERED$interaction_index),]
  
  # remove HLA's
  temp_df_ALL_raw_CAFs_rescaled_justCAFs_FILTERED <-
    temp_df_ALL_raw_CAFs_rescaled_justCAFs_export_FILTERED[!grepl('^area$|HLA', temp_df_ALL_raw_CAFs_rescaled_justCAFs_export_FILTERED$interaction_index),]
  
  
  # cell IDs for plotting
  temp_target_ids <- 
    unique(temp_df_ALL_raw_CAFs_rescaled_justCAFs_FILTERED$target_cluster)
  temp_target_ids <- 
    temp_target_ids[! temp_target_ids %in% c("CAF1A", "CAF1B", "CAF2A", "CAF2B")]
  
  # export
  temp_df_ALL_raw_CAFs_rescaled_justCAFs_FILTERED_export <- 
    temp_df_ALL_raw_CAFs_rescaled_justCAFs_FILTERED[,colnames(temp_df_ALL_raw_CAFs_rescaled_justCAFs_FILTERED) 
                                                    %in% c("interaction_index", 
                                                           "rescale",
                                                           "rescalejustCAFs",
                                                           "sending_cluster",
                                                           "target_cluster",
                                                           "product_of_original")]
  
  # for plot
  temp_interactions <- 
    temp_df_ALL_raw_CAFs_rescaled_justCAFs_FILTERED
  
}

# PLOT ACROSS ALL CELLTYPES
for(celltype in c("T_cells", "Cancer", "Endothelial", "Myeloid")){
  print(celltype)
  
  if(celltype == "T_cells"){
    temp_png_function <-
      function(x) {
        png(
          file = (x), 
          width = 12, 
          height = 18, 
          res = 600, 
          units = 'in'
        )
      }
    
    temp_pdf_function <-
      function(x) {
        pdf(
          file = (x),
          width = 12,
          height =18,
          useDingbats=F
        )
      }
    
  }
  if(celltype == "Cancer"){
    temp_png_function <-
      function(x) {
        png(
          file = (x), 
          width = 8, 
          height = 18, 
          res = 600, 
          units = 'in'
        )
      }
    
    temp_pdf_function <-
      function(x) {
        pdf(
          file = (x),
          width = 8,
          height =18,
          useDingbats=F
        )
      }
  }
  if(celltype == "Endothelial"){
    temp_png_function <-
      function(x) {
        png(
          file = (x), 
          width = 4, 
          height = 18, 
          res = 600, 
          units = 'in'
        )
      }
    
    temp_pdf_function <-
      function(x) {
        pdf(
          file = (x),
          width = 4,
          height =18,
          useDingbats=F
        )
      }
  }
  if(celltype == "Myeloid"){
    temp_png_function <-
      function(x) {
        png(
          file = (x), 
          width = 4, 
          height = 18, 
          res = 600, 
          units = 'in'
        )
      }
    
    temp_pdf_function <-
      function(x) {
        pdf(
          file = (x),
          width = 4,
          height =18,
          useDingbats=F
        )
      }
  }
  
  temp_interactions_rescale_subset <- NULL
  if(celltype == "T_cells"){
    temp_interactions_rescale_subset <- 
      rbind(temp_interactions_rescale_subset,
            temp_interactions[temp_interactions$target_cluster %in% paste0("CD4_T_Cells") ,],
            temp_interactions[temp_interactions$target_cluster %in% paste0("CD8_T_Cells") ,],
            temp_interactions[temp_interactions$target_cluster %in% paste0("T_Regs") ,],
            temp_interactions[temp_interactions$target_cluster %in% paste0("T_Cells_Cycling") ,])
  }
  if(celltype == "Cancer"){
    temp_interactions_rescale_subset <- 
      rbind(temp_interactions_rescale_subset,
            temp_interactions[temp_interactions$target_cluster %in% paste0("Epithelial_Basal") ,],
            temp_interactions[temp_interactions$target_cluster %in% paste0("Epithelial_Basal_Cycling") ,])
  }
  if(celltype == "Endothelial"){
    temp_interactions_rescale_subset <- 
      rbind(temp_interactions_rescale_subset,
            temp_interactions[temp_interactions$target_cluster %in% paste0("Endothelial") ,])
    
  }
  if(celltype == "Myeloid"){
    temp_interactions_rescale_subset <- 
      rbind(temp_interactions_rescale_subset,
            temp_interactions[temp_interactions$target_cluster %in% paste0("Myeloid") ,])
  }
  
  
  
  # keep only top 2 ligand - receptor pairs
  # temp_interactions_to_keep <-
  #   temp_interactions_rescale_subset %>% group_by(ligand) %>%
  #   top_n(2,
  #         product_of_original)
  temp_interactions_to_keep <- temp_interactions_rescale_subset
  
  temp_interactions_rescale_to_keep <-
    as.data.frame(temp_interactions_to_keep)
  
  temp_interactions_rescale_to_keep <-
    unique(temp_interactions_rescale_to_keep$interaction_index)
  
  temp_interactions_rescale_subset <-
    temp_interactions_rescale_subset[temp_interactions_rescale_subset$interaction_index %in% temp_interactions_rescale_to_keep,]
  
  # remove irrelevant columns
  temp_interactions_rescale_subset <- 
    temp_interactions_rescale_subset[,colnames(temp_interactions_rescale_subset) 
                                     %in% c("interaction_index",
                                            "cluster_index",
                                            "rescale",
                                            "rescalejustCAFs",
                                            "sending_cluster",
                                            "target_cluster",
                                            "product_of_original")]
  
  temp_interactions_rescale_subset <- arrange(temp_interactions_rescale_subset,
                                              # (interaction_index),
                                              desc(product_of_original))

  temp_interactions_rescale_subset <- 
    temp_interactions_rescale_subset[temp_interactions_rescale_subset$product_of_original > 0.01,]
  
  # export filtered
  temp_interactions_rescale_subset$newCAF <- NULL
  temp_interactions_rescale_subset$newCAF[ temp_interactions_rescale_subset$sending_cluster == "CAF1A" ] <- temp_CAF1A_short
  temp_interactions_rescale_subset$newCAF[ temp_interactions_rescale_subset$sending_cluster == "CAF1B" ] <- temp_CAF1B_short
  temp_interactions_rescale_subset$newCAF[ temp_interactions_rescale_subset$sending_cluster  == "CAF2A" ] <- temp_CAF2A_short
  temp_interactions_rescale_subset$newCAF[ temp_interactions_rescale_subset$sending_cluster  == "CAF2B" ] <- temp_CAF2B_short
  
  print("  # interactions post filtering")
  print(length(unique(temp_interactions_rescale_subset$interaction_index)))
  
  # filter for interactions for clustering
  temp_top <- 
    unique(temp_interactions_rescale_subset$interaction_index)[1:100]
  temp_interactions_top <- 
    temp_interactions_rescale_subset[temp_interactions_rescale_subset$interaction_index %in% temp_top,]
  
  
  # rescale by interaction
  temp_interactions_top <-
    transform(temp_interactions_top,
              rescale = ave(product_of_original,
                            interaction_index,
                            FUN = scale))
  
  temp_interactions.r <- temp_interactions_top[,colnames(temp_interactions_top) %in% c("interaction_index",
                                                                                       "rescale",
                                                                                       "cluster_index")]
  
  temp_interactions.r <- reshape(temp_interactions.r, 
                                 timevar="cluster_index", 
                                 idvar=c("interaction_index"), 
                                 direction="wide")
  
  temp_interactions.r[is.na(temp_interactions.r)] <- 0
  rownames(temp_interactions.r) <- temp_interactions.r$interaction_index
  temp_interactions.r <- temp_interactions.r[, ! colnames(temp_interactions.r) %in% c("interaction_index")]
  
  # cluster interactions
  temp_d <- dist(temp_interactions.r, 
                 method = "euclidean")
  
  temp_hc1 <- hclust(temp_d, 
                     method = "complete")
  
  temp_order_hclust <- temp_hc1$labels[temp_hc1$order]
  
  # order by hclust results
  temp_interactions_subset <- temp_interactions_top
  temp_interactions_subset$interaction_index <- factor(temp_interactions_subset$interaction_index,
                                                       levels = rev(temp_order_hclust))
  
  # new labels
  temp_CAF1A_short <- "myCAFs"
  temp_CAF1B_short <- "iCAFs"
  temp_CAF2A_short <- "dVDSCs"
  temp_CAF2B_short <- "imVDSCs"
  temp_combined_short <- c(temp_CAF1A_short,
                           temp_CAF1B_short,
                           temp_CAF2A_short,
                           temp_CAF2B_short)
  
  temp_interactions_subset$newCAF <- NULL
  temp_interactions_subset$newCAF[ temp_interactions_subset$sending_cluster == "CAF1A" ] <- temp_CAF1A_short
  temp_interactions_subset$newCAF[ temp_interactions_subset$sending_cluster == "CAF1B" ] <- temp_CAF1B_short
  temp_interactions_subset$newCAF[ temp_interactions_subset$sending_cluster  == "CAF2A" ] <- temp_CAF2A_short
  temp_interactions_subset$newCAF[ temp_interactions_subset$sending_cluster  == "CAF2B" ] <- temp_CAF2B_short
  
  temp_interactions_subset$newCAF <- factor(temp_interactions_subset$newCAF,
                                            levels=temp_combined_short)
  
  temp_interactions_subset$newtarget <- NULL
  
  # factorise
  if(celltype == "T_cells"){
    temp_interactions_subset$newtarget[ temp_interactions_subset$target_cluster == "CD4_T_Cells" ] <- "CD4+ T-Cells"
    temp_interactions_subset$newtarget[ temp_interactions_subset$target_cluster == "CD8_T_Cells" ] <- "CD8+ T-Cells"
    temp_interactions_subset$newtarget[ temp_interactions_subset$target_cluster  == "T_Regs" ] <- "T-Regs"
    temp_interactions_subset$newtarget[ temp_interactions_subset$target_cluster  == "T_Cells_Cycling" ] <- "CD8+ T-Cells Cycling"
    
    temp_interactions_subset$newtarget <- factor(temp_interactions_subset$newtarget,
                                                 levels=c("CD8+ T-Cells", "CD8+ T-Cells Cycling", "CD4+ T-Cells", "T-Regs"))
  }
  if(celltype == "Cancer"){
    temp_interactions_subset$newtarget[ temp_interactions_subset$target_cluster == "Epithelial_Basal" ] <- "Cancer"
    temp_interactions_subset$newtarget[ temp_interactions_subset$target_cluster == "Epithelial_Basal_Cycling" ] <- "Cancer Cycling"
    
    temp_interactions_subset$newtarget <- factor(temp_interactions_subset$newtarget,
                                                 levels=c("Cancer", "Cancer Cycling"))
  }
  if(celltype == "Endothelial"){
    temp_interactions_subset$newtarget[ temp_interactions_subset$target_cluster == "Endothelial" ] <- "Endothelial"
  }
  if(celltype == "Myeloid"){
    temp_interactions_subset$newtarget[ temp_interactions_subset$target_cluster == "Myeloid" ] <- "Myeloid"
  }
  
  # plot
  p1 <- ggplot(temp_interactions_subset) +
    geom_point(aes(x=newCAF,
                   y=interaction_index, 
                   size=rescale, 
                   colour=rescale)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text.y = element_text(size=7,
                                     angle = 0),
          strip.text.x = element_text(size = 16)
    ) +
    scale_size(range = c(1,7),
               guide = 'none') +
    scale_colour_gradientn(colours = c("#4575b4", "#fdae61", "#d73027"),
                           limits=c(min(temp_interactions_subset$rescale),
                                    max(temp_interactions_subset$rescale)), 
                           labels=scales::number_format(accuracy = 1,
                                                        decimal.mark = ',')) +
    facet_grid(~ newtarget) +
    labs(size=NULL, colour="Interaction\nStrength") + # /n new line
    xlab(label = " ") +
    ylab(label = " ")
  
  # write.csv(temp_interactions_subset, "temp.csv")
  
  temp_png_function(paste0("Output/Figures/05_",celltype,"_ALL_INTERACTIONS_COMBINED.png"))
  print(p1)
  dev.off() 

}







