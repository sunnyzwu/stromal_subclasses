# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 06D PYSCENIC VISUALISATION
# R 3.5.0
#
#
# 01 LOAD PACKAGES -----------------------------------------------------------

library(Seurat)
library(reshape2)
library(dplyr)
library(pheatmap)

# 02 COMMAND LINE ARGUMENTS  ------------------------------------------------------------

temp_start_time <- date()
temp_args <-
  commandArgs(trailingOnly = T)

# input path to seurat object
temp_seurat_object_path <- 
  temp_args[1]

# input path to seurat object
temp_aucell_path <- 
  temp_args[2]

# 03 SET UP AND FUNCTIONS ---------------------------------------------------------------

dir.create("output_visualisation")

# 04 LOAD DATA ---------------------------------------------------------------

seurat_10X <- 
  readRDS(temp_seurat_object_path)

temp_auc_df <- as.data.frame(read.csv(temp_aucell_path))

# 05 APPEND PYSCENIC DATA ----------------------------------------------------

temp_cells_AUC_matrix <- data.frame(row.names = names(temp_auc_df)[2:length(names(temp_auc_df))])

temp_auc_df_trimmed <- t(temp_auc_df[2:nrow(temp_auc_df),2:length(temp_auc_df)])
colnames(temp_auc_df_trimmed) <- temp_auc_df$Regulon[2:nrow(temp_auc_df)]

temp_cells_AUC_matrix <- cbind(temp_cells_AUC_matrix,
                               temp_auc_df_trimmed)

temp_cells_AUC_matrix_sorted <- 
  temp_cells_AUC_matrix[rownames(seurat_10X@meta.data),,drop=FALSE]

temp_cells_AUC_matrix_sorted <- 
  as.data.frame.matrix(temp_cells_AUC_matrix_sorted)

temp_first_geneset <- 
  (ncol(seurat_10X@meta.data))+1

seurat_10X <- AddMetaData(seurat_10X, 
                                   metadata = temp_cells_AUC_matrix_sorted)

temp_last_geneset <-
  (ncol(seurat_10X@meta.data))

temp_gene_set_names <-
  colnames(seurat_10X@meta.data[temp_first_geneset:temp_last_geneset])

temp_AUCell_sig_IDs <- temp_gene_set_names


# 06 pySCENIC HEATMAP AVERAGED --------------------------------

temp_pos <- c()
for(i in temp_AUCell_sig_IDs){
  
  if(stringr::str_detect(i,"[+]")) {
    temp_pos <- c(temp_pos,i)
  }
  
}

temp_data_frame_GRNs <- 
  as.data.frame(seurat_10X@meta.data[,temp_pos])

temp_data_frame_2 <- data.frame(barcode = row.names(temp_data_frame_GRNs),
                                cluster = seurat_10X@meta.data$stromal_celltype)

temp_data_frame_combined <- cbind(temp_data_frame_2,
                                  temp_data_frame_GRNs)

temp_data_frame_combined <- arrange(temp_data_frame_combined,
                                    cluster)

temp_data_frame_combined$barcode <- factor(temp_data_frame_combined$barcode,
                                           levels = (temp_data_frame_combined$barcode)[order(temp_data_frame_combined$cluster)])

temp_data_frame_combined <- arrange(temp_data_frame_combined,
                                    barcode)

temp_data_frame_combined.m <- melt(temp_data_frame_combined)

# statistical test
temp_data_frame_t.m.df <- 
  temp_data_frame_combined.m[, names(temp_data_frame_combined.m) %in% c("cluster","variable", "value")]

temp_GRN_sig <- NULL
for(GRN in unique(temp_data_frame_t.m.df$variable)) {
  #http://www.sthda.com/english/wiki/one-way-anova-test-in-r
  
  # print(GRN)
  temp_subset <- temp_data_frame_t.m.df[temp_data_frame_t.m.df$variable == GRN,]
  temp_subset <- temp_subset[, names(temp_subset) %in% c("cluster", "value")]
  
  temp.aov <- aov(value ~ cluster, 
                  data = temp_subset)
  
  temp.aov <- data.frame(unclass(summary(temp.aov)), check.names = FALSE, stringsAsFactors = FALSE)
  
  temp_pval <- temp.aov[1,5]
  # print(temp_pval)
  
  if(! temp_pval == "NaN") {
    if(temp_pval < 0.001) {
      temp_row <- c(GRN, temp_pval, T)
      temp_GRN_sig <- rbind(temp_GRN_sig,temp_row)
    }
    if(temp_pval > 0.001) {
      temp_row <- c(GRN, temp_pval, F)
      temp_GRN_sig <- rbind(temp_GRN_sig,temp_row)
    }}
  if( temp_pval == "NaN") {
    temp_row <- c(GRN, temp_pval, T)
    temp_GRN_sig <- rbind(temp_GRN_sig,temp_row)
  }
  
  
}

temp_GRN_sig <- 
  as.data.frame(temp_GRN_sig)
colnames(temp_GRN_sig) <- 
  c("GRN", "p_val", "sig_0.05")
temp_GRN_sig$p_val <- 
  as.numeric(as.character(temp_GRN_sig$p_val))

temp_GRN_sig <- arrange(temp_GRN_sig,
                        p_val)

temp_GRN_sig <- 
  temp_GRN_sig[temp_GRN_sig$sig_0.05 == "TRUE",]

temp_data_frame_combined.m <- 
  temp_data_frame_combined.m[temp_data_frame_combined.m$variable %in% temp_GRN_sig$GRN,]

# top x highest AUC values by average
temp_highest_expressed_cluster_average <-
  aggregate(.~variable,
            temp_data_frame_combined.m,
            mean)
temp_highest_expressed_cluster_average <- arrange(temp_highest_expressed_cluster_average,
                                                  -value)
clustercol <- T
topn <- 50
    
temp_highest_expressed_cluster_average_topn <- 
      (temp_highest_expressed_cluster_average %>% top_n(topn, value))

temp_top_50_TFs <- temp_highest_expressed_cluster_average_topn$variable
    
    temp_data_frame_combined.m_topn <- 
      temp_data_frame_combined.m[temp_data_frame_combined.m$variable %in% temp_highest_expressed_cluster_average_topn$variable,]
    
    # cluster average
    temp_data_frame_combined.m.agg <-
      aggregate(.~cluster+variable,
                temp_data_frame_combined.m_topn,
                mean)
    
    temp_data_frame_combined.m.agg.dcast <-
      dcast(data = temp_data_frame_combined.m.agg,
            formula = variable~cluster,
            fun.aggregate = sum,
            value.var = "value")
    
    
    #row names to variables
    rownames(temp_data_frame_combined.m.agg.dcast) <- 
      temp_data_frame_combined.m.agg.dcast$variable
    
    temp_data_frame_combined.m.agg.dcast <- 
      temp_data_frame_combined.m.agg.dcast[, ! colnames(temp_data_frame_combined.m.agg.dcast) %in% "variable"]
    
    temp_data_frame_hclust <- as.matrix(temp_data_frame_combined.m.agg.dcast)
    # temp_annotation <- data.frame(row.names = temp_data_frame_combined$barcode,
    #                               cluster = temp_data_frame_combined$cluster,
    #                               # sampleID = temp_data_frame_combined$sampleID,
    #                               subtype = temp_data_frame_combined$clinical_subtype)
    
    # plot
    temp_png_function <-
      function(x) {
        png(
          file = (x), 
          width = 10, 
          height = 15, 
          res = 600, 
          units = 'in'
        )
      }
    
    temp_pdf_function <-
      function(x) {
        pdf(
          file = (x), 
          width = 10, 
          height = 15, 
          useDingbats = F
        )
      }
    
    temp_png_function(paste0("output_visualisation/01_pheatmap_pySCENIC_clusteravg_",clustercol,"_top",topn,".png"))
    pheatmap(temp_data_frame_hclust, 
             color = c("#fff7ec",
                       "#fee8c8",
                       "#fdd49e",
                       "#fdbb84",
                       "#fc8d59",
                       "#ef6548",
                       "#d7301f",
                       "#990000"),
             cluster_cols = clustercol, 
             cluster_rows = T, 
             scale = "row", 
             fontsize_row = 15,
             show_colnames = T, 
             fontsize_col = 20,
             annotation_legend = T,
             cellheight=15, 
             cellwidth =30,
             gaps_col = NULL, 
             annotation_names_col = T, 
             angle_col = 20,
             treeheight_row = 50, 
             legend = T,
             border_color=FALSE)
    dev.off()

# 07 CORRELATION PLOTS -------------------------------------------------------

  temp_GRNS <- temp_top_50_TFs
  temp_combined_Rsquared <- NULL
  for (i in temp_GRNS) {
    print(i)
    
    temp_AUCEll_name <- i
    temp_gene_name <- gsub("([()])", "", temp_AUCEll_name)
    temp_gene_name <- gsub("([(-)])", "", temp_gene_name)
    temp_gene_name <- gsub("([-])", "", temp_gene_name)
    temp_gene_name <- gsub("([+])", "", temp_gene_name)
    
    # temp_Names <- c(paste0(temp_gene_name, "_GRN"),
    #                 paste0(temp_gene_name, "_Expression"))
    
    temp_AUCEll_data <-
      data.frame(
        AUC_value = seurat_10X@meta.data[, temp_AUCEll_name],
        Expression = seurat_10X@data[temp_gene_name, ]
      )
    
    temp_linearMod <- lm(AUC_value ~ Expression,
                         data = temp_AUCEll_data)
    temp <- summary(temp_linearMod)
    
    
    temp_label <- paste0(temp_gene_name)
    
    temp_Rsquared <- paste0("R = ", round(temp$adj.r.squared, digits=3))
    
    temp_ggplot <- 
      ggplot(temp_AUCEll_data, 
             aes(x=AUC_value, y=Expression)) +
      geom_point(size=2, shape=23) +
      geom_smooth(method="lm") + 
      labs(title = temp_label) +
      annotate("text", 
               label = temp_Rsquared,  x = -Inf, y = Inf, 
               hjust = 0, vjust = 1,parse = F,
               colour = "red",
               size=10) +
      theme(
        plot.title = element_text(color = '#666666',face = 'bold',size = 25
        ))
    
    print(temp_ggplot)
    dev.off()
    
    temp_df <- data.frame(GRN = temp_gene_name,
                          R_squared = temp$adj.r.squared)
    
    temp_combined_Rsquared <- rbind(temp_combined_Rsquared, temp_df)
    n <- paste0("temp_ggplot_",temp_gene_name)
    assign(n,temp_ggplot)
  }

# COMBINED PLOT
temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 24, 
      height = 30, 
      res = 300, 
      units = 'in'
    )
  }

temp_grid <- 
  plot_grid(plotlist=mget(paste0("temp_ggplot_",temp_combined_Rsquared$GRN)),
            labels = NULL,
            nrow = 10,
            ncol = 5)

temp_png_function("output_visualisation/02_COMBINED.png")
print(temp_grid)
dev.off()

