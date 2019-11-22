# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 05D GENE ONTOLOGY SCORING SCRIPT CLUSTER PROFILER HPC
# R 3.4.1
#
### inputs (in order); 
#     - FindAllMarkers Output (Differentially expressed genes)
# 
### REQUIRES R v.3.4.1
### QSUB ARGUMENTS
# qsub \
# -cwd \
# -pe smp 8 \
# -l mem_requested=10G \
# -P TumourProgression \
# -b y \
# -j y \
# -V \
# -N ${JOBNAME}\
# "${R} CMD BATCH \
# --no-save '--args \
# ${SEURATOBJECTPATH}' \
# ${SCRIPT}"
#
#
# 01 LOAD PACKAGES -----------------------------------------------------------

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# 02 COMMAND LINE ARGUMENTS  ------------------------------------------------------------

temp_start_time <- date()
temp_args <-
  commandArgs(trailingOnly = T)

# input path to seurat object
temp_DEG_path <- 
  temp_args[1]

# 03 SET UP AND FUNCTIONS ---------------------------------------------------------------

dir.create("Output")
dir.create("Output/Figures")
dir.create("Output/GO_enrichment")

# 04 LOAD OP 250 GENES FROM EACH CLUSTER -------

temp_findallmarkers <- read.csv(temp_DEG_path)

# filter top 250 for each cluster
temp_top250 <- NULL
for(i in unique(temp_findallmarkers$cluster)) {
  
  temp_df <- temp_findallmarkers[temp_findallmarkers$cluster == i,]
  temp_df <- temp_df[1:250,]
  temp_top250 <- rbind(temp_top250,temp_df)
  
}
table(temp_top250$cluster)


# 05 RUN GO ------------------------------------------------------------------

temp_clusterMarkers <- NULL
for (cluster in as.character(unique(temp_findallmarkers$cluster))){
  
  temp_signatureEntrezIDs <- bitr(temp_findallmarkers[temp_findallmarkers$cluster == cluster,]$gene,
                                  fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = "org.Hs.eg.db"
  )
  
  temp_cm <-
    data.frame(EntrezID = temp_signatureEntrezIDs$ENTREZID,
               clusterID = cluster)
  if (!is.null(temp_clusterMarkers))
  {
    temp_clusterMarkers <- rbind(temp_clusterMarkers, 
                                 temp_cm)
  } else {
    temp_clusterMarkers <- temp_cm
  }}

temp_go <- compareCluster(
  EntrezID ~ clusterID,
  data = temp_clusterMarkers,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "ALL",
  readable = T
)

temp_go_df <- as.data.frame(temp_go)

# df for plotting
temp_go_df$logp_val <- (-(log10(temp_go_df$p.adjust)))
write.csv(temp_go_df,
          paste0("Output/GO_enrichment/TOP250_ALLGENES_compareCluster_stromal.csv"), 
          row.names = F)

# 06 SIMPLIFIED DF --------------------------------------------------------------------

temp_ids <- unique(temp_go_df$clusterID)
for(i in temp_ids) {
  
  temp_df <- temp_go_df[temp_go_df$Cluster == i,]
  temp_df$newcol <- paste0(temp_df$Description, "_", temp_df$ID)
  
  temp_df <- temp_df[,c("Cluster","newcol","logp_val")]
  
  n <- paste0("temp_enrichGO_df_",
              i)
  assign(n,
         temp_df)
  
}

temp_merged_1 <- merge(get(paste0("temp_enrichGO_df_",temp_ids[1])), 
                       get(paste0("temp_enrichGO_df_",temp_ids[2])), 
                       by = "newcol", 
                       all = T)
temp_merged_2 <- merge(get(paste0("temp_enrichGO_df_",temp_ids[3])), 
                       get(paste0("temp_enrichGO_df_",temp_ids[4])), 
                       by = "newcol", 
                       all = T)
temp_merged_all <- merge(temp_merged_1, 
                         temp_merged_2, 
                         by = "newcol", 
                         all = T)

temp_merged_all <- temp_merged_all[, -grep("Cluster.", colnames(temp_merged_all))]
temp_merged_all[is.na(temp_merged_all)] <- 0
names(temp_merged_all) <- c("Description", temp_ids)

temp_merged_all$CAFDIFF <- (temp_merged_all[,temp_ids[1]] - temp_merged_all[,temp_ids[2]])
temp_merged_all$VDSCDIFF <- (temp_merged_all[,temp_ids[3]] - temp_merged_all[,temp_ids[4]])

write.csv(temp_merged_all,
          "Output/GO_enrichment/TOP250_ALLGENES_compareCluster_stromal_simplified.csv",
          row.names = F)












