# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 05C GO JOB SUBMISSION TO HPC SCRIPT
#
# GLOBS
source activate Renv
  ## R path (NOTE USE R V.3.4.1)
  R="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.4.1/bin/R"
  ## script path
  TEMPPWD=$(pwd)
  SCRIPT="${TEMPPWD}/05D_GO_HPC_R341.R"
  ## DEG list
  DEGPATH="/share/ScratchGeneral/sunwu/projects/Stromal_heterogeneity_paper/REPROCESSING_FROM_RAW/04_seurat_annotation/output/Output/Stromal_gene_signatures/01_FindAllMarkers_stromal_celltype.csv"
  # job name
  JOBNAME="GO_enrichment"

# DIRECTORIES
mkdir "output"
#
# SUBMIT JOB TO CLUSTER
    cd "output"
    qsub \
    -cwd \
    -pe smp 8 \
    -l mem_requested=10G \
    -P TumourProgression \
    -b y \
    -j y \
    -V \
    -N ${JOBNAME}\
    "${R} CMD BATCH \
    --no-save '--args \
    ${DEGPATH}' \
    ${SCRIPT}"

cd ${TEMPPWD}
