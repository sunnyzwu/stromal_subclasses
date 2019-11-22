# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 05A AUCELL JOB SUBMISSION TO HPC SCRIPT
#
# GLOBS
source activate Renv
  ## R path
  R="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.5.0/bin/R"
  ## script path
  TEMPPWD=$(pwd)
  SCRIPT="${TEMPPWD}/05B_AUCell_HPC.R"
  ## seurat objects path
  SEURATOBJECTPATH="/share/ScratchGeneral/sunwu/projects/Stromal_heterogeneity_paper/REPROCESSING_FROM_RAW/04_seurat_annotation/output/Output/Rdata/Rdata_Stromal_Immune.Rdata"
  ## genesets path
  GENESETS="${TEMPPWD}/gene_signatures.csv"
  # job name
  JOBNAME="AUCell_annotation"

# DIRECTORIES
mkdir "output"
#
# SUBMIT JOB TO CLUSTER
    cd "output"
    qsub \
    -cwd \
    -pe smp 16 \
    -l mem_requested=10G \
    -P TumourProgression \
    -b y \
    -j y \
    -V \
    -N ${JOBNAME}\
    "${R} CMD BATCH \
    --no-save '--args \
    ${SEURATOBJECTPATH} \
    ${GENESETS}' \
    ${SCRIPT}"

cd ${TEMPPWD}
