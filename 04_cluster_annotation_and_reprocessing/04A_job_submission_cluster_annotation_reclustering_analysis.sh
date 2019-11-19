# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 03A SEURAT V2 INTEGRATION JOB SUBMISSION TO HPC SCRIPT
#
# GLOBS
source activate Renv
  ## R path
  R="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.5.0/bin/R"
  ## species
  SPECIES='human'
  ## script path
  TEMPPWD=$(pwd)
  SEURATSCRIPT="${TEMPPWD}/03B_R_seurat_integration.R"
  ## seurat individual objects path
  SEURATOUTPUTPATH="/share/ScratchGeneral/sunwu/projects/Stromal_heterogeneity_paper/REPROCESSING_FROM_RAW/02_seurat_individual_processing/output/"
  # job name
  SEURATJOBNAME="seurat_integration"
  # project name
  PROJECTNAME="intregration"
#
# DIRECTORIES
mkdir "output"
#
# SUBMIT JOB TO CLUSTER
    cd "output"
    qsub \
    -cwd \
    -pe smp 16 \
    -l mem_requested=20G \
    -P TumourProgression \
    -b y \
    -j y \
    -V \
    -N ${SEURATJOBNAME}\
    "${R} CMD BATCH \
    --no-save '--args \
    ${PROJECTNAME} \
    ${SPECIES} \
    ${SEURATOUTPUTPATH}' \
    ${SEURATSCRIPT}"

cd ${TEMPPWD}
