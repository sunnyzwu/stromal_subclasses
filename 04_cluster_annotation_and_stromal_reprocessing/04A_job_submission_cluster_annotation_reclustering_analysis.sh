# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 04A SEURAT V2 ANNOTATION AND REPROCESSING JOB SUBMISSION TO HPC SCRIPT
#
# GLOBS
source activate Renv
  ## R path
  R="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.5.0/bin/R"
  ## species
  SPECIES='human'
  ## script path
  TEMPPWD=$(pwd)
  SEURATSCRIPT="${TEMPPWD}/04B_R_seurat_integration.R"
  ## seurat individual objects path
  SEURATOUTPUTPATH="/share/ScratchGeneral/sunwu/projects/Stromal_heterogeneity_paper/REPROCESSING_FROM_RAW/03_seurat_integration/output/"
  # job name
  SEURATJOBNAME="seurat_annotation"
  # project name
  PROJECTNAME="annotation"
#
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
    -N ${SEURATJOBNAME}\
    "${R} CMD BATCH \
    --no-save '--args \
    ${PROJECTNAME} \
    ${SPECIES} \
    ${SEURATOUTPUTPATH}' \
    ${SEURATSCRIPT}"

cd ${TEMPPWD}
