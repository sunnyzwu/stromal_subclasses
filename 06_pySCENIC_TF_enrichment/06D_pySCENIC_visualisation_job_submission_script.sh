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
  SCRIPT="${TEMPPWD}/06E_pySCENIC_visualisation.R"
  ## seurat objects path
  SEURATOBJECTPATH="/share/ScratchGeneral/sunwu/projects/Stromal_heterogeneity_paper/REPROCESSING_FROM_RAW/04_seurat_annotation/output/Output/Rdata/Rdata_all_stromal_dataset.Rdata"
  ## aucell path
  AUCELLPATH="/share/ScratchGeneral/sunwu/projects/Stromal_heterogeneity_paper/REPROCESSING_FROM_RAW/06_pySCENIC/output/aucell_matrix.csv"
  # job name
  JOBNAME="pySCENIC_input"

# DIRECTORIES
mkdir "output"
#
# SUBMIT JOB TO CLUSTER
    cd "output"
    qsub \
    -cwd \
    -pe smp 4 \
    -l mem_requested=10G \
    -P TumourProgression \
    -b y \
    -j y \
    -V \
    -N ${JOBNAME}\
    "${R} CMD BATCH \
    --no-save '--args \
    ${SEURATOBJECTPATH} \
    ${AUCELLPATH}' \
    ${SCRIPT}"

cd ${TEMPPWD}
