# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 02A SEURAT V2 JOB SUBMISSION TO HPC SCRIPT
#
# GLOBS
source activate Renv # for UMAP
R="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.5.0/bin/R"
SPECIES='human'
TEMPPWD=$(pwd)
SEURATSCRIPT="${TEMPPWD}/02B_R_seurat_processing_individual_objects.R"
CELLRANGEROUTPUTPATH="/share/ScratchGeneral/sunwu/projects/Stromal_heterogeneity_paper/REPROCESSING_FROM_RAW/01_cellranger_count/"
# SAMPLE IDS
SAMPLEIDS="CID44041_P1 CID44971_P2 CID44991_P3 CID4513_P4 CID4515_P5"
#
#
# DIRECTORIES
mkdir "output"
for SAMPLENAME in ${SAMPLEIDS}; do
  mkdir "output/${SAMPLENAME}"
done
#
# SUBMIT JOBS TO CLUSTER
for SAMPLENAME in ${SAMPLEIDS}; do
    echo ${SAMPLENAME}

    SEURATJOBNAME="seurat_${SAMPLENAME}"
    PROJECTNAME=${SAMPLENAME}
    MATRIXPATH="/share/ScratchGeneral/sunwu/projects/Stromal_heterogeneity_paper/REPROCESSING_FROM_RAW/01_cellranger_count/output/${SAMPLENAME}/count_${SAMPLENAME}_GRCh38/outs/raw_gene_bc_matrices/GRCh38/"

    cd "output/${SAMPLENAME}/"
    qsub \
    -cwd \
    -pe smp 4 \
    -l mem_requested=10G \
    -P TumourProgression \
    -b y \
    -j y \
    -V \
    -N $SEURATJOBNAME\
    "${R} CMD BATCH \
    --no-save '--args \
    ${PROJECTNAME} \
    ${SPECIES} \
    ${MATRIXPATH}' \
    ${SEURATSCRIPT}"

    cd ${TEMPPWD}
  done
