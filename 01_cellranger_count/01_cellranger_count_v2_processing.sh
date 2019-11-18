# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 01 CELLRANGER COUNT PROCESSING
#
# GLOBS
## CELLRANGER VARIABLES
CELLRANGER="/share/ClusterShare/software/contrib/CTP_single_cell/cellranger/cellranger-2.2.0/cellranger"
TRANSCRIPTOME="/share/ClusterShare/software/contrib/CTP_single_cell/cellranger/refdata/refdata-cellranger-GRCh38-1.2.0/"
FASTQPATH="/share/ScratchGeneral/sunwu/projects/Stromal_heterogeneity_paper/REPROCESSING_FROM_RAW/data/"
CELLSTOEXPECT=6000
#
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
TEMPPWD=$(pwd)
for SAMPLENAME in ${SAMPLEIDS}; do
  echo ${SAMPLENAME}

  COUNTJOBNAME="count_${SAMPLENAME}"
  OUTPUT_ID_STRING="count_${SAMPLENAME}_GRCh38"

  cd "output/${SAMPLENAME}"
      qsub \
      -cwd \
      -pe smp 24 \
      -l mem_requested=15G \
      -P TumourProgression \
      -b y \
      -j y \
      -V \
      -N ${COUNTJOBNAME} \
      "${CELLRANGER} count \
      --id=${OUTPUT_ID_STRING} \
      --sample=${SAMPLENAME} \
      --transcriptome=${TRANSCRIPTOME} \
      --fastqs=${FASTQPATH} \
      --expect-cells=${CELLSTOEXPECT}"

    cd ${TEMPPWD}

  done
