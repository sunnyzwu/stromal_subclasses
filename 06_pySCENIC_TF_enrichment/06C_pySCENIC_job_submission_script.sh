# STROMAL SUBCLASSES MANUSCRIPT
# author: Sunny Z. Wu
# dated: 20191115
#
# SCRIPT 06C PYSCENIC JOB SUBMISSION TO HPC SCRIPT
#
  export PATH=/share/ClusterShare/software/contrib/CTP_single_cell/tools/anaconda2/bin:$PATH
  source activate py36
#
# INPUTS
    PROJECTNAME="Stromal_subclasses"
  ## TF annotation table
    ANNFNAME="/share/ClusterShare/software/contrib/CTP_single_cell/SCENIC/human/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
  ## TF database
    DATABASE="/share/ClusterShare/software/contrib/CTP_single_cell/SCENIC/human/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
  ## TF names file
    TFFILEPATH="/share/ClusterShare/software/contrib/CTP_single_cell/SCENIC/human/human_TFs.txt"
  ## MATRIXNAME
    MATRIXNAME="./output/data/Stromal_normalized_expression_matrix_filtered.csv"

# STEP 1; GRNBOOST CO-EXPRESSION

PYSCENICJOBID1="output/PYSC1_${PROJECTNAME}_GRNBOOST_jobid.txt"
PYSCENICJOBNAME1="PYSC1_${PROJECTNAME}"

qsub \
-cwd \
-pe smp 16 \
-l mem_requested=20G \
-b y \
-j y \
-V \
-P TumourProgression \
-N $PYSCENICJOBNAME1 \
pyscenic \
grnboost ${MATRIXNAME} ${TFFILEPATH} \
-t \
-o output/grnboost_output.tsv \
--num_workers 8 \
> $PYSCENICJOBID1

# STEP 2; MOTIF ENRICHMENT

PYSCENICJOBID2="output/PYSC2_${PROJECTNAME}_MOTIFENRICH_jobid.txt"
PYSCENICJOBNAME2="PYSC2_${PROJECTNAME}"
JOBIDHOLD=$(cut -d ' ' -f 3 output/PYSC1_${PROJECTNAME}_GRNBOOST_jobid.txt)

qsub \
-cwd \
-pe smp 16 \
-l mem_requested=10G \
-b y \
-j y \
-V \
-P TumourProgression \
-hold_jid $JOBIDHOLD \
-N $PYSCENICJOBNAME2 \
pyscenic ctx \
--annotations_fname ${ANNFNAME} \
--expression_mtx_fname ${MATRIXNAME} \
-t \
--all_modules \
--num_workers 8 \
-o output/ctx_output.csv \
output/grnboost_output.tsv \
${DATABASE} \
> $PYSCENICJOBID2

# STEP 3; AUCELL

PYSCENICJOBID3="output/PYSC3_${PROJECTNAME}_AUCELL_jobid.txt"
PYSCENICJOBNAME3="PYSC3_${PROJECTNAME}"
JOBIDHOLD=$(cut -d ' ' -f 3 output/PYSC2_${PROJECTNAME}_MOTIFENRICH_jobid.txt)

qsub \
-cwd \
-pe smp 16 \
-l mem_requested=10G \
-b y \
-j y \
-P TumourProgression \
-V \
-hold_jid $JOBIDHOLD \
-N $PYSCENICJOBNAME3 \
pyscenic aucell \
${MATRIXNAME} \
output/ctx_output.csv \
-t \
--num_workers 8 \
-o output/aucell_matrix.csv
> $PYSCENICJOBID3

echo all jobs submitted !
