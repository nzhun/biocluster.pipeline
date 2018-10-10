#!/bin/bash
CONFIG_FILE="/share/shenlab/WES/PCGC_0706/src/config.sh"
source $CONFIG_FILE
source $REF_FILE

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Map FastQC
#WD=$PROJECT_DIR/QC/FASTQC
#mkdir -p $WD; cd $WD
#JOB_NUM=`wc -l $FASTQ_LIST|cut -d ' ' -f1`
#echo $JOB_NUM
#RUN_FASTQC=$(qsub -t 1-$JOB_NUM $FASTQC -i $FASTQ_LIST -r $REF_FILE)
#echo $RUN_FASTQC
#exit
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Map Fastq to Bam
#WD=$PROJECT_DIR/RAW_BAM
#mkdir -p $WD; cd $WD
#JOB_NUM=`wc -l $FASTQ_TABLE|cut -d ' ' -f1`
#RUN_BWA=$(qsub -t 1-$JOB_NUM $BWA -i $FASTQ_TABLE -r $REF_FILE)
#echo $RUN_BWA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RE-Map the Bam to Bam
#WD=$PROJECT_DIR/RAWBAM
#mkdir -p $WD; cd $WD
#JOB_NUM=`wc -l $OLD_BAM_LIST|cut -d ' ' -f1`
#RUN_REMAP=$(qsub -t 11-60 $REMAP -i $OLD_BAM_LIST -r $REF_FILE)
#echo $RUN_REMAP
#COLLECT_RESULT=$(qsub -hold_jid BamBWABam $COLLECT_BAM $CONFIG_FILE)
#echo $COLLECT_RESULT
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GATK Haplotype Caller
WD=$PROJECT_DIR/GVCF
mkdir -p $WD; cd $WD
#BAM_LIST=
JOB_NUM=`wc -l $BAM_LIST|cut -d ' ' -f1`
#RUN_HAPCALLER=$(qsub -t 1-$JOB_NUM -hold_jid Collect_BAMs $HAPCALLER -i $BAM_LIST -r $REF_FILE -t $TARGET_INTERVAL)
RUN_HAPCALLER=$(qsub -q shenlab.q -t 1-$JOB_NUM $HAPCALLER -i $BAM_LIST -r $REF_FILE -t $TARGET_INTERVAL)
echo $RUN_HAPCALLER
COLLECT_RESULT=$(qsub -hold_jid HaplotypeCaller_GVCFmode $COLLECT_GVCF $CONFIG_FILE)
echo $COLLECT_RESULT
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
#THIRD=$(qsub -t 1-10 -hold_jid job2.pbs job3.pbs)
#echo $THIRD
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

