#!/bin/bash

### Made by Do Hyeon Cha TNL@KAIST
### Inspired by IBM workflows (https://www.ibm.com/downloads/cas/ZJQD0QAL)
### Args from wrapper_aligner.py
### Scatter-and-gather approach for BQSR, to boost up speed

MEM=${1}
DIR=${2}
REF=${3}
dbSNP=${4}
Mills_1000G=${5}

sample=${6}
FQ_F=${7}
FQ_R=${8}
Cores=${9}
Queue=${10}

### Other parameters
bwaThreadMultiple=2
bwaThread=`expr ${Cores} \* ${bwaThreadMultiple}`
memPerCores=`expr 1024 \* ${MEM} / ${Cores}`

### Tools

tools=/home/public/tools
bwa=${tools}/bwa-0.7.17/bwa
samtools=${tools}/samtools/bin/samtools
gatk=${tools}/gatk-4.2.0.0/gatk
Picard=${tools}/picard-2.23.9/picard.jar
bedtools=${tools}/bedtools2/bin/bedtools


### Set working directories
FQDIR=${DIR}/00.fastq
BAMDIR=${DIR}/01.bam
JAVA_SPACE=${BAMDIR}/java_space
SPARK_SPACE=${BAMDIR}/spark_space
PBS_DIR=${DIR}/pbs_jobs
PBS_LOGS=${DIR}/pbs_stdouts

mkdir -p ${BAMDIR}
mkdir -p ${BAMDIR}/tmp
mkdir -p ${SPARK_SPACE}
mkdir -p ${JAVA_SPACE}
mkdir -p ${PBS_DIR}
mkdir -p ${PBS_LOGS}

### DATE

DATE=$(date '+%y%m%d.%H-%M-%S')

###########################################################################
########################### GATK Best Practices ###########################
###########################################################################

### 0. Prep for multi-core works
${gatk} --java-options "-Xmx${MEM}g  -Djava.io.tmpdir=${JAVA_SPACE}" SplitIntervals\
 -R ${REF}\
 -L ${REF}.picard.ACGT.interval_list\
 --scatter-count ${Cores}\
 -O ${BAMDIR}/tmp/${sample}_BQSR_chunks

if [ "${Cores}" -gt 1 ]; then
for itr in `seq -f %04g 1 $(expr ${Cores} - 1)`; do
endToDiscard=`grep -vP "^@" ${BAMDIR}/tmp/${sample}_BQSR_chunks/${itr}-scattered.interval_list | head -n 1 | awk '{print $1,$2}'`
python3 -c "cp='${endToDiscard}'.split(); wf=open('${BAMDIR}/tmp/${sample}_BQSR_chunks/${itr}-discardEnd.bed','w'); wf.write(cp[0]+'\t'+str(int(cp[1])-2)+'\t'+str(int(cp[1])-1)); wf.close()"
done
fi

### 1. BWA-MEM + sortedBAM
cat > ${PBS_DIR}/${sample}_01-1.align_sort.${DATE}.pbs << EOF
#PBS -N ${sample}_01-1.align_sort.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${sample}_01-1.align_sort.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=${Cores}

echo "QUERY:01.BWA+SamToSortedBam:Sample=${sample}"
${bwa} mem -t ${bwaThread} -Ma -R "@RG\tID:${sample}\tLB:${sample}\tSM:${sample}\tPL:ILLUMINA" ${REF} ${FQDIR}/${FQ_F} ${FQDIR}/${FQ_R} > ${BAMDIR}/tmp/${sample}.sam
${samtools} sort -@ ${Cores} -m ${memPerCores}M -T ${JAVA_SPACE} ${BAMDIR}/tmp/${sample}.sam >  ${BAMDIR}/tmp/${sample}.bam
${samtools} index ${BAMDIR}/tmp/${sample}.bam
sleep 10
exit 0

EOF
FIRST_JOB=`qsub ${PBS_DIR}/${sample}_01-1.align_sort.${DATE}.pbs`

### 2. MarkDuplicate Picard (Spark doesn't work well)
cat > ${PBS_DIR}/${sample}_01-2.mark_split.${DATE}.pbs << EOF
#PBS -N ${sample}_01-2.mark_split.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${sample}_01-2.mark_split.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

${gatk} --java-options "-Xmx${MEM}g -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${JAVA_SPACE}" MarkDuplicates\
 -I ${BAMDIR}/tmp/${sample}.bam\
 -O ${BAMDIR}/tmp/${sample}.marked.bam\
 -M ${BAMDIR}/tmp/${sample}.marked_dup_metrics.txt\
 --MAX_RECORDS_IN_RAM 5000000\
 -MAX_SEQS 5000000\
 --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500\
 -VALIDATION_STRINGENCY SILENT\
 -MAX_FILE_HANDLES 1000

${samtools} index ${BAMDIR}/tmp/${sample}.marked.bam


${gatk} --java-options "-Xmx${MEM}g  -Djava.io.tmpdir=${JAVA_SPACE}" SplitIntervals\
 -R ${REF}\
 -L ${REF}.picard.ACGT.interval_list\
 --scatter-count ${Cores}\
 -O ${BAMDIR}/tmp/${sample}_BQSR_chunks
sleep 10
exit 0

EOF
SECOND_JOB=`qsub -W depend=afterok:${FIRST_JOB} ${PBS_DIR}/${sample}_01-2.mark_split.${DATE}.pbs`

### 3. BQSR Scattered jobs
THIRD_MULTI_JOBS=""
allBQSR_Bams=""

for itr in `seq -f %04g 0 $(expr ${Cores} - 1)`; do
    allBQSR_Bams+="${BAMDIR}/tmp/${sample}_BQSR_chunks/${sample}.marked.${itr}_recal.bam "
done

cat > ${PBS_DIR}/${sample}_01-3.BQSR_chunk0000.${DATE}.pbs << EOF
#PBS -N ${sample}_01-3.BQSR_chunk0000.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${sample}_01-3.BQSR_chunk0000.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

echo "QUERY:03.BQSR_calculate_scattered:Sample=${sample}:Chunk:0000"
${gatk} --java-options "-Xmx${memPerCores}m -Djava.io.tmpdir=${JAVA_SPACE}" BaseRecalibrator\
 -R ${REF}\
 -I ${BAMDIR}/tmp/${sample}.marked.bam\
 -L ${BAMDIR}/tmp/${sample}_BQSR_chunks/0000-scattered.interval_list\
 --known-sites ${Mills_1000G}\
 --known-sites ${dbSNP}\
 -O ${BAMDIR}/tmp/${sample}_BQSR_chunks/${sample}.recal_data.0000.grp

echo "QUERY:04.BQSR_APPLY_scattered:Sample=${sample}:Chunk:0000"
${gatk} --java-options "-Xmx${memPerCores}m -Djava.io.tmpdir=${JAVA_SPACE}" ApplyBQSR\
 -R ${REF}\
 -I ${BAMDIR}/tmp/${sample}.marked.bam\
 -L ${BAMDIR}/tmp/${sample}_BQSR_chunks/0000-scattered.interval_list\
 --bqsr-recal-file ${BAMDIR}/tmp/${sample}_BQSR_chunks/${sample}.recal_data.0000.grp\
 -O ${BAMDIR}/tmp/${sample}_BQSR_chunks/${sample}.marked.0000_recal.bam
sleep 10
exit 0
EOF
THIRD_MULTI_JOBS+=":"`qsub -W depend=afterok:${SECOND_JOB} ${PBS_DIR}/${sample}_01-3.BQSR_chunk0000.${DATE}.pbs`

if [ "${Cores}" -gt 1 ]; then
for itr in `seq -f %04g 1 $(expr ${Cores} - 1)`; do

cat > ${PBS_DIR}/${sample}_01-3.BQSR_chunk${itr}.${DATE}.pbs << EOF
#PBS -N ${sample}_01-3.BQSR_chunk${itr}.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${sample}_01-3.BQSR_chunk${itr}.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

echo "QUERY:03.BQSR_calculate_scattered:Sample=${sample}:Chunk:${itr}"
${gatk} --java-options "-Xmx${memPerCores}m -Djava.io.tmpdir=${JAVA_SPACE}" BaseRecalibrator\
 -R ${REF}\
 -I ${BAMDIR}/tmp/${sample}.marked.bam\
 -L ${BAMDIR}/tmp/${sample}_BQSR_chunks/${itr}-scattered.interval_list\
 --known-sites ${Mills_1000G}\
 --known-sites ${dbSNP}\
 -O ${BAMDIR}/tmp/${sample}_BQSR_chunks/${sample}.recal_data.${itr}.grp

echo "QUERY:04.BQSR_APPLY_scattered:Sample=${sample}:Chunk:${itr}"
${gatk} --java-options "-Xmx${memPerCores}m -Djava.io.tmpdir=${JAVA_SPACE}" ApplyBQSR\
 -R ${REF}\
 -I ${BAMDIR}/tmp/${sample}.marked.bam\
 -L ${BAMDIR}/tmp/${sample}_BQSR_chunks/${itr}-scattered.interval_list\
 --bqsr-recal-file ${BAMDIR}/tmp/${sample}_BQSR_chunks/${sample}.recal_data.${itr}.grp\
 -O ${BAMDIR}/tmp/${sample}_BQSR_chunks/tmp_${sample}.marked.${itr}_recal.bam

${bedtools} intersect\
 -abam ${BAMDIR}/tmp/${sample}_BQSR_chunks/tmp_${sample}.marked.${itr}_recal.bam\
 -b ${BAMDIR}/tmp/${sample}_BQSR_chunks/${itr}-discardEnd.bed\
 -v > ${BAMDIR}/tmp/${sample}_BQSR_chunks/${sample}.marked.${itr}_recal.bam
${samtools} index ${BAMDIR}/tmp/${sample}_BQSR_chunks/${sample}.marked.${itr}_recal.bam

rm ${BAMDIR}/tmp/${sample}_BQSR_chunks/tmp_${sample}.marked.${itr}_recal.bam
rm ${BAMDIR}/tmp/${sample}_BQSR_chunks/tmp_${sample}.marked.${itr}_recal.bai

sleep 10
exit 0
EOF

THIRD_MULTI_JOBS+=":"`qsub -W depend=afterok:${SECOND_JOB} ${PBS_DIR}/${sample}_01-3.BQSR_chunk${itr}.${DATE}.pbs`
done
fi

### 3. Gathering BQSR scatters
cat > ${PBS_DIR}/${sample}_01-4.Gather_BQSR.${DATE}.pbs << EOF
#PBS -N ${sample}_01-4.Gather_BQSR.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${sample}_01-4.Gather_BQSR.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1
echo "QUERY:05.Gathering_BQSR_chunks:Sample=${sample}"
${samtools} cat -o ${BAMDIR}/${sample}.marked.recal.bam ${allBQSR_Bams}
${samtools} index ${BAMDIR}/${sample}.marked.recal.bam
EOF

FOURTH_JOB=`qsub -W depend=afterok${THIRD_MULTI_JOBS} ${PBS_DIR}/${sample}_01-4.Gather_BQSR.${DATE}.pbs`