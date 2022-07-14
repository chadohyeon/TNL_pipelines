#!/bin/bash

### Made by Do Hyeon Cha TNL@KAIST

### Args from wrapper_caller.py

MEM=${1}
DIR=${2}
REF=${3}
PON=${4}
bamSuffix=${5}
commaChroms=${6}
Cores=${7}
Queue=${8}
matched_sample=${9}
somatic_sample=${10}

### Other parameters
nativeHmmThread=4
memPerCores=`expr 1024 \* ${MEM} / ${Cores}`

### Tools
tools=/home/public/tools
scriptsPath=/home/public/scripts
gatk=${tools}/gatk-4.2.0.0/gatk
bcftools=${tools}/bcftools/bin/bcftools
tabix=${tools}/htslib/bin/tabix

# set working dirs and files
BAMDIR=${DIR}/01.bam
SOMDIR=${DIR}/04.somatic_vcf
JAVA_SPACE=${SOMDIR}/java_space

SOMpassDIR=${DIR}/05.PASS_somatic_vcf

PBS_DIR=${DIR}/pbs_jobs
PBS_LOGS=${DIR}/pbs_stdouts

mkdir -p ${SOMDIR}
mkdir -p ${SOMDIR}/mutect2_tmp
mkdir -p ${JAVA_SPACE}
mkdir -p ${SOMpassDIR}

mkdir -p ${PBS_DIR}
mkdir -p ${PBS_LOGS}

### DATE

DIR2REF=`python3 -c "print('/'.join('${REF}'.split('/')[:-2]))"`
GRCH=`python3 -c "print('${DIR2REF}'.split('/')[-1].split('_')[-1])"`

DATE=$(date '+%y%m%d.%H-%M-%S')

######################################################################
############################# Mutect2 ################################
######################################################################

### 01. Scattering Jobs
allRaw_vcfs=""
itrs=""
stats_to_merge=""
all_f1r2_input=""

for itr in `seq -f %04g 0 $(expr ${Cores} - 1)`; do
    itrs+="${itr} "
        allRaw_vcfs+="${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.mutect2.unfiltered.vcf.gz "
    stats_to_merge+="-stats ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.mutect2.unfiltered.vcf.gz.stats "
    all_f1r2_input+="-I ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}-f1r2.tar.gz "
done

cat > ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-1.mutect2.splitIntervals.${DATE}.pbs << EOF
#PBS -N ${somatic_sample}.pair_${matched_sample}_03-1.mutect2.splitIntervals.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${somatic_sample}.pair_${matched_sample}_03-1.mutect2.splitIntervals.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

python3 ${scriptsPath}/commaChromsToIntervalList.py\
 ${REF}.picard.ACGT.interval_list\
 ${commaChroms}\
 ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.interval_list

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" SplitIntervals\
 -R ${REF}\
 -L ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.interval_list\
 --scatter-count ${Cores}\
 -O ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}_chunks
sleep 10
exit 0
EOF

FIRST_JOB=`qsub ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-1.mutect2.splitIntervals.${DATE}.pbs`

### 02. Scattered variant calling
SECOND_MULTI_JOBS=""

for itr in ${itrs}; do
cat > ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-2.mutect2.chunk_${itr}.${DATE}.pbs << EOF
#PBS -N ${somatic_sample}.pair_${matched_sample}_03-2.mutect2.chunk_${itr}.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${somatic_sample}.pair_${matched_sample}_03-2.mutect2.chunk_${itr}.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

${gatk} --java-options "-Xmx${memPerCores}m -Djava.io.tmpdir=${JAVA_SPACE}" Mutect2\
 --native-pair-hmm-threads ${nativeHmmThread}\
 -R ${REF}\
 -I ${BAMDIR}/${matched_sample}${bamSuffix}\
 -normal ${rgPrefix}${matched_sample}\
 -I ${BAMDIR}/${somatic_sample}${bamSuffix}\
 -tumor ${rgPrefix}${somatic_sample}\
 --germline-resource  ${DIR2REF}/annot/af-only-gnomad.${GRCH}.vcf.gz\
 --panel-of-normals ${PON}\
 -L ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${itr}-scattered.interval_list\
 --f1r2-tar-gz ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}-f1r2.tar.gz\
 -O ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.mutect2.unfiltered.vcf.gz
sleep 10
exit 0
EOF

SECOND_MULTI_JOBS+=":"`qsub -W depend=afterok:${FIRST_JOB} ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-2.mutect2.chunk_${itr}.${DATE}.pbs `
done


### 03. Remainders
cat > ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-3.mutect2.merge.${DATE}.pbs << EOF
#PBS -N ${somatic_sample}.pair_${matched_sample}_03-3.mutect2.merge.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${somatic_sample}.pair_${matched_sample}_03-3.mutect2.merge.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

${bcftools} concat -a ${allRaw_vcfs}\
 -Oz -o ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.mutect2.unfiltered.vcf.gz

${tabix} -p vcf ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.mutect2.unfiltered.vcf.gz

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" MergeMutectStats\
 ${stats_to_merge}\
 -O ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.mutect2.unfiltered.vcf.gz.stats

## Jointly learning read-orientation model
${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}"  LearnReadOrientationModel\
 ${all_f1r2_input}\
 -O ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.read-orientation-model.tar.gz


## Estimate contamination and errors
${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" GetPileupSummaries\
 -I ${BAMDIR}/${somatic_sample}${bamSuffix}\
 -V ${DIR2REF}/annot/small_exac_common_3.${GRCH}.vcf.gz\
 -L ${DIR2REF}/annot/small_exac_common_3.${GRCH}.vcf.gz\
 -O ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.matchedSomatic.getpileupsummaries.table

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" GetPileupSummaries\
 -I ${BAMDIR}/${matched_sample}${bamSuffix}\
 -V ${DIR2REF}/annot/small_exac_common_3.${GRCH}.vcf.gz\
 -L ${DIR2REF}/annot/small_exac_common_3.${GRCH}.vcf.gz\
 -O ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.matchedNormal.getpileupsummaries.table

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" CalculateContamination\
 -I ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.matchedSomatic.getpileupsummaries.table\
 -matched ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.matchedNormal.getpileupsummaries.table\
 -O ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.calculatecontamination.table\
 --tumor-segmentation ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.segments.table

## FilterMutectCalls
${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" FilterMutectCalls\
 -R ${REF}\
 -V ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.mutect2.unfiltered.vcf.gz\
 --tumor-segmentation ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.segments.table\
 --contamination-table ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.calculatecontamination.table\
 --ob-priors ${SOMDIR}/mutect2_tmp/${somatic_sample}.pair_${matched_sample}.read-orientation-model.tar.gz\
 -O ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.mutect2.filtered.vcf.gz

## PASS calls

${bcftools} view -f PASS ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.mutect2.filtered.vcf.gz -Oz\
 -o ${SOMpassDIR}/${somatic_sample}.pair_${matched_sample}.PASS.mutect2.filtered.vcf.gz

${tabix} -p vcf ${SOMpassDIR}/${somatic_sample}.pair_${matched_sample}.PASS.mutect2.filtered.vcf.gz
EOF

THIRD_JOB=`qsub -W depend=afterok${SECOND_MULTI_JOBS} ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-3.mutect2.merge.${DATE}.pbs`
