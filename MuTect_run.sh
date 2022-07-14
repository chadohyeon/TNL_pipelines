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
memPerCores=`expr 1024 \* ${MEM} / ${Cores}`

### Tools
tools=/home/public/tools
scriptsPath=/home/public/scripts
mutect=${tools}/mutect-1.1.7/mutect-1.1.7.jar
gatk=${tools}/gatk-4.2.0.0/gatk
bcftools=${tools}/bcftools/bin/bcftools
bgzip=${tools}/htslib/bin/bgzip
tabix=${tools}/htslib/bin/tabix

### Set working directories
BAMDIR=${DIR}/01.bam
SOMDIR=${DIR}/04.somatic_vcf
JAVA_SPACE=${SOMDIR}/java_space

SOMpassDIR=${DIR}/05.PASS_somatic_vcf

PBS_DIR=${DIR}/pbs_jobs
PBS_LOGS=${DIR}/pbs_stdouts

mkdir -p ${SOMDIR}
mkdir -p ${SOMDIR}/mutect_tmp
mkdir -p ${JAVA_SPACE}
mkdir -p ${SOMpassDIR}

mkdir -p ${PBS_DIR}
mkdir -p ${PBS_LOGS}


### DATE

DIR2REF=`python3 -c "print('/'.join('${REF}'.split('/')[:-2]))"`
GRCH=`python3 -c "print('${DIR2REF}'.split('/')[-1].split('_')[-1])"`

DATE=$(date '+%y%m%d.%H-%M-%S')

java7=/usr/bin/java7

######################################################################
############################# MuTect #################################
######################################################################

### 01. Scattering Jobs
allRaw_vcfs=""
itrs=""

for itr in `seq -f %04g 0 $(expr ${Cores} - 1)`; do
    itrs+="${itr} "
    allRaw_vcfs+="${SOMDIR}/mutect_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.mutect.filtered.vcf.gz "
done

cat > ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-1.mutect.splitIntervals.${DATE}.pbs << EOF
#PBS -N ${somatic_sample}.pair_${matched_sample}_03-1.mutect.splitIntervals.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${somatic_sample}.pair_${matched_sample}_03-1.mutect.splitIntervals.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

python3 ${scriptsPath}/commaChromsToIntervalList.py\
 ${REF}.picard.ACGT.interval_list\
 ${commaChroms}\
 ${SOMDIR}/mutect_tmp/${somatic_sample}.pair_${matched_sample}.interval_list

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" SplitIntervals\
 -R ${REF}\
 -L ${SOMDIR}/mutect_tmp/${somatic_sample}.pair_${matched_sample}.interval_list\
 --scatter-count ${Cores}\
 -O ${SOMDIR}/mutect_tmp/${somatic_sample}.pair_${matched_sample}_chunks
sleep 10
exit 0
EOF

FIRST_JOB=`qsub ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-1.mutect.splitIntervals.${DATE}.pbs`

### 02. Scattered variant calling
SECOND_MULTI_JOBS=""

for itr in ${itrs}; do
cat > ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-2.mutect.chunk_${itr}.${DATE}.pbs << EOF
#PBS -N ${somatic_sample}.pair_${matched_sample}_03-2.mutect.chunk_${itr}.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${somatic_sample}.pair_${matched_sample}_03-2.mutect.chunk_${itr}.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

${java7} -Xmx${memPerCores}m -Djava.io.tmpdir=${JAVA_SPACE} -jar ${mutect} --analysis_type MuTect\
 --reference_sequence ${REF}\
 --cosmic ${DIR2REF}/annot/Cosmic.${GRCH}.MuTect.vcf.gz\
 --dbsnp ${DIR2REF}/annot/dbSNP154.${GRCH}.vcf.gz\
 --intervals ${SOMDIR}/mutect_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${itr}-scattered.interval_list\
 --input_file:tumor ${BAMDIR}/${somatic_sample}${bamSuffix}\
 --input_file:normal ${BAMDIR}/${matched_sample}${bamSuffix}\
 --vcf ${SOMDIR}/mutect_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.mutect.filtered.vcf

${bgzip} ${SOMDIR}/mutect_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.mutect.filtered.vcf
${tabix} ${SOMDIR}/mutect_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.mutect.filtered.vcf.gz
sleep 10
exit 0
EOF

SECOND_MULTI_JOBS+=":"`qsub -W depend=afterok:${FIRST_JOB} ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-2.mutect.chunk_${itr}.${DATE}.pbs `
done

### 03. Remainders
cat > ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-3.mutect.merge.${DATE}.pbs << EOF
#PBS -N ${somatic_sample}.pair_${matched_sample}_03-3.mutect.merge.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${somatic_sample}.pair_${matched_sample}_03-3.mutect.merge.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

${bcftools} concat -a ${allRaw_vcfs}\
 -Oz -o ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.mutect.filtered.vcf.gz

${tabix} -p vcf ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.mutect.filtered.vcf.gz

${bcftools} view -f PASS ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.mutect.filtered.vcf.gz -Oz \
        -o ${SOMpassDIR}/${somatic_sample}.pair_${matched_sample}.PASS.mutect.vcf.gz

${tabix} -p vcf ${SOMpassDIR}/${somatic_sample}.pair_${matched_sample}.PASS.mutect.vcf.gz
EOF
THIRD_JOB=`qsub -W depend=afterok${SECOND_MULTI_JOBS} ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-3.mutect.merge.${DATE}.pbs`
