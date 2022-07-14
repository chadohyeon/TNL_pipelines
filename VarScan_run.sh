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
varscan=${tools}/VarScan.v2.3.9/VarScan.v2.3.9.jar
gatk=${tools}/gatk-4.2.0.0/gatk
bgzip=${tools}/htslib/bin/bgzip
bcftools=${tools}/bcftools/bin/bcftools
tabix=${tools}/htslib/bin/tabix
samtools=${tools}/samtools/bin/samtools

### Set working directories
BAMDIR=${DIR}/01.bam/
SOMDIR=${DIR}/04.somatic_vcf
JAVA_SPACE=${SOMDIR}/java_space

SOMpassDIR=${DIR}/05.PASS_somatic_vcf

PBS_DIR=${DIR}/pbs_jobs
PBS_LOGS=${DIR}/pbs_stdouts

mkdir -p ${SOMDIR}
mkdir -p ${SOMDIR}/varscan_tmp
mkdir -p ${JAVA_SPACE}
mkdir -p ${SOMpassDIR}

mkdir -p ${PBS_DIR}
mkdir -p ${PBS_LOGS}

### DATE

DIR2REF=`python3 -c "print('/'.join('${REF}'.split('/')[:-2]))"`
GRCH=`python3 -c "print('${DIR2REF}'.split('/')[-1].split('_')[-1])"`

DATE=$(date '+%y%m%d.%H-%M-%S')

######################################################################
############################# VarScan ################################
######################################################################

### 01. Scattering Jobs
allRaw_snvs=""
allRaw_indels=""
itrs=""

for itr in `seq -f %04g 0 $(expr ${Cores} - 1)`; do
    itrs+="${itr} "    allRaw_snvs+="${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.somatic.snv.vcf.gz "
    allRaw_indels+="${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.somatic.indel.vcf.gz "
done

cat > ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-1.varscan.splitIntervals.${DATE}.pbs << EOF
#PBS -N ${somatic_sample}.pair_${matched_sample}_03-1.varscan.splitIntervals.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${somatic_sample}.pair_${matched_sample}_03-1.varscan.splitIntervals.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

python3 ${scriptsPath}/commaChromsToIntervalList.py\
 ${REF}.picard.ACGT.interval_list\
 ${commaChroms}\
 ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}.interval_list

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" SplitIntervals\
 -R ${REF}\
 -L ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}.interval_list\
 --scatter-count ${Cores}\
 -O ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks
 
python3 ${scriptsPath}/intervalListsToBed3.py\
 ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks
sleep 10
exit 0
EOF

FIRST_JOB=`qsub ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-1.varscan.splitIntervals.${DATE}.pbs`

### 02. Scattered variant calling
SECOND_MULTI_JOBS=""
for itr in ${itrs}; do
cat > ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-2.varscan.chunk_${itr}.${DATE}.pbs << EOF
#PBS -N ${somatic_sample}.pair_${matched_sample}_03-2.varscan.chunk_${itr}.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${somatic_sample}.pair_${matched_sample}_03-2.varscan.chunk_${itr}.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

rm -f ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${itr}-scattered.interval_list

${samtools} mpileup -q 1 -f ${REF}\
 -l ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${itr}-scattered.bed\
 ${BAMDIR}/${matched_sample}${bamSuffix} >\
 ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${matched_sample}.chunk_${itr}.mpileup
 
${samtools} mpileup -q 1 -f ${REF}\
 -l ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${itr}-scattered.bed\
 ${BAMDIR}/${somatic_sample}${bamSuffix} >\
 ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.chunk_${itr}.mpileup

java -Xmx${memPerCores}m -Djava.io.tmpdir=${JAVA_SPACE} -jar ${varscan} somatic\
 ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${matched_sample}.chunk_${itr}.mpileup\
 ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.chunk_${itr}.mpileup\
 ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}\
 --output-vcf 1

grep '#' ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.snp.vcf\
 > ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.somatic.snv.vcf
grep -v '#' ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.snp.vcf\
 | grep 'SOMATIC' >> ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.somatic.snv.vcf

${bgzip} ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.somatic.snv.vcf
${tabix} ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.somatic.snv.vcf.gz


grep '#' ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.indel.vcf\
 > ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.somatic.indel.vcf
grep -v '#' ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.indel.vcf\
 | grep 'SOMATIC' >> ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.somatic.indel.vcf

${bgzip} ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.somatic.indel.vcf
${tabix} ${SOMDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}_chunks/${somatic_sample}.pair_${matched_sample}.chunk_${itr}.somatic.indel.vcf.gz

sleep 10
exit 0
EOF
SECOND_MULTI_JOBS+=":"`qsub -W depend=afterok:${FIRST_JOB} ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-2.varscan.chunk_${itr}.${DATE}.pbs `
done

### 03. Remainders
cat > ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-3.varscan.merge.${DATE}.pbs << EOF
#PBS -N ${somatic_sample}.pair_${matched_sample}_03-3.varscan.merge.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${somatic_sample}.pair_${matched_sample}_03-3.varscan.merge.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

mkdir -p ${SOMpassDIR}/varscan_tmp
${bcftools} concat -a ${allRaw_snvs}\
 -Oz -o ${SOMpassDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}.PASS.SNV.varscan.vcf.gz
${tabix} -p vcf ${SOMpassDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}.PASS.SNV.varscan.vcf.gz

${bcftools} concat -a ${allRaw_indels}\
 -Oz -o ${SOMpassDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}.PASS.INDEL.varscan.vcf.gz
${tabix} -p vcf ${SOMpassDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}.PASS.INDEL.varscan.vcf.gz

${bcftools} concat -a ${SOMpassDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}.PASS.SNV.varscan.vcf.gz\
 ${SOMpassDIR}/varscan_tmp/${somatic_sample}.pair_${matched_sample}.PASS.INDEL.varscan.vcf.gz\
 -Oz -o ${SOMpassDIR}/${somatic_sample}.pair_${matched_sample}.PASS.MERGED.varscan.vcf.gz
 
${tabix} ${SOMpassDIR}/${somatic_sample}.pair_${matched_sample}.PASS.MERGED.varscan.vcf.gz
 
EOF
THIRD_JOB=`qsub -W depend=afterok${SECOND_MULTI_JOBS} ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-3.varscan.merge.${DATE}.pbs`

