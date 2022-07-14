#!/bin/bash

### Made by Do Hyeon Cha, TNL@KAIST

### Args from wrapper_caller.py

MEM=${1}
DIR=${2}
REF=${3}
PON=${4}
bamSuffix=${5}
commaChroms=${6}
Cores=${7}
Queue=${8}
Sample=${9}

### Other parameters
nativeHmmThread=4
memPerCores=`expr 1024 \* ${MEM} / ${Cores}`

### GATK
tools=/home/public/tools
scriptsPath=/home/public/scripts
gatk=${tools}/gatk-4.2.0.0/gatk
bcftools=${tools}/bcftools/bin/bcftools
tabix=${tools}/htslib/bin/tabix

### Set working directories
BAMDIR=${DIR}/01.bam
VCFDIR=${DIR}/02.germline_vcf
JAVA_SPACE=${VCFDIR}/java_space

VCFpassDIR=${DIR}/03.PASS_germline_vcf

PBS_DIR=${DIR}/pbs_jobs
PBS_LOGS=${DIR}/pbs_stdouts

mkdir -p ${VCFDIR}
mkdir -p ${VCFDIR}/haplotypecaller_tmp
mkdir -p ${VCFpassDIR}
mkdir -p ${JAVA_SPACE}

mkdir -p ${PBS_DIR}
mkdir -p ${PBS_LOGS}

### DATE

DIR2REF=`python3 -c "print('/'.join('${REF}'.split('/')[:-2]))"`
GRCH=`python3 -c "print('${DIR2REF}'.split('/')[-1].split('_')[-1])"`

DATE=$(date '+%y%m%d.%H-%M-%S')

######################################################################
########################## HaplotypeCaller ###########################
######################################################################

### 01. Scattering Jobs
allRaw_vcfs=""
itrs=""

for itr in `seq -f %04g 0 $(expr ${Cores} - 1)`; do
    itrs+="${itr} "
    allRaw_vcfs+="${VCFDIR}/haplotypecaller_tmp/${Sample}_chunks/${Sample}.chunk_${itr}.haplotypecaller.vcf.gz "
done

cat > ${PBS_DIR}/${Sample}_02-1.haplotypecaller.splitIntervals.${DATE}.pbs << EOF
#PBS -N ${Sample}_02-1.haplotypecaller.splitIntervals.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${Sample}_02-1.haplotypecaller.splitIntervals.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

python3 ${scriptsPath}/commaChromsToIntervalList.py\
 ${REF}.picard.ACGT.interval_list\
 ${commaChroms}\
 ${VCFDIR}/haplotypecaller_tmp/${Sample}.interval_list

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" SplitIntervals\
 -R ${REF}\
 -L ${VCFDIR}/haplotypecaller_tmp/${Sample}.interval_list\
 --scatter-count ${Cores}\
 -O ${VCFDIR}/haplotypecaller_tmp/${Sample}_chunks
sleep 10
exit 0
EOF

FIRST_JOB=`qsub ${PBS_DIR}/${Sample}_02-1.haplotypecaller.splitIntervals.${DATE}.pbs`


### 02. Scattered variant calling
SECOND_MULTI_JOBS=""

for itr in ${itrs}; do
cat > ${PBS_DIR}/${Sample}_02-2.haplotypecaller.chunk_${itr}.${DATE}.pbs << EOF
#PBS -N ${Sample}_02-2.haplotypecaller.chunk_${itr}.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${Sample}_02-2.haplotypecaller.chunk_${itr}.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

${gatk} --java-options "-Xmx${memPerCores}m -Djava.io.tmpdir=${JAVA_SPACE}" HaplotypeCaller\
 --native-pair-hmm-threads ${nativeHmmThread}\
 -R ${REF}\
 -L ${VCFDIR}/haplotypecaller_tmp/${Sample}_chunks/${itr}-scattered.interval_list\
 -I ${BAMDIR}/${Sample}${bamSuffix}\
 --dbsnp ${DIR2REF}/annot/dbSNP154.${GRCH}.vcf.gz\
 -O ${VCFDIR}/haplotypecaller_tmp/${Sample}_chunks/${Sample}.chunk_${itr}.haplotypecaller.vcf.gz
sleep 10
exit 0
EOF

SECOND_MULTI_JOBS+=":"`qsub -W depend=afterok:${FIRST_JOB} ${PBS_DIR}/${Sample}_02-2.haplotypecaller.chunk_${itr}.${DATE}.pbs`
done


### 03. Remainders
cat > ${PBS_DIR}/${Sample}_02-3.haplotypecaller.merge.${DATE}.pbs << EOF
#PBS -N ${Sample}_02-3.haplotypecaller.merge.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${Sample}_02-3.haplotypecaller.merge.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1


${bcftools} concat -a ${allRaw_vcfs}\
-Oz -o ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.vcf.gz

${tabix} -p vcf ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.vcf.gz


### 02. Extract the SNV from the call set
${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" SelectVariants\
 -R ${REF}\
 -V ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.vcf.gz\
 --select-type-to-include SNP\
 -O ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.SNV_select.vcf.gz

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" VariantFiltration\
 -R ${REF}\
 -V ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.SNV_select.vcf.gz\
 -filter "QD < 2.0" --filter-name "QD_under2"\
 -filter "QUAL < 30.0" --filter-name "QUAL_under30"\
 -filter "FS > 60.0" --filter-name "FS_over60"\
 -filter "MQ < 40.0" --filter-name "MQ_under40"\
 -filter "MQRankSum < -12.5" --filter-name "MQRS_under-12.5"\
 -filter "ReadPosRankSum < -8.0" --filter-name "RPRS_under-8"\
 -O ${VCFDIR}/${Sample}.haplotypecaller.SNV_select.filter_marked.vcf.gz

### 03. Extract the INDEL from the call set
${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" SelectVariants\
 -R ${REF}\
 -V ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.vcf.gz\
 --select-type-to-include INDEL\
 -O ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.INDEL_select.vcf.gz


${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" VariantFiltration\
 -R ${REF}\
 -V ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.INDEL_select.vcf.gz\
 -filter "QD < 2.0" --filter-name "QD_under2"\
 -filter "QUAL < 30.0" --filter-name "QUAL_under30"\
 -filter "FS > 200.0" --filter-name "FS_over200"\
 -filter "ReadPosRankSum < -8.0" --filter-name "RPRS_under-8"\
 -O ${VCFDIR}/${Sample}.haplotypecaller.INDEL_select.filter_marked.vcf.gz

### 04. Grep PASS-only rows and index again
${bcftools} view -f PASS\
 ${VCFDIR}/${Sample}.haplotypecaller.SNV_select.filter_marked.vcf.gz -Oz\
 -o ${VCFpassDIR}/${Sample}.PASS.SNV.haplotypecaller.vcf.gz
${tabix} -p vcf ${VCFpassDIR}/${Sample}.PASS.SNV.haplotypecaller.vcf.gz

${bcftools} view -f PASS\
 ${VCFDIR}/${Sample}.haplotypecaller.INDEL_select.filter_marked.vcf.gz -Oz\
 -o ${VCFpassDIR}/${Sample}.PASS.INDEL.haplotypecaller.vcf.gz
${tabix} -p vcf ${VCFpassDIR}/${Sample}.PASS.INDEL.haplotypecaller.vcf.gz

### 05. Merge SNV and INDEL passed calls
${bcftools} concat -a ${VCFDIR}/${Sample}.haplotypecaller.SNV_select.filter_marked.vcf.gz\
 ${VCFDIR}/${Sample}.haplotypecaller.INDEL_select.filter_marked.vcf.gz\
 -Oz -o ${VCFDIR}/${Sample}.haplotypecaller.MERGED.filter_marked.vcf.gz
${tabix} -p vcf ${VCFDIR}/${Sample}.haplotypecaller.MERGED.filter_marked.vcf.gz


${bcftools} concat -a ${VCFpassDIR}/${Sample}.PASS.SNV.haplotypecaller.vcf.gz\
 ${VCFpassDIR}/${Sample}.PASS.INDEL.haplotypecaller.vcf.gz\
 -Oz -o ${VCFpassDIR}/${Sample}.PASS.MERGED.haplotypecaller.vcf.gz

${tabix} -p vcf ${VCFpassDIR}/${Sample}.PASS.MERGED.haplotypecaller.vcf.gz

mkdir -p ${VCFpassDIR}/haplotypecaller_tmp
mv ${VCFpassDIR}/${Sample}.PASS.SNV.haplotypecaller.vcf.gz* ${VCFpassDIR}/haplotypecaller_tmp
mv ${VCFpassDIR}/${Sample}.PASS.INDEL.haplotypecaller.vcf.gz* ${VCFpassDIR}/haplotypecaller_tmp
EOF

THIRD_JOB=`qsub -W depend=afterok${SECOND_MULTI_JOBS} ${PBS_DIR}/${Sample}_02-3.haplotypecaller.merge.${DATE}.pbs`
