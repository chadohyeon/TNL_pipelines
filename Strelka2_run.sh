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
strelkaPath=${tools}/strelka-2.9.10/bin
gatk=${tools}/gatk-4.2.0.0/gatk
bgzip=${tools}/htslib/bin/bgzip
bcftools=${tools}/bcftools/bin/bcftools
tabix=${tools}/htslib/bin/tabix

### Set working directories
BAMDIR=${DIR}/01.bam
SOMDIR=${DIR}/04.somatic_vcf
SOMpassDIR=${DIR}/05.PASS_somatic_vcf

PBS_DIR=${DIR}/pbs_jobs
PBS_LOGS=${DIR}/pbs_stdouts

STRELKATMP=${SOMDIR}/strelka2_tmp
DIR2REF=`python3 -c "print('/'.join('${REF}'.split('/')[:-2]))"`
GRCH=`python3 -c "print('${DIR2REF}'.split('/')[-1].split('_')[-1])"`

mkdir -p ${SOMDIR}
mkdir -p ${STRELKATMP}

mkdir -p ${SOMpassDIR}

mkdir -p ${PBS_DIR}
mkdir -p ${PBS_LOGS}

### DATE

DATE=$(date '+%y%m%d.%H-%M-%S')

######################################################################
############################# Strelka2 ###############################
######################################################################

### 01. Scattering Jobs
allRaw_snvs=""
allRaw_indels=""
chroms=$(echo ${commaChroms} | sed "s/,/ /g")

for chrom in ${chroms}; do
    allRaw_snvs+="${STRELKATMP}/${somatic_sample}.pair_${matched_sample}_chunks/${chrom}/results/variants/somatic.snvs.vcf.gz "
    allRaw_indels+="${STRELKATMP}/${somatic_sample}.pair_${matched_sample}_chunks/${chrom}/results/variants/somatic.indels.vcf.gz "
done

### 01. Scattered variant calling
FIRST_MULTI_JOBS=""

for chrom in ${chroms}; do
cat > ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-2.strelka2.chunk_${chrom}.${DATE}.pbs << EOF
#PBS -N ${somatic_sample}.pair_${matched_sample}_03-2.strelka2.chunk_${chrom}.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${somatic_sample}.pair_${matched_sample}_03-2.strelka2.chunk_${chrom}.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1


python2 ${strelkaPath}/configureStrelkaSomaticWorkflow.py\
 --referenceFasta ${REF}\
 --normalBam ${BAMDIR}/${matched_sample}${bamSuffix}\
 --tumorBam ${BAMDIR}/${somatic_sample}${bamSuffix}\
 --callRegions ${DIR2REF}/annot/perChromBED/perChromBED.${chrom}.bed.gz\
 --runDir ${STRELKATMP}/${somatic_sample}.pair_${matched_sample}_chunks/${chrom}

python2 ${STRELKATMP}/${somatic_sample}.pair_${matched_sample}_chunks/${chrom}/runWorkflow.py -m local -j 1

sleep 10
exit 0
EOF

FIRST_MULTI_JOBS+=":"`qsub ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-2.strelka2.chunk_${chrom}.${DATE}.pbs `
done

### 03. Remainders
cat > ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-3.strelka2.merge.${DATE}.pbs << EOF
#PBS -N ${somatic_sample}.pair_${matched_sample}_03-3.strelka2.merge.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${somatic_sample}.pair_${matched_sample}_03-3.strelka2.merge.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

${bcftools} concat -a ${allRaw_snvs}\
 -Oz -o ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.SNV.strelka2.vcf.gz
${tabix} -p vcf ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.SNV.strelka2.vcf.gz

${bcftools} concat -a ${allRaw_indels}\
 -Oz -o ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.INDEL.strelka2.vcf.gz
${tabix} -p vcf ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.INDEL.strelka2.vcf.gz

${bcftools} concat -a ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.SNV.strelka2.vcf.gz\
 ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.INDEL.strelka2.vcf.gz\
 -Oz -o ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.MERGED.strelka2.vcf.gz
${tabix} -p vcf ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.MERGED.strelka2.vcf.gz

mkdir -p ${SOMpassDIR}/strelka2_tmp
${bcftools} view -f PASS ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.SNV.strelka2.vcf.gz\
 -Oz -o ${SOMpassDIR}/strelka2_tmp/${somatic_sample}.pair_${matched_sample}.PASS.SNV.strelka2.vcf.gz
${tabix} -p vcf ${SOMpassDIR}/strelka2_tmp/${somatic_sample}.pair_${matched_sample}.PASS.SNV.strelka2.vcf.gz

${bcftools} view -f PASS ${SOMDIR}/${somatic_sample}.pair_${matched_sample}.INDEL.strelka2.vcf.gz\
 -Oz -o ${SOMpassDIR}/strelka2_tmp/${somatic_sample}.pair_${matched_sample}.PASS.INDEL.strelka2.vcf.gz
${tabix} -p vcf ${SOMpassDIR}/strelka2_tmp/${somatic_sample}.pair_${matched_sample}.PASS.INDEL.strelka2.vcf.gz

${bcftools} concat -a ${SOMpassDIR}/strelka2_tmp/${somatic_sample}.pair_${matched_sample}.PASS.SNV.strelka2.vcf.gz\
 ${SOMpassDIR}/strelka2_tmp/${somatic_sample}.pair_${matched_sample}.PASS.INDEL.strelka2.vcf.gz\
 -Oz -o ${SOMpassDIR}/${somatic_sample}.pair_${matched_sample}.PASS.MERGED.strelka2.vcf.gz
 
${tabix} ${SOMpassDIR}/${somatic_sample}.pair_${matched_sample}.PASS.MERGED.strelka2.vcf.gz
 
EOF
SECOND_JOB=`qsub -W depend=afterok${FIRST_MULTI_JOBS} ${PBS_DIR}/${somatic_sample}.pair_${matched_sample}_03-3.strelka2.merge.${DATE}.pbs`
