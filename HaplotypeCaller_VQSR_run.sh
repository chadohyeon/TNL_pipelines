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

EOF

THIRD_JOB=`qsub -W depend=afterok${SECOND_MULTI_JOBS} ${PBS_DIR}/${Sample}_02-3.haplotypecaller.merge.${DATE}.pbs`

### 04. VQSR-SNV
cat > ${PBS_DIR}/${Sample}_02-4.haplotypecaller.VQSR_SNV.${DATE}.pbs << EOF
#PBS -N ${Sample}_02-4.haplotypecaller.VQSR_SNV.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${Sample}_02-4.haplotypecaller.VQSR_SNV.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" SelectVariants\
 -R ${REF}\
 -V ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.vcf.gz\
 --select-type-to-include SNP\
 -O ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.SNV_select.vcf.gz

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" VariantRecalibrator\
 -V ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.SNV_select.vcf.gz\
 --trust-all-polymorphic\
 -mode SNP\
 -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0\
 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR\
 --max-gaussians 6\
 --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${DIR2REF}/annot/hapmap_3.3.${GRCH}.vcf.gz\
 --resource:omni,known=false,training=true,truth=true,prior=12.0 ${DIR2REF}/annot/1000G_omni2.5.${GRCH}.vcf.gz\
 --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${DIR2REF}/annot/1000G_phase1.snps.high_confidence.${GRCH}.vcf.gz\
 --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 ${DIR2REF}/annot/dbSNP154.${GRCH}.vcf.gz\
 --tranches-file ${VCFDIR}/haplotypecaller_tmp/${Sample}.recalibrate_SNP.tranches\
 -O ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.SNV_select.vqsr.recal

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" ApplyVQSR\
 -V ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.SNV_select.vcf.gz\
 --recal-file  ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.SNV_select.vqsr.recal\
 --tranches-file ${VCFDIR}/haplotypecaller_tmp/${Sample}.recalibrate_SNP.tranches\
 --truth-sensitivity-filter-level 99.5\
 --create-output-variant-index true\
 -mode SNP\
 -O ${VCFDIR}/${Sample}.haplotypecaller.SNV_select.VQSR_applied.vcf.gz\

EOF

FOURTH_JOB=`qsub -W depend=afterok:${THIRD_JOB} ${PBS_DIR}/${Sample}_02-4.haplotypecaller.VQSR_SNV.${DATE}.pbs`

### 05. VQSR-INDEL
cat > ${PBS_DIR}/${Sample}_02-5.haplotypecaller.VQSR_INDEL.${DATE}.pbs << EOF
#PBS -N ${Sample}_02-5.haplotypecaller.VQSR_INDEL.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${Sample}_02-5.haplotypecaller.VQSR_INDEL.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" SelectVariants\
 -R ${REF}\
 -V ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.vcf.gz\
 --select-type-to-include INDEL\
 -O ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.INDEL_select.vcf.gz

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" VariantRecalibrator\
 -V ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.INDEL_select.vcf.gz\
 --trust-all-polymorphic\
 -mode INDEL\
 -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0\
 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR\
 --max-gaussians 4\
 --resource:mills,known=false,training=true,truth=true,prior=15.0 ${DIR2REF}/annot/Mills_and_1000G_gold_standard.indels.${GRCH}.vcf.gz\
 --resource:axiomPoly,known=false,training=true,truth=false,prior=12.0 ${DIR2REF}/annot/Axiom_Exome_Plus.genotypes.all_populations.poly.${GRCH}.vcf.gz\
 --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 ${DIR2REF}/annot/dbSNP154.${GRCH}.vcf.gz\
 --tranches-file ${VCFDIR}/haplotypecaller_tmp/${Sample}.recalibrate_INDEL.tranches\
 -O ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.INDEL_select.vqsr.recal

${gatk} --java-options "-Xmx${MEM}g -Djava.io.tmpdir=${JAVA_SPACE}" ApplyVQSR\
 -V ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.INDEL_select.vcf.gz\
 --recal-file  ${VCFDIR}/haplotypecaller_tmp/${Sample}.haplotypecaller.INDEL_select.vqsr.recal\
 --tranches-file ${VCFDIR}/haplotypecaller_tmp/${Sample}.recalibrate_INDEL.tranches\
 --truth-sensitivity-filter-level 99.0\
 --create-output-variant-index true\
 -mode INDEL\
 -O ${VCFDIR}/${Sample}.haplotypecaller.INDEL_select.VQSR_applied.vcf.gz\

EOF

FIFTHJOB=`qsub -W depend=afterok:${THIRD_JOB} ${PBS_DIR}/${Sample}_02-5.haplotypecaller.VQSR_INDEL.${DATE}.pbs`

### 06. Merge VQSR SNV+INDEL
cat > ${PBS_DIR}/${Sample}_02-6.haplotypecaller.mergeVQSR.${DATE}.pbs << EOF
#PBS -N ${Sample}_02-6.haplotypecaller.mergeVQSR.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${Sample}_02-6.haplotypecaller.mergeVQSR.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=1

### Grep PASS-only rows and index again
${bcftools} view -f PASS\
 ${VCFDIR}/${Sample}.haplotypecaller.SNV_select.VQSR_applied.vcf.gz -Oz\
 -o ${VCFpassDIR}/${Sample}.VQSR_PASS.SNV.haplotypecaller.vcf.gz
${tabix} -p vcf ${VCFpassDIR}/${Sample}.VQSR_PASS.SNV.haplotypecaller.vcf.gz

${bcftools} view -f PASS\
 ${VCFDIR}/${Sample}.haplotypecaller.INDEL_select.VQSR_applied.vcf.gz -Oz\
 -o ${VCFpassDIR}/${Sample}.VQSR_PASS.INDEL.haplotypecaller.vcf.gz
${tabix} -p vcf ${VCFpassDIR}/${Sample}.VQSR_PASS.INDEL.haplotypecaller.vcf.gz

### Merge SNV and INDEL passed calls
${bcftools} concat -a ${VCFDIR}/${Sample}.haplotypecaller.SNV_select.VQSR_applied.vcf.gz\
 ${VCFDIR}/${Sample}.haplotypecaller.INDEL_select.VQSR_applied.vcf.gz\
 -Oz -o ${VCFDIR}/${Sample}.haplotypecaller.MERGED.VQSR_applied.vcf.gz
${tabix} -p vcf ${VCFDIR}/${Sample}.haplotypecaller.MERGED.VQSR_applied.vcf.gz


${bcftools} concat -a ${VCFpassDIR}/${Sample}.VQSR_PASS.SNV.haplotypecaller.vcf.gz\
 ${VCFpassDIR}/${Sample}.VQSR_PASS.INDEL.haplotypecaller.vcf.gz\
 -Oz -o ${VCFpassDIR}/${Sample}.VQSR_PASS.MERGED.haplotypecaller.vcf.gz

${tabix} -p vcf ${VCFpassDIR}/${Sample}.VQSR_PASS.MERGED.haplotypecaller.vcf.gz

mkdir -p ${VCFpassDIR}/haplotypecaller_tmp
mv ${VCFpassDIR}/${Sample}.VQSR_PASS.SNV.haplotypecaller.vcf.gz* ${VCFpassDIR}/haplotypecaller_tmp
mv ${VCFpassDIR}/${Sample}.VQSR_PASS.INDEL.haplotypecaller.vcf.gz* ${VCFpassDIR}/haplotypecaller_tmp

EOF

SIXTH_JOB=`qsub -W depend=afterok:${FOURTH_JOB}:${FIFTH_JOB} ${PBS_DIR}/${Sample}_02-6.haplotypecaller.mergeVQSR.${DATE}.pbs`
