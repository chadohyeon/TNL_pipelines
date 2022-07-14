#!/bin/bash

### Made by Do Hyeon Cha TNL@KAIST

### Args from wrapper_trimGalore.py

DIR=${1}
PLATFORM=${2}
sample=${3}
Cores=${4}
Queue=${5}
FQ_F=${6}
FQ_R=${7}

### Tools

tools=/home/public/tools
trim_galore=${tools}/TrimGalore-0.6.6/trim_galore



### Set working directories
FQDIR=${DIR}/00.fastq
TMPDIR=${FQDIR}/trimTmp
PBS_DIR=${DIR}/pbs_jobs
PBS_LOGS=${DIR}/pbs_stdouts

mkdir -p ${TMPDIR}
mkdir -p ${PBS_DIR}
mkdir -p ${PBS_LOGS}

### DATE
DATE=$(date '+%y%m%d.%H-%M-%S')

###########################################################################
###########################   Run TrimGalore!   ###########################
###########################################################################

cat > ${PBS_DIR}/${sample}_00.trimGalore.${DATE}.pbs << EOF
#PBS -N ${sample}_00.trimGalore.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${sample}_00.trimGalore.${DATE}.stdout.txt
#PBS -l nodes=1:ppn=${Cores}

echo "QUERY:trimGalore:Sample=${sample}"
${trim_galore} -o ${FQDIR} --${PLATFORM} --paired --basename ${sample} --cores ${Cores} ${FQDIR}/${FQ_F} ${FQDIR}/${FQ_R}

mv ${FQDIR}/${FQ_F} ${TMPDIR}/${FQ_F}
mv ${FQDIR}/${FQ_R} ${TMPDIR}/${FQ_R}
mv ${FQDIR}/${sample}_val_1.fq.gz ${FQDIR}/${sample}_1.fq.gz
mv ${FQDIR}/${sample}_val_2.fq.gz ${FQDIR}/${sample}_2.fq.gz 
mv ${FQDIR}/*trimming_report.txt ${TMPDIR}
sleep 10
exit 0
EOF

FIRSTJOB=`qsub ${PBS_DIR}/${sample}_00.trimGalore.${DATE}.pbs`
