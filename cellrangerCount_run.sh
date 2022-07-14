#!/bin/bash

### Name of cellranger output foler
DIR=${1}
ID=${2}
samples=${3}
transcriptome=${4}
expectedCells=${5}
cores=${6}
mem=${7}
Queue=${8}
singleNucleus=${9}


### cellRanger
cellRanger=/home/public/tools/cellranger-6.0.0/bin/cellranger

### Set working directories
FQDIR=${DIR}/00.fastq/${ID}
countDIR=${DIR}/01.cellrangerCount

PBS_DIR=${DIR}/pbs_jobs
PBS_LOGS=${DIR}/pbs_stdouts

mkdir -p ${countDIR}
mkdir -p ${PBS_DIR}
mkdir -p ${PBS_LOGS}

### DATE

DATE=$(date '+%y%m%d.%H-%M-%S')

######################################################################
########################## cellranger-count ##########################
######################################################################

cat > ${PBS_DIR}/${ID}_cellrangerCount.${DATE}.pbs << EOF
#PBS -N ${ID}_cellrangerCount.${DATE}
#PBS -j oe
#PBS -q ${Queue}
#PBS -o ${PBS_LOGS}/${ID}_cellrangerCount.${DATE}.pbs.stdout.txt
#PBS -l nodes=1:ppn=${cores}

cd ${countDIR}

${cellRanger} count --id=${ID}\
 --transcriptome=${transcriptome}\
 --fastqs=${FQDIR}\
 --sample=${samples}\
 --expect-cells=${expectedCells}\
 --localcores=${cores}\
 --localmem=${mem}\
 ${singleNucleus}

sleep 10
exit 0
EOF

FIRST_JOB=`qsub ${PBS_DIR}/${ID}_cellrangerCount.${DATE}.pbs`
