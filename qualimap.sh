#!/bin/bash

bamdir=${PWD}/01.bam
outdir=${PWD}/qualimap_results

mkdir -p ${outdir}

bams=$(ls ${bamdir} | grep .marked.recal.bam)

#unset DISPLAY

for bam in ${bams}; do
/home/public/tools/qualimap_v2.2.1/qualimap bamqc \
	-bam ${bamdir}/${bam} \
	--java-mem-size=10G -c \
	-outdir ${outdir} \
	-outfile ${bam}.pdf \
	-outformat PDF \
	-nt 6
done
