#!/bin/bash

path=${1}

cd ${path}

phy="${path}cellTypes_dna.phylip"

mkdir -p ${path}raxml
mkdir -p ${path}raxml/iter100

#cd raxml/iter100

raxmlHPC -s ${phy} -m GTRGAMMA -# 100 -p 12345 -n outgroup -o Blood -w ${path}raxml/iter100
#cp RAxML_bestTree.outgroup ..

#cd ..
raxmlHPC -s ${phy} -m GTRGAMMA -p 12345 -n anc -f A -t ${path}raxml/iter100/RAxML_bestTree.outgroup -w ${path}raxml




