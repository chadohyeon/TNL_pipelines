#!/home/dcha/python3

## Written by Do Hyeon Cha, TNL@KAIST
## Wrapper for cellrangerCount_run.sh
## Usage
## wrapper_cellrangerCount.py -s human -c 10 -m 20 -n y
   # 10 cores will be utilized for one-sample-job
   # 20g would be the maximum memory for one-sample-job
   # Do not exceed total available threads


import os, sys, argparse
from itertools import product
import re

class cellRangerCount():
    def __init__(self, script, DIR, transcriptome, expectedCells, cores, mem, Queue, singleNucleus):
        self.script=script
        self.DIR=DIR
        self.transcriptome=transcriptome
        self.expectedCells=expectedCells
        self.cores=cores
        self.mem=mem
        self.singleNucleus=singleNucleus
        self.Queue=Queue
        
    def sample_define(self):
        sampleDic={ids:"" for ids in os.listdir(self.DIR+"/00.fastq")}
        for ids in sampleDic.keys():
            fqs=[]
            for j in os.listdir(self.DIR+"/00.fastq/"+ids):
                if j.endswith(".fq.gz") or j.endswith(".fastq.gz") or j.endswith(".fq") or j.endswith(".fastq"):
                  fqs.append(j.split("_")[0])  
            fqs=list(set(fqs))
            #fqs=list(set([j.split("_")[0] for j in os.listdir(self.DIR+"/00.fastq/"+ids) if j.endswith("fq.gz")]))
            sampleDic[ids]=",".join(fqs)
        return sampleDic

    def run(self):
        self.sampleDic=self.sample_define()
        for ids in self.sampleDic.keys():
            sampleNames=self.sampleDic[ids]
            self.singleRun(ids, sampleNames)

    def singleRun(self, ids, sampleNames):                               
            commandLi = ["bash", self.script, self.DIR, ids, sampleNames, self.transcriptome,self.expectedCells, self.cores, self.mem, self.Queue, self.singleNucleus]
            command = " ".join(commandLi)
            print(command)
            os.system(command)

def main():
    script="/home/public/scripts/cellrangerCount_run.sh"

    humanTranscriptome="/home/public/ref/00_Ref_GRCh38/refdata-gex-GRCh38-2020-A"
    mouseTranscriptome="/home/public/ref/00_Ref_mm10/refdata-gex-mm10-2020-A"
    
    parser=argparse.ArgumentParser()
    parser.add_argument('-s', '--species', required=True, help="human or mouse", type=str, metavar='')
    parser.add_argument('-n', '--nucleus', required=False, help="single-nucleus seqeuncing (y/n)", type=str, default="n", metavar='')
    parser.add_argument('-c', '--cores',required=False, type=str, help='Multiple cores for parallelization', default='1', metavar='')
    parser.add_argument('-m', '--mem', required=False, type=str, help='Max GB memory per sample', default='30', metavar='')
    parser.add_argument('-e', '--expectCells', required=False, type=str, help='expect-cells option for cellranger', default='10000', metavar='')
    parser.add_argument('-d', '--directory', required=False, type=str, help="Default workdir will be PWD", default=os.getcwd(), metavar='') ## Provide it when needed
    parser.add_argument('-q', '--queue', required=False, type=str, help="PBS queue", default='long', metavar='') 

    args=parser.parse_args()
    #parser.print_help()

    ### Species detection
    if args.species=="human":
        transcriptome=humanTranscriptome
    elif args.species=="mouse":
        transcriptome=mouseTranscriptome
    else:
        print("Please provide right species information")
        return None
    
    ### Single-cell or Single-nucleus RNA sequencing
    if args.nucleus=="n":
        singleNucleus=""
    elif args.nucleus=="y":
        singleNucleus="--include-introns"
    else:
        print("Please provide right scRNAseq/snRNAseq information")
        return None                                        
                                                
    
    cellRangerCountInstance=cellRangerCount(script, args.directory, transcriptome, args.expectCells, args.cores, args.mem, args.queue, singleNucleus)
    cellRangerCountInstance.run()
                                                
if __name__ == "__main__":
    main()
