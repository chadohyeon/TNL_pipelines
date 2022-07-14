#!/home/dcha/python3

## Written by Do Hyeon Cha, TNL@KAIST
## Wrapper for TrimGalore_single.sh, which is master shell script of a PBS job
## Usage
## wrapper_trimgalore.py -c 8
   # 8 cores will be utilized for one-sample-job
   # Jobs for multiple samples in FASTQ directory would be submitted.
    
import os, sys, math, time, argparse
#from multiprocessing import Process
from itertools import product
import re

def sample_define(direc):
    sample_dict={}
    _rawfiles=[i for i in os.listdir(direc)]
    rawfiles=[]
    for f in _rawfiles:
        for suf in [".fq", ".fastq", ".fq.gz", ".fastq.gz"]:
            if f.endswith(suf):
                rawfiles.append(f)
    for i in rawfiles:
        rule=fq_nomen(i)
        x=i.split(rule)[0]
        if x[-1] in [".", "-", "_"]:
            rule=x[-1]+rule
            y=x[:-1]
        else: y=x
        if y not in sample_dict.keys():
            sample_dict[y]=[rule]
        else: sample_dict[y].append(rule)
    for j in sample_dict.keys():
        sample_dict[j]=sorted(sample_dict[j])
    sample_dict = {k:sample_dict[k] for k in sample_dict.keys() if len(sample_dict[k])==2}
    #sample_dict = {k:sample_dict[k] for k in sample_dict.keys() if not k.endswith("_val")}
    return sample_dict

def fq_nomen(fn):
    pe_rules=['R[1-2]', '[1-2]']
    fq_rules=[".fq", ".fastq","_001.fq", "_001.fastq"]
    rules=[k[0]+k[1] for k in list(product(*[pe_rules,fq_rules]))]
    p=list(map(lambda i:re.compile(i), rules))
    for i in p:
        if len(i.findall(fn))!=0:
            tmpName=i.findall(fn)[0]
            if fn.endswith(".gz"):
                tmpName+=".gz"
            return(tmpName)


def main():
    script="/home/public/scripts/TrimGalore_run.sh"
    refPath="/home/public/ref/"

    parser=argparse.ArgumentParser()
    parser.add_argument('-c', '--core',required=False, type=str, help='cores', default='1', metavar='')
    parser.add_argument('-s', '--sequencing', required=False, type=str, help='Sequencing platform', default='illumina', metavar='')
    parser.add_argument('-d', '--directory', required=False, type=str, help="Default workdir will be PWD", default=os.getcwd(), metavar='') ## Provide it when needed
    parser.add_argument('-q', '--queue', required=False, type=str, help="PBS queue", default='long', metavar='') 
    args=parser.parse_args()
    #parser.print_help()

    try:
        os.listdir("{0}/00.fastq".format(args.directory))
    except:
        print("fastq (gzipped or not) should be in 00.fastq")
        return None
  
    sample_dict=sample_define("{0}/00.fastq".format(args.directory))
    if sample_dict=={}:
        print("fastq (gzipped or not) not found")
        return None
    
    for k in sample_dict.keys():
        print("Sample {0} are being trim-galored".format(k))
        fq1 = k + sample_dict[k][0]
        fq2 = k + sample_dict[k][1]
        os.system(" ".join(["bash", script, args.directory,args.sequencing ,k, args.core, args.queue, fq1, fq2]))
        
if __name__ == "__main__":
    main()
