#!/home/dcha/python3

## Written by Do Hyeon Cha, TNL@KAIST
## Wrapper for Aligner_single.sh, which is master shell script of every PBS jobs
## Usage
## wrapper_align.py -r GRCh38 -c 18 -m 30
   # 18 cores will be utilized for one-sample-job
   # 30g = 30g would be the maximum memory for one-sample-job
   # Do not exceed total available threads


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
    return sample_dict

def fq_nomen(fn):
    pe_rules=['R[1-2]', '[1-2]']
    fq_rules=[".fq", ".fastq", "_001.fq", "_001.fastq"]
    rules=[k[0]+k[1] for k in list(product(*[pe_rules,fq_rules]))]
    p=list(map(lambda i:re.compile(i), rules))
    for i in p:
        if len(i.findall(fn))!=0:
            tmpName = i.findall(fn)[0]
            if fn.endswith(".gz"):
                tmpName+=".gz"
            return(tmpName)


def multi(script, cores, mem, direc, fa, annot, sample_dict, queue):
    for sampleName in sample_dict.keys():
        single(script, sampleName, mem, cores, direc, fa, annot, sample_dict, queue)   
    #sample_heads=list(sample_dict.keys())
    #size=float(size)
    #num_iter=math.ceil(float(len(sample_heads))/size)
    #for i in range(num_iter):
    #    start=int(size)*i
    #    end=int(size)*(i+1)
    #    if end<=len(sample_heads): batch=sample_heads[start:end]
    #    else: batch=sample_heads[start:]

    #    procs=[Process(target=single, args=(script, s,mem,cores,direc,fa,annot,sample_dict)) for s in sample_heads]
    #    for p in procs: p.start()
    #    for p in procs: p.join()

def single(script, sampleName, mem, cores, direc, fa, annot, sample_dict, queue):
    try:
        fqf=sampleName+sample_dict[sampleName][0]
        fqr=sampleName+sample_dict[sampleName][1]
        commandLi = ["bash", script, mem, direc, fa, annot['dbSNP'], annot['Mills_1000G'], sampleName, fqf, fqr, cores, queue]
        command = " ".join(commandLi)
        print(command)
        os.system(command)
    except:
        print("An error has occured")

def annotDic(dbSNP, Mills_1000G):
    annotHash={}
    annotHash["dbSNP"]=dbSNP
    annotHash["Mills_1000G"]=Mills_1000G

    return annotHash

def main():
    refPath="/home/public/ref/"

    fasta38="Homo_sapiens_assembly38.fa"
    fasta37="Homo.Sapiens.GRCh37.75.primary_assembly.full.fa"

    annot38=annotDic(dbSNP="dbSNP154.GRCh38.vcf.gz", Mills_1000G="Mills_and_1000G_gold_standard.indels.GRCh38.vcf.gz")
    annot37=annotDic(dbSNP="dbSNP154.GRCh37.vcf.gz", Mills_1000G="Mills_and_1000G_gold_standard.indels.GRCh37.vcf.gz")
    
    parser=argparse.ArgumentParser()
    parser.add_argument('-r', '--refFa', required=True, help="GRCh38 or GRCh37", type=str, metavar='')
    parser.add_argument('-c', '--cores',required=False, type=str, help='Multiple cores for parallelization', default='1', metavar='')
    parser.add_argument('-m', '--mem', required=False, type=str, help='Max GB memory per sample', default='30', metavar='')
    parser.add_argument('-d', '--directory', required=False, type=str, help="Default workdir will be PWD", default=os.getcwd(), metavar='') ## Provide it when needed
    parser.add_argument('-q', '--queue', required=False, type=str, help="PBS queue", default='long', metavar='') 
    parser.add_argument('-s', '--spark', required=False, type=str, help="Using Spark", default='y', metavar='') 

    args=parser.parse_args()
    #parser.print_help()

    print("START: "+time.strftime('%c', time.localtime(time.time())))
    
    fastaPath=refPath+"00_Ref_{0}/{0}_fasta/".format(args.refFa)
    annotPath=refPath+"00_Ref_{0}/annot/".format(args.refFa)
    
    if args.spark == "y":
        script="/home/public/scripts/AlignerSpark_run.sh"
    elif args.spark == "n":
        script="/home/public/scripts/Aligner_run.sh"
    else:
        print("Please provide right --spark information")
        return(None)
    
    ## Reference fasta and common annotation sources
    if args.refFa == "GRCh38":
        fa=fastaPath+fasta38
        annot={k:annotPath+annot38[k] for k in annot38.keys()}
    elif args.refFa == "GRCh37":
        fa=fastaPath+fasta37
        annot={k:annotPath+annot37[k] for k in annot37.keys()}
    else:
        print("\nReference genome is not in GRCh38 or GRCh37\n")
        return None

    sample_dict=sample_define("{0}/00.fastq".format(args.directory))
    multi(script, args.cores, args.mem, args.directory, fa, annot, sample_dict, args.queue)

if __name__ == "__main__":
    main()
