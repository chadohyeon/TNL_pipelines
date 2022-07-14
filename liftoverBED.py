#!/usr/bin/python3

import os, sys, argparse

def replaceAll(string, pat1, pat2):
    while pat1 in string:
        string=string.replace(pat1, pat2)
    return(string)

def mergeMultipleBed(bedLi, wfn):
    for b in bedLi:
        os.system("mv {0} {0}.unsorted".format(b))
        os.system("sort -k1,1 -k2,2n {0}.unsorted > {0}".format(b))
        os.system("rm {0}.unsorted".format(b))
    os.system("bedops --everything {0} | sort -k1,1 -k2,2n > union.bed".format(" ".join(bedLi)))
    os.system("bedops --partition {0} | sort -k1,1 -k2,2n > partition.bed".format(" ".join(bedLi)))
    os.system("bedmap --echo --echo-map-id --delim '\t' partition.bed union.bed | sort -k1,1V -k2,2n > {0}".format(wfn))
    
def bedReduction(dupBedFn):
    tmpSortedBedFn=dupBedFn.replace(".bed", ".tmpSorted.bed")
    tmpSortedSpaceFilledBedFn=dupBedFn.replace(".bed", ".tmpSorted.spaceFilled.bed")
    finalBedSpaceFilledFn=dupBedFn.replace(".bed", ".final.spaceFilled.bed")
    finalBedFn=dupBedFn.replace(".bed", ".trimmed.bed")
    os.system("sort -k1,1 -k2,2n {0} > {1}".format(dupBedFn,tmpSortedBedFn))
    with open(tmpSortedBedFn) as rf:
        with open(tmpSortedSpaceFilledBedFn, "w") as wf:
            for k in rf:
                k=replaceAll(k, " ", "%%")
                wf.write(k)
    annotLi=[]
    with open(tmpSortedSpaceFilledBedFn) as bed:
        for line in bed:
            annotLi=list(set(annotLi+[line.rstrip().split("\t")[-1]]))
    for itr in range(len(annotLi)):
        ann=annotLi[itr]
        os.system('grep -P "\t{0}$" {1} > {2}_chunk.{1}'.format(ann, tmpSortedSpaceFilledBedFn, itr))
    tmpBedLiJoined=["{0}_chunk.{1}".format(itr, tmpSortedSpaceFilledBedFn) for itr in range(len(annotLi))]
    mergeMultipleBed(tmpBedLiJoined,finalBedSpaceFilledFn)
    with open(finalBedSpaceFilledFn) as rf:
        with open(finalBedFn, "w") as wf:
            for k in rf:
                k=replaceAll(k, "%%", " ")
                wf.write(k)
    
    for fn in [tmpSortedBedFn, tmpSortedSpaceFilledBedFn, finalBedSpaceFilledFn]+tmpBedLiJoined:
        os.system("rm -f {0}".format(fn))
    return(finalBedFn)

def liftingOver(fn, reduce):
    direc="/home/public/tools/liftover/"
    liftover=direc+"liftOver"
    chain19_38=direc+"hg19ToHg38.over.chain"
    chain38_19=direc+"hg38ToHg19.over.chain"
    contigs37=[str(k) for k in range(1,23)] + ["X","Y","MT"]
    contigs38=["chr"+str(k) for k in range(1,23)] + ["chrX", "chrY", "chrM"]
    dic37to38={contigs37[itr]:contigs38[itr] for itr in range(len(contigs37))}

    if "hg19" in fn:
        unsortedFn=fn.replace("hg19","hg38")+".unsorted"
        chain=chain19_38
    elif "GRCh37" in fn:
        unsortedFn=fn.replace("GRCh37", "GRCh38")+".unsorted"
        chain=chain19_38
    elif "hg38" in fn:
        unsortedFn=fn.replace("hg38","hg19")+".unsorted"
        chain=chain38_19
    elif "GRCh38" in fn:
        unsortedFn=fn.replace("GRCh38", "GRCh37")+".unsorted"
        chain=chain38_19
    else:
        print("please refer to files with 'hg19','GRCh37','hg38','GRCh38' in name")
        return(None)
    wfn=unsortedFn.replace(".unsorted","")
    chromAdd=False
    with open(fn) as bed:
        first = bed.readline().rstrip().split("\t")
        second = bed.readline().rstrip().split("\t")
        if second[0] in contigs37:
            chromAdd=True
    intermFn=fn
    
    firstlineHeader=False

    if chromAdd:
        intermFn=fn.replace(".bed", ".chromAdd.bed")
        wf=open(intermFn, "w")
        with open(fn) as bed:
            first = bed.readline().rstrip().split("\t")
            if not first[0] in contigs37:
                firstlineHeader=True
        with open(fn) as bed:
            first=bed.readline().rstrip().split("\t")
            if not firstlineHeader:
                first[0]=dic37to38[first[0]]
                wf.write("\t".join(first)+"\n")
            for k in bed:
                line=k.rstrip().split("\t")
                line[0]=dic37to38[line[0]]
                wf.write("\t".join(line)+"\n")


    command1="mkdir {0}".format(os.getcwd()+"/remainders")
    command2="{0} {1} {2} {3} {4}".format(liftover, intermFn, chain, unsortedFn, os.getcwd()+"/remainders/"+fn+".remainder")
    command3="sort -k1,1V -k2,2n {0} > {1}".format(unsortedFn, wfn)
    command4="rm -f {0}".format(unsortedFn)
    command5="rm -rf {0}".format(os.getcwd()+"/remainders/")
    
    os.system(command1)
    os.system(command2)
    os.system(command3)
    os.system(command4)
    os.system(command5)
    
    if fn!=intermFn:
        command6="rm "+intermFn
        os.system(command6)
    
    if reduce=="y":
        finalFn=bedReduction(wfn)
        os.system("rm -f {0}".format(wfn))
        os.system("mv {0} {1}".format(finalFn, wfn))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help="Input BED", type=str, metavar='')
    parser.add_argument('-r', '--reduce',required=False, type=str, help='Reducing BED when overlapped regions exist(y/n)', default='n', metavar='')

    args=parser.parse_args()
    parser.print_help()
    
    reduce=args.reduce.lower()
    if not reduce in ["y", "n"]:
        print("Please assign right -r/--reduce option (y/n)")
        return(None)
        
    liftingOver(args.input, reduce)

if __name__ == "__main__":
    main()
