#!/home/dcha/python3

### Written by Do Hyeon Cha, TNL@KAIST
### Wrapper for Germline, Somatic variant callers
###
### USAGE
### wrapper_caller.py -v HaplotypeCaller -r GRCh38 -c 10 -g ABC_blood,ABC_PFC,ABC_Cbll
### wrapper_caller.py -v HaplotypeCaller -r GRCh38 -c 10 -d auto -m 30 -g ABC_blood,ABC_PFC
### wrapper_caller.py -v Mutect2 -r GRCh37 -c 10 -d mito -g ABC_blood -s ABC_PFC,ABC_Cbll
### wrapper_caller.py -v MuTect -r GRCh38 -c 15 -g ABC_blood -s ABC_PFC,ABC_Cbll
### wrapper_caller.py -v MuTect-frac0 -r GRCh38 -c 12 -g ABC_blood -s ABC_PFC,ABC_Cbll
### wrapper_caller.py -v Strelka2 -r GRCh38 -c 13 -g ABC_blood -s ABC_PFC,ABC_Cbll
### wrapper_caller.py -v VarScan -r GRCh38 -c 10 -g ABC_blood -s ABC_PFC,ABC_Cbll

import os, argparse, glob
from multiprocessing import Process

class caller():
    def __init__(self, singleRunSh, bamSuffix, fa, annot, chrExist, sample, chrom, mem, directory, cores, queue):
        #self.per=per
        self.singleRunSh=singleRunSh

        self.bamSuffix=bamSuffix
        self.fa=fa

        self.PON=annot["PON"]

        self.chrExist=chrExist
        self.sample=sample
        self.chrom=chrom
        self.mem=mem
        self.directory=directory
        self.cores=cores
        self.queue=queue

    def singleRun(self, commaChroms):
        comLi = ["bash", self.singleRunSh, self.mem, self.directory, self.fa, self.PON, self.bamSuffix, commaChroms, self.cores, self.queue]
        comLi += self.sample
        command=" ".join(comLi)
        print(command)
        os.system(command)

    def whole(self):
        chroms=[str(i) for i in range(1,23)]
        if self.chrom=="sex":
            chroms+=['X', 'Y']
        elif self.chrom=="auto":
            pass
        elif self.chrom=="mito":
            chroms+=['X','Y','MT']
        else:	
            chroms+=['X','Y']
        
        if self.chrExist=='y':
            chroms2=['chrM' if "MT" in i else 'chr'+i for i in chroms]
        elif self.chrExist=='n':
            chroms2=[i for i in chroms] # hard copy
	
        print("\nSample: {0}".format(self.sample))
        print("Assigned contigs: {0}\n".format(",".join(chroms2)))

        commaChroms=",".join(chroms2)
        self.singleRun(commaChroms)

def checkSamples(sampleLi, bamSuffix):
    endsWithSuffix=list(set([g.endswith(bamSuffix) for g in sampleLi]))
    if endsWithSuffix==[True]:
        if  list(set([True if glob.glob("01.bam/"+g) else False for g in sampleLi])) == [True]:
            sampleLi=[g.replace(bamSuffix, "") for g in sampleLi]
        elif list(set([True if glob.glob(g) else False for g in sampleLi])) == [True]:
            sampleLi=[g.replace(bamSuffix, "").split("01.bam/")[1] for g in sampleLi]
        else:
            return(None)
    elif endsWithSuffix==[False]:
        if list(set([True if glob.glob("01.bam/"+g+bamSuffix) else False for g in sampleLi])) == [True]:
            sampleLi=[g for g in sampleLi]
        elif list(set([True if glob.glob(g+bamSuffix) else False for g in sampleLi])) == [True]:
            sampleLi=[g.split("01.bam/")[1] for g in sampleLi]
        else:
            return(None)
    else: 
        return(None)
    return(sampleLi)

def main():
    somatic_callers=["Mutect2", "MuTect", "MuTect-frac0", "Strelka2", "VarScan"]
    germline_callers=["HaplotypeCaller", "HaplotypeCaller_hardfilter"]
    scriptPath="/home/public/scripts/"
    refPath="/home/public/ref/"
    fasta38="Homo_sapiens_assembly38.fa"
    fasta37="Homo.Sapiens.GRCh37.75.primary_assembly.full.fa"
    fasta37_chr="Homo.Sapiens.GRCh37.75.primary_assembly.full.chr.fa"
    
    ### Define annotation dictionaries
    annot38={"PON":"1000g_pon.GRCh38.WGS.vcf.gz"}
    annot37={"PON":"1000g_pon.GRCh37.WGS.vcf.gz"}
    
    parser=argparse.ArgumentParser()
    parser.add_argument('-v', '--variantCaller', type=str, required=True, help="(Germline) HaplotypeCaller HaplotypeCaller_VQSR (Somatic) Mutect2 MuTect MuTect-frac0 Strelka2 VarScan", metavar='')
    parser.add_argument('-r', '--refFa', type=str, required=True, help='GRCh38 or GRCh37', metavar='') ## GRCh37, GRCh38
    parser.add_argument('-c', '--cores', required=False, type=str, help='multiple cores for parallelization', default='1', metavar='') ## auto, sex (default), mito
    parser.add_argument('-d', '--chrom', required=False, type=str, help='chromosome types', default='sex', metavar='') ## auto, sex (default), mito
    parser.add_argument('-e','--chrExist', required=False, type=str, help='[chr] in contig or not', default='.', metavar='') ## y: chr1,...chrM / n: 1, ..., MT
    parser.add_argument('-m','--mem', required=False, type=str, help='Max memory per sample', default='40', metavar='') ## Per-chromosome calling 
    parser.add_argument('-g', '--germline', required=False, type=str, help='Germline samples, by default, 01.bam/[sample].marked.bqsr.recal.bam', default=".", metavar='') ## NA12375_blood,NA12333_blood, ...
    parser.add_argument('-s', '--somatic', required=False, type=str, help='Should match the order with germline samples', default=".", metavar='') ## NA12375_GBM,NA_12333_GBM, ...
    parser.add_argument('-p', '--pon', required=False, type=str, help="PON for Mutect2", default=".", metavar='')   ## If exist, please provide context-specific PON.vcf(.gz) with absolute path
    parser.add_argument('-w', '--directory', required=False, type=str, help="Default workdir will be PWD", default=os.getcwd(), metavar='') ## Provide it when needed
    parser.add_argument('-b', '--bamSuffix', required=False, type=str, help="BAM file suffix", default='.marked.recal.bam', metavar='') ## Provide it when needed
    parser.add_argument('-q', '--queue', required=False, type=str, help="PBS queue", default='long', metavar='')
    args=parser.parse_args()
    #parser.print_help()

    ## Germline samples (automatically or assigned externally)
    if args.germline==".":
        germline_samples=[i.split(args.bamSuffix)[0] for i in os.listdir(args.directory+"/01.bam") if args.bamSuffix in i if not "bai" in i]
    else:
        germline_samples=checkSamples([k.strip() for k in args.germline.split(",")], args.bamSuffix)
        if not germline_samples:
            print("Germline samples were not properly assigned")
            return(None)

        
    ## Somatic samples (should assigned externally)
    if args.somatic!=".":
        somatic_samples=checkSamples([k.strip() for k in args.somatic.split(",")], args.bamSuffix)
        if not somatic_samples:
            print("Somatic samples were not properly assigned")
            return(None)

    ## Variant callers, germline (single) or somatic (matched) and final tuplized samples
    samples=None
    if args.variantCaller in germline_callers:
        samples=[[i] for i in germline_samples]
    elif args.variantCaller in somatic_callers:
        if args.somatic==".":
            print("\nPlease provide germline and somatic samples\n")
            return None
        elif len(germline_samples)==1:
            germline_samples=[germline_samples[0] for k in somatic_samples]
            samples=[[germline_samples[itr], somatic_samples[itr]] for itr in range(len(germline_samples))]
        elif len(germline_samples)!=len(somatic_samples):
            print("\nPlease match germline and somatic samples\n")
            return None
        else:
            samples=[[germline_samples[itr], somatic_samples[itr]] for itr in range(len(germline_samples))]
    else:
        print("\nGiven variant caller is not available, please check the spell\n")
        return(None)
    
    print("{0} has been selected".format(args.variantCaller))
    if args.variantCaller == "HaplotypeCaller":
        print("VQSR will be implemented (Not hard filtering)")
    
    if samples==None:
        print("\nPlease assign samples properly\n")
        return None

    print("\nSamples in order")
    print(samples)

    ## GRCh38 - chr1,...,chrM, while GRCh37: 1,...,MT as default. You can add 'chr' on GRCh37 too.
    if args.chrExist==".":
        if args.refFa=="GRCh38":
            chrExist="y"
        elif args.refFa=="GRCh37":
            chrExist="n"
    elif args.chrExist in ["y", "n"]:
        chrExist=args.chrExist
    else:
        chrExist="y"

    fastaPath=refPath+"00_Ref_{0}/{0}_fasta/".format(args.refFa)
    annotPath=refPath+"00_Ref_{0}/annot/".format(args.refFa)
    
    ## Reference fasta and common annotation sources
    if args.refFa == "GRCh38":
        fa=fastaPath+fasta38
        annot={k:annotPath+annot38[k] for k in annot38.keys()}
        null_PON=annotPath+"null.GRCh38.pon.vcf.gz"
    elif args.refFa == "GRCh37":
        if chrExist == "y":
            fa=fastaPath+fasta37_chr
            annot={k:annotPath+annot37_chr[k] for k in annot37_chr.keys()}
            null_PON=annotPath+"null.GRCh37.pon.chr.vcf.gz"
        elif chrExist =="n":
            fa=fastaPath+fasta37
            annot={k:annotPath+annot37[k] for k in annot37.keys()}
            null_PON=annotPath+"null.GRCh37.pon.vcf.gz"
    else:
        print("\nReference genome is not in GRCh38 or GRCh37\n")
        return None
    
    ## PON for Mutect2 or others, when assigned externally
    if args.pon == ".":
        pass
    elif args.pon == "n":
        annot["PON"]=null_PON
    elif args.pon.startswith("/"):
        annot["PON"]=args.pon
    else:
        annot["PON"]=args.directory+"/"+args.pon


    #per=scriptPath+args.variantCaller+"_chrom.sh"
    singleRunSh=scriptPath+args.variantCaller+"_run.sh"
    #print(samples)
    for sample in samples:
        callerInstance=caller(singleRunSh, args.bamSuffix, fa, annot, chrExist, sample, args.chrom, args.mem, args.directory, args.cores, args.queue)
        callerInstance.whole()

if __name__ == "__main__":
    main()
