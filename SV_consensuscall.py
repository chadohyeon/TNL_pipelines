#!/home/dcha/python3
import os
import time
import numpy as np
import argparse
from multiprocessing import Process


refDIR = "/home/public/ref"
pwd = os.getcwd()
class SV_call():
    def __init__(self, normalbam, tumorbam, refFa, core, DIR, bindDIR):
        self.normalbam = normalbam
        self.tumorbam = tumorbam
        self.refFa = refFa
        if self.refFa == "GRCh38":
            refFasta = "/00_Ref_GRCh38/GRCh38_fasta/Homo_sapiens_assembly38.fa"
        elif self.refFa == "GRCh37":
            refFasta = "/00_Ref_GRCh37/GRCh37_fasta/Homo.Sapiens.GRCh37.75.primary_assembly.full.fa"
        else:
            print("please check reference fasta file. only GRCh38, GRCh37 supports")
        self.refFasta = refFasta
        self.DIR = DIR
        self.bindDIR = bindDIR
        self.core = core
        #tumor = self.tumorbam.replace(".marked.recal.bam", "").split("/")[-1]
        #tu_no = tumor+"_"+"somatic"
        #self.tu_no = tu_no
    
    def samplename(self):
        sn = self.tumorbam.replace(".marked.recal.bam", "").split("/")[-1]
        sn += "_"+"somatic"
        return(sn)
    
    def manta(self, DIR, bindDIR, refFa, core, tumorbam, normalbam):
        core_manta = int((self.core-4)/2)
        if core_manta < 1:
            print("Manta needs at least one core. This pipeline needs at least 6 cores")
            return(None)
        os.system("mkdir {0}/mantaresult".format(self.DIR))
        _cml = ["singularity exec --cleanenv --bind {0}:/mnt".format(self.bindDIR),
              "/home/public/singularityImages/strelka2-manta_latest.sif",
              "configManta.py",
              "--normalBam {0} --tumorBam {1}".format(self.normalbam, self.tumorbam),
              "--referenceFasta /mnt{0}".format(self.refFasta),
              "--runDir {0}/mantaresult".format(self.DIR)]
        cml = " ".join(_cml)
        os.system(cml)
        _cml = ["singularity exec --cleanenv --bind {0}:/mnt".format(self.bindDIR),
              "/home/public/singularityImages/strelka2-manta_latest.sif",
                "{0}/mantaresult/runWorkflow.py -j {1} --quiet".format(self.DIR, core_manta)]
        cml = " ".join(_cml)
        print("running manta")
        start = time.time()
        os.system(cml)
        print(cml)
        print("manta finished"+ " it takes {0} minutes".format((time.time() - start)//60))
        
    def delly(self, DIR, refFa, tumorbam, normalbam  ):
        os.system("mkdir {0}/dellyresult".format(self.DIR))
        #_cml = [delly call -x /home/jsh/jobs/pcawg/delly/excludeTemplates/human.hg38.excl.tsv 
        #        -o /home/jsh/jobs/pcawg/dellyresult/ALS187.delly.bcf 
        #        -g /home/jsh/jobs/ref/Homo_sapiens_assembly38.fa 
        #        /home/jsh/jobs/ALS187/A187_brain.marked.recal.bam /home/neurologistkim/ALS/Sanger_filters/01.bam/02421_val.marked.recal.bam]
        _cml = ["delly call -x /home/public/ref/00_Ref_{0}/annot/human.{0}.excl.tsv".format(self.refFa),
               "-o {0}/dellyresult/{1}.delly.bcf".format(self.DIR, self.samplename()),
               "-g /home/public/ref{0}".format(self.refFasta),
               "{0} {1}".format(self.tumorbam, self.normalbam)]
        cml = " ".join(_cml)
        print("running delly")
        start = time.time()
        print(cml)
        os.system(cml)
        sampletsv = [i.rstrip() for i in list(os.popen("bcftools query -l {0}/dellyresult/{1}.delly.bcf".format(self.DIR, self.samplename())))]
        with open("{0}/samples.tsv".format(self.DIR), "w") as wf:
            wf.write(sampletsv[0] + "\t" + "tumor\n")
            wf.write(sampletsv[1] + "\t" + "control")
        _cml2 = ["delly filter -f somatic -o {0}/dellyresult/{1}.delly.filtered.bcf".format(self.DIR, self.samplename()),
               "-s {0}/samples.tsv {0}/dellyresult/{1}.delly.bcf".format(self.DIR, self.samplename())]
        cml2 = " ".join(_cml2)
        os.system(cml2)
        print("delly finished" + " it takes {0} minutes".format((time.time() - start)//60))
        
    def gridss2(self, DIR, bindDIR, refFa, core, tumorbam, normalbam):
        core_gridss = int((self.core-4)/2)
        if core_gridss < 1:
            print("gridss needs at least one core. This pipeline needs at least 6 cores")
            return(None)
        os.system("mkdir {0}/gridssresult".format(self.DIR))
        #singularity run --bind /home/neurologistkim:/mnt /home/jsh/jobs/pcawg/gridss2.sif gridss 
        #    -r /home/jsh/jobs/ref/Homo_sapiens_assembly38.fa 
        #    -o /home/jsh/sandbox/ALS187.gridss.vcf 
        #    -w /home/jsh/sandbox /home/jsh/jobs/ALS187/A187_brain.marked.recal.bam /mnt/ALS/Sanger_filters/01.bam/02421_val.marked.recal.bam
        _cml = ["singularity run --bind {0}:/mnt".format(self.bindDIR),
               "/home/public/singularityImages/gridss2.sif gridss",
               "-r /mnt{0}".format(self.refFasta),
                "-t {0}".format(core_gridss),
               "-o {0}/gridssresult/{1}.gridss.vcf".format(self.DIR, self.samplename()),
               "-w {0}".format(self.DIR),
                "{0} {1}".format(self.normalbam, self.tumorbam)]
        
        cml = " ".join(_cml)
        print("running gridss")
        start = time.time()
        os.system(cml)
        print("gridss finished"+ " it takes {0} minutes".format((time.time() - start)//60))
        
    def lumpy(self, DIR, tumorbam, normalbam):
        os.system("mkdir {0}/lumpyresult".format(self.DIR))
        tu_sn = self.tumorbam[:-4]
        no_sn = self.normalbam[:-4]
        _cml = ["singularity exec /home/public/singularityImages/lumpy.sif lumpyexpress",
               "-B {0},{1}".format(self.tumorbam, self.normalbam),
               "-S {0}.splitters.sorted.bam,{1}.splitters.sorted.bam".format(tu_sn, no_sn),
               "-D {0}.discordants.sorted.bam,{1}.discordants.sorted.bam".format(tu_sn, no_sn),
               "-o {0}/lumpyresult/{1}.lumpy.vcf".format(self.DIR, self.samplename())]
        cml = " ".join(_cml)
        print("running lumpmy")
        start = time.time()
        os.system(cml)
        print("lumpy finished"+ " it takes {0} minutes".format((time.time() - start)//60))
        
    def discordant(self, bam):
        sn = bam[:-4]
        cml = "samtools view -b -F 1294 {0} > {1}.discordants.unsorted.bam".format(bam, sn)
        os.system(cml)
        cml2 = "samtools sort {0}.discordants.unsorted.bam -o {0}.discordants.sorted.bam".format(sn)
        os.system(cml2)
        cml3 = "samtools index {0}.discordants.sorted.bam".format(sn)
        os.system(cml3)
        
    def splitread(self, bam):
        sn = bam[:-4]
        _cml = ["samtools view -h {0}".format(bam),
              "| extractSplitReads_BwaMem -i stdin | samtools view -Sb - > {0}.splitters.unsorted.bam".format(sn)]
        cml = " ".join(_cml)
        os.system(cml)
        cml2 = "samtools sort {0}.splitters.unsorted.bam -o {0}.splitters.sorted.bam".format(sn)
        os.system(cml2)
        cml3 = "samtools index {0}.splitters.sorted.bam".format(sn)
        os.system(cml3)
        
    def run(self):
        manta_process = Process(target = self.manta, args=(self.DIR, self.bindDIR, self.refFa, self.core, self.tumorbam, self.normalbam))
        delly_process = Process(target = self.delly, args=(self.DIR, self.refFa, self.tumorbam, self.normalbam ))
        gridss_process = Process(target = self.gridss2, args=(self.DIR, self.bindDIR, self.refFa, self.core, self.tumorbam, self.normalbam))
        discordant_process_tumor = Process(target = self.discordant, args=[self.tumorbam])
        discordant_process_normal = Process(target = self.discordant, args=[self.normalbam])
        splitread_process_tumor = Process(target = self.splitread, args=[self.tumorbam])
        splitread_process_normal = Process(target = self.splitread, args=[self.normalbam])
        lumpy_process = Process(target = self.lumpy, args=(self.DIR, self.tumorbam, self.normalbam))
        
        print("run manta for {0} as tumor sample and {1} as normal sample in {2}".format(self.tumorbam, self.normalbam, self.DIR))
        manta_process.start()
        
        print("run delly for {0} as tumor sample and {1} as normal sample in {2}".format(self.tumorbam, self.normalbam, self.DIR))
        delly_process.start()
        
        print("run gridss for {0} as tumor sample and {1} as normal sample in {2}".format(self.tumorbam, self.normalbam, self.DIR))
        gridss_process.start()
        
        print("make discordant paried-end alignments for lumpy")
        discordant_process_normal.start()
        discordant_process_tumor.start()

        print("make split-read alignments for lumpy")
        splitread_process_tumor.start()
        splitread_process_normal.start()
        
        manta_process.join()
        delly_process.join()
        gridss_process.join()
        discordant_process_normal.join()
        discordant_process_tumor.join()
        splitread_process_tumor.join()
        splitread_process_normal.join()
        lumpy_process.start()
        lumpy_process.join()

    def lumpy_preprocess(self):
        discordant_process_tumor = Process(target = self.discordant, args=[self.tumorbam])
        discordant_process_normal = Process(target = self.discordant, args=[self.normalbam])
        splitread_process_tumor = Process(target = self.splitread, args=[self.tumorbam])
        splitread_process_normal = Process(target = self.splitread, args=[self.normalbam])
        
        print("make discordant paried-end alignments for lumpy")
        discordant_process_normal.start()
        discordant_process_tumor.start()

        print("make split-read alignments for lumpy")
        splitread_process_tumor.start()
        splitread_process_normal.start()
        
        discordant_process_normal.join()
        discordant_process_tumor.join()
        splitread_process_tumor.join()
        splitread_process_normal.join()
    
def survivor(DIR, samplename):
    
    os.system("mkdir {0}/SURVIVOR".format(DIR))
    os.system("cp {0}/mantaresult/results/variants/somaticSV.vcf.gz {0}/SURVIVOR/{1}.manta.vcf.gz".format(DIR, samplename))
    os.system("gunzip {0}/SURVIVOR/{1}.manta.vcf.gz".format(DIR, samplename))
#     dellyresult_list = os.popen("ls {0}/dellyresult/*.delly.filtered.bcf".format(DIR))
#     dellyresult_file = [k.rstrip() for k in list(dellyresult_list)]
#     dellyresult_file = "".join(dellyresult_file).replace(".bcf", "")
    dellyresult_file = os.listdir(DIR + "/dellyresult/")[0]
    os.system("bcftools view {0}/dellyresult/{1}.bcf > {0}/dellyresult/{1}.vcf".format(DIR, dellyresult_file.replace(".bcf", "")))
    os.system("cp {0}/dellyresult/*.vcf {0}/SURVIVOR/".format(DIR))
    os.system("cp {0}/gridssresult/*.vcf {0}/SURVIVOR/".format(DIR))
    os.system("cp {0}/lumpyresult/*.vcf {0}/SURVIVOR/".format(DIR))
    cml = "ls *.vcf > {0}/{1}.sample_files".format(DIR, samplename)
    with open(DIR + "/{0}.sample_files".format(samplename), "w") as wf:
        for vcf in list(os.popen("ls {0}/SURVIVOR/*.vcf".format(DIR))):
            wf.write(vcf)
#     os.system("bash {0}/sample_files.sh".format(DIR))
    os.system("SURVIVOR merge {0}/{1}.sample_files 1000 2 1 1 0 30 {0}/SURVIVOR/{1}.merged.vcf".format(DIR, samplename))
    
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-n', '--normalbam', required=True, type=str, help='input control bam file', metavar='')
    parser.add_argument('-t', '--tumorbam', required=True, type=str, help='input tumor bam file', metavar='')
    parser.add_argument('-r', '--refFa', required=True, type=str, help='GRCh38 or GRCh37', metavar='') ## GRCh37, GRCh38
    parser.add_argument('-c', '--core', required=False, type=str, help='Multiple cores to use more than 6 needed', default="6", metavar='')
    parser.add_argument('-w', '--directory', required=False, type=str, help="Default workdir will be PWD", default=os.getcwd(), metavar='') ## Provide it when needed
    parser.add_argument('-b', '--bindDIR', required=False, type=str, help="Default bind directory", default="/home/public/ref", metavar='')
    args=parser.parse_args()
    parser.print_help()
        
    instance = SV_call(args.normalbam, args.tumorbam, args.refFa, int(args.core), args.directory, args.bindDIR)
    instance.run()
    survivor(args.directory, instance.samplename())
    
if __name__=="__main__":
    main()