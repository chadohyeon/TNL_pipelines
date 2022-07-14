#!/home/dcha/python3

### Ellis, et al. Nat Protoc (2020)
### Calculating ASRD/CLPM by mimicking CaVEMAN-style VCF, then utilizing SangerLCMfiltering(https://github.com/MathijsSanders/SangerLCMFiltering)
### Written by DH Cha @ TNL-KAIST
### pigz, bcftools should be on $PATH
### Using CGP-WGS and SangerLCMFiltering singularity image
### lcm_sanger_ellis.py -r GRCh38 -v Mutect2 -i 05.somatic_VCF/abcd.cvcv.mutect2.vcf  -g 01.bam/abcd.BLOOD.recal.bam -s 01.bam/abcd.BRAIN.recal.bam -c 20 
### lcm_sanger_ellis.py -r GRCh38 -v Strelka2 -i 05.somatic_VCF/abcd.cvcv.strelka2.vcf  -g 01.bam/abcd.BLOOD.recal.bam -s 01.bam/abcd.BRAIN.recal.bam -c 20 

import os, sys, argparse, glob
import subprocess
import multiprocessing
import decimal,math
import re
import time

class sangerLCM_ellis_filter():
    def __init__(self, DIR, refFa, refFasta, germline, somatic, variantCaller, vcf, mem, clpm, asmd, cores):
        self.DIR=DIR
        self.refFa=refFa        
        if self.refFa=="GRCh38":
            self.refFaSingularity="GRCh38"
        elif self.refFa=="GRCh37":
            self.refFaSingularity="GRCh37d5"
            
        self.refFasta=refFasta
        self.germline=germline
        self.somatic=somatic
        
        self.variantCaller=variantCaller
        self.vcf=vcf
        self.falseCavemanMeta="/home/public/ref/00_Ref_{0}/annot/caveman/falseCavemanMeta.txt".format(self.refFa)
        with open(self.falseCavemanMeta) as rf:
            self.falseCavemanMetaLineNum=len(rf.readlines())
        self.falseFlagCavemanMeta="/home/public/ref/00_Ref_{0}/annot/caveman/falseFlagCavemanMeta.txt".format(self.refFa)
        with open(self.falseFlagCavemanMeta) as rf:
            self.falseFlagCavemanMetaLineNum=len(rf.readlines())
        self.mem=mem
        
        self.clpm=clpm
        self.asmd=asmd
        self.cores=cores


    def splitVCF(self, vcfFn, splitNum):
        metaTxts,vcf_header,vcf_body = self.decomposeVCF(vcfFn)
        chunkSize=max(1, int(len(vcf_body)/(splitNum)))
        chunks = [vcf_body[x:x+chunkSize] for x in range(0, len(vcf_body), chunkSize)]

        splittedVcfFn=[vcfFn.replace(".vcf", ".chunk_{0}.vcf".format(str(itr).zfill(4))) for itr in range(len(chunks))]
        for itr in range(len(chunks)):
            with open(splittedVcfFn[itr], "w") as chunkVcf:
                chunkVcf.write(metaTxts)
                chunkVcf.write("\t".join(vcf_header)+"\n")
                for line in chunks[itr]:
                    chunkVcf.write("\t".join(line)+"\n")
        return(splittedVcfFn)
    
    def snvOnly(self, vcfFn):
        metaTxt, vcf_header, vcf_body = self.decomposeVCF(vcfFn)
    
        wfn=vcfFn.replace(".vcf", ".snvs.vcf")
        wfn2=vcfFn.replace(".vcf", ".indels_multi.vcf")
        
        wf=open(wfn, "w")
        wf2=open(wfn2, "w")
        
        wf.write(metaTxt)
        wf2.write(metaTxt)
        
        wf.write("\t".join(vcf_header)+"\n")
        wf2.write("\t".join(vcf_header)+"\n")
        
        for l in vcf_body:
            if len(l[vcf_header.index("REF")])==1 and len(l[vcf_header.index("ALT")])==1:
                wf.write("\t".join(l)+"\n")
            else:
                wf2.write("\t".join(l)+"\n")
                
        wf.close()
        wf2.close()
        return(wfn)
    
    def decomposeVCF(self,vcfFn):
        metaTxts=""
        vcf_body=[]
        with open(vcfFn) as vcf:
            for line in vcf:
                if line.startswith("##"):
                    metaTxts+=line
                elif line.startswith("#CHROM"):
                    vcf_header=line.rstrip().split("\t")
                else:
                    vcf_body+=[line.rstrip().split("\t")]
        return(metaTxts,vcf_header,vcf_body)
    
    def flagByLine(self, vcfFn, queue=None):
        metaTxts,vcf_header,vcf_body = self.decomposeVCF(vcfFn)
        folderFn=vcfFn.replace(".vcf", "")
        flagFn=vcfFn.replace(".vcf",".flag.vcf")
        with open(flagFn, "w") as wf:
            with open(self.falseFlagCavemanMeta) as rf:
                for l in rf:
                    wf.write(l)
        preFolderFn="{0}/preFlag".format(folderFn)
        postFolderFn="{0}/postFlag".format(folderFn)
        
        os.system("mkdir -p {0}".format(folderFn))
        os.system("mkdir -p {0}".format(preFolderFn))
        os.system("mkdir -p {0}".format(postFolderFn))
        
        lineNum=0
        for l in vcf_body:
            lineNum+=1
            oneLineFn="{0}/lineNumber{1}.vcf".format(preFolderFn, str(lineNum).zfill(6))
            oneLineOutFn="{0}/lineNumber{1}.vcf".format(postFolderFn, str(lineNum).zfill(6))
            with open(oneLineFn, "w") as oneLine:
                oneLine.write(metaTxts)
                oneLine.write("\t".join(vcf_header)+"\n")
                oneLine.write("\t".join(l)+"\n")
            self.flagging(oneLineFn, oneLineOutFn)
        lineDic={k:len(open(postFolderFn+"/"+k).readlines()) for k in os.listdir(postFolderFn) if k.endswith(".vcf")}

        for k in sorted(lineDic.keys()):
            if lineDic[k]<=self.falseFlagCavemanMetaLineNum:
                beforeCLPM=os.popen('tail -n 1 {0}'.format(preFolderFn+"/"+k)).read()
                x=beforeCLPM.rstrip().split("\t")
                infos=x[vcf_header.index("INFO")].split(";")
                infos+=["ASRD=.","CLPM=.","ASMD=."]
                x[vcf_header.index("INFO")]=";".join(infos)
                with open(flagFn, "a") as wf:
                    wf.write("\t".join(x)+"\n")
            elif lineDic[k]==self.falseFlagCavemanMetaLineNum+1:
                os.system('tail -n 1 {0} >> {1}'.format(postFolderFn+"/"+k, flagFn))
            else: print(".")
        
        if queue:
            queue.put(flagFn)
            return None
        else:
            return(flagFn)

    
    def flagging(self, rfn, wfn):
        comLi=["singularity exec --bind /home/public/ref/00_Ref_{0}/:/mnt".format(self.refFa),
               "/home/public/singularityImages/cgpwgs_2.1.1.sif /opt/wtsi-cgp/bin/cgpFlagCaVEMan.pl",
               "-i {0} -o {1}".format(rfn, wfn),
               "-m {0} -n {1}".format(self.somatic, self.germline),
               "-b /mnt/annot/caveman/flagging -umv /mnt/annot/caveman/flagging -g /mnt/annot/caveman/flagging/empty.bed.gz",
               "-ref /mnt/{0}_fasta/{1}.fai -s Human -sa {2} -t genomic".format(self.refFa, self.refFasta, self.refFaSingularity)]
        command=" ".join(comLi)

        recursionCutoff=2
        my_timeout=20
        p = subprocess.Popen(command.split(" "))
        
        try:
            p.wait(my_timeout)
            return(wfn)
        
        except subprocess.TimeoutExpired:
            p.kill()
            print("Perhabs High-depth region neglected: {0}".format(rfn))
            return(wfn)
    
    def Sanger_preselect(self, rfn, queue=None):
        wfn = rfn.replace(".vcf", ".preselect.vcf")
        comLi=["singularity run --bind /home/public/ref/00_Ref_{0}/:/mnt".format(self.refFa),
               "--app preselect",
               "/home/public/singularityImages/SangerLCMFilteringSingularity_latest.sif",
               "-a {0} -c {1}".format(self.asmd, self.clpm),
               "-v {0} > {1}".format(rfn, wfn)]
        command=" ".join(comLi)
        os.system(command)
        
        if queue:
            queue.put(wfn)
            return(None)
        else:
            return(wfn)

    def Sanger_imtateANNOVAR(self, rfn, queue=None):
        wfn = rfn.replace(".vcf", ".ANNOVARformat.txt")
        comLi=["singularity run --bind /home/public/ref/00_Ref_{0}/:/mnt".format(self.refFa),
               "--app imitateANNOVAR",
               "/home/public/singularityImages/SangerLCMFilteringSingularity_latest.sif",
               "-v {0} > {1}".format(rfn, wfn)]
        command=" ".join(comLi)
        os.system(command)
        
        if queue:
            queue.put(wfn)
            return(None)
        else:
            return(wfn)
        
    def Sanger_annotateBAMstat(self, rfn, queue=None):
        wfn = rfn.replace(".txt", ".annotateBAMstat.txt")
        
        comLi=["singularity run --bind /home/public/ref/00_Ref_{0}/:/mnt".format(self.refFa),
               "--app annotateBAMStatistics",
               "/home/public/singularityImages/SangerLCMFilteringSingularity_latest.sif",
               "-a {0}".format(rfn),
               "-b {0}".format(self.somatic),
               "-t {0}".format(self.cores),
               " > {0}".format(wfn)]
        command=" ".join(comLi)
        os.system(command)
        
        if queue:
            queue.put(wfn)
            return(None)
        else:
            return(wfn)
    
    def Sanger_additionalBAMstat(self, rfn, queue=None):
        wfn = rfn.replace(".txt", ".additional.txt")
        
        comLi=["singularity run --bind /home/public/ref/00_Ref_{0}/:/mnt".format(self.refFa),
               "--app additionalBAMStatistics",
               "/home/public/singularityImages/SangerLCMFilteringSingularity_latest.sif",
               "-a {0}".format(rfn),
               "-b {0}".format(self.somatic),
               "-t {0}".format(self.cores),
               "-c {0}".format(self.mem),
               "-r /mnt/{0}_fasta/{1}".format(self.refFa, self.refFasta),
               "-s /mnt/annot/dbSNP151_common.{0}.vcf.gz".format(self.refFa),
               " > {0}".format(wfn)]
        command=" ".join(comLi)
        print(command)
        os.system(command)
        
        if queue:
            queue.put(wfn)
            return(None)
        else:
            return(wfn)
        
    def Sanger_cruciformFilter(self, vcfFn, rfn, queue=None):
        wVcfPrefix = vcfFn.replace(".vcf", ".cruciformfilter")
        outDir="/".join(wVcfPrefix.split("/")[:-1])
        preFix=wVcfPrefix.split("/")[-1]

        comLi=["singularity run --bind /home/public/ref/00_Ref_{0}/".format(self.refFa),
               "--app filtering",
               "/home/public/singularityImages/SangerLCMFilteringSingularity_latest.sif",
               "-a {0}".format(rfn),
               "-v {0}".format(vcfFn),
               "-o {0}".format(outDir),
               "-p {0}".format(preFix)]
        command=" ".join(comLi)
        os.system(command)
        
        res=[wVcfPrefix+"_passed.vcf", wVcfPrefix+"_filtered.vcf"]
        if queue:
            queue.put(res)
            return(None)
        else:
            return(res)

    def run(self):
        if self.vcf.endswith(".vcf.gz"):
            print("decompressing .vcf.gz > .vcf")
            os.system("pigz -p {0} -k -d {1}".format(self.cores, self.vcf))
            self.vcf=self.vcf.replace(".vcf.gz", ".vcf")
            
        print("### Selecting SNVs only (indels are to be separate)")
        snvVcf=self.snvOnly(self.vcf)
        if self.variantCaller.lower() in ["mutect2", "mutect"]:
            false_caveVcf=self.mutect2ToCaveman(snvVcf, self.falseCavemanMeta)
        elif self.variantCaller.lower() in ["strelka2", "strelka"]:
            false_caveVcf=self.strelka2ToCaveman(snvVcf, self.falseCavemanMeta)
        else:
            print("[ERROR] Variant Caller should be one of ['Mutect2', 'Strelka2']")
            return(None)
            
        effCores=max(1, self.cores-1)
        
        print("{0} processes will be used".format(effCores))
        splittedVcfFn=self.splitVCF(false_caveVcf, effCores)
        
        print("### Line-by-line CLPM, ASRD calculation with WTSI-CGP Singularity")
        
        flaggingQ=multiprocessing.Queue()
        procsflaggingQ=[multiprocessing.Process(target=self.flagByLine, args=(vcfFn,flaggingQ)) for vcfFn in splittedVcfFn]
        for p in procsflaggingQ: p.start()
        for p in procsflaggingQ: p.join()
        chunkedFlagVcfs=sorted([flaggingQ.get() for p in procsflaggingQ])
        
        firstChunkedFlagVcfs=chunkedFlagVcfs[0].split(".")
        firstChunkedFlagVcfs.remove("chunk_0000")
        flagVcfName=".".join(firstChunkedFlagVcfs)

        if len(chunkedFlagVcfs)>1:
            os.system("bcftools concat {0} > {1}".format(" ".join(chunkedFlagVcfs), flagVcfName))
        else:
            os.system("cp {0} {1}".format(chunkedFlagVcfs[0], flagVcfName))
        
        
        os.system("rm -rf {0}.chunk*".format(false_caveVcf.replace(".vcf","")))
        
        os.system("rm {0}".format(false_caveVcf))
        
        
        snv_metaTxt, snv_vcf_header, snv_vcf_body=self.decomposeVCF(snvVcf)
        cave_metaTxt, cave_vcf_header, cave_vcf_body=self.decomposeVCF(flagVcfName)
        
        mergedSnvVCF=snvVcf.replace(".snvs.vcf", ".snvs.ASRD_CLPM.vcf")
        if len(snv_vcf_body)==len(cave_vcf_body):
            print("\nMerging {0}_SNV and CaVEMan post-processing filters\n".format(self.variantCaller))
            with open(mergedSnvVCF, "w") as wf:
                wf.write(snv_metaTxt)
                wf.write('##INFO=<ID=ASRD,Number=1,Type=Float,Description="A soft flag median (read length adjusted) alignment score of reads showing the variant allele">\n')
                wf.write('##INFO=<ID=CLPM,Number=1,Type=Float,Description="A soft flag median number of soft clipped bases in variant supporting reads">\n')
                wf.write('##INFO=<ID=ASMD,Number=1,Type=Float,Description="A soft flag median alignement score of reads showing the variant allele">\n')
                wf.write("\t".join(snv_vcf_header)+"\n")
                for itr in range(len(snv_vcf_body)):
                    x=snv_vcf_body[itr]
                    x_info=x[snv_vcf_header.index("INFO")].split(";")
                    
                    y=cave_vcf_body[itr]
                    y_info=y[cave_vcf_header.index("INFO")].split(";")
                    x_info+=y_info[-3:]
                    x[snv_vcf_header.index("INFO")]=";".join(x_info)
                    wf.write("\t".join(x)+"\n")
        
            os.system("rm {0}".format(flagVcfName))
        else:
            print("Flagged and original VCF differs... Something went wrong")
            return(None)
       
        print("### Sanger LCM Filter (https://github.com/MathijsSanders/SangerLCMFiltering), with Singularity image\n")
        print("### Sanger LCM Filter: 01. Preselect VCF\n")
        preselectVCF=self.Sanger_preselect(mergedSnvVCF)
        print("### Sanger LCM Filter: 02. Immitate ANNOVAR\n")
        imtateANNOVARtxt=self.Sanger_imtateANNOVAR(preselectVCF)
        print("### Sanger LCM Filter: 03. Annotate Somatic BAM statistics\n")
        annotateBAMstatTxt=self.Sanger_annotateBAMstat(imtateANNOVARtxt)
        print("### Sanger LCM Filter: 04. Additional Somatic BAM statistics\n")
        additionalBAMstatTxt=self.Sanger_additionalBAMstat(annotateBAMstatTxt)
        print("### Sanger LCM Filter: 05. Filtering VCF with annotated LCM parameters\n")
        sangerFilteredPASS, sangerFilteredRemainder=self.Sanger_cruciformFilter(mergedSnvVCF, additionalBAMstatTxt)
        print("### Sanger LCM Filter: DONE!")
        
        os.system("rm {0}".format(preselectVCF))
        os.system("rm {0}".format(imtateANNOVARtxt))
        os.system("rm {0}".format(annotateBAMstatTxt))
        os.system("rm {0}".format(sangerFilteredRemainder))
        os.system("rm {0}".format(mergedSnvVCF))
        
    def mutect2ToCaveman(self, vcfFn, falseCavemanMeta, randomness=False):
        metaTxt, vcf_header, vcf_body = self.decomposeVCF(vcfFn)
        falseMetaTxt, _, _ = self.decomposeVCF(falseCavemanMeta)
        wfn=vcfFn.replace(".vcf", ".caveman.vcf")
        vcf_header[-2]="NORMAL"
        vcf_header[-1]="TUMOUR"

        def toDecimal(num):
            return(('%.1E' % decimal.Decimal(str(float(num)))).replace("E","e"))
        with open(wfn, "w") as wf:
            wf.write(falseMetaTxt)
            wf.write("\t".join(vcf_header)+"\n")
            for l in vcf_body:
                formats=l[vcf_header.index("FORMAT")].split(":")
                ref=l[vcf_header.index("REF")]
                alt=l[vcf_header.index("ALT")]
                dp=0
                l[vcf_header.index("FORMAT")]="GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM"
                for status in ["NORMAL", "TUMOUR"]:
                    gts=["FA", "FC", "FG", "FT", "RA", "RC", "RG", "RT"]
                    gtDic={gt:"0" for gt in gts}
                    fR, rR, fA, rA = l[vcf_header.index(status)].split(":")[formats.index("SB")].split(",")
                    gtDic["F"+ref], gtDic["R"+ref], gtDic["F"+alt], gtDic["R"+alt] = fR, rR, fA, rA
                    dp=dp+int(fA)+int(rA)+int(fR)+int(rR)
                    try:
                        pm=toDecimal((int(fA)+int(rA))/(int(fA)+int(rA)+int(fR)+int(rR)))
                    except:
                        pm=toDecimal(0)
                    if status=="NORMAL":
                        l[vcf_header.index("NORMAL")]=":".join(["0|0"]+[gtDic[gt] for gt in gts]+[pm])
                    elif status=="TUMOUR":
                        l[vcf_header.index("TUMOUR")]=":".join(["0|1"]+[gtDic[gt] for gt in gts]+[pm])
                        try:
                            vaf=(int(fA)+int(rA))/(int(fA)+int(rA)+int(fR)+int(rR))*10
                        except:
                            vaf=0

                rou=int(round(vaf, 0))
                cei=math.ceil(vaf)
                flo=math.floor(vaf)

                if rou==0:
                    gt1num=1
                    gt2num=0
                else:
                    gt1num=rou
                    if cei==rou:
                        gt2num=flo
                    else:
                        gt2num=cei

                TG="".join([ref,ref,"/"]+sorted([alt for k in range(gt1num)]+[ref for k in range(10-gt1num)]))
                SG="".join([ref,ref,"/"]+sorted([alt for k in range(gt2num)]+[ref for k in range(10-gt2num)]))
                mp, gp, tp, sp = "7.0e-1", "5.0e-2", "4.5e-1", "2.2e-1"
                if randomness:
                        gp=np.random.rand()*0.1
                        mp=0.4+np.random.rand()*0.6

                        tp=(0.5+np.random.rand()*0.5)*mp
                        sp=min((np.random.rand()*0.5)*mp, (mp-tp)*(0.8+np.random.rand()*0.2))
                mp, gp, tp, sp = toDecimal(mp), toDecimal(gp), toDecimal(tp), toDecimal(sp)  
                l[vcf_header.index("INFO")]="DP={0};MP={3};GP={4};TG={1};TP={5};SG={2};SP={6}".format(dp, TG, SG, mp, gp, tp ,sp)
                wf.write("\t".join(l)+"\n")

        return(wfn)
    
    def strelka2ToCaveman(self, vcfFn, falseCavemanMeta, randomness=False):
        metaTxt, vcf_header, vcf_body = self.decomposeVCF(vcfFn)
        falseMetaTxt, _, _ = self.decomposeVCF(falseCavemanMeta)
        wfn=vcfFn.replace(".vcf", ".caveman.vcf")
        vcf_header[-1]="TUMOUR"

        def toDecimal(num):
            return(('%.1E' % decimal.Decimal(str(float(num)))).replace("E","e"))
        with open(wfn, "w") as wf:
            wf.write(falseMetaTxt)
            wf.write("\t".join(vcf_header)+"\n")
            for l in vcf_body:
                formats=l[vcf_header.index("FORMAT")].split(":")
                ref=l[vcf_header.index("REF")]
                alt=l[vcf_header.index("ALT")]
                dp=0
                l[vcf_header.index("FORMAT")]="GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM"
                for status in ["NORMAL", "TUMOUR"]:
                    gts=["FA", "FC", "FG", "FT", "RA", "RC", "RG", "RT"]
                    gtDic={gt:"0" for gt in gts}
                    gtDic["FA"]= l[vcf_header.index(status)].split(":")[formats.index("AU")].split(",")[0]
                    gtDic["FC"]= l[vcf_header.index(status)].split(":")[formats.index("CU")].split(",")[0]
                    gtDic["FG"]= l[vcf_header.index(status)].split(":")[formats.index("GU")].split(",")[0]
                    gtDic["FT"]= l[vcf_header.index(status)].split(":")[formats.index("TU")].split(",")[0]
                    dp=dp+int(gtDic["FA"])+int(gtDic["FC"])+int(gtDic["FG"])+int(gtDic["FT"])
                    try:
                        pm=toDecimal(int(gtDic["F"+alt])/(int(gtDic["FA"])+int(gtDic["FC"])+int(gtDic["FG"])+int(gtDic["FT"])))
                    except:
                        pm=toDecimal(0)
                    if status=="NORMAL":
                        l[vcf_header.index("NORMAL")]=":".join(["0|0"]+[gtDic[gt] for gt in gts]+[pm])
                    elif status=="TUMOUR":
                        l[vcf_header.index("TUMOUR")]=":".join(["0|1"]+[gtDic[gt] for gt in gts]+[pm])
                        try:
                            vaf=int(gtDic["F"+alt])/(int(gtDic["FA"])+int(gtDic["FC"])+int(gtDic["FG"])+int(gtDic["FT"]))*10
                        except:
                            vaf=0
                rou=int(round(vaf, 0))
                cei=math.ceil(vaf)
                flo=math.floor(vaf)

                if rou==0:
                    gt1num=1
                    gt2num=0
                else:
                    gt1num=rou
                    if cei==rou:
                        gt2num=flo
                    else:
                        gt2num=cei

                TG="".join([ref,ref,"/"]+sorted([alt for k in range(gt1num)]+[ref for k in range(10-gt1num)]))
                SG="".join([ref,ref,"/"]+sorted([alt for k in range(gt2num)]+[ref for k in range(10-gt2num)]))
                mp, gp, tp, sp = "7.0e-1", "5.0e-2", "4.5e-1", "2.2e-1"
                if randomness:
                        gp=np.random.rand()*0.1
                        mp=0.4+np.random.rand()*0.6

                        tp=(0.5+np.random.rand()*0.5)*mp
                        sp=min((np.random.rand()*0.5)*mp, (mp-tp)*(0.8+np.random.rand()*0.2))
                mp, gp, tp, sp = toDecimal(mp), toDecimal(gp), toDecimal(tp), toDecimal(sp)  
                l[vcf_header.index("INFO")]="DP={0};MP={3};GP={4};TG={1};TP={5};SG={2};SP={6}".format(dp, TG, SG, mp, gp, tp ,sp)
                wf.write("\t".join(l)+"\n")

        return(wfn)
def checkBam(bam, bamSuffix):
    endsWithSuffix=bam.endswith(bamSuffix)
    if endsWithSuffix:
        if glob.glob("01.bam/"+bam):
            bam="01.bam/"+bam
        elif glob.glob(bam):
            bam=bam
        else:
            return(None)
    elif bam.endswith(".bam"):
        if glob.glob("01.bam/"+bam):
            bam="01.bam/"+bam
        elif glob.glob(bam):
            bam=bam
        else:
            return(None)
    else:
        if glob.glob("01.bam/"+bam+bamSuffix):
            bam="01.bam/"+bam+bamSuffix
        elif glob.glob(bam+bamSuffix):
            bam=bam+bamSuffix
        else:
            return(None)
    return(bam)   

def main():    
    parser=argparse.ArgumentParser()
    parser.add_argument('-r', '--refFa', required=True, type=str, help='GRCh38 or GRCh37', metavar='') ## GRCh37, GRCh38
    parser.add_argument('-v', '--variantCaller', required=True, type=str, help='Mutect2 or Strelka2', metavar='') ## Mutect2, Strelka2
    parser.add_argument('-i', '--input', required=True, type=str, help='Input VCF', default=".", metavar='') ## NA12375_blood,NA12333_blood, ...
    parser.add_argument('-g', '--germline', required=True, type=str, help='Germline sample',metavar='') ## NA12375_blood,NA12333_blood, ...
    parser.add_argument('-s', '--somatic', required=True, type=str, help='Somatic sample', metavar='') ## NA12375_GBM,NA_12333_GBM, ...
    
    parser.add_argument('-c', '--cores', required=False, type=str, help='Multiple cores to use', default="1", metavar='')
    parser.add_argument('-m', '--mem', required=False, type=str, help='Java Xmx heap size', default="10", metavar='')

    parser.add_argument('--clpm', required=False, type=str, help='CLPM score threshold (Default: 0)', default="0", metavar='')
    parser.add_argument('--asmd', required=False, type=str, help='ASMD score threshold (Default: 140)', default="140", metavar='')

    parser.add_argument('-b', '--bamSuffix', required=False, type=str, help="BAM file suffix", default='.marked.recal.bam', metavar='') ## Provide it when needed
    parser.add_argument('-w', '--directory', required=False, type=str, help="Default workdir will be PWD", default=os.getcwd(), metavar='') ## Provide it when needed
    args=parser.parse_args()
    #parser.print_help()
    
    # Germline samples (automatically or assigned externally)
    germline_bam=checkBam(args.germline, args.bamSuffix)
    if not germline_bam:
        print("Germline samples were not properly assigned")
        return(None)
    if glob.glob(germline_bam.replace(".bam", ".bai")):
        os.system("mv {0}.bai {0}.bam.bai".format(germline_bam.replace(".bam", "")))
        
    ## Somatic samples (should assigned externally)
    somatic_bam=checkBam(args.somatic, args.bamSuffix)
    if not somatic_bam:
        print("Somatic samples were not properly assigned")
        return(None)
    
    if glob.glob(somatic_bam.replace(".bam", ".bai")):
    #if os.popen("ls {0}".format(somatic_bam.replace(".bam", ".bai"))).read():
        os.system("mv {0}.bai {0}.bam.bai".format(somatic_bam.replace(".bam", "")))
    
    print("Germline: {0}, Somatic: {1} detected".format(germline_bam, somatic_bam))
    
    refFa=args.refFa
    
    if refFa=="GRCh38":
        refFasta="Homo_sapiens_assembly38.fa"
    elif refFa=="GRCh37":
        refFasta="Homo.Sapiens.GRCh37.75.primary_assembly.full.fa"
    else:
        print("Assembly {0} doesn't exist".format(refFa))
        return(None)
    
    instance=sangerLCM_ellis_filter(args.directory, refFa, refFasta, germline_bam, somatic_bam, args.variantCaller, args.input, args.mem, args.clpm, args.asmd, int(args.cores))
    instance.run()                                          
if __name__ == "__main__":
    main()
