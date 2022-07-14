#!/home/dcha/python3

## Written by Do Hyeon Cha, TNL@KAIST
## Wrapper and customized scripts utilize bcftools / VEP / pigz, ... 
## Custom annotation can be added when properly made as BED, bgzipped and tabix-indexed
##
## <Things to do first>
## Source /etc/profile
## conda activate
##
## fullAnnotVCF.py -r GRCh38 -c 20 -i abcd.vcf (Only Final TSV; Standard usage) <- Automatically detect the caller used by parsing metatexts

## fullAnnotVCF.py -r GRCh38 -c 20 -i abcd.vcf -m n <- Only annotating dbNSFP+SpliceAI+LOFTEE
## fullAnnotVCF.py -r GRCh38 -c 20 -i abcd.vcf -m g <- general (Annotating dbNSFP+SpliceAI+LOFTEE and numerous custom annotations)
## fullAnnotVCF.py -r GRCh38 -c 20 -i abcd.vcf -m b <- Brain-specific annotation mode (general+brain)
## fullAnnotVCF.py -r GRCh38 -c 20 -i abcd.vcf -m c <- cancer-specific annotation mode (general+cancer)
## fullAnnotVCF.py -r GRCh38 -c 20 -i abcd.vcf -m a <- Annotating all possible sources (general+brain+cancer)

## fullAnnotVCF.py -r GRCh38 -c 20 -i abcd.vcf -f y (Final TSV + Final VCF; slower)
## fullAnnotVCF.py -r GRCh38 -c 20 -v Mutect2 -i abcd.vcf (Mutect2 as variant caller)

    
import sys, os
import subprocess
import multiprocessing
import argparse
import pandas as pd
import re
from collections import OrderedDict

def bcftoolsAnnot(vcf, bcftools, annotDBdir, customDict, queue=None):
    tmpVcf=vcf.replace(".vcf", ".bcftools_tmp.vcf")
    annotVcf=vcf.replace(".vcf", ".bcftools.vcf")
    
    os.system("cp {0} {1}".format(vcf, tmpVcf))
    for k,v in customDict["multiple"].items():
        command="{0} annotate --force -a {1} -C {2} -h {3} {4} > {5}".format(bcftools, annotDBdir+"/"+v["file"], annotDBdir+"/"+v["columnsFn"], annotDBdir+"/"+v["infoFn"], tmpVcf, annotVcf)
        with open(annotDBdir+"/"+v["infoFn"]) as rf:
            print("".join(rf.readlines()))
        os.system(command)
        os.system("mv {0} {1}".format(annotVcf, tmpVcf))
    
    for colName,fn in customDict["single"].items():
        falseInfo='##INFO=<ID={0},Number=.,Type=String,Description="{0}">'.format(colName)
        with open(tmpVcf+"_tmpFalseInfo", "w") as wf:
            wf.write(falseInfo)
        print(falseInfo)
        command="{0} annotate --force -a {1} -c CHROM,FROM,TO,{2} -h {3}_tmpFalseInfo {3} > {4}".format(bcftools, annotDBdir+"/"+fn, colName, tmpVcf, annotVcf)
        os.system(command)
        os.system("mv {0} {1}".format(annotVcf, tmpVcf))
        os.system("rm {0}_tmpFalseInfo".format(tmpVcf))
        
    for colName,fn in customDict["allelic"].items():
        falseInfo='##INFO=<ID={0},Number=.,Type=String,Description="{0}">'.format(colName)
        with open(tmpVcf+"_tmpFalseInfo", "w") as wf:
            wf.write(falseInfo)
        print(falseInfo)
        command="{0} annotate --force -a {1} -c CHROM,FROM,-,REF,ALT,{2} -h {3}_tmpFalseInfo {3} > {4}".format(bcftools, annotDBdir+"/"+fn, colName, tmpVcf, annotVcf)
        os.system(command)
        os.system("mv {0} {1}".format(annotVcf, tmpVcf))
        os.system("rm {0}_tmpFalseInfo".format(tmpVcf))
    
    os.system("mv {0} {1}".format(tmpVcf, annotVcf))
    
    if queue:
        queue.put(annotVcf)
        return None
    else:
        return(annotVcf)

def vepAnnot(vcf, outFormat, vepPath, assembly, ref_genome, annotDBdir, suffix="Annot", customDict=None, queue=None):
    ### PERL5LIB setup
    perl5lib=os.environ['PERL5LIB']
    if os.environ['PERL5LIB'].endswith(":"):
        os.environ['PERL5LIB']+="/home/public/tools/ensembl-vep/.vep/Plugins/loftee_{0}:".format(assembly)
    else:
        os.environ['PERL5LIB']+=":/home/public/tools/ensembl-vep/.vep/Plugins/loftee_{0}:".format(assembly)
    
    vcfPrefix=vcf.split(".vcf")[0]
    #cmd=["{0}/vep --dir  {0}/.vep --merged --assembly {1} --hgvs --fasta {2}".format(vepPath, assembly, ref_genome),
    cmd=["{0}/vep --dir {0}/.vep --cache --merged --assembly {1} --use_given_ref".format(vepPath, assembly),
        "--format vcf -i {0} -o {1}.{2}.{3} --{3}".format(vcf, vcfPrefix, suffix, outFormat),
        "--no_stats --force_overwrite --offline", 
        "--canonical --appris --protein --biotype --uniprot --tsl --ccds",
        "--pick_allele --pick_order canonical, appris, tsl, biotype, ccds, qrank, length",
        "--variant_class --pubmed",
        "--dir_plugins {0}/.vep/Plugins".format(vepPath)]
    
    ### dbNSFP v4.2
    cmd+=["--plugin dbNSFP,{0}/dbNSFP4.2a_{1}.gz,ALL".format(annotDBdir, assembly)]

    ### 5'-UTR annotator
    cmd+=["--plugin UTRannotator,{0}/.vep/Plugins/UTRannotator/uORF_starts_ends_{1}_PUBLIC.txt".format(vepPath, assembly)]
    
    ### SpliceAI
    cmd+=["--plugin SpliceAI,snv={0}/spliceai_scores.raw.snv.{1}.vcf.gz,indel={0}/spliceai_scores.raw.indel.{1}.vcf.gz,cutoff=0.5".format(annotDBdir, assembly)]
    
    ### LOFTEE
    if assembly=="GRCh38":
        cmd+=["--plugin LoF,loftee_path:{0}/.vep/Plugins/loftee_GRCh38,human_ancestor_fa:{1}/human_ancestor.fa.gz,conservation_file:{1}/loftee_conservation.sql,gerp_bigwig:{1}/gerp_conservation_scores.homo_sapiens.GRCh38.bw".format(vepPath, annotDBdir)]
    elif assembly=="GRCh37":
        cmd+=["--plugin LoF,loftee_path:{0}/.vep/Plugins/loftee_GRCh37,human_ancestor_fa:{1}/human_ancestor.fa.gz,conservation_file:{1}/loftee_conservation.sql".format(vepPath, annotDBdir)]

    os.system(" ".join(cmd))
    os.environ['PERL5LIB']=perl5lib
    output="{0}.{1}.{2}".format(vcfPrefix,suffix,outFormat)
    if queue:
        queue.put(output)
        return None
    else:
        return(output)

def splitVCF(vcfFn, splitNum):
    metaTxts,vcf_header,vcf_body = decomposeVCF(vcfFn)
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
            
def decomposeVCF(vcfFn):
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

def vcfAddMoreInfo(vcfFn, caller):
    wholeLineAnnotDic={"QUAL":[], "FILTER":[], "INFO":[], "FORMAT":[]}
    _, vcf_header, vcf_body=decomposeVCF(vcfFn)
    
    preCallerInfoCols=[]
    preCallerFormatCols=[]
    
    for line in vcf_body:
        #print(line)
        wholeLineAnnotDic["QUAL"]+=[{"QUAL":line[vcf_header.index("QUAL")]}]
        wholeLineAnnotDic["FILTER"]+=[{"FILTER":line[vcf_header.index("FILTER")]}]
        
        vcfLineInfo=line[vcf_header.index("INFO")]
        
        localCallerInfoDic={}
        for j in vcfLineInfo.split(";"):
            if "=" in j:
                localCallerInfoDic[j.split("=")[0]]=j.split("=")[1]
            else:
                localCallerInfoDic[j]="YES"
        wholeLineAnnotDic["INFO"]+=[localCallerInfoDic]

        preCallerInfoCols+=localCallerInfoDic.keys()
        preCallerInfoCols=list(set(preCallerInfoCols))
        
        sampleNames=vcf_header[vcf_header.index("FORMAT")+1:]
        vcfLineFormatHeader=line[vcf_header.index("FORMAT")].split(":")
        vcfLineFormatContentsLi=[string.split(":") for string in line[vcf_header.index("FORMAT")+1:]]
        vcfLineFormatContentsDic={}
        
        vcfLineFormatContentsDic={}
        for itr1 in range(len(sampleNames)):
            sampleN=sampleNames[itr1]
            localFormatLi=vcfLineFormatContentsLi[itr1]
            for itr2 in range(len(vcfLineFormatHeader)):
                if vcfLineFormatHeader[itr2]=="GT":
                    localFormatLi[itr2]="|".join(localFormatLi[itr2].split("/"))
                vcfLineFormatContentsDic[vcfLineFormatHeader[itr2]+":"+sampleN]=localFormatLi[itr2]

        
        for sN in sampleNames:
            if caller=="HaplotypeCaller":
                perSampleAD=[float(readDepth) for readDepth in vcfLineFormatContentsDic["AD:"+sN].split(",")]
                perSampleTotalDepth=float(vcfLineFormatContentsDic["DP:"+sN])
                if perSampleTotalDepth==0:
                    perSampleVAF=",".join(["NA" for altD in perSampleAD[1:]])
                else:
                    perSampleVAF=",".join([str(round(altD/perSampleTotalDepth,3)) for altD in perSampleAD[1:]])
                
                vcfLineFormatContentsDic["VAF:"+sN]=perSampleVAF
            elif caller=="Mutect2":
                vafEstim=vcfLineFormatContentsDic["AF:"+sN]
                forRoundUp=max(len(k.split(".")[1]) for k in vafEstim.split(","))
                perSampleAD=[float(readDepth) for readDepth in vcfLineFormatContentsDic["AD:"+sN].split(",")]
                perSampleTotalDepth=float(vcfLineFormatContentsDic["DP:"+sN])
                if perSampleTotalDepth==0:
                    perSampleVAF=",".join(["NA" for altD in perSampleAD[1:]])
                else:
                    perSampleVAF=",".join([str(round(altD/perSampleTotalDepth,forRoundUp)) for altD in perSampleAD[1:]])
                
                vcfLineFormatContentsDic["VAF_estimated:"+sN]=vafEstim
                vcfLineFormatContentsDic["VAF_calculated:"+sN]=perSampleVAF
            
            elif caller=="Strelka2":
                if "TAR" in vcfLineFormatHeader:
                    tar=int(vcfLineFormatContentsDic["TAR:"+sN].split(",")[0])
                    tir=int(vcfLineFormatContentsDic["TIR:"+sN].split(",")[0])
                    totalReads=tar+tir
                    altReads=tir
                else:
                    altSeq=line[vcf_header.index("ALT")]
                    totalReads=int(vcfLineFormatContentsDic["DP:"+sN])-int(vcfLineFormatContentsDic["FDP:"+sN])
                    altReads=int(vcfLineFormatContentsDic[altSeq+"U:"+sN].split(",")[0])
                
                if totalReads==0:
                    perSampleVAF="NA"
                else:
                    perSampleVAF=str(round(float(altReads/totalReads),3))
                vcfLineFormatContentsDic["VAF:"+sN]=perSampleVAF  
            
            elif caller=="VarScan":
                perSampleVAF=str(round(float(vcfLineFormatContentsDic["FREQ:"+sN].split("%")[0])/100,3))
                vcfLineFormatContentsDic["VAF:"+sN]=perSampleVAF
            elif caller=="MuTect":
                vcfLineFormatContentsDic["VAF:"+sN]=vcfLineFormatContentsDic["FA:"+sN]

        
        wholeLineAnnotDic["FORMAT"]+=[vcfLineFormatContentsDic]
        
        preCallerFormatCols+=vcfLineFormatContentsDic.keys()
        preCallerFormatCols=list(set(preCallerFormatCols))
        
    def columnMatcher(wholeDic, query, colList):
        for localDic in wholeDic[query]:
            for col in colList:
                if not col in localDic.keys():
                    localDic[col]="-"
        return wholeDic
        
    wholeLineAnnotDic=columnMatcher(wholeLineAnnotDic, "INFO", preCallerInfoCols)
    wholeLineAnnotDic=columnMatcher(wholeLineAnnotDic, "FORMAT", preCallerFormatCols)
    
    vcfLineNum=len(vcf_body)
    
    wholeLineAnnotDicList=[]
    for itr in range(vcfLineNum):
        reshapedDic={}
        for k1 in wholeLineAnnotDic.keys():
            localDic=wholeLineAnnotDic[k1][itr]
            for k2 in sorted(list(localDic.keys())):
                reshapedDic[k2]=localDic[k2]
        wholeLineAnnotDicList+=[reshapedDic]
    return(wholeLineAnnotDicList)

def callerParse(fn):
    callers=["Mutect2", "MuTect", "VarScan", "Strelka2", "HaplotypeCaller"]
    for caller in callers:
        if caller in fn or caller.lower() in fn:
            return(caller)
    metatxt, _, _ = decomposeVCF(fn)
    sources=[line.split("##source=")[1] for line in metatxt.split("\n") if line.startswith("##source=")]
    
    for caller in callers:
        if caller in sources or caller.lower() in sources:
            return(caller)
        if re.findall("[a-zA-Z]+", caller)[0] in sources or re.findall("[a-zA-Z]+", caller)[0] in sources:
            return(caller)
    return("")

def vcfcol2tsv(tempTab, finalVcf, callerArgs, queue=None):
    trimmedTab=tempTab.replace(".trimmed_tmp.tsv",".trimmed.tsv")
    wlad=vcfAddMoreInfo(finalVcf, callerArgs)
    print(len(wlad))
    with open(tempTab) as _tab:
        tabHeader=_tab.readline().rstrip().split("\t")
        tabBody=[line.rstrip().split("\t") for line in _tab.readlines()]
    with open(trimmedTab, "w") as wf:
        header2add=list(wlad[0].keys())
        breakPoint=tabHeader.index("Existing_variation")
        tabHeader_pre=tabHeader[:breakPoint]
        tabHeader_post=tabHeader[breakPoint:]
        newHeader=tabHeader_pre+header2add+tabHeader_post
            
        wf.write("\t".join(newHeader)+"\n")
        for itr in range(len(tabBody)):
            line2add=[wlad[itr][col] for col in header2add]
            lineBody_pre=tabBody[itr][:breakPoint]
            lineBody_post=tabBody[itr][breakPoint:]
            newLine=lineBody_pre+line2add+lineBody_post
            if "PUBMED" in newHeader:
                newLine[newHeader.index("PUBMED")]="_".join(newLine[newHeader.index("PUBMED")].split(","))
            wf.write("\t".join(newLine)+"\n")
    
    missenseTier(trimmedTab) ### Tier classification
    
    if queue:
        queue.put(trimmedTab)
        return None
    else:
        return(trimmedTab)
    
def multiAllelicGTcorrection(colonJoinedFormat, formatHeader):
    colonJoinedFormatSplitted=colonJoinedFormat.split(":")
    gt=colonJoinedFormatSplitted[formatHeader.index("GT")]
    if gt=="1/0":
        colonJoinedFormatSplitted[formatHeader.index("GT")]="0/1"
    elif len(gt.split("/"))>2 and list(set(gt.split("/")))==["0","1"]:
        colonJoinedFormatSplitted[formatHeader.index("GT")]="0/1"
    colonJoinedFormat=":".join(colonJoinedFormatSplitted)
    return(colonJoinedFormat)

def splitMultiAllelic(vcf_input, bcftools):
    vcfCheckedTmpFn=vcf_input.replace(".vcf", ".allele_tmp.vcf")
    vcfCheckedFn=vcf_input.replace(".vcf", ".allele.vcf")
    
    print("Splitting multi-allelic sites")
    command1="{0} norm -m -any {1} > {2}".format(bcftools, vcf_input, vcfCheckedTmpFn)
    os.system(command1)
    metaTxts, vcf_header, vcf_body=decomposeVCF(vcfCheckedTmpFn)
    metaTxts+='##INFO=<ID=bcftools_multiAllelic,Number=1,Type=String,Description="Multiallelic locus examined by BCFtools">\n'
    print("Counting multi-allelic sites")
    chromPosDic={}
    for line in vcf_body:
        chrom=line[vcf_header.index("#CHROM")]
        pos=line[vcf_header.index("POS")]
        chromPos=chrom+":"+pos
        if not chromPos in chromPosDic.keys():
            chromPosDic[chromPos]=1
        else:
            chromPosDic[chromPos]+=1
    multiAllelicSites=[k for k in chromPosDic.keys() if chromPosDic[k]>1]
    
    with open(vcfCheckedFn, "w") as vcfChecked:
        vcfChecked.write(metaTxts)
        vcfChecked.write("\t".join(vcf_header)+"\n")
        for line in vcf_body:
            chrom=line[vcf_header.index("#CHROM")]
            pos=line[vcf_header.index("POS")]
            chromPos=chrom+":"+pos
            if chromPos in multiAllelicSites:
                line[vcf_header.index("INFO")]="bcftools_multiAllelic=multi_allelic;"+line[vcf_header.index("INFO")]
                formatHeader=line[vcf_header.index("FORMAT")].split(":")
                for itr in range(len(vcf_header)-vcf_header.index("FORMAT")-1):
                    line[vcf_header.index("FORMAT")+itr+1]=multiAllelicGTcorrection(line[vcf_header.index("FORMAT")+itr+1], formatHeader)
            vcfChecked.write("\t".join(line)+"\n")
    print("Multi-allelic sites were marked")

    command2="rm {0}".format(vcfCheckedTmpFn)
    os.system(command2)
    return(vcfCheckedFn)

def mergingTsvs(tsvs):
    firstTsvName=tsvs[0].split(".")
    firstTsvName.remove("chunk_0000")
    mergedTsvName=".".join(firstTsvName)
    if len(tsvs)>1:
        tsvsDF=[pd.read_csv(t, sep='\t') for t in tsvs]
        concatTsvDF=pd.concat(tsvsDF, axis=0, ignore_index=True)
        concatTsvDF.to_csv(mergedTsvName, index=False, sep="\t")
    else:
        os.system("cp {0} {1}".format(tsvs[0], mergedTsvName))
    return mergedTsvName

def missenseTier(tsvFn):
    print("Classifying missense variants into damaging TIERs\n")
    aaDic={'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K', 'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N', 'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', 'Ter': '*', 'Ala': 'A', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}
    comSep=["VEP_canonical","Ensembl_transcriptid", "Ensembl_proteinid", "HGVSc_VEP", "HGVSp_VEP", "Interpro_domain"]
    comSep+=["SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "MPC_score", "FATHMM_pred", "MetaRNN_pred", "PROVEAN_pred"]
    oneDam=["CADD_phred", "CADD_raw_rankscore", "REVEL_score", "MetaSVM_pred", "MetaLR_pred", "BayesDel_noAF_pred", "Eigen-raw_coding_rankscore", "Eigen-PC-raw_coding_rankscore", "M-CAP_score", "DANN_rankscore", "LRT_pred", "fathmm-MKL_coding_pred", "VEST4_rankscore", "PrimateAI_pred"]
    mutTas=["MutationTaster_AAE", "MutationTaster_pred"]
    tmpFn=tsvFn+"_tmp"
    os.system("mv {0} {1}".format(tsvFn,tmpFn))
    canonicalFn=tsvFn.replace(".tsv", ".canonical.tsv")
    trashFn=tsvFn.replace(".tsv", ".dontUseThisRemainders.tsv")
    
    wf1=open(tsvFn, "w")
    wf2=open(canonicalFn, "w")
    wf3=open(trashFn, "w")
    with open(tmpFn) as tmp:
        header=tmp.readline().rstrip().split("\t")
        newHeader="\t".join(header+["MissenseTier"])+"\n"
        wf1.write(newHeader)
        wf2.write(newHeader)
        wf3.write(newHeader)
        lineCount=0
        for line in tmp:
            lineCount+=1
            x=line.rstrip().split("\t")
            level="-"
            damDic={}
            if x[header.index("Consequence")].startswith("missense_variant") and x[header.index("BIOTYPE")]=="protein_coding":
                canon=x[header.index("VEP_canonical")].split(",")
                if canon.count("YES")==1:
                    pickN=canon.index("YES")
                else:
                    ensts=x[header.index("Ensembl_transcriptid")].split(",")
                    try:
                        enstPick=x[header.index("HGVSc")].split(":")[0].split(".")[0]
                        pickN=ensts.index(enstPick)
                    except:
                        if canon.count("YES")>1:
                            pickN=canon.index("YES")
                        else:
                            wf3.write("\t".join(x+[level])+"\n")
                            continue
                for c in comSep:
                    damDic[c]=x[header.index(c)].split(",")[pickN]
                for c in oneDam:
                    damDic[c]=x[header.index(c)]
                    
                hgvsp_vep=damDic["HGVSp_VEP"]
                if "p." in hgvsp_vep:
                    aac3=damDic["HGVSp_VEP"].split("p.")[1]
                    aas_12=[aaDic[aa] for aa in re.findall("[a-zA-Z]+", aac3)]
                    aas_pos=re.findall("[0-9]+", aac3)[0]
                    aac1=aas_12[0]+aas_pos+aas_12[1]
                else:
                    aac1="-"
                
                mutTasAAE=x[header.index("MutationTaster_AAE")].split(",")
                if aac1 in mutTasAAE:
                    mutTasIndex=mutTasAAE.index(aac1)
                    damDic["MutationTaster_pred"]=x[header.index("MutationTaster_pred")].split(",")[mutTasIndex]
                    damDic["MutationTaster_AAE"]=mutTasAAE[mutTasIndex]

                else:
                    damDic["MutationTaster_pred"]="-"
                    damDic["MutationTaster_AAE"]="-"
                try:
                    if damDic["MutationTaster_pred"] in ["D","A"] and damDic["LRT_pred"]=="D" and damDic["SIFT_pred"]=="D" and damDic["Polyphen2_HDIV_pred"]=="D" and damDic["Polyphen2_HVAR_pred"]=="D":
                        level="2"
                        try:
                            if float(damDic["M-CAP_score"])>0.025 and damDic["MetaLR_pred"]=="D" and damDic["MetaSVM_pred"]=="D" and damDic["PROVEAN_pred"]=="D" and damDic["fathmm-MKL_coding_pred"]=="D" and damDic["FATHMM_pred"]=="D":
                                level="3"
                                try:
                                    if float(damDic["VEST4_rankscore"])>0.9 and float(damDic["CADD_raw_rankscore"])>0.9 and float(damDic["DANN_rankscore"])>0.9 and float(damDic["Eigen-raw_coding_rankscore"])>0.9 and float(damDic["Eigen-PC-raw_coding_rankscore"])>0.9:
                                        level="4"
                                except:pass
                        except:pass
                    elif damDic["MutationTaster_pred"] not in ["D","A"] and damDic["LRT_pred"]!="D" and damDic["SIFT_pred"]!="D" and damDic["Polyphen2_HDIV_pred"]!="D" and damDic["Polyphen2_HVAR_pred"]!="D":
                        level="0"
                    else:
                        level="1"
                    
                except:pass
                
            wf1.write("\t".join(x+[level])+"\n")

            wf2_toWrite=[]
            for h in header:
                if h in damDic.keys():
                    wf2_toWrite+=[damDic[h]]
                else:
                    wf2_toWrite+=[x[header.index(h)]]
            wf2.write("\t".join(wf2_toWrite+[level])+"\n")
    wf1.close()
    wf2.close()
    wf3.close()
    
def run(vcf_input, vepPath, ref_genome, assembly, annotDBdir, customAnnot, sortConfig, bcftools, vcfProd, cores=1, callerArgs=""):
    if vcf_input.endswith(".vcf.gz"):
        print("decompressing .vcf.gz > .vcf")
        os.system("pigz -p {0} -k -d {1}".format(cores, vcf_input))
        vcf_input=vcf_input.replace(".vcf.gz", ".vcf")
    
    if not callerArgs:
        callerArgs=callerParse(vcf_input)
    
    if callerArgs:
        print("\n{0} was detected as the variation caller\n".format(callerArgs))
    else:
        print("No variation caller was detected\n")    
    
    ### Get default INFO annotations 
    defMeta=[]
    with open(vcf_input) as initVcf:
        for line in initVcf:
            if line.startswith("#"):
                defMeta.append(line)
            else:
                break
    
    defColnames=[j.split("ID=")[1].split(",")[0] for j in [k for k in defMeta if k.startswith("##INFO")]]
    defFormat=["QUAL", "FILTER"]+[j.split("ID=")[1].split(",")[0] for j in [k for k in defMeta if k.startswith("##FORMAT")]]+["VAF", "VAF_estimated", "VAF_calculated"]
    sampleNames=[k for k in defMeta if k.startswith("#CHROM")][0].rstrip().split("\t")[9:]
    for f in defFormat:
        for sn in sampleNames:
            defColnames.append(f+":"+sn)    
    
    ### Multi-allelicity Checking
    effCores=max(1, cores-1)
    print("{0} processes will be used".format(effCores))
    
    splittedVcfFn=splitVCF(vcf_input, effCores)
    splittedVcfCheckedFn=[splitMultiAllelic(vcf_indiv, bcftools) for vcf_indiv in splittedVcfFn]
    
    print("\n1. VCFs were splitted, while multi-allelicity resolved\n".format(effCores))

    ### BCFtools annotation of gene-level (dbNSFP_gene_4.2 + COSMIC_CMC + GTEx_v8) and custom ones.
    bcfVcfQ=multiprocessing.Queue()
    procsBcftoolsVcf=[multiprocessing.Process(target=bcftoolsAnnot, args=(vcfFn, bcftools, annotDBdir, customAnnot,bcfVcfQ)) for vcfFn in splittedVcfCheckedFn]
    for p in procsBcftoolsVcf: p.start()
    for p in procsBcftoolsVcf: p.join()
    finalBcftoolsVcfs=sorted([bcfVcfQ.get() for p in procsBcftoolsVcf])
        
    print("\n2. bcftools annotate: gene-level and custom annotations DONE\n".format(effCores))

    
    ### VEP annotations (TAB)
    vepTabQ=multiprocessing.Queue()
    procsVepTab=[multiprocessing.Process(target=vepAnnot, args=(vcfFn, "tab", vepPath, assembly, ref_genome, annotDBdir, "vep",customAnnot,vepTabQ)) for vcfFn in finalBcftoolsVcfs]
    for p in procsVepTab: p.start()
    for p in procsVepTab: p.join()
    finalVepTabs=sorted([vepTabQ.get() for p in procsVepTab])
    vepTrimmedTsvs=[t.replace(".tab",".trimmed_tmp.tsv") for t in finalVepTabs]
    for itr in range(len(finalVepTabs)):
        os.system("cat {0} | grep -v '##' > {1}".format(finalVepTabs[itr], vepTrimmedTsvs[itr]))
        
    ### VEP annotations (VCF)
    if vcfProd:
        vepVcfQ=multiprocessing.Queue()
        procsVepVcf=[multiprocessing.Process(target=vepAnnot, args=(vcfFn, "vcf", vepPath, assembly, ref_genome, annotDBdir, "vep",customAnnot,vepVcfQ)) for vcfFn in finalBcftoolsVcfs]
        for p in procsVepVcf: p.start()
        for p in procsVepVcf: p.join()
        finalVepVcfs=sorted([vepVcfQ.get() for p in procsVepVcf])

        firstVepVcfName=finalVepVcfs[0].split(".")
        firstVepVcfName.remove("chunk_0000")
        mergedVepVcfName=".".join(firstVepVcfName)

        if len(finalVepVcfs)>1:
            os.system("{0} concat {1} > {2}".format(bcftools, " ".join(finalVepVcfs), mergedVepVcfName))
        else:
            os.system("cp {0} {1}".format(finalVepVcfs[0], mergedBcftoolsVcfName))
        os.system("rm {0}".format(" ".join(finalVepVcfs)))
        
    print("\n2. VEP annotations with dbNSFP, SpliceAI, LOFTEE DONE\n".format(effCores))    
    
    print("Merging annotated VCF to annotated TSV\n")
    
    vcfCol2TsvQ=multiprocessing.Queue()
    procsVcfCol2Tsv=[multiprocessing.Process(target=vcfcol2tsv, args=(vepTrimmedTsvs[itr], finalBcftoolsVcfs[itr], callerArgs, vcfCol2TsvQ)) for itr in range(len(finalBcftoolsVcfs))]
    for p in procsVcfCol2Tsv: p.start()
    for p in procsVcfCol2Tsv: p.join()
        
    finalVepTsvNames=sorted([vcfCol2TsvQ.get() for p in procsVcfCol2Tsv])
    canonicalTsvNames=[tsv.replace(".tsv", ".canonical.tsv") for tsv in finalVepTsvNames]
    trashTsvNames=[tsv.replace(".tsv", ".dontUseThisRemainders.tsv") for tsv in finalVepTsvNames]

    preTierTsvNames=[tsv+"_tmp" for tsv in finalVepTsvNames]
    
    finalVepTsvName=mergingTsvs(finalVepTsvNames)
    canonicalTsvName=mergingTsvs(canonicalTsvNames)
    trashTsvName=mergingTsvs(trashTsvNames)
    
    selectedVepTsvName=finalVepTsvName.replace(".trimmed",".selected")
    selectedCanonicalTsvName=canonicalTsvName.replace(".trimmed",".selected")
    selectedTrashTsvName=trashTsvName.replace(".trimmed",".selected")

    with open(annotDBdir+"/"+sortConfig) as sortConfigFile:
        sortConfigCols=[line.rstrip() for line in sortConfigFile.readlines()]+defColnames
    with open(finalVepTsvName) as finalVepTsv:
        existingColNames=finalVepTsv.readline().rstrip().split("\t")

    selectedConfigCols=[c for c in sortConfigCols if c in existingColNames]
       
    finalVepTsvDf=pd.read_csv(finalVepTsvName, sep='\t')
    selectedVepTsvDf=finalVepTsvDf[selectedConfigCols]
    selectedVepTsvDf.to_csv(selectedVepTsvName, index=False, sep="\t")

    canonicalTsvDf=pd.read_csv(canonicalTsvName, sep='\t')
    selectedCanonicalTsvDf=canonicalTsvDf[selectedConfigCols]
    selectedCanonicalTsvDf.to_csv(selectedCanonicalTsvName, index=False, sep="\t")

    trashTsvDf=pd.read_csv(trashTsvName, sep='\t')
    selectedTrashTsvDf=trashTsvDf[selectedConfigCols]
    selectedTrashTsvDf.to_csv(selectedTrashTsvName, index=False, sep="\t")
    
    os.system("rm {0}".format(" ".join(splittedVcfFn)))
    os.system("rm {0}".format(" ".join(splittedVcfCheckedFn)))
    os.system("rm {0}".format(" ".join(finalVepTabs)))
    os.system("rm {0}".format(" ".join(finalBcftoolsVcfs)))
    os.system("rm {0}".format(" ".join(vepTrimmedTsvs)))
    
    os.system("rm {0}".format(" ".join(finalVepTsvNames)))
    os.system("rm {0}".format(" ".join(canonicalTsvNames)))
    os.system("rm {0}".format(" ".join(trashTsvNames)))
    os.system("rm {0}".format(" ".join(preTierTsvNames)))
    
    os.system("rm {0}".format(finalVepTsvName))
    os.system("rm {0}".format(canonicalTsvName))
    os.system("rm {0}".format(trashTsvName))
    print("Finished")
    
    
def run2(vcf_input, vepPath, ref_genome, assembly, annotDBdir, customAnnot, sortConfig, bcftools, vcfProd, cores=1, callerArgs=""):
    if vcf_input.endswith(".vcf.gz"):
        print("decompressing .vcf.gz > .vcf")
        os.system("pigz -p {0} -k -d {1}".format(cores, vcf_input))
        vcf_input=vcf_input.replace(".vcf.gz", ".vcf")
    
    if not callerArgs:
        callerArgs=callerParse(vcf_input)
    
    if callerArgs:
        print("\n{0} was detected as the variation caller\n".format(callerArgs))
    else:
        print("No variation caller was detected\n")    
    effCores=max(1, cores-1)
    splittedVcfFn=splitVCF(vcf_input, effCores)
    vepTrimmedTsvs=[i.replace(".vcf", ".allele.bcftools.vep.trimmed_tmp.tsv") for i in splittedVcfFn]
    finalBcftoolsVcfs=[i.replace(".vcf", ".allele.bcftools.vcf") for i in splittedVcfFn]
    
    print("Merging annotated VCF to annotated TSV\n")
   
    
    
    vcfCol2TsvQ=multiprocessing.Queue()
    procsVcfCol2Tsv=[multiprocessing.Process(target=vcfcol2tsv, args=(vepTrimmedTsvs[itr], finalBcftoolsVcfs[itr], callerArgs, vcfCol2TsvQ)) for itr in range(len(finalBcftoolsVcfs))]
    for p in procsVcfCol2Tsv: p.start()
    for p in procsVcfCol2Tsv: p.join()
        
    finalVepTsvNames=sorted([vcfCol2TsvQ.get() for p in procsVcfCol2Tsv])
    canonicalTsvNames=[tsv.replace(".tsv", ".canonical.tsv") for tsv in finalVepTsvNames]
    trashTsvNames=[tsv.replace(".tsv", ".dontUseThisRemainders.tsv") for tsv in finalVepTsvNames]

    preTierTsvNames=[tsv+"_tmp" for tsv in finalVepTsvNames]
    
    finalVepTsvName=mergingTsvs(finalVepTsvNames)
    canonicalTsvName=mergingTsvs(canonicalTsvNames)
    trashTsvName=mergingTsvs(trashTsvNames)
    
    selectedVepTsvName=finalVepTsvName.replace(".trimmed",".selected")
    selectedCanonicalTsvName=canonicalTsvName.replace(".trimmed",".selected")
    selectedTrashTsvName=trashTsvName.replace(".trimmed",".selected")

    with open(annotDBdir+"/"+sortConfig) as sortConfigFile:
        sortConfigCols=[line.rstrip() for line in sortConfigFile.readlines()]+defColnames
    with open(finalVepTsvName) as finalVepTsv:
        existingColNames=finalVepTsv.readline().rstrip().split("\t")

    selectedConfigCols=[c for c in sortConfigCols if c in existingColNames]
       
    finalVepTsvDf=pd.read_csv(finalVepTsvName, sep='\t')
    selectedVepTsvDf=finalVepTsvDf[selectedConfigCols]
    selectedVepTsvDf.to_csv(selectedVepTsvName, index=False, sep="\t")

    canonicalTsvDf=pd.read_csv(canonicalTsvName, sep='\t')
    selectedCanonicalTsvDf=canonicalTsvDf[selectedConfigCols]
    selectedCanonicalTsvDf.to_csv(selectedCanonicalTsvName, index=False, sep="\t")

    trashTsvDf=pd.read_csv(trashTsvName, sep='\t')
    selectedTrashTsvDf=trashTsvDf[selectedConfigCols]
    selectedTrashTsvDf.to_csv(selectedTrashTsvName, index=False, sep="\t")
    
    os.system("rm {0}".format(" ".join(splittedVcfFn)))
    os.system("rm {0}".format(" ".join(splittedVcfCheckedFn)))
    os.system("rm {0}".format(" ".join(finalVepTabs)))
    os.system("rm {0}".format(" ".join(finalBcftoolsVcfs)))
    os.system("rm {0}".format(" ".join(vepTrimmedTsvs)))
    
    os.system("rm {0}".format(" ".join(finalVepTsvNames)))
    os.system("rm {0}".format(" ".join(canonicalTsvNames)))
    os.system("rm {0}".format(" ".join(trashTsvNames)))
    os.system("rm {0}".format(" ".join(preTierTsvNames)))
    
    os.system("rm {0}".format(finalVepTsvName))
    os.system("rm {0}".format(canonicalTsvName))
    os.system("rm {0}".format(trashTsvName))
    print("Finished")
    
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='Input vcf(.gz)', type=str, metavar='')
    parser.add_argument('-r', '--refFa', required=True, help='GRCh38 or GRCh37', type=str, metavar='') ## GRCh37, GRCh38
    parser.add_argument('-v', '--variantCaller', required=False, default='',  help="HaplotypeCaller, Mutect2, MuTect, VarScan, Strelka2" ,type=str, metavar='')
    parser.add_argument('-c', '--cores', required=False, type=str, help='multiple cores for parallelization', default='1', metavar='')
    parser.add_argument('-f', '--vcf', required=False, type=str, help='Final VCF output (y/n) default="n"; only TSV', default='n', metavar='')
    parser.add_argument('-m', '--mode', required=False, type=str, help='Custom annotation source mode: a(all=b+c+g), b(brain), c(cancer), g(general), n(none)', default='a', metavar='')

    args=parser.parse_args()
        
    refFa=args.refFa
    
    refPath="/home/public/ref/"
    fasta38="Homo_sapiens_assembly38.fa"
    #fasta37="Homo.Sapiens.GRCh37.75.primary_assembly.full.fa"
    fasta37="Homo.Sapiens.GRCh37.75.primary_assembly.full.vepAnnot.fa.gz"
    
    fastaPath=refPath+"00_Ref_{0}/{0}_fasta/".format(refFa)
    annotDBdir=refPath+"00_Ref_{0}/annot".format(refFa)

    if args.vcf == "y":
        vcfProd=True
    else:
        vcfProd=False

    if refFa == "GRCh38":
        ref_genome=fastaPath+fasta38
    elif refFa == "GRCh37":
        ref_genome=fastaPath+fasta37
        
    vepPath="/home/public/tools/ensembl-vep"
    bcftools="/home/public/tools/bcftools/bin/bcftools"
    
    sortConfig="columnSort.config"
    
    if not args.mode in ["a", "b", "c", "g", "n"]:
        print("--mode should be [a,b,c,g,n]")
        mod=a
    else:
        mod=args.mode
    
    ######### CUSTOM Annotations
    ### multipleCustomAnnot: Multiple custom annot with dictionary structures as below
    ### singleCustomAnnot: Singleton annotations (overlapping region)
    ### allelicCustomAnnot: Singleton Allelic annotations (exact variation)

    multipleCustomAnnot, singleCustomAnnot, allelicCustomAnnot={}, {}, {}
    
    if mod in ["a", "b", "c", "g", "n"]:
        multipleCustomAnnot["geneLevel"]={"file":"genes.dbNSFP4.2_gtexV8_cosmic.{0}.bed.gz".format(refFa),\
                                      "columnsFn":"genes.dbNSFP.columns.txt",\
                                     "infoFn":"genes.dbNSFP.vcfheaders.txt"}    
    if mod in ["a", "b", "c", "g"]:
        
        singleCustomAnnot["Super_enhancers[dbSUPER]"]="dbSUPER_{0}.gnomAD_SV.bed.gz".format(refFa)
        singleCustomAnnot["DnaseI_HS[ENCODE.V3]"]="encode.DHS_annot.V3.gnomAD_SV.{0}.bed.gz".format(refFa)
        singleCustomAnnot["TF_BindingSites[ENCODE.V3]"]="encode.TFBS_annot.V3.gnomAD_SV.{0}.trimmed.bed.gz".format(refFa)
        
        singleCustomAnnot["EnhancerAtlas_2.0"]="human_enhancerAtlas.{0}.bed.gz".format(refFa)
        singleCustomAnnot["VISTA_validated_enhancer"]="vista_enhancer.{0}.bed.gz".format(refFa)
        singleCustomAnnot["TAD_boundaries[IMR90]"]="TAD_boundaries.IMR90.gnomAD_SV.{0}.trimmed.bed.gz".format(refFa)
        
        singleCustomAnnot["Human_accelerated_region"]="Human_Accelerated_Region.{0}.sorted.bed.gz".format(refFa)
        singleCustomAnnot["Ultraconserved_non-coding"]="UCNE.{0}.bed.gz".format(refFa)
        
        allelicCustomAnnot["GWAS_Catalog_SNPs[allelic]"]="gwasCatalog.{0}.allelic.bed.gz".format(refFa)
        allelicCustomAnnot["Splicing_dbscSNV[allelic]"]="dbscSNV1.1_splicing.{0}.bed.gz".format(refFa)
    
    if mod in ["a", "c"]:
        multipleCustomAnnot["COSMIC_variants"]={"file":"cmc_variants.{0}.sorted.bed.gz".format(refFa),\
                                  "columnsFn":"cmc_variants.columns.txt",\
                                 "infoFn":"cmc_variants.vcfheaders.txt"}
    if mod in ["a", "b"]:
        singleCustomAnnot["ATAC-seq_Peaks[BOCA_Brain]"]="brain_ATAC-seq_consensusPeaks_boca.{0}.bed.gz".format(refFa)
        singleCustomAnnot["Bivalent_enhancer[RoadmapEpigenomics_Brain]"]="brain.bivalent_enhancer.chromHMM.{0}.overlapped.bed.gz".format(refFa)
        singleCustomAnnot["Enhancer[RoadmapEpigenomics_Brain]"]="brain.enhancer.chromHMM.{0}.overlapped.bed.gz".format(refFa)
        singleCustomAnnot["Genic_enhancer[RoadmapEpigenomics_Brain]"]="brain.genic_enhancer.chromHMM.{0}.overlapped.bed.gz".format(refFa)
        singleCustomAnnot["Repressed_polycomb[RoadmapEpigenomics_Brain]"]="brain.repressed_polycomb.chromHMM.{0}.overlapped.bed.gz".format(refFa)

        singleCustomAnnot["TAR[PsychENCODE]"]="psychEncode_active_TAR.{0}.bed.gz".format(refFa)
        singleCustomAnnot["H3K27ac_peaks[PsychENCODE]"]="psychEncode_H3K27ac_peaks.{0}.bed.gz".format(refFa)
        singleCustomAnnot["Enhancer_PFC[PsychENCODE]"]="psychEncode_PFC_enhancer.{0}.bed.gz".format(refFa)
        singleCustomAnnot["TAD_boundaries_DLPFC[PsychENCODE]"]="psychEncode_TAD_Boundary_DLPFC.{0}.trimmed.bed.gz".format(refFa)
        
        allelicCustomAnnot["cis-eQTL[PsychENCODE_allelic]"]="psychEncode_cis-eQTL.FDR_sig.{0}.allelic.bed.gz".format(refFa)
    

    
    ### 



    customAnnot={"multiple":multipleCustomAnnot, "single":singleCustomAnnot, "allelic": allelicCustomAnnot}

    run(args.input, vepPath, ref_genome, refFa, annotDBdir, customAnnot, sortConfig, bcftools, vcfProd, int(args.cores), args.variantCaller)
    #testRun(args.input, vepPath, ref_genome, refFa, annotDBdir, customAnnot, sortConfig, bcftools, vcfProd, int(args.cores), args.variantCaller)

    #onlyMerge(args.input, vepPath, ref_genome, refFa, annotDBdir, customAnnot, sortConfig, bcftools, vcfProd, int(args.cores), args.variantCaller)
