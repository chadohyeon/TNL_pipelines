#!/usr/bin/python3

### Updated 2021.04.13.

import os, sys
import argparse
import pickle

config = "/home/public/ref/metaTexts/mouseHumanOrthologs.txt"
humanGencodeGenes = "/home/public/ref/metaTexts/gencode.v35.annotation_geneListOnly.txt"
mouseGencodeGenes = "/home/public/ref/metaTexts/gencode.m23_geneListOnly.txt"

def buildDB():
    hom2mus={}
    mus2hom={}	
	
    geneDic={}
    uniqDic={}
    with open(config) as _rf:
        rf = [k.strip().split("\t")[0:4] for k in _rf.readlines()]
    for k in rf:
        if k[1]=="mouse, laboratory":
            k[1]="mouse"

    ### Matching by MGI numbers, which are shared across orthologs
    for k in rf:
        if not k[0] in geneDic.keys():
            geneDic[k[0]]=[[k[1],k[3]]]
        else:
            geneDic[k[0]]+=[[k[1],k[3]]]
    for k in geneDic.keys():
        unique=list(set([j[0] for j in geneDic[k]]))
        if sorted(unique) == ["human", "mouse"]:
            uniqDic[k] = geneDic[k]
    for k in rf:
        if "human" in k:
            homGene=k[3]
            num=k[0]
            if not num in uniqDic.keys():
                continue
            hom2mus[homGene]=[j[1] for j in uniqDic[num] if j[0]=="mouse"]
        elif "mouse" in k:
            musGene=k[3]
            num=k[0]
            if not num in uniqDic.keys():
                continue
            mus2hom[musGene]=[j[1] for j in uniqDic[num] if j[0]=="human"]
  

    ### Matching by gene nomenclature structure
    with open(mouseGencodeGenes) as _mgg:
        mgg=_mgg.readlines()
    with open(humanGencodeGenes) as _hgg:
        hgg=_hgg.readlines()

    mouseAllGenes=[k.rstrip() for k in mgg]
    humanAllGenes=[k.rstrip() for k in hgg]
    
    for mg in mouseAllGenes:
        if not mg in mus2hom.keys():
            if mg.upper() in humanAllGenes:
                mus2hom[mg]=[mg.upper()]
    for hg in humanAllGenes:
        if not hg in hom2mus.keys():
            hgSplit=list(hg)
            hg2mg="".join([hgSplit[0]]+[k.lower() for k in hgSplit[1:]])
            if hg2mg in mouseAllGenes:
                hom2mus[hg]=[hg2mg]
    with open('/home/public/ref/metaTexts/hom2mus.pickle', 'wb') as wf:
        pickle.dump(hom2mus, wf)
    with open('/home/public/ref/metaTexts/mus2hom.pickle', 'wb') as wf:
        pickle.dump(mus2hom, wf)

    return hom2mus,mus2hom

def query(fn, species="human", ncol=2, load=False):
    if load:
        with open('/home/public/ref/metaTexts/hom2mus.pickle', 'rb') as rf:
            hom2mus=pickle.load(rf)
        with open('/home/public/ref/metaTexts/mus2hom.pickle', 'rb') as rf:
            mus2hom=pickle.load(rf)
    else:
        hom2mus, mus2hom = buildDB()
    with open(fn) as _rf:
        qGenes=[g.strip() for g in _rf.readlines()]
    wfn = toWrite(fn, qGenes, hom2mus, mus2hom, species, ncol)
    return wfn

def toWrite(fn, qGenes, hom2mus, mus2hom, species, ncol):
    if species == "human":
        header = "human\tmouse\n"
        dic = hom2mus
        prefix="hom2mus."

    elif species == "mouse":
        header = "mouse\thuman\n"
        dic = mus2hom
        prefix="mus2hom."
   
    wfn=prefix+fn
    wf = open(wfn,'w')
    if ncol==2:
        wf.write(header)
    for g in qGenes:
        if not g in dic.keys():
            rgs=["NA"]
        else:
            rgs=dic[g]
        if len(rgs)==1:
            if ncol==2:
                wf.write(g+"\t"+rgs[0]+"\n")
            else:
                if not rgs[0]=="NA":
                    wf.write(rgs[0]+"\n")
        else:
            for rg in rgs:
                if ncol==2:
                    wf.write(g+"\t"+rg+"\n")
                else:
                    if not rg=="NA":
                        wf.write(rg+"\n")
    wf.close()

    return wfn

def countMultimap(fn):
    suffix="."+fn.split(".")[-1]
    wfn = fn.replace(suffix, ".count"+suffix)
    wf = open(wfn, 'w')
    with open(fn) as rf:
        header=rf.readline()
        body=rf.readlines()
    
    wf.write(header.replace("\n", "\tdescription\n"))
    [col1,col2] = header.rstrip().split("\t")
    li1 = [k.rstrip().split("\t")[0] for k in body]
    li2 = [k.rstrip().split("\t")[1] for k in body]

    dic1={}
    dic2={}
    for l in li1:
        if l not in dic1.keys():
            dic1[l]=1
        else:   dic1[l]+=1
    for l in li2:
        if l not in dic2.keys():
            dic2[l]=1
        else:   dic2[l]+=1

    for k in body:
        msg=[]
        [v1, v2] = k.rstrip().split("\t")
        if dic1[v1]>1:
            msg+=["multi_{0}_genes".format(col1)]
        if v2!="NA":
            if dic2[v2]>1:
                msg+=["multi_{0}_genes".format(col2)]
        
        if msg==[]:
            add="."
        else:
            add=",".join(msg)
        wf.write(k.rstrip()+"\t"+add+"\n")
    wf.close()



def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help="gene lists", type=str, metavar='')
    parser.add_argument('-s', '--species',required=True, type=str, help='human/mouse', metavar='')
    parser.add_argument('-c', '--colNum', required=False, type=str, help='output column numbers', default='2', metavar='')
    args=parser.parse_args()
    parser.print_help()

    wfn=query(args.input, args.species, int(args.colNum), load=True)
    if int(args.colNum)==2:
        countMultimap(wfn)
        os.system("rm {0}".format(wfn))

if __name__ == "__main__":
    #buildDB()
    main()
