#!/home/dcha/python3

import os,sys
import re
from ete3 import Tree
import copy
import math

def muts_dic_gen(mutsFn):
    muts_dic={"Blood":[]}
    with open(mutsFn) as muts:
        for k in muts:
            x=k.rstrip().split("\t")
            ct=x[0]
            mutsPerCell=x[1].split(",")
            muts_dic[ct]=mutsPerCell
    return(muts_dic)

def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return(text)

def nhx_conversion(path, treeFn, muts_dic):
    with open(treeFn) as _tree:
        treeString=_tree.readline().rstrip()
    specDic={}
    wholeDic={}
    wholeDic["ROOT"]=[]
    specDic["ROOT"]=[]
    origName=["ROOT"]
    unionMutsLi=[]
    for k in muts_dic.keys():
        unionMutsLi+=muts_dic[k]
    unionMuts=set(unionMutsLi)

    t=Tree(treeString,format=1)
    wholeDic, specDic, origName=traversing(t, origName, specDic, wholeDic, muts_dic, unionMuts)
  
    #print(origName)
    print(specDic.keys())
    print([len(specDic[k]) for k in specDic.keys()])
    print([len(wholeDic[k]) for k in specDic.keys()])

    with open(treeFn+".nhx", "w") as wf:
        _nhx=copy.copy(treeString)
        rep=dict(zip(origName, specDic.keys()))
        rep[",Blood"]=""
        rep["ROOT"]="Blood:0.000"
        nhx=replace_all(_nhx, rep)
        for k in specDic.keys():
            num=len(specDic[k])
            nhx_format="{0}:{1}[&&NHX:mut={2}]".format(k,round(math.log10(num+1),3), num)
            print(nhx_format)
            nhx = nhx.replace(k, nhx_format)
        print(nhx)
        wf.write(nhx)

    with open(path+"branch_mutations.txt", "w") as wf:
        for k in specDic.keys():
            toWrite=k+"\t"+",".join(specDic[k])+"\n"
            wf.write(toWrite)

def nodeRenamer(t):
    if re.findall('^[0-9]+$', t.name):
        return("seq"+t.name)
    else:   return(t.name)

def traversing(t,origName, specDic, wholeDic, muts_dic, unionMuts): ## Recurrent function
    root=nodeRenamer(t)

    for c in t.children:
        cn=nodeRenamer(c)
        if cn!="Blood":
            origName.append(c.name)
            leaves=[leaf.name for leaf in c]
            mut_set=copy.copy(unionMuts)

            for l in leaves:
                mut_set=mut_set.intersection(set(muts_dic[l]))
            mut_li=list(mut_set)
            wholeDic[cn]=mut_li
            specDic[cn]=[k for k in mut_li if k not in wholeDic[root]]

            wholeDic, specDic, origName =traversing(c,origName, specDic, wholeDic, muts_dic, unionMuts)
    return(wholeDic, specDic, origName)

def branchVCF(path):

    branchFn=path+"branch_mutations.txt"
    vcfFn=path+"vcfMerge.rm.vcf"
    idFn=path+"vcfMerge.rm.identifiers.txt"

    brdic={}
    with open(branchFn) as brf:
        _=brf.readline()
        br=[k.rstrip().split("\t") for k in brf.readlines()]
    brdic={k[0]:k[1].split(",") for k in br}
    headerLi=[]
    vcfBody=[]
    
    with open(idFn) as idf:
        idents=[l.rstrip() for l in idf.readlines()]
    identsDic={idents[itr]:itr for itr in range(len(idents))}
    print(identsDic)

    with open(vcfFn) as vcf:
        for l in vcf:
            if l.startswith("#"):
                headerLi+=[l]
            else:
                vcfBody+=[l]

    header="".join(headerLi)

    for k in brdic.keys():
        wfn1=path+"{0}.tmpBranchMuts.vcf".format(k)
        with open(wfn1, "w") as wf:
            for mut in brdic[k]:
                try:
                    wf.write(vcfBody[identsDic[mut]])
                except: pass
    
        wfn2=path+"{0}.branchMuts.vcf".format(k)
        with open(wfn2, "w") as wf:
            wf.write(header)

        os.system("sort -k1,1V -k2,2n {0}{1}.tmpBranchMuts.vcf >> {0}{1}.branchMuts.vcf".format(path, k))
        os.system("rm -f {0}{1}.tmpBranchMuts.vcf".format(path, k))
    
def main(path, treeFn, mutsFn):
    muts_dic=muts_dic_gen(mutsFn)
    nhx_conversion(path, treeFn, muts_dic)
    
if __name__ == "__main__":
    path=sys.argv[1]
    treeFn=path+"raxml/RAxML_nodeLabelledRootedTree.anc"
    mutsFn=path+"cellTypes_mutations.txt"
    os.system(" ".join(["bash",
            "/home/jsh/02.Tools/cellTypeGenome_codes/run_raxml_rooted.sh",
            path]))
    main(path, treeFn, mutsFn)
    branchVCF(path)
