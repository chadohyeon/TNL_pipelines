#!/home/dcha/python3
import os,sys
def main(refInterval, commaChroms, outFile):
    chroms=commaChroms.split(",")
    with open(refInterval) as _ref:
        ref=_ref.readlines()
    with open(outFile, "w") as wf:
        for l in ref:
            if l.startswith("@"):   wf.write(l)
            else:
                if l.split("\t")[0] in chroms:
                    wf.write(l)

if __name__=="__main__":
    refInterval=sys.argv[1]
    commaChroms=sys.argv[2]
    outFile=sys.argv[3]
    main(refInterval, commaChroms, outFile)

