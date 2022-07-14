#!/home/dcha/python3
import os,sys

def main(DIR):
    intervalLists=[k for k in os.listdir(DIR) if k.endswith(".interval_list")]
    for fn in intervalLists:
        wfn=fn.replace(".interval_list", ".bed")
        with open (DIR+"/"+fn) as _rf:
            rf=["\t".join(k.split("\t")[:3]) for k in _rf.readlines() if not k.startswith("@")]
        with open(DIR+"/"+wfn, "w") as wf:
            for k in rf:
                wf.write(k+"\n")
    return None

if __name__ == "__main__":
    DIR=sys.argv[1]
    main(DIR)
