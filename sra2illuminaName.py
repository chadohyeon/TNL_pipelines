#!/home/dcha/python3
import os, sys

def main():
    fqgzs=[k for k in os.listdir(os.getcwd()) if ".fq" in k]

    for k in fqgzs:
        newPrefix=k.split("_")[0]+"_S1_L001"
        newSuffix="_R"+k.split("_")[1].replace(".fq","_001.fastq")
        command="mv {0} {1}{2}".format(k, newPrefix,newSuffix)
        print(command)
        os.system(command)
if __name__ == "__main__":
    main()