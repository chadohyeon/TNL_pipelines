#!/home/dcha/python3

import os, argparse

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-f', '--file', required=False, help="scriptFiles", type=str, default="STDIN", metavar='')
    parser.add_argument('-l', '--line', required=False, help='line Command with quote, e.g., "samtools view -h ABC.bam"', type=str, default="Null", metavar='')
    parser.add_argument('-c', '--cores', required=False, type=str, help='multiple cores for parallelization', default='1', metavar='')
    parser.add_argument('-q', '--queue', required=False, type=str, help='PBS queue', default='long', metavar='')
    args=parser.parse_args()
    
    if args.file=="STDIN" and args.line=="Null":
        print("Please specify file with '-f' or quoted command line with '-l'")
        return None
    
    date=os.popen('echo `date "+%y%m%d.%H-%M-%S"`').read().strip()
    os.system("mkdir -p ./pbs_jobs && mkdir -p ./pbs_stdouts")

    if args.file=="STDIN":
        fn=args.file+"."+date
    else:
        fn=args.file.split("/")[-1]+"."+date
    
    with open("./pbs_jobs/{0}.pbs".format(fn), "w") as pbs:
        to_write=["#PBS -N {0}".format(fn),
                 "#PBS -j oe",
                 "#PBS -q {0}".format(args.queue),
                 "#PBS -o {0}/pbs_stdouts/{1}.stdout.txt".format(os.getcwd(),fn),
                 "#PBS -l nodes=1:ppn={0}".format(args.cores),
                 "source /etc/environment",
                 "source /etc/profile",
                 "source ~/.bashrc",
                 "conda activate",
                 "cd ${PBS_O_WORKDIR}"]
        if args.file!="STDIN":
            os.system("chmod 755 {0}".format(args.file))
            to_write+=["./{0}".format(args.file)]
        else:
            to_write+=[args.line]
        pbs.write("\n".join(to_write))
    toQsub="qsub {0}/pbs_jobs/{1}.pbs".format(os.getcwd(), fn)
    print("Job submitted: "+toQsub)
    print("Stdout+Stderr: {0}/pbs_stdouts/{1}.stdout.txt".format(os.getcwd(),fn))
    os.system(toQsub)
    
if __name__ == "__main__":
    main()

