#!/usr/bin/env python2

import sys
import os
import threading
import re
from subprocess import call
import argparse
import gzip

#generic worker thread function
class workerThread(threading.Thread):
    def __init__(self, threadID, target, *args):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self._target = target
        self._args = args
        threading.Thread.__init__(self)

    def run(self): 
        self._target(*self._args)

def run_bwa(ref, fastqs, outdir, sname, thread, memory):
    outname = outdir + sname
    print(outname)
    print("Checking for ref index")
    exts = [".sa",".amb",".ann",".pac",".bwt"]
    indexPresent = True
    for i in exts:
        if not os.path.exists(ref + i):
            indexPresent = False
            print "Could not find " + ref + i + ", building BWA index from scratch. This could take > 60 minutes"
            break

    if not indexPresent:
        cmd = "bwa index " + ref
        call(cmd, shell=True)


    print("Performing alignment and sorting")
    cmd = "{{ bwa mem -t {} {} {} | samtools view -Shu - | samtools sort -m {}G -@{} \
        -o {}.cs.bam -; }} 2>{}_aln_stage.stderr".format(thread, ref, fastqs, memory, thread, outname, outname)

    print(cmd)
    call(cmd, shell=True)
    print("Performing duplicate removal & indexing")
    cmd_list = ["samtools", "rmdup", "-s", "{}.cs.bam".format(outname), "{}.cs.rmdup.bam".format(outname)]
    print(" ".join(cmd_list))
    call(cmd_list)
    print("Running samtools index")
    cmd_list = ["samtools", "index", "{}.cs.rmdup.bam".format(outname)]
    print(" ".join(cmd_list))
    call(cmd_list)
    print("Removing temp BAM")
    cmd = "rm {}.cs.bam".format(outname)
    call(cmd, shell=True)
    return outname + ".cs.rmdup.bam"

### MAIN ###
if __name__ == '__main__':
    # Parses the command line arguments
    parser = argparse.ArgumentParser(description="Run bwa genome mapping for AmpliconArchitect")
    parser.add_argument("-o", "--output_directory",
                        help="output directory names (will create if not already created)")
    parser.add_argument("-s", "--sample_name",
                        required=True,
                        help="sample name")
    parser.add_argument("-t", "--thread",
                        type=int,
                        required=True,
                        help="Number of threads to use in BWA and CNV calling")
    parser.add_argument("-m", "--memory",
                        type=int,
                        default=60,
                        help="memory to use in BWA")
    parser.add_argument("--ref",
                       choices=["hg19","GRCh37","GRCh38"],
                       default="GRCh38",
                       help="Reference genome version.")
    parser.add_argument("--fastqs",
                       nargs=2,
                       help="Fastq files (r1.fq r2.fq)")

    args = parser.parse_args()

    #Check if AA_REPO set, print error and quit if not
    try:
        AA_REPO = os.environ['AA_DATA_REPO'] + "/"

    except KeyError:
        sys.stderr.write("AA_DATA_REPO bash variable not found. AmpliconArchitect may not be properly installed.\n")
        sys.exit(1)

    #Paths of all the repo files needed
    refFnames = {"hg19":"hg19full.fa","GRCh37":"human_g1k_v37.fasta","GRCh38":"hg38full.fa"}
    gdir = AA_REPO + args.ref + "/"
    ref = gdir + refFnames[args.ref]
    ref_genome_size_file = gdir + args.ref + "_noAlt.fa.fai"

    #set an output directory if user did not specify
    if not args.output_directory:
        args.output_directory = os.getcwd()

    #make the output directory location if it does not exist
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    sname =  args.sample_name
    outdir = args.output_directory + "/"

    print("Running bwa mapping on sample: " + sname)

    #Check if Fastqs provided
    memory = int(args.memory * 0.75 / args.thread)
    if args.fastqs:
        #Run BWA
        fastqs = " ".join(args.fastqs)
        print("Running bwa on " + fastqs)
        run_bwa(ref, fastqs, outdir, sname, args.thread, memory)
    print("Completed\n")
