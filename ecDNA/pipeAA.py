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

def run_bwa(ref, fastqs, outdir, sname, thread, memory, usingDeprecatedSamtools = False):
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
    if usingDeprecatedSamtools:
         cmd = "{{ bwa mem -t {} {} {} | samtools view -Shu - | samtools sort -m {}G -@{} \
             - {}.cs; }} 2>{}_aln_stage.stderr".format(thread, ref, fastqs, memory, thread, outname, outname)
    else:
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

def run_cnvkit(thread, outdir, bamfile, vcf=None):
    #CNVkit cmd-line args
    # -m wgs: wgs data
    # -y: assume chrY present
    # -n: create flat reference (cnv baseline)
    # -p: number of threads
    # -f: reference genome fasta
    bamBase = os.path.splitext(os.path.basename(bamfile))[0]
    ckRef = AA_REPO + args.ref + "/" + args.ref + "_cnvkit_filtered_ref.cnn"
    cmd = "cnvkit.py batch -m wgs -y -r {} -p {} -d {} {}".format(ckRef, thread, outdir, bamfile)
    print(cmd)
    call(cmd,shell=True)

    cnrFile = outdir + bamBase + ".cnr"
    cnsFile = outdir + bamBase + ".cns"

    cmd = "cnvkit.py segment {} -p {} -o {}".format(cnrFile, thread, cnsFile)
    print(cmd)
    call(cmd,shell=True)

#Read the CNVkit .cns files
def convert_cnvkit_cnv_to_seeds(cnvkit_output_directory,bam):
    base = os.path.splitext(os.path.basename(bam))[0]
    with open(cnvkit_output_directory + base + ".cns") as infile, open(cnvkit_output_directory + base + "_CNV_GAIN.bed",'w') as outfile:
        head = infile.next().rstrip().rsplit("\t")
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            s,e = int(fields[1]), int(fields[2])
            cn_r = float(fields[4])
            cn = 2**(cn_r + 1)
            if cn >= args.cngain and e - s >= args.cnsize_min:
                outline = "\t".join(fields[0:3] + ["CNVkit",str(cn)]) + "\n"
                outfile.write(outline)

    return cnvkit_output_directory + base + "_CNV_GAIN.bed"

def run_amplified_intervals(CNV_seeds_filename, sorted_bam, output_directory, sname, cngain, cnsize_min):
    print "Running amplified_intervals"
    AA_seeds_filename = "{}_AA_CNV_SEEDS".format(output_directory + sname)
    cmd = "python {}/amplified_intervals.py --ref {} --bed {} --bam {} --gain {} \
        --cnsize_min {} --out {}".format(AA_SRC, args.ref, CNV_seeds_filename, sorted_bam, str(cngain), str(cnsize_min), AA_seeds_filename)
    print cmd
    call(cmd,shell=True)

    return AA_seeds_filename + ".bed"

def run_AA(amplified_interval_bed, sorted_bam, AA_outdir, sname, downsample, ref):
    print("Running AA with default arguments (& downsample " + str(downsample) + "). To change settings run AA separately.")
    cmd = "python {}/AmpliconArchitect.py --ref {} --downsample {} --bed {} \
        --bam {} --out {}/{}".format(AA_SRC, ref, str(downsample), amplified_interval_bed, sorted_bam, AA_outdir, sname)
    print cmd
    call(cmd,shell=True)

def get_ref_sizes(ref_genome_size_file):
    chr_sizes = {}
    with open(ref_genome_size_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            chr_sizes[fields[0]] = str(int(fields[1])-1)

    return chr_sizes

def get_ref_centromeres(ref_name):
    centromere_dict = {} 
    fnameD = {"GRCh38":"GRCh38_centromere.bed", "GRCh37":"human_g1k_v37_centromere.bed", "hg19":"hg19_centromere.bed"}
    with open(AA_REPO + ref_name + "/" + fnameD[ref_name]) as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            if fields[0] not in centromere_dict:
                centromere_dict[fields[0]] = (fields[1],fields[2])

            else:
                pmin = min(int(centromere_dict[fields[0]][0]),int(fields[1]))
                pmax = max(int(centromere_dict[fields[0]][1]),int(fields[2]))
                #pad with 20kb
                centromere_dict[fields[0]] = (str(pmin-20000),str(pmax+20000))

    return centromere_dict

### MAIN ###
if __name__ == '__main__':
    # Parses the command line arguments
    parser = argparse.ArgumentParser(description="Modified pipeline for AmpliconArchitect (derived from PrepareAA.py, supports cnvkit.py only)")
    parser.add_argument("-o", "--output_directory",
                        help="output directory names (will create if not already created)")
    parser.add_argument("-s", "--sample_name",
                        required=True,
                        help="sample name")
    parser.add_argument("-t","--thread",
                        type=int,
                        required=True,
                        help="Number of threads to use in BWA and CNV calling")
    parser.add_argument("-m", "--memory",
                        type=int,
                        default=60,
                        help="memory to use in BWA")
    parser.add_argument("--run_AA",
                        help="Run AA after all files prepared. Default off.", action='store_true')
    parser.add_argument("--ref",
                        help="Reference genome version.",
                        choices=["hg19","GRCh37","GRCh38"],default="hg19")
    parser.add_argument("--cngain",
                        type=float,
                        default=4.999999,
                        help="CN gain threshold to consider for AA seeding")
    parser.add_argument("--cnsize_min",
                        type=int,
                        default=50000,
                        help="CN interval size (in bp) to consider for AA seeding")
    parser.add_argument("--downsample",
                        type=float,
                        help="AA downsample argument (see AA documentation)",default=5)
    parser.add_argument("--skip_mapping",
                        action='store_true',
                        default=False,
                        help="Skip genome mapping (bwa)")
    parser.add_argument("--skip_cnv",
                        action='store_true',
                        default=False,
                        help="Skip cnv calling (cnvkit.py)")
    parser.add_argument("--cnv_bed",
                        default="",
                        help="BED file of CNV changes. Fields in the bed file should be: chr start end name cngain")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sorted_bam",
                        help= "Sorted BAM file (aligned to an AA-supported reference.)")
    group.add_argument("--fastqs",
                        nargs=2,
                        help="Fastq files (r1.fq r2.fq)")

    args = parser.parse_args()

    #Todo: Implement support for non-hg19. Canvas is not well behaved and has documented bugs regarding reference genome selection.
    #Todo: Implement support for different pipeline tools to be used instead.

    #Check if AA_REPO set, print error and quit if not
    try:
        AA_REPO = os.environ['AA_DATA_REPO'] + "/"

    except KeyError:
        sys.stderr.write("AA_DATA_REPO bash variable not found. AmpliconArchitect may not be properly installed.\n")
        sys.exit(1)
    
    try:
        try:
            AA_SRC = os.environ['AA_SRC']
        except KeyError as e:
            AA_SRC = list(filter(lambda x:re.search(r'AmpliconArchitect\/src', x), os.environ['PATH'].split(":")))[0]
    except KeyError:
        sys.stderr.write("AA_SRC bash variable not found. AmpliconArchitect may not be properly installed.\n")
        sys.exit(1)

    runCNV = "CNVkit"

    #Paths of all the repo files needed
    refFnames = {"hg19":"hg19full.fa", "GRCh37":"human_g1k_v37.fasta", "GRCh38":"hg38full.fa"}
    gdir = AA_REPO + args.ref + "/"
    ref = gdir + refFnames[args.ref]
    ref_genome_size_file = gdir + args.ref + "_noAlt.fa.fai"

    #set an output directory if user did not specify
    if not args.output_directory:
        args.output_directory = os.getcwd()

    #make the output directory location if it does not exist
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    if args.cnv_bed and not os.path.isfile(args.cnv_bed):
        sys.stderr.write("Specified CNV bed file does not exist: " + args.cnv_bed + "\n")
        sys.exit(1)

    sname =  args.sample_name
    outdir = args.output_directory + "/"

    print("Running PrepareAA on sample: " + sname)

     #Check if Fastqs provided
    if args.fastqs:
        #Run BWA
        fastqs = " ".join(args.fastqs)
        memory = int(args.memory * 0.9 / args.thread)
        print("Running pipeline on " + fastqs)
        if args.skip_mapping:
            print("Skip the mapping step.")
            args.sorted_bam = outdir + sname + ".cs.rmdup.bam"
        else:
            args.sorted_bam = run_bwa(ref,fastqs, outdir, sname, args.thread, memory, False)

    if not os.path.isfile(args.sorted_bam + ".bai"):
        print(args.sorted_bam + ".bai not found, calling samtools index")
        call(["samtools","index", args.sorted_bam])
        print("Finished indexing")

    centromere_dict = get_ref_centromeres(args.ref)

    #chunk the genome by chr
    chr_sizes = get_ref_sizes(ref_genome_size_file)
    regions = []
    for key,value in chr_sizes.iteritems():
        try:
            cent_tup = centromere_dict[key]
            regions.append((key,"0-" + cent_tup[0],"p"))
            regions.append((key,cent_tup[1] + "-" + value,"q"))
        except KeyError:
            regions.append((key,"0-" + value,""))
    
    #coordinate CNV calling
    if bool(args.skip_cnv) is False:
        if bool(args.cnv_bed) is False:
            cnvkit_output_directory = args.output_directory + "/cnvkit_output/"
            if not os.path.exists(cnvkit_output_directory):
                    os.mkdir(cnvkit_output_directory)
            print('Calling CNV with CNVkit.py')
            run_cnvkit(args.thread, cnvkit_output_directory, args.sorted_bam)
            args.cnv_bed = convert_cnvkit_cnv_to_seeds(cnvkit_output_directory, args.sorted_bam)
        
        #Run amplified_intervals.py
        print('Determining amplified intervals!')
        amplified_interval_bed = run_amplified_intervals(args.cnv_bed, args.sorted_bam, outdir, sname,args.cngain, args.cnsize_min)
    else:
        amplified_interval_bed = "{}_AA_CNV_SEEDS.bed".format(outdir + sname)

    #Run AA
    if args.run_AA:
        AA_outdir = outdir + "/" + sname + "_AA_results"
        if not os.path.exists(AA_outdir):
            os.mkdir(AA_outdir)
        print('Running AmpliconArchitect!')
        run_AA(amplified_interval_bed, args.sorted_bam, AA_outdir, sname, args.downsample, args.ref)
        
    print("Completed\n")
