#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "K.R.Chow"

import sys
import os
import re
import argparse
import datetime
import shutil
from multiprocessing import Pool
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-cpu', action='store', type=int,
                    default=10, help='cores for multiprocessing')
parser.add_argument('-chrsize', action='store', type=str,
                    help='chromosome size file (eg. hg38.chrom.sizes)')
parser.add_argument('-dtype', action='store', type=str,
                    default='narrowPeak', help='The --input-file-type \
                    (narrowPeak|broadPeak|bed|gff)')
parser.add_argument('-ip', nargs='+', type=str,
                    help='The ip files')
parser.add_argument('-input', nargs='+', type=str,
                    help='The input file strings')
parser.add_argument('-gsize', action='store', type=str,
                    default='hs', help='genome size')
parser.add_argument('-format', action='store', type=str,
                    help='input format (BED, BAM...)')
parser.add_argument('-shift', action='store', type=int,
                    default=0, help='shift value')
parser.add_argument('-extsize', action='store', type=int,
                    default=147, help='extsize value')
parser.add_argument('-name', nargs='+', type=str,
                    help='The prefix name')
parser.add_argument('-output', nargs='+', type=str,
                    help='The output result directories')
parser.add_argument('-pval', action='store', type=float,
                    default=1e-2, help='The p-value cutoff')
parser.add_argument('-other', action='store', type=str,
                    help='The other commands, \
                    "," separated (eg.--broad,--broad)')
parser.add_argument('-rank', action='store', type=str,
                    default='signal.value', help='Which column to use to rank peaks \
                    (idr, signal.value|p.value|q.value)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()


def runMacs2(ip, control, name, other, output,
             extsize, gsize, shift, pval,
             dFormat, dtype, rank, chrsize, sortMem):
    # calling peaks
    logFile = '{output}/../{name}.macs2AndTrack.log 2>&1'.format(**vars())
    macs2Command = 'macs2 callpeak -t {ip} \
        -c {control} -f {dFormat} -g {gsize} \
        -p {pval} --nomodel --shift {shift} \
        --extsize {extsize} -n {name} --outdir {output} \
        {other} > {logFile} 2>&1'.format(**vars())
    subprocess.run(macs2Command, shell=True)

    # delete treat_pileup.bdg and control_lambda.bdg 
    # for pseudo-replicates
    # keep bdg files for pooled and true replicates
    # for macs2-bdfdiff analysis
    if bool(re.search(r'(pr[12])|(pr[12]_pooled)', name)):
        prefix = '{output}/{name}'.format(**vars())
        os.remove('{0}_treat_pileup.bdg'.format(prefix))
        os.remove('{0}_control_lambda.bdg'.format(prefix))

    # sort -S {sortMem}% peaks by signal.value | -log10(p-value) | -log10(p-value)
    peakFile = os.path.join(output, name + '_peaks.{0}'.format(dtype))
    broadFlag = False
    if 'broadPeak' in dtype:
        broadFlag = True
    if rank == 'signal.value':
        sortCommand = 'sort -S {sortMem}% -k7,7n {peakFile} -o {peakFile}'.format(**vars())
    elif rank == 'p.value':
        sortCommand = 'sort -S {sortMem}% -k8,8nr {peakFile} -o {peakFile}'.format(**vars())
    else:
        sortCommand = 'sort -S {sortMem}% -k9,9nr {peakFile} -o {peakFile}'.format(**vars())
    subprocess.run(sortCommand, shell=True)

    # generate signal track
    # skip generate signal track: broad peak or pseudo-replicates
    if broadFlag is False and bool(re.search(r'(pr[12])|(pr[12]_pooled)', name)) is False:
        prefix = '{output}/{name}'.format(**vars())
        signalCommand = 'macs2 bdgcmp -t {prefix}_treat_pileup.bdg \
            -c {prefix}_control_lambda.bdg \
            --outdir {output} -o {name}_FE.bdg -m FE \
            >> {logFile} 2>&1'.format(**vars())
        subprocess.run(signalCommand, shell=True)

        # Remove coordinates outside chromosome sizes (stupid MACS2 bug)
        removeCommand = "slopBed -i {prefix}_FE.bdg -g {chrsize} \
            -b 0 | awk '{{if ($3 != -1) print $0}}' | bedClip stdin {chrsize} \
            {prefix}.fc.signal.bedgraph >> {logFile} 2>&1".format(**vars())
        subprocess.run(removeCommand, shell=True)
        os.remove('{prefix}_FE.bdg'.format(**vars()))

        # Convert bedgraph to bigwig, sort with -S N%, to reduce memory
        sortCommand = 'sort -S {sortMem}% -k1,1 -k2,2n {prefix}.fc.signal.bedgraph \
            -o {prefix}.fc.signal.bedgraph'.format(**vars())
        subprocess.run(sortCommand, shell=True)
        mergeCommand = 'bedtools merge -i {prefix}.fc.signal.bedgraph -c 4 -d -1 -o max \
            > {prefix}.fc.signal.sorted.bedgraph'.format(**vars())
        subprocess.run(mergeCommand, shell=True)
        sortCommand = 'sort -S {sortMem}% -k1,1 -k2,2n {prefix}.fc.signal.sorted.bedgraph \
            -o {prefix}.fc.signal.sorted.bedgraph'.format(**vars())
        subprocess.run(sortCommand, shell=True)
        convertCommand = "bedGraphToBigWig {prefix}.fc.signal.sorted.bedgraph \
            {chrsize} {prefix}.fc.signal.bw >> {logFile} 2>&1".format(**vars())
        subprocess.run(convertCommand, shell=True)
        os.remove('{0}.fc.signal.bedgraph'.format(prefix))
        os.remove('{0}.fc.signal.sorted.bedgraph'.format(prefix))

        # generate count signal track
        # postive strand
        countCommand = "zcat -f {ip} | sort -S {sortMem}% -k1,1 -k2,2n | \
            bedtools genomecov -5 -bg -strand + -g {chrsize} \
            -i stdin > {output}/TMP.POS.BED".format(**vars())
        subprocess.run(countCommand, shell=True)
        countCommand = "bedGraphToBigWig {output}/TMP.POS.BED \
            {chrsize} {output}/{name}.positive.bw >> \
            {logFile} 2>&1".format(**vars())
        subprocess.run(countCommand, shell=True)
        # negative strand
        countCommand = "zcat -f {ip} | sort -S {sortMem}% -k1,1 -k2,2n | \
            bedtools genomecov -5 -bg -strand - -g {chrsize} \
            -i stdin > {output}/TMP.POS.BED".format(**vars())
        subprocess.run(countCommand, shell=True)
        countCommand = "bedGraphToBigWig {output}/TMP.POS.BED {chrsize} \
            {output}/{name}.negative.bw >> {logFile} 2>&1".format(**vars())
        subprocess.run(countCommand, shell=True)
        os.remove('{output}/TMP.POS.BED'.format(**vars()))
    return 1


# public variables
ipList = args.ip
inputList = args.input
nameList = args.name
otherList = args.other.split(',')
chrsize = args.chrsize
outputList = args.output
extsize = args.extsize
gsize = args.gsize
shift = args.shift
pval = args.pval
dFormat = args.format
dtype = args.dtype
rank = args.rank

# judge elements in same length
if len(ipList) == len(inputList) and len(ipList) == len(nameList):
    if len(ipList) == len(otherList):
        if len(ipList) == len(outputList):
            pass
        else:
            sys.exit('not equal parameters')
    else:
        sys.exit('not equal parameters')
else:
    sys.exit('not equal parameters')

starttime = datetime.datetime.now()
sys.stderr.write("Parallel peak-calling with macs2!\n")

# multiple call macs2
sortMem = int(95 / args.cpu)
pool = Pool(processes=args.cpu)
for i in range(len(ipList)):
    ip = ipList[i]
    control = inputList[i]
    name = nameList[i]
    other = otherList[i]
    output = outputList[i]
    # try to make dirs
    try:
        shutil.rmtree(output)
    except FileNotFoundError:
        pass
    os.makedirs(output, exist_ok=True)
    # runMacs2(ip, control, name, other,
    #    output, extsize, gsize, shift, pval, dFormat, dtype, rank, chrsize)
    result = pool.apply_async(runMacs2, args=(ip, control, name, other,
                                              output, extsize, gsize,
                                              shift, pval, dFormat,
                                              dtype, rank, chrsize, sortMem))
pool.close()
pool.join()

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("All jobs done!")
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
