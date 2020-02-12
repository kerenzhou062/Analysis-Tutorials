#!/usr/bin/env python3
'''
Created on Oct 31, 2019
@author: Keren Zhou
@lab: Jianjun Chen Lab, City of Hope
'''
import sys
import argparse
import re
from shutil import which
import subprocess
from collections import defaultdict

# usage: scanMotif.py -input peaks.bed -format bed12 \
# -output motif.bed -fasta genome.fa -motif RRACH -tag 3

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-input', action='store', type=str, required=True,
                    help='input bed')
parser.add_argument('-fasta', action='store', type=str, required=True,
                    help='the genome fasta')
parser.add_argument('-format', action='store', type=str, required=True,
                    default='bed12', help='bed format')
parser.add_argument('-motif', action='store', type=str, required=True,
                    help='Motif for scan (eg. RRACH)')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='the outut result')
parser.add_argument('-keep', action='store_true',
                    default=False, help='keep duplicate output records')
parser.add_argument('-tag', action='store', type=int, required=True,
                    default=0,
                    help='tagged position in motif \
                    (eg. 3 for tagging A in RRACH)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()


def decodeMotif(motif):
    tempSeq = motif.upper()
    motifList = [tempSeq]
    for i in range(len(tempSeq)):
        codeList = IUPACcodeDict[tempSeq[i]]
        tempList = list()
        for seq in motifList:
            for code in codeList:
                seqSplit = list(seq)
                seqSplit[i] = code
                tempList.append(''.join(seqSplit))
        motifList = tempList
    return sorted(motifList)


# public arguments
IUPACcodeDict = {
    'A': ('A'),
    'C': ('C'),
    'U': ('T'),
    'T': ('T'),
    'G': ('G'),
    'R': ('A', 'G'),
    'Y': ('C', 'T'),
    'S': ('G', 'C'),
    'W': ('A', 'T'),
    'K': ('G', 'T'),
    'M': ('A', 'C'),
    'B': ('C', 'G', 'T'),
    'D': ('A', 'G', 'T'),
    'H': ('A', 'C', 'T'),
    'V': ('A', 'C', 'G'),
    'N': ('A', 'C', 'G', 'T')
}

fasta = args.fasta
bed = args.input
bedFormat = args.format
motifLength = len(args.motif)

# check bedtools

if bool(which('bedtools')) is False:
    print('Error: bedtools missing in path!')
    sys.exit()

if args.tag > motifLength and args.tag < 0:
    print('Error: -tag larger than -motif!')
    sys.exit()

# program start
motifCompileList = list(map(lambda x: re.compile(
    r'{0}'.format(x)), decodeMotif(args.motif)))
posMapDict = defaultdict(dict)
chromMapDict = defaultdict(dict)
headSkip = re.compile(r'^#')
with open(bed, 'r') as f:
    for line in f:
        if headSkip.match(line):
            continue
        row = line.strip().split('\t')
        # ENSG00000233750.3::chr1:133406-133466(+)
        chrom = row[0]
        start = int(row[1])
        end = row[2]
        name = row[3]
        strand = row[5]
        key = '{0}::{1}:{2}-{3}({4})'.format(name, chrom, start, end, strand)
        chromMapDict[key]['name'] = name
        chromMapDict[key]['chrom'] = chrom
        chromMapDict[key]['start'] = str(start)
        chromMapDict[key]['end'] = end
        chromMapDict[key]['strand'] = strand
        if bedFormat == 'bed12':
            blockLengthList = row[10].split(',')
            blockStartsList = row[11].split(',')
            # remove the last ''
            if bool(blockLengthList[-1]) is False:
                __ = blockLengthList.pop()
            if bool(blockStartsList[-1]) is False:
                __ = blockStartsList.pop()
            # skip non-valid line
            if len(blockLengthList) != len(blockStartsList):
                continue
            blockNum = len(blockLengthList)
            blockLengthList = list(map(lambda x: int(x), blockLengthList))
            blockStartsList = list(map(lambda x: int(x), blockStartsList))
            faStart = 0
            faEnd = 0
            for i in range(blockNum):
                blockStart = start + blockStartsList[i]
                faStart = faEnd
                faEnd = faEnd + blockLengthList[i]
                for j in range(faStart, faEnd):
                    posMapDict[key][j] = blockStart + j - faStart
        else:
            blockLength = int(end) - start
            for i in range(blockLength):
                posMapDict[key][i] = start + i

# get fasta by using
split = '-split'
if bedFormat != 'bed12':
    split = ''
getfastaCommand = 'bedtools getfasta -fi {fasta} \
    -bed {bed} -name+ -tab -s {split}'.format(**vars())
bedFasta = bytes.decode(subprocess.check_output(
    getfastaCommand, shell=True)).split('\n')

tagMapDict = dict()
tagMapDict['+'] = args.tag - 1
tagMapDict['-'] = motifLength - args.tag

recordDict = defaultdict(dict)
with open(args.output, 'w') as out:
    row = ['#chrom', 'tagStart', 'tagEnd', 'motif-id', '1-based position',
           'strand', 'motifChrom', 'motifStart', 'motifEnd', 'motifSeq',
           'motifSeqMask']
    out.write('\t'.join(map(str, row)) + '\n')
    for line in bedFasta:
        if bool(line) is False:
            continue
        key, seq = line.strip().split('\t')
        seqUpper = seq.upper()
        seqLength = len(seq)
        startList = []
        for motifCompile in motifCompileList:
            for match in motifCompile.finditer(seqUpper):
                startList.append(match.start())
        if len(startList):
            name = chromMapDict[key]['name']
            chrom = chromMapDict[key]['chrom']
            start = chromMapDict[key]['start']
            end = chromMapDict[key]['end']
            strand = chromMapDict[key]['strand']
            for i in range(len(startList)):
                matchEnd = startList[i] + motifLength
                startFaPos = startList[i]
                endFaPos = startList[i] + motifLength - 1
                motifGStart = posMapDict[key][startFaPos]
                motifGEnd = posMapDict[key][endFaPos] + 1
                tagStartFaPos = startFaPos + args.tag - 1
                if strand == '-':
                    tagStartFaPos = startFaPos + args.tag - 1
                    tagStartFaPos = seqLength - 1 - tagStartFaPos
                    startFaPos = seqLength - 1 - startFaPos
                    endFaPos = seqLength - 1 - endFaPos
                    motifGStart = posMapDict[key][endFaPos]
                    motifGEnd = posMapDict[key][startFaPos] + 1
                motifSeq = seqUpper[startList[i]:matchEnd]
                motifSeqMask = seq[startList[i]:matchEnd]
                ## output
                if args.tag == 0:
                    motifName = '|'.join(
                        [name, chrom, start, end, strand, str(i + 1)])
                    row = [chrom, motifGStart, motifGEnd, motifName, '0',
                           strand, chrom, start, end,
                           motifSeq, motifSeqMask]
                else:
                    ## determin tag position
                    tagStart = posMapDict[key][tagStartFaPos]
                    tagEnd = tagStart + 1
                    ## record coordinates
                    record = '\t'.join([chrom, str(tagStart), str(tagEnd), strand])
                    if record not in recordDict:
                        recordDict[record] = 1
                    else:
                        if args.keep is False:
                            continue
                    tagName = '|'.join(
                        [name, chrom, start, end, strand, str(i + 1)])
                    row = [chrom, tagStart, tagEnd, tagName, tagEnd,
                           strand, chrom, motifGStart, motifGEnd,
                           motifSeq, motifSeqMask]
                out.write('\t'.join(map(str, row)) + '\n')
