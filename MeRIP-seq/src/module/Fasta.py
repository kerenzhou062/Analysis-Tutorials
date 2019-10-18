#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##########################################
#     handle the fasta file              #
#          2018.1.20                    #
##########################################
__author__ = "K.R.Chow"
__version__ = "v1.1"

import sys
import os
from collections import defaultdict


#############################################################################
# argument examples:
#   SeqDict(fastaFile)
#############################################################################

def SeqDict(file):
    fastaDict = defaultdict(dict)
    tempDict = defaultdict(dict)
    with open(file, 'r', errors='ignore') as f:
        index = 0
        for line in f.readlines():
            if line[0] == '>':
                chromosome = line.rstrip('\n')[1:]
                tempDict[index]['seq'] = list()
                tempDict[index]['chr'] = chromosome
                index += 1
            else:
                tempDict[index-1]['seq'].append(line.rstrip('\n'))
    for x in range(0, index):
        fastaDict[tempDict[x]['chr']] = ''.join(tempDict[x]['seq'])
    return fastaDict

if __name__ == '__main__':
    print(SeqDict("/data/zhoukr/phasiRNA/genome/RAP-DB/fasta/IRGSP_1_0_genome.fasta")['chr01'][0:10])
