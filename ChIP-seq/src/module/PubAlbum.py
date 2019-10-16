'''
Created on Sep 29, 2012
@author: Keren Zhou
@lab: Jianjun Chen Lab, City of Hope
'''

import os
import sys
from collections import defaultdict
import subprocess

class Anno:
    '''Class to map annotation file'''
    publicPath = '/home/kzhou/zhoukr/public/genome'
    def __init__(self, path=publicPath):
        self.path = path
    
    def pathjoin(self, file):
        return os.path.realpath(os.path.join(self.path, file))
    
    def blacklist(self, key, ):
        #ftype: gtf|gff3
        pubDict = defaultdict(dict)
        pubDict['hg38'] = self.pathjoin('annotation/hg38/hg38-blacklist.v2.bed')
        pubDict['hg19'] = self.pathjoin('annotation/hg19/hg19-blacklist.v2.bed')
        pubDict['mm10'] = self.pathjoin('annotation/mm10/mm10-blacklist.v2.bed')
        try:
            file = os.path.realpath(pubDict[key])
        except KeyError:
            file = os.path.realpath(key)
        return file
    
    def fasta(self, key, ):
        #ftype: gtf|gff3
        pubDict = defaultdict(dict)
        pubDict['hg38'] = self.pathjoin('fasta/hg38/hg38.fa')
        pubDict['hg19'] = self.pathjoin('fasta/hg19/hg19.fa')
        pubDict['mm10'] = self.pathjoin('fasta/mm10/hg19.fa')
        try:
            file = os.path.realpath(pubDict[key])
        except KeyError:
            file = os.path.realpath(key)
        return file
    
    def gtf(self, key):
        pubDict = defaultdict(dict)
        pubDict['hg38v31'] = self.pathjoin('annotation/hg38/v31/gencode.v31.annotation.gtf')
        try:
            file = os.path.realpath(pubDict[key])
        except KeyError:
            file = os.path.realpath(key)
        return file
    
    def sqlite(self, key, ftype='gtf'):
        #ftype: gtf|gff3
        pubDict = defaultdict(dict)
        pubDict['hg38v31']['gtf'] = self.pathjoin('annotation/hg38/v31/gencode.v31.annotation.gtf.db')
        pubDict['hg38v31']['gff3'] = self.pathjoin('annotation/hg38/v31/gencode.v31.annotation.gff3.db')
        try:
            file = os.path.realpath(pubDict[key][ftype])
        except KeyError:
            file = os.path.realpath(key)
        return file
    
    def gindex(self, key):
        pubDict = defaultdict(dict)
        pubDict['hg38'] = self.pathjoin('index/bowtie2/hg38')
        if key in pubDict:
            return os.path.realpath(pubDict[key])
        else:
            return os.path.realpath(key)

class GeneClass(object):
    '''return tuple'''
    def __init__(self):
        self.source = 'gencode'
    
    def filter(self, genetype, extend=[]):
        genetypes = list()
        genetypeDict = dict()
        genetypeDict['protein_coding'] = ["protein_coding", "TR_V_gene", "TR_D_gene",
            "TR_J_gene", "TR_C_gene", "IG_LV_gene",
            "IG_V_gene", "IG_J_gene", "IG_C_gene", "IG_D_gene"]
        genetypeDict['lncRNA'] = ["processed_transcript", "lincRNA",
            "3prime_overlapping_ncrna", "3prime_overlapping_ncRNA",
            "antisense", "non_coding", "sense_intronic", 
            "sense_overlapping", "TEC", "known_ncrna", "macro_lncRNA",
            "bidirectional_promoter_lncrna", "bidirectional_promoter_lncRNA"]
        genetypeDict['sncRNA'] = ["snRNA", "snoRNA", "rRNA", "Mt_tRNA",
            "Mt_rRNA", "misc_RNA", "miRNA", "ribozyme", "sRNA",
            "scRNA", "scaRNA", "vaultRNA"]
        genetypeDict['pseudogene'] = ["transcribed_unprocessed_pseudogene",
            "translated_processed_pseudogene", "unprocessed_pseudogene",
            "processed_pseudogene", "transcribed_processed_pseudogene",
            "unitary_pseudogene", "pseudogene", "polymorphic_pseudogene",
            "transcribed_unitary_pseudogene", "TR_V_pseudogene",
            "TR_J_pseudogene", "IG_V_pseudogene", "IG_C_pseudogene",
            "IG_D_pseudogene", "IG_pseudogene", "IG_J_pseudogene"]
        
        try:
            genetypes = genetypeDict[genetype]
        except KeyError:
            genetypes = [genetype]
        if extend:
            genetypes.extend(extend)
        return tuple(genetypes)
