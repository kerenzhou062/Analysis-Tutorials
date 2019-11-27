#!/usr/bin/env python3
import os
import sys
import argparse
import gffutils
from PubAlbum import Anno

#usage: buildSqlite.py or buildSqlite.py <fastq dir>

parser = argparse.ArgumentParser()
parser.add_argument('-input', action='store', type=str,
                     help='input annotation file (gtf|gff3)')
parser.add_argument('-output', action='store', type=str,
                    help='sqlite3-file based databse file')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

print("Start building sqlite...")
sqliteDb = gffutils.create_db(args.input, dbfn=args.output, force=True, 
    keep_order=True, merge_strategy='merge', sort_attribute_values=True)

print("Jobs finished!")
