#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict
import subprocess

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input', action='store', type=str, required=True,
                    help='the root path of input directory')
parser.add_argument('--output', action='store', type=str, required=True,
                    help='the root path of output directory')
parser.add_argument('--final', action='store', type=str, required=True,
                    help='the final name of the folder in the path')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

inputDir = os.path.realpath(args.input)
outputDir = os.path.realpath(args.output)
command = 'find {0} -type d'.format(inputDir)
folderList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')

for folder in folderList:
    baseName = os.path.basename(folder)
    if baseName == args.final:
        outputFolder = folder.replace(inputDir, outputDir)
        command = 'mkdir -p {0}'.format(outputFolder)
        subprocess.run(command, shell=True)
        command = 'cp -r {0}/* {1}'.format(folder, outputFolder)
        subprocess.run(command, shell=True)
