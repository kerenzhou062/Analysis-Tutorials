#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
import plotly
import plotly.express as px
import time

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', action='store', type=str, required=True,
                    help='input matrix file')
parser.add_argument('-l', '--label', action='store', type=str, required=True,
                    help='name of label column')
parser.add_argument('-o', '--output', action='store', type=str, required=True,
                    help='output pdf')
parser.add_argument('-t', '--title', action='store', type=str,
                    default="Pie Chart",
                    help='title of the pie chart')
parser.add_argument('-v', '--value', action='store', type=str, required=True,
                    help='name of value column (count)')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

#plotly.io.orca.config.executable = '/home/kzhou/.conda/envs/py3/bin/orca'
plotly.io.orca.ensure_server()
time.sleep(20)
#
df = pd.read_csv(args.input, sep='\t')

fig = px.pie(df, values=args.value, names=args.label, title=args.title)

plotly.io.write_image(fig, args.output)

#fig.write_image(args.output)
