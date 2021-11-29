#!/bin/bash
#SBATCH --job-name=hic2post    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org      # Where to send mail
#SBATCH -n 4                         # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p gpu                        # default queue is all if you don't specify
#SBATCH --mem=100G                      # Amount of memory in GB
#SBATCH --gres=gpu:4
#SBATCH --core-spec=0
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=hic2post.log   # Standard output and error log

##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########
# Juicer postprocessing script.
# This will run the major post-processing on the HiC file, including finding
# loops with HiCCUPS, finding motifs of these loops with MotifFinder, and
# finding contact domains with Arrowhead.
# Juicer version 1.5

## Read arguments

function showHelp {
  echo -e "This script is designed to be run on a slurm server with GPUs
  Usage: sbatch -p gpu -o <log> -n <int> --gres \"gpu:4\" juicer2hic.sh -d <distance> -t 4 --jtools <path> --knorm <Normalizations> --f <fdr>
  --loops <merge thresholds> --peak <peak_width> -o <output directory> --resolution <resolution> --window <window width>"
  echo -e "options:
    -b|--bed: path to a bed directory [RAD21, SMC3, and CTCF], finding DNA motifs for loops (https://github.com/aidenlab/juicer/wiki/MotifFinder) <str>
    -c|--hic: the input .hic file <str>
    -d|--distance: Distances used for merging nearby pixels to a centroid (comma-separated) [20000,20000,30000,40000,50000] <str>
    -e|--prefix: prefix of output files [hic2post] <str>
    -f|--fdr: FDR values actually corresponding to max_q_val (comma-separated) [0.1,0.1,0.1,0.1,0.1] <str>
    -g|--genome: genome ID or genome size file (chrom.size) <str>
    -k|--knorm: Normalizations (case sensitive) that can be selected (NONE/VC/VC_SQRT/KR) [KR] <str>
    -j|--jtools: path to juicer_tools <str>
    -l|--loop: Thresholds for merging loop lists of different resolutions (comma-separated). [0.02,1,1.5,1.75,2] <str>
    -o|--output: the output directory ['./hic2post'] <str>
    -p|--peak: Peak width used for finding enriched pixels in HiCCUPS (comma-separated) [4,3,2,1,1] <str>
    -r|--resolution: Resolution(s) for which HiCCUPS will be run. (comma-separated) [5000,10000,15000,20000,25000] <str>
    -s|--size: Maximum size of the submatrix within the chromosome passed on to GPU (Must be an even number) [3000] <int>
    -t|--threads: number of threads when running BWA alignment (limited to the number of GPU cards, 4 in GEMINI) [4] <int>
    -w|--window: Window width used for finding enriched pixels in HiCCUPS (comma-separated) <str>
    -h|--help: print this help and exit"
    exit 2
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

## Read arguments
TEMP=`getopt -o b:c:d:e:f:g:k:j:o:p:r:s:t:w:h, --long bed:,prefix:,hic:,distance:,fdr:,genome:,knorm:, \
  --long jtools:,output:,peak:,resolution:,size:,threads:,window:, \
  --long help, \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# default arguments
threads=4
distance="20000,20000,30000,40000,50000"
fdr="0.1,0.1,0.1,0.1,0.1"
knorm="KR"
loop="0.02,1,1.5,1.75,2"
peak="4,3,2,1,1"
resolution="5000,10000,15000,20000,25000"
window="7,6,5,4,3"
size=3000
output="./hic2post"
prefix="hic2post"

while true; do
  case "$1" in
    -h |--help) showHelp; shift ;;
    -bed |--bed) bed="$2"; shift 2 ;;
    -c |--hic) hicfile="$2"; shift 2 ;;
    -d |--distance) distance="$2"; shift 2 ;;
    -e |--prefix) prefix="$2"; shift 2 ;;
    -f |--fdr) fdr="$2"; shift 2 ;;
    -g |--genome) genome="$2"; shift 2 ;;
    -k |--knorm) knorm="$2"; shift 2 ;;
    -j |--jtools) juicer_tools="$2"; shift 2 ;;
    -l |--loop) loop="$2"; shift 2 ;;
    -o |--output) output="$2"; shift 2 ;;
    -p |--peak) peak="$2"; shift 2 ;;
    -r |--resolution) size="$2"; shift 2 ;;
    -s |--size) size="$2"; shift 2 ;;
    -t |--threads) threads="$2"; shift 2 ;;
    -w |--window) window="$2"; shift 2 ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

## print input arguments
echo "hicfile=$hicfile"
echo "distance=$distance"
echo "fdr=$fdr"
echo "knorm=$knorm"
echo "juicer_tools=$juicer_tools"
echo "motif=$motif"
echo "output=$output"
echo "peak=$peak"
echo "size=$size"
echo "threads=$threads"
echo "window=$window"
echo ""

start=$(date +%s.%N)
## check arguments
if [[ ! -f "$juicer_tools" ]]; then
  echo "-j|--jtools Wrong! No juicer_tools found!";
  exit 1
else
  jtoolsName=$(basename "$juicer_tools")
  if [ "$jtoolsName" != "juicer_tools" ]; then
    echo "-j|--jtools Wrong! No juicer_tools found!";
    exit 1
  fi
fi

if [[ ! -f $hicfile ]]; then
  echo "-c|hic Wrong! No .hic file found!";
  exit 1
fi

if [[ ! -d $output ]]; then
  mkdir -p $output
else
  rm -rf $output
  mkdir -p $output
fi

## loading modules
module load Java/1.8.0_162
module load CUDA/8.0.61_375.26-GCC-5.4.0-2.26

## variables
arrowhead_dir="$output/arrowhead"
loops_dir="$output/loops"
apa_dir="$output/apa"
motif_dir="$output/motif"

#arrowhead
echo -e "Creating arrowhead directory...\n"
mkdir -p "$arrowhead_dir"

echo -e "Starting to generate ARROWHEAD...\n"
echo -e "arrowhead for 5k"
${juicer_tools} arrowhead ${hicfile} -k $knorm -m $size --threads $threads -r 5000 "${arrowhead_dir}"
mv "${arrowhead_dir}/5000_blocks.bedpe" "${arrowhead_dir}/${prefix}.blocks.${knorm}.5k.bedpe"
##bedpe to bed
bedpe2bed.py -i "${arrowhead_dir}/${prefix}.blocks.${knorm}.5k.bedpe" --addchr \
  --prefix "${prefix}_${knorm}_5k_blocks" --output "${arrowhead_dir}/${prefix}.blocks.${knorm}.5k.bed"

echo -e "arrowhead for 10k"
${juicer_tools} arrowhead ${hicfile} -k $knorm -m $size --threads $threads -r 10000 "${arrowhead_dir}"
mv "${arrowhead_dir}/10000_blocks.bedpe" "${arrowhead_dir}/${prefix}.blocks.${knorm}.10k.bedpe"
##bedpe to bed
bedpe2bed.py -i "${arrowhead_dir}/${prefix}.blocks.${knorm}.10k.bedpe" --addchr \
  --prefix "${prefix}_${knorm}_10k_blocks" --output "${arrowhead_dir}/${prefix}.blocks.${knorm}.10k.bed"

if [ $? -ne 0 ]; then
    echo "***! Problem while running Arrowhead";
    exit 1
fi
echo -e "ARROWHEAD done!\n"

loops_bedpe="${loops_dir}/${prefix}.mloops.bedpe"

echo -e "Creating loop directory...\n"
mkdir -p "$loops_dir"

echo -e "Starting to generate HiCCUPS..\n"
if hash nvcc 2>/dev/null 
then
    ${juicer_tools} hiccups --threads $threads -m $size -r "$resolution" -k "$knorm" -f "$fdr" \
        -p "$peak" -i "$window" -t "$loop" -d "$distance" --ignore-sparsity "$hicfile" "$loops_dir"
    if [ $? -ne 0 ]; then
    echo "***! Problem while running HiCCUPS";
    exit 1
    fi
    mv "${loops_dir}/merged_loops.bedpe" "$loops_bedpe"
    ##bedpe to bed
    bedpe2bed.py -i "${loops_bedpe}" --addchr \
      --prefix "${prefix}_merge_loops" --output "${loops_dir}/${prefix}.mloops.bed"
else 
    echo "GPUs are not installed so HiCCUPs cannot be run";
fi

if [ -f "$loops_bedpe" ]
then
    echo -e "Generating APA:\n"
    mkdir -p $apa_dir

    ${juicer_tools} apa "$hicfile" "$loops_bedpe" "$apa_dir"
    ## Check that bed folder exists    
    if [[ -z $bed ]] || [[ -d $genome ]]; then
      echo "***! Can't find bed folder or chrom.size file";
      echo "***! Not running motif finder";
    else
      if [ ! -e "${bed}" ] | [ ! -f "${genome}" ]; then
      echo "***! Can't find bed folder or chrom.size file";
      echo "***! Not running motif finder";
      else
      echo -e "\nMOTIF FINDER:\n"
      motif_loop="${motif_loop}/${prefix}.motif_loops.txt"
      ${juicer_tools} motifs "${genome}" "${bed}" "${motif_loop}"
      fi
      echo -e "\n(-: Feature annotation successfully completed (-:"
    fi

else
  # if loop lists do not exist but Juicer Tools didn't return an error, likely 
  # too sparse
    echo -e "\n(-: Postprocessing successfully completed, maps too sparse to annotate or GPUs unavailable (-:"
fi

## end time
end=$(date +%s.%N)
dt=$(echo "$end - $start" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

LC_NUMERIC=C printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds
