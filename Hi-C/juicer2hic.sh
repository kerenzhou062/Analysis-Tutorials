#!/bin/bash
#SBATCH --job-name=juicer2hic    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kzhou@coh.org     # Where to send mail 
#SBATCH -n 15                            # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p all                        # default queue is all if you don't specify
#SBATCH --mem=150G                      # Amount of memory in GB
#SBATCH --time=120:10:00               # Time limit hrs:min:sec
#SBATCH --output=juicer2hic.log               # Time limit hrs:min:sec

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.

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
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##########
#
# Single CPU version of Juicer.
#
# Alignment script. Sets the reference genome and genome ID based on the input
# arguments (default human, none). Optional arguments are description for 
# stats file, stage to relaunch at, paths to various files if 
# needed, path to scripts directory, and the top-level directory (default 
# current directory). In lieu of setting the genome ID, you can instead set the
# reference sequence and the chrom.sizes file path, but the directory 
# containing the reference sequence must also contain the BWA index files.
#
# Aligns the fastq files, sorts, and merges.
#
# If all is successful, takes the final merged file, removes name duplicates,
# removes PCR duplicates, and creates the hic job and stats job.  Final
# product will be hic file and stats file in the aligned directory.
#                                                                       
# [output]/fastq  - Should contain the fastq files. This code assumes that
#                   there is an "_R" in the appropriate files, i.e. *_R*.fastq
# From the top-level directory, the following two directories are created:
#                                                                              
# [output]/splits  - Where to write the scratch split files (fastq files and
#                    intermediate SAM files). This can be deleted after 
#                    execution.
# [output]/aligned - Where to write the final output files.
#
# The following globals should be set correctly before proceeding:
#
# read1str  - portion of fastq filename that indicates this is the "read 1"
#             file; used to loop over only the read 1 and within that loop,
#             also align read 2 and merge.  If this is not set correctly,
#             script will not work. The error will often manifest itself
#             through a "*" in the name because the wildcard was not able to
#             match any files with the read1str.   
#set -e ## This is causing problems; need better error detection
# Juicer version 1.6
shopt -s extglob
export LC_ALL=C

juicer_version="1.6" 
### LOAD BWA AND SAMTOOLS


# fastq files should look like filename_R1.fastq and filename_R2.fastq 
# if your fastq files look different, change this value
read1str="_R1" 
read2str="_R2" 

## Default options, overridden by command line arguments

# Juicer directory, contains scripts/, references/, and restriction_sites/
# can also be set in options via -D
juicerdir="/opt/juicer"
# top level directory, can also be set in options
output=$(pwd)
# restriction enzyme, can also be set in options
site="MboI"
# genome ID, default to human, can also be set in options
genomeid="hg19"
# description, default empty
about=""
# do not include fragment delimited maps by default
nofrag=1
# use wobble for dedupping by default (not just exact matches)
justexact=0

resolution="5000,10000,15000,20000,25000,50000,100000,250000,500000,1000000,2500000"
## stage
## "merge" when alignment has finished but the merged_sort file has not yet been created.
## "dedup" when the files have been merged into merged_sort but merged_nodups has not yet been created.
## "final" when the reads have been deduped into merged_nodups but the final stats and hic files have not yet been created.
## "postproc" when the hic files have been created and only postprocessing feature annotation remains to be completed.
## "early" for an early exit, before the final creation of the stats and hic files

function showHelp {
  echo -e "Usage: juicer2hic.sh -d <Juicer directory> -t 15 --refseq <BWA index> --site <MboI(ENCODE)> --gsize <chrom.sizes>
  --fastq1 <read1> --fastq2 <read2> -o <output directory> --gsize <chrom.sizes> --sitefile <restriction site file>"
  echo -e "options:
    -a|--about: the description of experiment, enclosed in single quotes <str>
    -b|--ligation: use this string when counting ligation junctions <str>
    -d|--juicerdir: set the Juicer directory, which should have scripts/ references/ and restriction_sites/ underneath it (default ${juicerdir}) <str>
    -e|--prefix: prefix of output .hic files <str>
    -f|--nofrag: include fragment-delimited maps in hic file creation <bool>
    -g|--genomeid: must be defined in the script, e.g. 'hg19' or 'mm10' (default '$genomeid'); alternatively, it can be defined using the -z command <str>
    -j|--justexact: just exact duplicates excluded at dedupping step <bool>
    -l|--gsize: path to chrom.sizes file <str>
    -n|--stage: must be one of 'merge', 'dedup', 'final', 'postproc', or 'early' (https://github.com/aidenlab/juicer/wiki/Usage) <str>
    -o|--output: the output directory (default './') <str>
    -p|--fastq1: input fastq-read1, separated by ',' (*_R1.fastq.gz) <str>
    -q|--fastq2: input fastq-read2, separated by ',' (*_R2.fastq.gz) <str>
    -r|--resolution: only calculate specific resolutions for .hic, comma-separated list
                     <5000,10000,15000,20000,25000,50000,100000,250000,500000,1000000,2500000> <str>
    -s|--site: enzyme used in the Hi-C experiments, e.g.  'HindIII' or 'MboI' (default '$site') <str>
    -t|--threads: number of threads when running BWA alignment [10] <int>
    -y|--sitefile: the path for restriction site file (locations of restriction sites in genome; generate_site_positions.py) <str>
    -z|--refseq: the path for reference sequence file, BWA index files must be in same directory <str>
    -h|--help: print this help and exit <bool>
    --delete: delete output directory before running <bool>"
    exit 2
}

# show help when no arguments input
if [[ $# == 0 ]]; then
  showHelp
  exit 2
fi

## Read arguments
TEMP=`getopt -o a:b:d:e:fg:jl:n:o:p:q:r:s:t:y:z:h, --long fastq1:,fastq2:,genomeid:,output:,site:,shortreadend:,about:, \
  --long prefix:,gsize:,sitefile:,refseq:,stage:,juicedir:,nofrag,ligation:,threads:,resolution:, \
  --long help,delete,justexact, \
  -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

#default arguments
prefix="juicer2hic"
delete=false

while true; do
  case "$1" in
    -h | --help) showHelp; shift ;;
    -a | --about) about="$2"; shift 2 ;;
    -b | --ligation) ligation="$2"; shift 2 ;;
    -d | --juicerdir) juicerdir="$2"; shift 2 ;;
    -e | --prefix) prefix="$2"; shift 2 ;;
    -f | --nofrag) nofrag=0; shift ;;
    -g | --genomeid) genomeid="$2"; shift 2 ;;
    -j | --justexact) justexact=1; shift ;;
    -l | --gsize) gsize="$2"; shift 2 ;;
    -n | --stage) stage="$2"; shift 2 ;;
    -o | --output) output="$2"; shift 2 ;;
    -p | --fastq1) fastq1="$2"; shift 2 ;;
    -q | --fastq2) fastq2="$2"; shift 2 ;;
    -r | --resolution) resolution="$2"; shift 2 ;;
    -s | --site) site="$2"; shift 2 ;;
    -t | --threads) threads="$2"; shift 2 ;;
    -y | --sitefile) sitefile="$2"; shift 2 ;;
    -z | --refseq) refseq="$2"; shift 2 ;;
    --help) delete=true; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

##print input arguments
echo "about=$about"
echo "delete=$delete"
echo "ligation=$ligation"
echo "juicerdir=$juicerdir"
echo "prefix=$prefix"
echo "nofrag=$nofrag"
echo "genomeid=$genomeid"
echo "justexact=$justexact"
echo "gsize=$gsize"
echo "stage=$stage"
echo "output=$output"
echo "fastq1=$fastq1"
echo "fastq2=$fastq2"
echo "resolution=$resolution"
echo "site=$site"
echo "threads=$threads"
echo "sitefile=$sitefile"
echo "refseq=$refseq"
echo ""

if [ ! -z "$stage" ]
then
    case $stage in
        merge) merge=1 ;;
        dedup) dedup=1 ;;
        early) earlyexit=1 ;;
        final) final=1 ;;
        chimeric) chimeric=1 ;;
        postproc) postproc=1 ;; 
        *) showHelp;
        exit 1
    esac
fi
  
## Set reference sequence based on genome ID
if [ -z "$refseq" ]
then 
    case $genomeid in
    mm9) refseq="${juicerdir}/references/Mus_musculus_assembly9_norandom.fasta";;
    mm10) refseq="${juicerdir}/references/Mus_musculus_assembly10.fasta";;
    hg38) refseq="${juicerdir}/references/Homo_sapiens_assembly38.fasta";;
    hg19) refseq="${juicerdir}/references/Homo_sapiens_assembly19.fasta";;
    *)  echo "$usageHelp"
            echo "$genomeHelp"
            exit 1
    esac
else
    # Reference sequence passed in, so gsize must be set for the .hic file
    # to be properly created
    if [ -z "$gsize" ]
    then
        echo "***! You must define a chrom.sizes file via the \"-p\" flag that delineates the lengths of the chromosomes in the genome at $refseq; you may use \"-p hg19\" or other standard genomes";
        exit 1;
    fi
fi

## Check that refseq exists 
if [ ! -e "$refseq" ]; then
    echo "***! Reference sequence $refseq does not exist";
    exit 1;
fi

## Check that index for refseq exists
if [ ! -e "${refseq}.bwt" ]; then
    echo "***! Reference sequence $refseq does not appear to have been indexed. Please run bwa index on this file before running juicer.";
    exit 1;
fi

## Set ligation junction based on restriction enzyme
if [ -z "$ligation" ]
then
    case $site in
    HindIII) ligation="AAGCTAGCTT";;
    DpnII) ligation="GATCGATC";;
    MboI) ligation="GATCGATC";;
    NcoI) ligation="CCATGCATGG";;
        Arima) ligation="'(GAATAATC|GAATACTC|GAATAGTC|GAATATTC|GAATGATC|GACTAATC|GACTACTC|GACTAGTC|GACTATTC|GACTGATC|GAGTAATC|GAGTACTC|GAGTAGTC|GAGTATTC|GAGTGATC|GATCAATC|GATCACTC|GATCAGTC|GATCATTC|GATCGATC|GATTAATC|GATTACTC|GATTAGTC|GATTATTC|GATTGATC)'" ;;
    none) ligation="XXXX";;
    *)  ligation="XXXX"
        echo "$site not listed as recognized enzyme. Using $sitefile as site file"
        echo "Ligation junction is undefined"
    esac
fi

## If DNAse-type experiment, no fragment maps
if [ "$site" == "none" ]
then
    nofrag=1;
fi

if [ -z "$sitefile" ]
then
    sitefile="${juicerdir}/restriction_sites/${genomeid}_${site}.txt"
fi

## Check that site file exists, needed for fragment number for merged_nodups
if [[ ! -e "$sitefile" ]] && [[ "$site" != "none" ]] &&  [[ ! "$sitefile" =~ "none" ]]
then
    echo "***! $sitefile does not exist. It must be created before running this script."
    exit 1
elif [[ "$site" != "none" ]] && [[ ! "$sitefile" =~ "none" ]]
then
    echo  "Using $sitefile as site file"
fi

## Set threads for sending appropriate parameters to cluster and string for BWA call
if [ ! -z "$threads" ]
then
    threadstring="-t $threads"
else
    threads=10
    #threads="$(getconf _NPROCESSORS_ONLN)"
    threadstring="-t $threads"
fi

sortThead=$( expr $threads - 2 )

if [ "$sortThead" -lt "1" ]; then
  sortThead=1
fi

if [ -z $fastq1 ]; then
  echo "--fastq1: Wrong! Reads NOT found!"
  exit 2
fi

if [ -z $fastq2 ]; then
  echo "--fastq2: Wrong! Reads NOT found!"
  exit 2
fi

if [[ ! -d "$output" ]]; then
  mkdir -p $output
fi

## record time
start=$(date +%s.%N)
output=$(realpath "$output")
## Directories to be created and regex strings for listing files
splitdir=${output}"/splits"
donesplitdir=${output}"/done_splits"
outputdir=${output}"/aligned"
tmpdir=${output}"/HIC_tmp"

## Check that fastq directory exists and has proper fastq files
if [ ! -d "$splitdir" ]; then
  echo "Directory \"$splitdir\" does not exist."
  echo "Create \"$splitdir\" and put fastq files to be aligned there."
  mkdir -p "${splitdir}"
else
  if [[ "$delete" = "true" ]]; then
    echo "Initialization: Delete all files in ${output}!"
    rm -rf "${output}"
    echo "Create folders!"
    mkdir -p "${splitdir}"
  fi
fi

## link fastqs to $splitdir
fastqfile=${splitdir}/*.fastq*
rm -f ${fastqfile}
if [[ $fastq1 =~ '*' ]] && [[ $fastq2 =~ '*' ]]; then
  fastq1=$(for i in "$fastq1";do echo $i; done)
  fastq2=$(for j in "$fastq2";do echo $j; done)

  IFS_OLD=$IFS
  IFS=' ' read -d '' -ra fqArr1 < <(printf '%s\0' "$fastq1")
  IFS=' ' read -d '' -ra fqArr2 < <(printf '%s\0' "$fastq2")
  IFS=$IFS_OLD
else
  IFS_OLD=$IFS
  IFS=',' read -d '' -ra fqArr1 < <(printf '%s\0' "$fastq1")
  IFS=',' read -d '' -ra fqArr2 < <(printf '%s\0' "$fastq2")
  IFS=$IFS_OLD
fi

for i in "${!fqArr1[@]}"
do
  fq1=$(realpath ${fqArr1[$i]})
  fq2=$(realpath ${fqArr2[$i]})
  if [ ! -f $fq1 ]; then
    echo "--fastq1: Wrong! Reads '$fq1' NOT found!"
    exit 2
  fi
  
  if [ ! -f $fq2 ]; then
    echo "--fastq2: Wrong! Reads '$fq1' NOT found!"
    exit 2
  fi

  fqName=$(echo "${fq1}" | perl -pe 's/.fastq(.gz)?$|.fq(.gz)?$//' | perl -pe 's/_R1$|_r1$|_1$//')
  fqName=$(echo $(basename "$fqName"))
  if [[ -z $gzipped ]]; then
    if [ "${fq1: -3}" == ".gz" ]; then
      appendix=".gz"
      gzipped=1
    else
      appendix=""
      gzipped=0
    fi
  fi
  fq1link=${splitdir}/${fqName}_R1.fastq${appendix}
  fq2link=${splitdir}/${fqName}_R2.fastq${appendix}
  ln -s ${fq1} ${fq1link}
  ln -s ${fq2} ${fq2link}
done

## Create output directory, only if not in merge, dedup, final, or postproc stages
if [[ -d "$outputdir" && -z "$final" && -z "$merge" && -z "$dedup" && -z "$postproc" ]] 
then
    echo "***! Move or remove directory \"$outputdir\" before proceeding."
    echo "***! Type \"juicer.sh -h \" for help"
    exit 1            
else
    if [[ -z "$final" && -z "$dedup" && -z "$merge" && -z "$postproc" ]]; then
        mkdir "$outputdir" || { echo "***! Unable to create ${outputdir}, check permissions." ; exit 1; } 
    fi
fi

## Create split directory
if [ -d "$splitdir" ]; then
    splitdirexists=1
else
    mkdir "$splitdir" || { echo "***! Unable to create ${splitdir}, check permissions." ; exit 1; }
fi

## Create temporary directory, used for sort later
if [ ! -d "$tmpdir" ] && [ -z "$final" ] && [ -z "$dedup" ] && [ -z "$postproc" ]; then
    mkdir "$tmpdir"
    chmod 777 "$tmpdir"
fi

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline

if [ $gzipped == "1" ]
then
    read1="${splitdir}"/*${read1str}*.fastq.gz
else
    read1="${splitdir}"/*${read1str}*.fastq
fi

echo "All read-1 fastqs"
ls $read1

headfile=${outputdir}/header
date > $headfile
# Experiment description
if [ -n "${about}" ]
then
    echo -ne 'Experiment description: ${about}; ' >> $headfile
else
    echo -ne 'Experiment description: ' >> $headfile
fi

# Get version numbers of all software
echo -ne "Juicer version $juicer_version;" >> $headfile
bwa 2>&1 | awk '$1=="Version:"{printf(" BWA %s; ", $2)}' >> $headfile 
echo -ne "$threads threads; " >> $headfile
java -version 2>&1 | awk 'NR==1{printf("%s; ", $0);}' >> $headfile 
${juicerdir}/common/juicer_tools -V 2>&1 | awk '$1=="Juicer" && $2=="Tools"{printf("%s; ", $0);}' >> $headfile
echo "$0 $@" >> $headfile

## ALIGN FASTQ AS SINGLE END, SORT BY READNAME, HANDLE CHIMERIC READS

## Not in merge, dedup, final, or postproc stage, i.e. need to align files.
if [ -z $merge ] && [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    echo -e "(-: Aligning files matching $fastqfile\n to genome $genomeid with site file $sitefile"
    if [ ! $splitdirexists ]
    then
        echo "(-: Created $splitdir and $outputdir."
        filename=$(basename "$i")
        filename=${filename%.*}
    #modified ln -s ${fastqfile} ${splitdir}/.
    #mv ${fastqfile} ${splitdir}/.
    else
        echo -e "---  Using already created files in $splitdir\n"
    fi

    ## Loop over all read1 fastq files and create jobs for aligning read1,
    ## aligning read2, and merging the two. Keep track of merge names for final
    ## merge. When merge jobs successfully finish, can launch final merge job.
    for i in ${read1}
    do
        ext=${i#*$read1str}
        name=${i%$read1str*}
        # these names have to be right or it'll break                     
        name1=${name}${read1str}
        name2=${name}${read2str}
        jname=$(basename $name)${ext}
        usegzip=0
        if [ ${ext: -3} == ".gz" ]
        then
            usegzip=1
        fi
        source ${juicerdir}/common/countligations.sh
        if [ -z "$chimeric" ]
            then
            # Align fastq pair
                echo "bwa mem -SP5M $threadstring $refseq $name1$ext $name2$ext > $name$ext.sam"
    
                bwa mem -SP5M $threadstring $refseq $name1$ext $name2$ext > $name$ext.sam
                if [ $? -ne 0 ]
                then
                    echo "***! Alignment of $name1$ext $name2$ext failed."
                    exit 1
                else                                                            
            echo "(-:  Align of $name$ext.sam done successfully"
                fi                                    
            fi                                                              
        ## delete fastq files
        rm -f "$fastqfile"
        # call chimeric_blacklist.awk to deal with chimeric reads; 
        # sorted file is sorted by read name at this point
        touch $name${ext}_abnorm.sam $name${ext}_unmapped.sam  
        awk -v "fname1"=$name${ext}_norm.txt -v "fname2"=$name${ext}_abnorm.sam -v "fname3"=$name${ext}_unmapped.sam -f ${juicerdir}/common/chimeric_blacklist.awk $name$ext.sam
        if [ $? -ne 0 ]
        then
                echo "***! Failure during chimera handling of $name${ext}"
                exit 1
        fi
            # if any normal reads were written, find what fragment they correspond to 
            # and store that
        if [ -e "$name${ext}_norm.txt" ] && [ "$site" != "none" ] && [ -e "$sitefile" ]
        then
                ${juicerdir}/common/fragment.pl $name${ext}_norm.txt $name${ext}.frag.txt $sitefile                                                                
        elif [ "$site" == "none" ] || [ "$nofrag" -eq 1 ] 
        then
                awk '{printf("%s %s %s %d %s %s %s %d", $1, $2, $3, 0, $4, $5, $6, 1); for (i=7; i<=NF; i++) {printf(" %s",$i);}printf("\n");}' $name${ext}_norm.txt > $name${ext}.frag.txt
        else                                                                    
                echo "***! No $name${ext}_norm.txt file created"
                exit 1
        fi                                                                      
        if [ $? -ne 0 ]
        then
                echo "***! Failure during fragment assignment of $name${ext}"
                exit 1
        fi                              
            # sort by chromosome, fragment, strand, and position                    
        sort --parallel $sortThead -T $tmpdir -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n $name${ext}.frag.txt > $name${ext}.sort.txt
        if [ $? -ne 0 ]
        then
                echo "***! Failure during sort of $name${ext}"
                exit 1
        else
                rm $name${ext}_norm.txt $name${ext}.frag.txt
        fi
    done
fi

#MERGE SORTED AND ALIGNED FILES
if [ -z $final ] && [ -z $dedup ] && [ -z $postproc ]
then
    if [ -d $donesplitdir ]
    then
        mv "$donesplitdir"/* $splitdir/.
    fi
    echo "Merging sorted alignments..."
    if ! sort --parallel $sortThead -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n "$splitdir"/*.sort.txt  > $outputdir/merged_sort.txt
    then 
        echo "***! Some problems occurred somewhere in creating sorted align files."
        exit 1
    else
        echo "(-: Finished sorting all sorted files into a single merge."
        rm -r ${tmpdir}
    fi
fi

#REMOVE DUPLICATES
if [ -z $final ] && [ -z $postproc ]
then
    touch ${outputdir}/dups.txt
    touch ${outputdir}/optdups.txt
    touch ${outputdir}/merged_nodups.txt
    echo "Duplicating sorted files ..."
    if [ "$justexact" -eq 1 ]
    then
    awk -f ${juicerdir}/common/dups.awk -v name=${outputdir}/ -v nowobble=1 ${outputdir}/merged_sort.txt
    else
    awk -f ${juicerdir}/common/dups.awk -v name=${outputdir}/ ${outputdir}/merged_sort.txt
    fi
    # for consistency with cluster naming in split_rmdups
    mv ${outputdir}/optdups.txt ${outputdir}/opt_dups.txt 
fi

#STATISTICS

#Skip if post-processing only is required
if [ -z $postproc ]
then        
    echo "Starting post-processing..."
    export _JAVA_OPTIONS=-Xmx16384m
    export LC_ALL=en_US.UTF-8 
    tail -n1 $headfile | awk '{printf"%-1000s\n", $0}' > $outputdir/inter.txt;
    cat "$splitdir"/*.res.txt | awk -f ${juicerdir}/common/stats_sub.awk >> $outputdir/inter.txt
    ${juicerdir}/common/juicer_tools LibraryComplexity $outputdir inter.txt >> $outputdir/inter.txt
    cp $outputdir/inter.txt $outputdir/inter_30.txt 
    ${juicerdir}/common/statistics.pl -s $sitefile -l $ligation -o $outputdir/inter.txt -q 1 $outputdir/merged_nodups.txt 
    cat "$splitdir"/*_abnorm.sam > $outputdir/abnormal.sam
    cat "$splitdir"/*_unmapped.sam > $outputdir/unmapped.sam
    awk -f ${juicerdir}/common/collisions.awk $outputdir/abnormal.sam > $outputdir/collisions.txt
    # Collisions dedupping: two pass algorithm, ideally would make one pass
    gawk -v fname=$outputdir/collisions.txt -f ${juicerdir}/common/collisions_dedup_rearrange_cols.awk $outputdir/collisions.txt | \
      sort -k3,3n -k4,4n -k10,10n -k11,11n -k17,17n -k18,18n -k24,24n -k25,25n -k31,31n -k32,32n | awk -v name=$outputdir/ -f ${juicerdir}/common/collisions_dups.awk
fi

if [ -z "$gsize" ]
then
    #If no path to genome is give, use genome ID as default.
    gsize=$genomeid
fi

#CREATE HIC FILES
# if early exit, we stop here, once the statistics are calculated
if [ -z "$earlyexit" ]
then
    #Skip if post-processing only is required
    if [ -z $postproc ]
    then        
        echo "Starting to call .hic..."
        echo "Running with resolution: $resolution."
        mkdir -p $tmpdir
        if [ "$nofrag" -eq 1 ]
        then 
            ${juicerdir}/common/juicer_tools pre -r "$resolution" -t $tmpdir -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/"${prefix}".inter.hic $gsize
        else 
            ${juicerdir}/common/juicer_tools pre -r "$resolution" -t $tmpdir -f $sitefile -s $outputdir/inter.txt -g $outputdir/inter_hists.m -q 1 $outputdir/merged_nodups.txt $outputdir/"${prefix}".inter.hic $gsize 
        fi 
        ${juicerdir}/common/statistics.pl -s $sitefile -l $ligation -o $outputdir/inter_30.txt -q 30 $outputdir/merged_nodups.txt
        if [ "$nofrag" -eq 1 ]
        then 
            ${juicerdir}/common/juicer_tools pre -r "$resolution" -t $tmpdir -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/"${prefix}".inter30.hic $gsize 
        else 
            ${juicerdir}/common/juicer_tools pre -r "$resolution" -t $tmpdir -f $sitefile -s $outputdir/inter_30.txt -g $outputdir/inter_30_hists.m -q 30 $outputdir/merged_nodups.txt $outputdir/"${prefix}".inter30.hic $gsize
        fi
    fi
    ## POSTPROCESSING
    #${juicerdir}/common/juicer_postprocessing.sh -j ${juicerdir}/common/juicer_tools -i ${outputdir}/inter_30.hic -m ${juicerdir}/references/motif -g ${genomeid}
fi

#CHECK THAT PIPELINE WAS SUCCESSFUL
export early=$earlyexit
export splitdir=$splitdir
source ${juicerdir}/common/check.sh

## delete all mapped sam, txt files
echo "Delete *.sam and non-used txt files except for inter*.txt"
## delete all mapped sam and txt files
rm -rf $tmpdir
rm -f $splitdir/*.sam
rm -f $splitdir/*.sort.txt
rm -f $outputdir/*.sam
find $outputdir -type f -name "*.txt" | grep -v "inter" | xargs -I {} rm -f {}

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
