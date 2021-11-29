# rMATS-turbo

rMATS-turbo is the C/Cython version of rMATS (refer to http://rnaseq-mats.sourceforge.net/user_guide.htm). The major difference between rMATS-turbo and rMATS is speed and space usage. The speed of rMATS-turbo is about 100 times faster and the output file is about 1000 times smaller than rMATS. These advantages make analysis and storage of large-scale dataset easy and convenient.

|  | Counting part | Statistical part |
|---------------|----------------|--------------|
| Speed (rMATS-turbo vs rMATS) | 20~100 times faster (one thread) | 300 times faster (6 threads) |
| Storage usage (rMATS-turbo vs rMATS) | 1000 times smaller (PC3E & DS689 dataset) | - |

Important note: rMATS-turbo is multi-threaded while rMATS is not.

Four pre-built binaries are avaliable:

### rMATS-turbo for Linux (Tested on CentOS 6 and Ubuntu 14):
built with Python 2.7.2 on CentOS 6.8.

- rMATS-turbo-Linux-UCS2
- rMATS-turbo-Linux-UCS4

### rMATS-turbo for Mac (Tested on Mac OS X Yosemite 10.10.5):
At the time of compiling this binary, the default C++ compiler on OS X doesn't support OpenMP. In order to enjoy multi-threading feature, GCC 5 is needed.

- rMATS-turbo-Mac-UCS2
- rMATS-turbo-Mac-UCS4

## Prerequisites

- Linux Distribution (CentOS or Ubuntu) or Mac OS X
- Python 2.7.2+
- Numpy
- blas
- lapack
- gsl
- gfortran (Fortran 77 library needed)
- Samtools (Optional)

Installation procedure used in our testing:

- For CentOS 6:
	- pip install numpy
	- yum install lapack-devel blas-devel
	- yum install gsl-devel.x86_64
	- yum install gcc-gfortran

- For Ubuntu 14:
	- pip install numpy
	- sudo apt-get install libblas-dev liblapack-dev
	- sudo apt-get install libgsl0ldbl
	- sudo apt-get install gfortran

- For Mac OS X Yosemite 10.10.5 (Using Homebrew for package management):
	- brew install gcc@5
	- brew install gsl
	- pip install numpy

## Usage

    usage: python rmats.py [options] arg1 arg2
	
    optional arguments:
		  -h, --help            show this help message and exit
		  --version             Version.
		  --gtf GTF             An annotation of genes and transcripts in GTF format.
		  --b1 B1               BAM configuration file.
		  --b2 B2               BAM configuration file.
		  --s1 S1               FASTQ configuration file.
		  --s2 S2               FASTQ configuration file.
		  --od OD               output folder.
		  -t {paired,single}    readtype, single or paired.
		  --libType {fr-unstranded,fr-firststrand,fr-secondstrand}
		                        Library type. Default is unstranded (fr-unstranded).
		                        Use fr-firststrand or fr-secondstrand for strand-
		                        specific data.
		  --readLength READLENGTH
		                        The length of each read.
		  --anchorLength ANCHORLENGTH
		                        The anchor length. (default is 1.)
		  --tophatAnchor TOPHATANCHOR
		                        The "anchor length" or "overhang length" used in the
		                        aligner. At least “anchor length” NT must be
		                        mapped to each end of a given junction. The default is
		                        6. (This parameter applies only if using fastq).
		  --bi BINDEX           The folder name of the STAR binary indexes (i.e., the
		                        name of the folder that contains SA file). For
		                        example, use ~/STARindex/hg19 for hg19. (Only if using
		                        fastq)
		  --nthread NTHREAD     The number of thread. The optimal number of thread
		                        should be equal to the number of CPU core.
		  --tstat TSTAT         the number of thread for statistical model.
		  --cstat CSTAT         The cutoff splicing difference. The cutoff used in the
		                        null hypothesis test for differential splicing. The
		                        default is 0.0001 for 0.01% difference. Valid: 0 ≤
		                        cutoff < 1.
		  --statoff             Turn statistical analysis off.

## Output

Example using fastq.
`
$cat s1.txt:
231ESRP.25K.rep-1.R1.fastq:231ESRP.25K.rep-1.R2.fastq,231ESRP.25K.rep-2.R1.fastq:231ESRP.25K.rep-2.R2.fastq
$cat s2.txt:
231EV.25K.rep-1.R1.fastq:231EV.25K.rep-1.R2.fastq,231EV.25K.rep-2.R1.fastq:231EV.25K.rep-2.R2.fastq
$
`

`
python rMATS-turbo-xxx-UCSx/rmats.py --s1 s1.txt --s2 s2.txt --gtf gtf/Homo_sapiens.Ensembl.GRCh37.72.gtf --bi ~/STARindex/hg19 --od out_test -t paired --nthread 6 --readLength 50 --tophatAnchor 8 --cstat 0.0001 --tstat 6
`

Example using bam.
`
$cat b1.txt:
231ESRP.25K.rep-1.bam,231ESRP.25K.rep-2.bam
$cat b2.txt:
231EV.25K.rep-1.bam,231EV.25K.rep-2.bam
$
`

`
python rMATS-turbo-xxx-UCSx/rmats.py --b1 b1.txt --b2 b2.txt -gtf gtf/Homo_sapiens.Ensembl.GRCh37.75.gtf -od bam_test -t paired --readLength 50 --cstat 0.0001 --libType fr-unstranded
`

All output files are in outputFolder:

- AS_Event.MATS.JC.txt evaluates splicing with only reads that span splicing junctions
	- IC_SAMPLE_1: inclusion counts for SAMPLE_1, replicates are separated by comma
	- SC_SAMPLE_1: skipping counts for SAMPLE_1, replicates are separated by comma
	- IC_SAMPLE_2: inclusion counts for SAMPLE_2, replicates are separated by comma
	- SC_SAMPLE_2: skipping counts for SAMPLE_2, replicates are separated by comma

- AS_Event.MATS.JCEC.txt evaluates splicing with reads that span splicing junctions and reads on target (striped regions on home page figure)
	- IC_SAMPLE_1: inclusion counts for SAMPLE_1, replicates are separated by comma
	- SC_SAMPLE_1: skipping counts for SAMPLE_1, replicates are separated by comma
	- IC_SAMPLE_2: inclusion counts for SAMPLE_2, replicates are separated by comma
	- SC_SAMPLE_2: skipping counts for SAMPLE_2, replicates are separated by comma

- Important columns contained in output files above
	- IncFormLen: length of inclusion form, used for normalization
	- SkipFormLen: length of skipping form, used for normalization
	- P-Value: (The meaning of p value???)
	- FDR: (The meaning of FDR???)
	- IncLevel1: inclusion level for SAMPLE_1 replicates (comma separated) calculated from normalized counts
	- IncLevel2: inclusion level for SAMPLE_2 replicates (comma separated) calculated from normalized counts
	- IncLevelDifference: average(IncLevel1) - average(IncLevel2)

- fromGTF.AS_Event.txt: all possible alternative splicing (AS) events derived from GTF and RNA.

- JC.raw.input.AS_Event.txt evaluates splicing with only reads that span splicing junctions
	- IJC_SAMPLE_1: inclusion junction counts for SAMPLE_1, replicates are separated by comma
	- SJC_SAMPLE_1: skipping junction counts for SAMPLE_1, replicates are separated by comma
	- IJC_SAMPLE_2: inclusion junction counts for SAMPLE_2, replicates are separated by comma
	- SJC_SAMPLE_2: skipping junction counts for SAMPLE_2, replicates are separated by comma
	- IncFormLen: length of inclusion form, used for normalization
	- SkipFormLen: length of skipping form, used for normalization

- JCEC.raw.input.AS_Event.txt evaluates splicing with reads that span splicing junctions and reads on target (striped regions on home page figure)
	- IC_SAMPLE_1: inclusion counts for SAMPLE_1, replicates are separated by comma
	- SC_SAMPLE_1: skipping counts for SAMPLE_1, replicates are separated by comma
	- IC_SAMPLE_2: inclusion counts for SAMPLE_2, replicates are separated by comma
	- SC_SAMPLE_2: skipping counts for SAMPLE_2, replicates are separated by comma
	- IncFormLen: length of inclusion form, used for normalization
	- SkipFormLen: length of skipping form, used for normalization

- bamX_Y STAR mapping result.

## Examples

Suppose we have 4 samples, and we create b1.txt and b2.txt to record these samples.

    $cat b1.txt:
    /xxx/xxx/1.bam,/xxx/xxx/2.bam

    $cat b2.txt:
    /xxx/xxx/3.bam,/xxx/xxx/4.bam

    $python rmats.py --b1 b1.txt --b2 b2.txt --gtf ./human.gtf --od output -t paired --nthread 4 --readLength 76 --anchorLength 1 --tstat 4

## Which version to use

rMATS-turbo was built with two different settings of Python interpreter. In order to know which rMATS-turbo version you should use, you need to check which Unicode type your Python is built with. Open python console and type in:

When built with --enable-unicode=ucs4:

```py
>>> import sys
>>> print sys.maxunicode
1114111
```

Then, download rMATS-turbo-UCS4.tar.gz.

When built with --enable-unicode=ucs2:

```py
>>> import sys
>>> print sys.maxunicode
65535
```

Download rMATS-turbo-UCS2.tar.gz instead.


## FAQ

**ImportError: /xxx/xxx/rmatspipeline.so: undefined symbol: PyUnicodeUCS2_DecodeUTF8**

If you have this error, you are using a different Python interpreter (Unicode characters are stored as UCS-4.) for running your code than the one used to compile rMATS-turbo (Unicode characters are stored as UCS-2.). Here are some common fixes:

- Recompile Python with:

    ./configure --enable-unicode=ucs2

- Use another rMATS-turbo build which is compiled with UCS4. This UCS4 version will always be released along with the UCS2 version.

**ImportError: /xxx/xxx/rmatspipeline.so: undefined symbol: PyUnicodeUCS4_DecodeUTF8**

- Recompile Python with UCS4 or download rMATS-turbo UCS2 version.

**OSError: Permission denied**

Unix and Unix-like systems will not execute a program unless it is marked with permission to execute.

rMATS-turbo has two binary file, rMATS_C/rMATSexe and rmatspipeline.so. In order to run/import these binary programs, you need to grant execute permission to them.

```sh
chmod +x rMATS_C/rMATSexe
```