#!/usr/bin/env python
# -*- coding: utf-8 -*-
##
# @file lite2.py
# @brief 
# @author Zhijie Xie
# @version 1.0.0
# @date 2015-11-27


import os
import shutil
import tempfile
import commands
import argparse
from subprocess import call
from rmatspipeline import run_pipe


VERSION = 'v3.1.0'
USAGE = '''usage: %(prog)s [options] arg1 arg2'''


def doSTARMapping(args): ## do STAR mapping
    fastqs = [args.s1, args.s2,]
    bams = [[], [],]

    for i in range(len(fastqs)):
        if fastqs[i] != '':
            sample = [pair.split(':') for pair in fastqs[i].split(',')]
            print "mapping the first sample"
            for rr, pair in enumerate(sample):
                map_folder = os.path.join(args.tmp, 'bam%d_%d' % (i+1, rr+1));

                if os.path.exists(map_folder):
                    if os.path.isdir(map_folder):
                        os.rmdir(map_folder)
                    else:
                        os.unlink(map_folder)

                os.makedirs(map_folder)
                cmd = 'STAR --chimSegmentMin 2 --outFilterMismatchNmax 3 --alignEndsType EndToEnd --runThreadN 4 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate ';
                cmd += '--alignSJDBoverhangMin ' + str(args.tophatAnchor) + ' --alignIntronMax 299999 --genomeDir ' + args.bIndex + ' --sjdbGTFfile ' + args.gtf; 
                cmd += ' --outFileNamePrefix ' + map_folder + '/ --readFilesIn ';
                cmd += ' '.join(pair)
                status,output = commands.getstatusoutput(cmd);
                print "mapping sample_%d, %s is done with status %s" % (i, ' '.join(pair), status)
                if (int(status)!=0): ## it did not go well
                    print "error in mapping sample_%d, %s: %s" % (i, ' '.join(pair),status)
                    print "error detail: %s" % output
                    raise Exception();
                print output
                bams[i].append(os.path.join(map_folder, 'Aligned.sortedByCoord.out.bam'))

    return ','.join(bams[0]), ','.join(bams[1])
##### end of doSTARMapping ####


def get_args():
    """TODO: Docstring for get_args.
    :returns: TODO

    """
    parser = argparse.ArgumentParser(usage=USAGE)

    group1 = parser.add_mutually_exclusive_group()
    group2 = parser.add_mutually_exclusive_group()

    parser.add_argument('--version', action='version',
                        help='Version.', version=VERSION)
    parser.add_argument('--gtf', action='store',
                        help='An annotation of genes and transcripts in GTF format.', dest='gtf')

    group1.add_argument('--b1', action='store', default='',
                        help='BAM configuration file.', dest='b1')
    group2.add_argument('--b2', action='store', default='',
                        help='BAM configuration file.', dest='b2')
    group1.add_argument('--s1', action='store', default='',
                        help='FASTQ configuration file.', dest='s1')
    group2.add_argument('--s2', action='store', default='',
                        help='FASTQ configuration file.', dest='s2')

    parser.add_argument('--od', action='store',
                        help='output folder of post step.', dest='od')
    parser.add_argument('-t', action='store', default='paired',
                        choices=['paired', 'single',],
                        help='readtype, single or paired.', dest='readtype')
    parser.add_argument('--libType', action='store', default='fr-unstranded',
                        choices=['fr-unstranded', 'fr-firststrand',
                                 'fr-secondstrand',],
                        help='Library type. Default is unstranded (fr-unstranded). Use fr-firststrand or fr-secondstrand for strand-specific data.', dest='dt')
    parser.add_argument('--readLength', action='store', type=int, default=0,
                        help='The length of each read.', dest='readLength')
    parser.add_argument('--anchorLength', action='store', type=int, default=1,
                        help='The anchor length. (default is 1.)', dest='anchorLength')
    parser.add_argument('--tophatAnchor', action='store', type=int, default=6,
                        help='The "anchor length" or "overhang length" used in the aligner. At least “anchor length” NT must be mapped to each end of a given junction. The default is 6. (This parameter applies only if using fastq).', dest='tophatAnchor')
    parser.add_argument('--bi', action='store', default='',
                        help='The folder name of the STAR binary indexes (i.e., the name of the folder that contains SA file). For example, use ~/STARindex/hg19 for hg19. (Only if using fastq)', dest='bIndex')
    parser.add_argument('--nthread', action='store', type=int, default=1,
                        help='The number of thread. The optimal number of thread should be equal to the number of CPU core.', dest='nthread')

    parser.add_argument('--tstat', action='store', type=int, default=1,
                        help='the number of thread for statistical model.', dest='tstat')
    parser.add_argument('--cstat', action='store', type=float, default=0.0001,
                        help='The cutoff splicing difference. The cutoff used in the null hypothesis test for differential splicing. The default is 0.0001 for 0.01%% difference. Valid: 0 ≤ cutoff < 1.', dest='cstat')

    parser.add_argument('--statoff', action='store_false',
                        help='Turn statistical analysis off.', dest='stat')

    args = parser.parse_args()

    args.tmp = tempfile.mkdtemp()
    if args.b1 == '' and args.b2 == '' and args.s1 == '' and args.s2 == '':
        print 'ERROR: BAM/FASTQ required. Please check b1, b2, s1 and s2.'
        exit(0)
    if args.gtf == '' or args.od == '':
        print 'ERROR: GTF file and output folder required. Please check --gtf and --od.'
        exit(0)
    if (args.s1 != '' or args.s2 != '') and args.bIndex == '':
        print 'ERROR: STAR binary indexes required. Please check --bi.'
        exit(0)

    if len(args.b1) > 0:
        with open(args.b1, 'r') as fp:
            args.b1 = fp.read().strip(' ,\n')
    if len(args.b2) > 0:
        with open(args.b2, 'r') as fp:
            args.b2 = fp.read().strip(' ,\n')
    if len(args.s1) > 0:
        with open(args.s1, 'r') as fp:
            args.s1 = fp.read().strip(' ,\n')
    if len(args.s2) > 0:
        with open(args.s2, 'r') as fp:
            args.s2 = fp.read().strip(' ,\n')

    if args.stat and (len(args.b1) * len(args.b2) == 0 and\
                      len(args.s1) * len(args.s2) == 0):
        print 'ERROR: while performing statistical analysis, user should provide two groups of samples. Please check b1,b2 or s1,s2.'
        exit(0)

    if args.b1 == '' and args.b2 == '' and (args.s1 != '' or args.s2 != ''):
        args.b1, args.b2 = doSTARMapping(args)

    args.bams = ','.join([args.b1, args.b2]).strip(',')
    args.junctionLength = 2 * (args.readLength - args.anchorLength)

    dt_map = {'fr-unstranded':0, 'fr-firststrand':1, 'fr-secondstrand':2}
    args.dt = dt_map[args.dt]

    return args


def check_integrity(bams, od):
    """TODO: Docstring for check_integrity.
    :returns: TODO

    """
    vbams = bams.split(',')
    flags = [0 for i in range(len(vbams))]
    all_files = [os.path.join(root, name) for root, dirs, files in os.walk(od)
                 for name in files if name.endswith('.rmats')]

    for name in all_files:
        with open(name, 'r') as fp:
            bams = fp.readline().strip().split(',')
            for i in range(len(bams)):
                idx = vbams.index(bams[i])
                flags[idx] += 1

    flags = [ele != 1 for ele in flags]
    if any(flags):
        print 'WARNING: the number of retrieved BAM files is different from the number of inputed BAM files.'
    else:
        print 'Ok.'


def validate_countfile(fn):
    """TODO: Docstring for validate_countfile.
    :returns: TODO

    """
    with open(fn, 'r+') as fp:
        data = fp.readlines()
        header = data[0]
        rest = [header,]

        for line in data[1:]:
            eles = line.split('\t')
            if len(eles) != 7:
                print 'ERROR: currpted read count file: %s' % (fn)
                exit(0)
            sum_ic_1 = 0
            sum_sc_1 = 0
            sum_ic_2 = 0
            sum_sc_2 = 0
            incv1 = map(int, eles[1].split(','))
            skpv1 = map(int, eles[2].split(','))
            incv2 = map(int, eles[3].split(','))
            skpv2 = map(int, eles[4].split(','))
            inc_len = int(eles[5])
            skp_len = int(eles[6])

            for i in range(len(incv1)):
                sum_ic_1 += incv1[i]
                sum_sc_1 += skpv1[i]
            for i in range(len(incv2)):
                sum_ic_2 += incv2[i]
                sum_sc_2 += skpv2[i]

            if (sum_ic_1 + sum_sc_1 > 0 and\
                    sum_ic_2 + sum_sc_2 > 0 and\
                    (sum_ic_1 != 0 or sum_ic_2 != 0) and\
                    (sum_sc_1 != 0 or sum_sc_2 != 0) and\
                    inc_len != 0 and skp_len != 0):
                rest.append(line)

        fp.seek(0)
        fp.writelines(rest)
        fp.truncate()

    return


def run_stat(istat, tstat, counttype, ase, cstat, od, tmp, stat):
    """TODO: Docstring for run_stat.
    :returns: TODO

    """
    validate_countfile(istat)

    efn = '%s/fromGTF.%s.txt' % (od, ase)
    sec_tmp = os.path.join(tmp, '%s_%s' % (counttype, ase))
    if os.path.exists(sec_tmp):
        if os.path.isdir(sec_tmp):
            os.rmdir(sec_tmp)
        else:
            os.unlink(sec_tmp)
    os.mkdir(sec_tmp)
    ostat = os.path.join(sec_tmp, 'rMATS_result_%s.txt' % ('%s'))

    FNULL = open(os.devnull, 'w')
    resfp = open(ostat % (''), 'w')
    finfp = os.path.join(od, '%s.MATS.%s.txt' % (ase, counttype))

    root_dir = os.path.abspath(os.path.dirname(__file__))
    rmats_c = os.path.join(root_dir, 'rMATS_C/rMATSexe')
    pas_out = os.path.join(root_dir, 'rMATS_P/paste.py')
    inc_lvl = os.path.join(root_dir, 'rMATS_P/inclusion_level.py')
    fdr_cal = os.path.join(root_dir, 'rMATS_P/FDR.py')
    join_2f = os.path.join(root_dir, 'rMATS_P/joinFiles.py')

    if stat:
        call(map(str, [rmats_c, '-i', istat, '-t', tstat, '-o', ostat % ('P-V'), '-c', cstat,]), stdout=FNULL)
        call(['python', pas_out, '-i', istat, '--o1', ostat % ('ID'), '--o2', ostat % ('INP'),], stdout=FNULL)
        call(['python', inc_lvl, ostat % ('INP'), ostat % ('I-L'),], stdout=FNULL)
        call(['python', fdr_cal, ostat % ('P-V'), ostat % ('FDR'),], stdout=FNULL)
        call(['paste', ostat % ('FDR'), ostat % ('I-L'),], stdout=resfp)
        call(['python', join_2f, efn, resfp.name, '0', '0', finfp,], stdout=FNULL)
    else:
        finfp = open(finfp, 'w')
        call(['python', pas_out, '-i', istat, '--o1', ostat % ('ID'), '--o2', ostat % ('INP'),], stdout=FNULL)
        call(['python', inc_lvl, ostat % ('INP'), ostat % ('I-L'),], stdout=FNULL)
        call(['paste', istat, ostat % ('I-L'),], stdout=finfp)
        finfp.close()

    FNULL.close()
    resfp.close()

    return


def main():
    """TODO: Docstring for main.
    :returns: TODO

    """
    args = get_args()

    if not os.path.exists(args.od) or not os.path.isdir(args.od):
        os.mkdir(args.od)

    run_pipe(args)

    jc_it = os.path.join(args.od, 'JC.raw.input.%s.txt')
    jcec_it = os.path.join(args.od, 'JCEC.raw.input.%s.txt')

    print 'Running the statistical part.'
    run_stat(jc_it % ('SE'), args.tstat, 'JC', 'SE', args.cstat, args.od, args.tmp, args.stat)
    run_stat(jcec_it % ('SE'), args.tstat, 'JCEC', 'SE', args.cstat, args.od, args.tmp, args.stat)
    run_stat(jc_it % ('MXE'), args.tstat, 'JC', 'MXE', args.cstat, args.od, args.tmp, args.stat)
    run_stat(jcec_it % ('MXE'), args.tstat, 'JCEC', 'MXE', args.cstat, args.od, args.tmp, args.stat)
    run_stat(jc_it % ('A3SS'), args.tstat, 'JC', 'A3SS', args.cstat, args.od, args.tmp, args.stat)
    run_stat(jcec_it % ('A3SS'), args.tstat, 'JCEC', 'A3SS', args.cstat, args.od, args.tmp, args.stat)
    run_stat(jc_it % ('A5SS'), args.tstat, 'JC', 'A5SS', args.cstat, args.od, args.tmp, args.stat)
    run_stat(jcec_it % ('A5SS'), args.tstat, 'JCEC', 'A5SS', args.cstat, args.od, args.tmp, args.stat)
    run_stat(jc_it % ('RI'), args.tstat, 'JC', 'RI', args.cstat, args.od, args.tmp, args.stat)
    run_stat(jcec_it % ('RI'), args.tstat, 'JCEC', 'RI', args.cstat, args.od, args.tmp, args.stat)
    print 'The statistical part is done.'

    shutil.rmtree(args.tmp)
    print 'Done.'

    return


if __name__ == "__main__":
    main()
