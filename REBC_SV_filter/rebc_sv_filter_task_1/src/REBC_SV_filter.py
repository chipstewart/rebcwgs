_author__ = 'stewart'
import sys
import argparse
import os
import subprocess
import numpy as np
import pandas as pd
from pandas import DataFrame


if not (sys.version_info[0] == 2  and sys.version_info[1] in [7]):
    raise "Must use Python 2.7.x"

def parseOptions():
    description = '''
    Parse Breakpointer somatic.details.txt file to extract supporting read pairs from tumor bam
    write read table [read_id, chr1 , pos1, str1, chr2, pos2, str2, ...] to output file
    '''

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-d', '--BP_input_file', metavar='BP_input_file', type=str, help='REBC BP details.txt.')
    parser.add_argument('-i', '--id', metavar='id', type=str, help='sample id.',default='')
    parser.add_argument('-o', '--outdir', metavar='outdir', type=str, help='output area', default='.')
    parser.add_argument('-s', '--STUB', metavar='STUB', type=str, help='output file name stub', default='')
    parser.add_argument('-t', '--TALT_thresh', metavar='TALT_thresh', type=str, help='min tumor alt read count threshold.')
    parser.add_argument('-n', '--NALT_thresh', metavar='NALT_thresh', type=str, help='max normal read count threshold.')
    parser.add_argument('-v', '--VAF_thresh', metavar='VAF_thresh', type=str, help='min variant allele fraction threshold.')
    parser.add_argument('-a', '--NALG_thresh', metavar='NALG_thresh', type=str, help='min Algorithm count threshold.')
    parser.add_argument('-b', '--BlackListFile', metavar='BlackListFile', type=str, help='Black list bed file.')

    args = parser.parse_args()

    return args

def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0


if __name__ == "__main__":

    args = parseOptions()
    pair_id = args.id
    BP_input_file = args.BP_input_file
    outdir = args.outdir
    STUB = args.STUB
    if (len(STUB)>0)&~('.' in STUB):
        STUB='.'+STUB
    TALT_thresh = int(args.TALT_thresh)
    NALT_thresh = int(args.NALT_thresh)
    VAF_thresh = float(args.VAF_thresh)
    NALG_thresh = int(args.NALG_thresh)
    BlackListFile = args.BlackListFile

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outFile  = outdir +'/'+pair_id+STUB+'_filtered_SV.tsv'
    #outFP_all = file(outFile_all,'wt')

    BP=pd.read_csv(BP_input_file, sep="\t", index_col=None,low_memory=False,dtype=str) #dRangerDetails(BP_input_file)

    if BP.count==0:
        print 'empty file'
        BP.to_csv(outFile, sep="\t", index=None)
        sys.exit(0)

    BP['pos1'] = BP['pos1'].apply(int)
    BP['pos2'] = BP['pos2'].apply(int)
    BP['VCF_TALT'] = BP['VCF_TALT'].apply(int)
    BP['VCF_TREF'] = BP['VCF_TREF'].apply(int)
    BP['VCF_NALT'] = BP['VCF_NALT'].apply(int)
    BP['VAF'] = BP['VCF_TALT'].apply(float) / ( BP['VCF_TALT']+ BP['VCF_TREF'])
    
    BP['NALG'] = 0*BP['pos1']
    for field in BP.columns:
        if field == 'dRanger':
            BP['NALG'] = BP['NALG'] + (BP['dRanger'].apply(int)>0)*1
 
        if field == 'pcawg_snowman':
            BP['NALG'] = BP['NALG'] + (BP['pcawg_snowman'].apply(int)>0)*1
    
        if field == 'SvABA':
            BP['NALG'] = BP['NALG'] + (BP['SvABA'].apply(int)>0)*1
 
        if field == 'Manta':
            BP['NALG'] = BP['NALG'] + (BP['Manta'].apply(int)>0)*1
    	   
    	        
#    def vaf(a,r):
#        return float(a)/(a + r + 1e-10)
#
#    BP['VAF'] = BP[['VCF_TALT', 'VCF_TREF']].apply(lambda x: vaf(*x), axis=1)
#
    BP0=BP.copy()
    print ('input BP file: %d\n'% BP.individual.count())


    BP=BP[BP['VCF_TALT']>=TALT_thresh]
    print ('BP VCF_TALT>=%d:\t %d\n ' % (TALT_thresh, BP.individual.count()))
    if BP.individual.count()==0:
        BP.to_csv(outFile, sep="\t", index=None)
        sys.exit(0)

    BP=BP[BP['VCF_NALT']<=NALT_thresh]
    print ('BP VCF_NALT<=%d:\t %d\n ' % (NALT_thresh, BP.individual.count()))
    if BP.individual.count()==0:
        BP.to_csv(outFile, sep="\t", index=None)
        sys.exit(0)


    BL=pd.read_csv(BlackListFile, sep="\t", header=None,  names=['chr', 'pos1', 'pos2', 'label'],skip_blank_lines=True,comment='t')

    BP.to_csv('tmp.tsv', sep="\t", index=None)

    #outFP_all.write('%s\n'%dR.headline)

    with open(BlackListFile)as bf:
        for Bline in bf:
            B = Bline.strip().split()
            if "track" in B:
                continue
            #print(Bline.strip())
            chr=B[0].replace('chr','')
            if chr in 'X':
                chr='23'
            if chr in 'Y':
                chr='24'
            pos1=int(B[1])
            pos2=int(B[2])
            BP=BP[((BP.chr1)!=chr)|((BP['chr1']==chr)&((BP['pos1']<pos1)|(BP['pos1']>pos2)))]
            BP=BP[((BP.chr2)!=chr)|((BP['chr2']==chr)&((BP['pos2']<pos1)|(BP['pos2']>pos2)))]
            #print(BP.individual.count())

    print ('BP after Blacklist:\t %d\n ' % (BP.individual.count()))

    BP=BP[BP['NALG']>=NALG_thresh]
    print ('BP NALG<=%d:\t %d\n ' % (NALG_thresh, BP.individual.count()))
    if BP.individual.count()==0:
        BP.to_csv(outFile, sep="\t", index=None)
        sys.exit(0)

    BP=BP[BP['VAF']>=VAF_thresh]
    print ('BP VAF>=%.3f:\t %d\n ' % (VAF_thresh, BP.individual.count()))
    if BP.individual.count()==0:
        BP.to_csv(outFile, sep="\t", index=None)
        sys.exit(0)

    BP.to_csv(outFile, sep="\t", index=None)
    print('\ndone')

