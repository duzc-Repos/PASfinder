import sys

HELP = """
This script to prepare annotation for pipeline.
This script takes an annotaion file in GENCODE or ENSEMBL <*.gtf> format and outputs a modified <*.bed> file.
Note:
    The input file must contain "transcript", "CDS" and "UTR" features at the third coloumn.
    For annotation files without "UTR" annotated, we regarded the last exon/CDS as the UTR. (example. yeast)
"""

import argparse
import os

parser = argparse.ArgumentParser(description=HELP, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", "--input_file", type=str, required=True,
                   help="an annotation file in <*.gtf> format.")
parser.add_argument("-o", "--output_file", type=str, required = True,
                   help="the name of output file.")
parser.add_argument("-e", "--last_cds", type=bool, default=False,
                   help="regard the last exon(cds) as 3'UTR. (For species without UTR annotated like yeast.)")
args = parser.parse_args()

trans = {}
for line in open(args.input_file, 'r'):
    if line.startswith('#'):
        continue

    s = line.strip().split('\t')
    infor = { i.strip().split(' ')[0]:i.strip().split(' ')[-1].strip('"') for i in s[-1].split(';') if i != ''}
    if s[2].upper().find('UTR') >= 0:
        s[2] = 'UTR'
    if s[2].upper() == "TRANSCRIPT":
        trans[ infor['transcript_id'] ] = {'chrom':s[0], 'region':s[3:5], 'strand':s[6],
                                           'anno':[ infor[i] for i in ['gene_id', 'transcript_id', 'gene_name', 'gene_type', 'gene_biotype'] if infor.get(i, '') != '' ],
                                           'CDS':[], 'UTR':[]}
    elif s[2].upper() in ['CDS', 'UTR']:
        if infor['transcript_id'] not in trans:
            print("ERROR: the 'transcript' feature of %s don't annotate before CDS or UTR, please check your <.gtf> file."%(s[2]) )
            exit()

        trans[ infor['transcript_id'] ][s[2].upper()].append([int(s[3]), int(s[4])])
print(len(trans))


out = open('primary_transcript.bed', 'w')
for tran in trans:
    if len(trans[tran]['CDS']) == 0:
        continue

    cds_sort = sorted(trans[tran]['CDS'])
    if args.last_cds:
        utr = [trans[tran]['CDS'][-1]] if trans[tran]['strand'] == '+' else [trans[tran]['CDS'][0]]
        out.write( '\t'.join( [ trans[tran]['chrom'],                                                      #1. chrom
                                str(int(trans[tran]['region'][0])-1)+'\t'+trans[tran]['region'][1],        #2. start, end
                                '|'.join(trans[tran]['anno']) + '\t0',                                     #3. name, col5
                                trans[tran]['strand'],                                                     #4. strand
                                '|'.join([ ','.join(map(str, i)) for i in utr])                            #5. UTR region
                              ]) + '\n')
        continue

    if trans[tran]['strand'] == '+':
        utr = [ [i[0]-1, i[1]] for i in trans[tran]['UTR'] if i[0] >= cds_sort[-1][-1] ]
    elif trans[tran]['strand'] == '-':
        utr = [ [i[0]-1, i[1]] for i in trans[tran]['UTR'] if i[-1] <= cds_sort[0][0] ]
    utr = sorted(utr)
    if len(utr) == 0:
        continue
    out.write( '\t'.join( [ trans[tran]['chrom'],                                                      #1. chrom
                            str(int(trans[tran]['region'][0])-1)+'\t'+trans[tran]['region'][1],        #2. start, end
                            '|'.join(trans[tran]['anno']) + '\t0',                                     #3. name, col5
                            trans[tran]['strand'],                                                     #4. strand
                            '|'.join([ ','.join(map(str, i)) for i in utr])                            #5. UTR region
                          ]) + '\n')
out.close()

cmd = 'less -S primary_transcript.bed |sortBed |bedtools cluster -s -i - >primary_transcript.cluster.bed'
os.system(cmd)

trans = {}
for line in open('primary_transcript.cluster.bed'):
    s = line.strip().split('\t')
    length = int(s[2])-int(s[1])
    if not s[-1] in trans:
        trans[ s[-1] ] = [length, '\t'.join(s[:-1])]
    elif length > trans[s[-1]][0]:
        trans[ s[-1] ] = [length, '\t'.join(s[:-1])]
print(len(trans))

cmd = 'rm primary_transcript.bed primary_transcript.cluster.bed'
os.system(cmd)

with open(args.output_file, 'w') as out:
    for line in trans:
        out.write(trans[line][1]+'\n')

