import os
import argparse

HELP="""
###############################################################################
-------------------------------------------------------------------------------
This script is to help deal with internal primed events.
-------------------------------------------------------------------------------
Most of the 3'-enriched RNA-seq approach are based on oligo d(T) primer to
enrich the sequence at mRNA 3' end, so if there a A-rich region in the middle
of the mRNA, it also has a potential to be enriched and cause false
positive of cleavage site identification. Therefore, we should filter these
events before cleavage site identification.

#################################
# The workflow is as follow:
#################################
    a. extract upstream 50bp sequence and find polyA signal(PAS).
    b. extract downstream 20bp sequence and find A rich region. (6 consecutive
        As or 7 As in 10bp window)
###############################################################################
"""

def check_requirement(args):
    if int(os.popen("bedtools --version").read().split(' ')[-1][1]) != 2:
        print('Require software: bedtools v2.*')
        exit()
    if not os.path.exists(args.out_path):
        os.makedirs(args.out_path)
    pass

def get_sequence(seq_direct, strand, genome, region, basecomplement):
    cmd = "samtools faidx %s %s"%(genome, region)
    seq = os.popen(cmd).read().strip().split('\n')
    seq = ''.join(seq[1:])
    if seq_direct == 5:
        seq = seq if strand == '+' else ''.join([basecomplement.get(i, 'N') for i in seq[::-1]])
    elif seq_direct == 3:
        seq = ''.join([basecomplement.get(i, 'N') for i in seq[::-1]]) if strand == '+' else seq
    return seq

def mark_internal_primed_event(sequence):
    flag = False
    if sequence.upper().find('AAAAAA') >= 0:
        flag = True
    else:
        for i in range(len(sequence)-10):
            if sequence[i:i+10].upper().count('A') >= 7:
                flag = True
                break
    return str(flag)

def mark_polyA_signal(sequence, polyA_signal):
    signal = False
    for i in polyA_signal:
        if sequence.upper().find(i) >= 0:
            signal = i
            break
    return str(signal)

parser = argparse.ArgumentParser(description=HELP, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-m", "--method", type=str, default='heuristic', choices=['heuristic', 'cleanUpdTseq'],
                    help="The method used to filter the internal primed events.")
parser.add_argument("-i", "--input_file", type=str, nargs="+",required=True,
                    help="a list of <*.sort.bam> file.")
parser.add_argument("-g", "--genome", type=str, required=True,
                    help="The <genome.fasta> file")
parser.add_argument("-o", "--out_path", type=str, default="./",
                    help="The output path.")
parser.add_argument("-d", "--downstream", type = int, default=20,
                    help="Extract #N bp window after cleavage site.")
parser.add_argument("-u", "--upstream", type = int, default=50,
                    help="Extract #N bp window before cleavage site.")
args = parser.parse_args()
check_requirement(args)

#basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
polyA_signal = ["AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AAGAAA", "AATATA", "AATACA", "CATAAA", "GATAAA", "AATGAA",
                "TTTAAA", "ACTAAA", "AATAGA" ]

for i in args.input_file:
    print('\n###############%s'%(i))
    prefix = os.path.split(i)[-1].split('.')[0]
    if args.method == 'heuristic':
        #1. get upstream and downstream region <*.bed>
        up_region_file = os.path.join(args.out_path, prefix+'.up_region.tab')
        down_region_file = os.path.join(args.out_path, prefix+'.down_region.tab')
        cmd = r"""less -S %s |awk 'BEGIN{OFS="\t"}{if(($3-%d) >= 0 && $6=="+") print $1, $2-%d, $2, $4, $5, $6; else if($6=="-") print $1, $2, $2+%d, $4, $5, $6}' |fastaFromBed -s -name -tab -fi %s -bed - >%s """%(i, args.upstream, args.upstream, args.upstream, args.genome, up_region_file)
        print(cmd)
        os.system(cmd)
        cmd = r"""less -S %s |awk 'BEGIN{OFS="\t"}{if($6=="+") print $1, $2, $2+%d, $4, $5, $6; else if(($2-%d)>=0 && $6=="-") print $1, $2-%d, $2, $4, $5, $6}' |fastaFromBed -s -name -tab -fi %s -bed - >%s """%(i, args.downstream, args.downstream, args.downstream, args.genome, down_region_file)
        print(cmd)
        os.system(cmd)

        #2. get internal primed events flag and polyA signal flag
        flag_internal = {}
        for line in open(down_region_file):
            s = line.strip().split('\t')
            key = s[0].split('::')[0]
            flag_internal[key] = mark_internal_primed_event(s[-1])
        flag_polyA_signal = {}
        for line in open(up_region_file):
            s = line.strip().split('\t')
            key = s[0].split('::')[0]
            flag_polyA_signal[key] = mark_polyA_signal(s[-1], polyA_signal)
            #print(key, polyA_signal[key])

        #3. annotate internal primed events and polyA signal
        outname = os.path.join(args.out_path, prefix+'.heuristic.bed')
        out = open(outname, 'w')
        for line in open(i, 'r'):
            s = line.strip().split('\t')
            if flag_internal.get(s[3], "True") == "True" and flag_polyA_signal.get(s[3], "False") == "False":
                continue
            out.writelines(line)
        out.close()

        #3.4 remove internal produced files
        os.system('rm %s %s'%(up_region_file, down_region_file))
    else:
        #1. get upstream and downstream region <*.bed>
        up_region_file = os.path.join(args.out_path, prefix+'.up_region.tab')
        down_region_file = os.path.join(args.out_path, prefix+'.down_region.tab')
        cmd = r"""less -S %s |awk 'BEGIN{OFS="\t"}{if(($3-%d) >= 0 && $6=="+") print $1, $2-%d, $2, $4, $5, $6; else if($6=="-") print $1, $2, $2+%d, $4, $5, $6}' |fastaFromBed -s -name -tab -fi %s -bed - >%s """%(i,40, 40, 40, args.genome, up_region_file)
        print(cmd)
        os.system(cmd)
        cmd = r"""less -S %s |awk 'BEGIN{OFS="\t"}{if($6=="+") print $1, $2, $2+%d, $4, $5, $6; else if(($2-%d)>=0 && $6=="-") print $1, $2-%d, $2, $4, $5, $6}' |fastaFromBed -s -name -tab -fi %s -bed - >%s """%(i, 30, 30, 30, args.genome, down_region_file)
        print(cmd)
        os.system(cmd)

        #2. merge into one file
        up = { line.strip().split('\t')[0].split('::')[0]:line.strip().split('\t')[-1] for line in open(up_region_file)}
        down = { line.strip().split('\t')[0].split('::')[0]:line.strip().split('\t')[-1] for line in open(down_region_file)}
        add_seq_bed = os.path.join(args.out_path, prefix+'.addseq.bed')
        with open(add_seq_bed, 'w') as out:
            for line in open(i):
                key = line.strip().split('\t')[3]
                if (key in up) and (key in down):
                    out.writelines(line.strip() + '\t' + up[key] + '\t' + down[key] + '\n')

        #3. cleanUpdTSeq
        non_internal_list = os.path.join(args.out_path, prefix+'.cleanUpdTSeq.list')
        cmd = 'Rscript cleanUpdTSeq.r %s %s'%(add_seq_bed, non_internal_list)
        print(cmd)
        os.system(cmd)

        flag_key = { line.strip():1 for line in open(non_internal_list, 'r') }
        outname = os.path.join(args.out_path, prefix+'.cleanUpdTSeq.bed')
        with open(outname, 'w') as out:
            for line in open(i):
                key = line.strip().split('\t')[3]
                if flag_key.get(key, 0) == 1:
                    out.writelines(line)

        os.system('rm %s %s %s %s'%(up_region_file, down_region_file, add_seq_bed, non_internal_list) )

