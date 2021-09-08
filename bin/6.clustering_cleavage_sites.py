import argparse
import os

HELP="""
"""

parser = argparse.ArgumentParser(description=HELP, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-d", "--distance", type=int, default=7,
                    help="The distance between cleavage sites that should be considered as the same site.")
parser.add_argument("-i", "--input_file", type=str, nargs="+",required=True,
                    help="a list of <*.bed> file whose 5th column is the height.")
parser.add_argument("-o", "--out_path", type=str, default="./",
                    help="The output path.")
parser.add_argument("-p", "--p_value", type=float, default=0.05,
                    help="The threshold of significance, should be within (0, 1)")
args = parser.parse_args()

if not os.path.exists(args.out_path):
    os.makedirs(args.out_path)


for i in args.input_file:
    prefix = os.path.split(i)[-1].split('.')[0]
    cmd  = r"""less -S %s |awk '{if($7<%f) print $0}' |bedtools cluster -s -d %d -i /dev/stdin"""%(i, args.p_value, args.distance)
    #print(cmd)
    #1. cluster
    cluster_stat = {}
    gene_stat = {}
    for line in os.popen(cmd).read().strip().split('\n'):
        s = line.strip().split('\t')
        key = 'cluster_' + s[-1]
        if not key in cluster_stat:
            cluster_stat[key] = {'transcript':[], 'sites':[], 'strand':s[5], 'chrom':s[0]}
        cluster_stat[key]['transcript'].append(s[7])
        cluster_stat[key]['sites'].append( [ int(s[1]), int(s[4]), float(s[6]) ] ) #pos, height, pval
        if not s[7] in gene_stat:
            gene_stat[s[7]] = '\t'.join(s[7:10])

    #2. output
    out_name = os.path.join(args.out_path, prefix+'.cluster.bed')
    out = open(out_name, 'w')
    for key in cluster_stat:
        #2.1. multiple transcripts were clustered into same cluster.
        if len(set(cluster_stat[key]['transcript'])) > 1:
            cluster_split = {}
            for tran, site in zip(cluster_stat[key]['transcript'], cluster_stat[key]['sites']):
                if tran not in cluster_split:
                    cluster_split[tran] = []
                cluster_split[tran].append(site)
            index = 1
            for tran in cluster_split:
                positions = [ j[0] for j in cluster_split[tran] ]
                reads = [ j[1] for j in cluster_split[tran]]
                cluster_part = [cluster_stat[key]['chrom'], str(min(positions)), str(max(positions)), key+'.%d'%index, str(sum(reads)), cluster_stat[key]['strand']]
                peak_part = cluster_split[tran][0]
                for j in cluster_split[tran]:
                    if j[1] >= peak_part[1]:
                        peak_part = j[:]
                if peak_part[1] != max( [ j[1] for j in cluster_split[tran] ]):
                    print('Error: choose wrong peak position.')
                    exit()
                out.write( '\t'.join(cluster_part) +'\t'+ '\t'.join(map(str, peak_part)) +'\t'+ gene_stat[tran] +'\n' )
        #2.2 
        else:
            positions = [ j[0] for j in cluster_stat[key]['sites'] ]
            reads = [ j[1] for j in cluster_stat[key]['sites'] ]
            cluster_part = [cluster_stat[key]['chrom'], str(min(positions)), str(max(positions)), key, str(sum(reads)), cluster_stat[key]['strand']]
            peak_part = cluster_stat[key]['sites'][0]
            for j in cluster_stat[key]['sites']:
                if j[1] >= peak_part[1]:
                    peak_part = j[:]
            if peak_part[1] != max( [ j[1] for j in cluster_stat[key]['sites'] ]):
                print('Error: choose wrong peak position.')
                exit()
            out.write( '\t'.join(cluster_part) +'\t'+ '\t'.join(map(str, peak_part)) +'\t'+ gene_stat[cluster_stat[key]['transcript'][0]] +'\n' )
    out.close()




