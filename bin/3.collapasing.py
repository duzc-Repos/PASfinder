import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--seq_direct", type=int, default=5, choices=[5, 3],
                   help="Specify sequencing direction based on mRNA.")
parser.add_argument("-i", "--input_file", type=str, nargs="+",required=True,
                    help="A list of <*.sort.bam> file.")
parser.add_argument("-l", "--polyA_length", type=int, default=8,
                    help="The min length(>=1) of marked polyA.")
parser.add_argument("-o", "--out_path", type=str, default="./",
                    help="The output path.")
parser.add_argument("-gs", "--genome_size", type=str, required=True,
                   help="The genome size file.")
args = parser.parse_args()

if not os.path.exists(args.out_path):
    os.makedirs(args.out_path)

chrom_size = { line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in open(args.genome_size)}
for i in args.input_file:
    #0. bam2bed
    prefix = os.path.split(i)[-1].split('.')[0]
    bedname = os.path.join(args.out_path, prefix+'.bed')
    cmd = 'bedtools bamtobed -i %s >%s'%(i, bedname)
    print(cmd)
    os.system(cmd)

    #1. estimate polyA length
    filter_polyA_length_name = os.path.join(args.out_path, prefix+'.filter_polyA_length.bed')
    if args.seq_direct == 5:
        with open(filter_polyA_length_name, 'w') as out:
            for line in open(bedname):
                s = line.strip().split('\t')
                if s[0] not in chrom_size:
                    continue
                if s[0].upper() in ['CHRM', 'M', 'MT']:
                    continue
                if s[3].find('_pA=') >= 0 and int(s[3].split('_pA=')[-1]) >= args.polyA_length:
                    out.writelines(line)
    else:
        with open(filter_polyA_length_name, 'w') as out:
            for line in open(bedname):
                s = line.strip().split('\t')
                if s[0] not in chrom_size:
                    continue
                if s[0].upper() in ['CHRM', 'M', 'MT']:
                    continue
                out.writelines(line)

    #2. generate putative C/S site files
    cs_file_name = os.path.join(args.out_path, prefix+'.pCS.bed')
    if args.seq_direct == 5:
        cmd = r"""bedtools genomecov -strand + -dz -3 -i %s -g %s |awk 'BEGIN{OFS="\t"}{print $1, $2, $2+1, "pCS_plus_"NR, $3, "+"}' >%s """%(filter_polyA_length_name, args.genome_size, cs_file_name)
        print(cmd)
        os.system(cmd)
        cmd = r"""bedtools genomecov -strand - -dz -3 -i %s -g %s |awk 'BEGIN{OFS="\t"}{print $1, $2, $2+1, "pCS_minus_"NR, $3, "-"}' >>%s """%(filter_polyA_length_name, args.genome_size, cs_file_name)
        print(cmd)
        os.system(cmd)
    else:
        cmd = r"""bedtools genomecov -strand + -dz -5 -i %s -g %s |awk 'BEGIN{OFS="\t"}{print $1, $2, $2+1, "pCS_plus_"NR, $3, "-"}' >%s """%(filter_polyA_length_name, args.genome_size, cs_file_name)
        print(cmd)
        os.system(cmd)
        cmd = r"""bedtools genomecov -strand - -dz -5 -i %s -g %s |awk 'BEGIN{OFS="\t"}{print $1, $2, $2+1, "pCS_minus_"NR, $3, "+"}' >>%s """%(filter_polyA_length_name, args.genome_size, cs_file_name)
        print(cmd)
        os.system(cmd)

    #3. 
    os.remove(filter_polyA_length_name)

