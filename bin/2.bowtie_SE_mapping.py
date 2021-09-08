import os
import argparse

parser = argparse.ArgumentParser(description="Align reads to genome using bowtie.\n Requie: bowtie, samtools (1.6) ")
parser.add_argument("-in", "--index", type=str, required=True,
                   help="The path to bowtie index")
parser.add_argument("-p", "--core", type=int, default=4,
                   help="number of alignment threads to launch, detail in BOWTIE.")
parser.add_argument("-i", "--input_file", type=str, required=True, nargs="+",
                   help="A list of input file in <*.fastq> or <*.fq.gz> format.")
parser.add_argument("-o", "--output_path", type=str, default='mapping',
                   help="he path of output files.")

args = parser.parse_args()
if not os.path.exists(args.output_path):
    os.makedirs(args.output_path)
if not os.popen("which bowtie"):
    print("Require: bowtie")
    exit()
if not os.popen("which samtools"):
    print("Reqire: samtools, v1.6 is recommanded.")
    exit()

for i in args.input_file:
    prefix = os.path.split(i)[-1].split('.')[0]
    out_name = os.path.join(args.output_path, prefix )

    cmd = 'less -S %s |bowtie -p %d -n 2 -l 25 -m 1 --best --strata \
    --chunkmbs 240 %s - -S /dev/stdout 2>%s.log |samtools sort -@ %d -o \
    %s.sort.bam'%(i, args.core, args.index, out_name, args.core, out_name)
    print(cmd)
    os.system(cmd)


