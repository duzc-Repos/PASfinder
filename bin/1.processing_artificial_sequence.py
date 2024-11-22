import argparse
import os
import gzip

HELP="""
###############################################################################
-------------------------------------------------------------------------------
This is a script to process the artificial sequence at both ends of reads.
-------------------------------------------------------------------------------

In general, 3'-enriched RNA-seq can be classified two catagories:
    1. sequencing from 5' to 3' (mRNA);
    2. sequencing from 3' to 5' (mRNA).

The artificial sequence can be divided into two catagories:
    1. expected:
        a. UMI: to quantify mRNA and remove PCR duplicate more precisely.
        b. barcode: to classify which sample the read is from.
        c. linker: some known sequence to verify the correct sequencing start
                   or end.
    2. unexpected:
        a. adaptor: the insert fragment is too short
        b. A rich at 3' end of reads: the fragment is so short that the read is readthrough.
        c. no insert fragmrnt: the primer self cicurlization.
Therefore, we should mark expected artificial sequence and trim unexpexted sequence.

#################################
# The workflow is as follow:
#################################
    a. mark 5prime sequence.                       (optional)
    b. trim 3' adaptor.                            (required)
    c. remove high density A/T reads.              (required)
    d. trim polyA from 3' or trim polyT from 5'.   (required)
    e. mark polyA length.                          (5' only)

###############################################################################

"""

def check_requirement(args):
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    if not os.popen('which cutadapt'):
        exit('Require: cutadapt')
    iupac = ['A', 'C', 'T', 'U', 'G', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
    for i, letter in enumerate(args.adaptor_sequence):
        if letter == '^' and i != 0:
            print('Require adaptor: "^" must be at the beginning of the adaptor.')
            exit()
        elif letter == '$' and i != (len(adaptor) - 1):
            print('Require adaptor: "$" must be at the end of the adaptor.')
            exit()
        elif letter == '.' and ( i == 0 or i == (len(adaptor) - 1) ):
            print('Require adaptor: "." must be at the middle of the adaptor.')
            exit()
        elif letter.upper() not in iupac:
            print('Require adaptor: letter must be in %s.'%(','.join(iupac)) )
            exit()

def mark_5prime_sequence(input_file, tag_length, outdir):
    """
    Func: trim the artificial sequence like UMI, barcode, linker at 5' end and mark the sequence.
    Parm input_file: str, the raw <.fq> file.
    Parm tag_length: int, the length of the tag, usually the first #N base.
    Parm outdir: str, the output path.
    return: the file after trim and mark.
    """

    def process_compress_file(infile, tag_length, outdir, out_name):
        out = open(out_name, 'w')
        file = gzip.open(infile, 'r')
        for line in file:
            line = str(line, encoding = "utf-8")
            if line.startswith('@'):
                name = line.strip().split(' ')[0]
                seq = str(file.readline(), encoding = "utf-8")
                name2 = str(file.readline(), encoding = "utf-8")
                qual = str(file.readline(), encoding = "utf-8")
            name = name +"_5endSeq=%s\n"%(seq[:tag_length])
            out.writelines([name, seq[tag_length:], "+\n", qual[tag_length:]])
        out.close()
        return out

    def process_plain_file(infile, tag_length, outdir, out_name):
        out = open(out_name, 'w')
        file = open(infile, 'r')
        for line in file:
            if line.startswith('@'):
                name = line.strip().split(' ')[0]
                seq = file.readline()
                name2 = file.readline()
                qual = file.readline()
            name = name +"_5endSeq=%s\n"%(seq[:tag_length])
            out.writelines([name, seq[tag_length:], "+\n", qual[tag_length:]])
        out.close()
        return out

    out_name = os.path.join(outdir, "%s.5prime.fq"%(os.path.split(input_file)[-1]))
    if input_file.endswith('.gz'):
        out = process_compress_file(input_file, tag_length, outdir, out_name)
    else:
        out = process_plain_file(input_file, tag_length, outdir, out_name)
    return out_name

def cut_adaptor(infile, adaptor_sequence, m, outdir, flag):
    """
    Func: trim adaptor at reads 3' end.
    Parm infile: str, input file in <*.fq> format.
    Parm adaptor_sequence: str, the potential adaptor sequence.
    Parm m: int, the min reads length after trim. detail in cutadapt
    Parm outdir: str, the output path.
    Parm flag: int, mark sequence direction.
    return: the file name after trim.
    """
    out = os.path.join(outdir, "%s.rmAdaptor.fq"%(os.path.split(infile)[-1]))
    if flag == 5:
        cmd = "cutadapt -q 30,0 -a %s -m %d --max-n 0.25 --trim-n --quiet -o %s %s"%(adaptor_sequence, m, out, infile)
    elif flag == 3:
        cmd = "cutadapt -q 0,30 -a %s -m %d  -e 0.2 --max-n 0.25 --trim-n --quiet -o %s %s"%(adaptor_sequence, m, out, infile)
    print("     ", cmd)
    os.system(cmd)
    return out


def remove_low_complexity_read(input_file, cutoff, outdir):
    """
    Func: remove low complex reads.
    Parm input_file: str, input_file in <*.fq> format.
    Parm cutoff: float, the cutoff of A/T/C/G in a reads.
    Parm outdir: str, the output path.
    return: the file name after remove.
    """

    out_name = os.path.join(outdir, "%s.lowComplexity.fq"%(os.path.split(input_file)[-1]))
    out = open(out_name, 'w')
    file = open(input_file, 'r')
    for line in file:
        if line.startswith('@'):
            name = line.strip().split(' ')[0]
            seq = file.readline().strip()
            name2 = file.readline().strip()
            qual = file.readline().strip()
        ratios = [ seq.count(i)/len(seq) for i in ['A', 'T', 'C', 'G'] ]
        if max(ratios) < cutoff:
            out.writelines([name+'\n', seq+'\n', '+\n', qual+'\n'])
    out.close()
    return out_name

def trim_polyA_or_polyT_tail(input_file, m, p, outdir, flag, base, overlap, flag_discard):
    """
    Func: trim poly(A) or poly(T) reads.
    Parm input_file: str, input file in <*.fq> format.
    Parm m: int, the min reads length after trim. detail in cutadapt.
    Parm p: int, The length of poly(A) or poly(T) sequence used for cutadapt.
    Parm outdir: str, the output path.
    Parm flag: int, mark sequence direction.
    Parm base: str, A or T base.
    Parm overlap: int, the min overlap between reads and adaptor.
    Parm flag_discard: [0, 1], whether retain reads without polyT.
    return: the file name after trim.
    """

    out_name = os.path.join(outdir, "%s.trimPoly%s.fq"%(os.path.split(input_file)[-1], base))
    if flag == 5:
        cmd = 'cutadapt -a A{%d} -m %s -n 2 -e 0.1 --quiet -O %d -o %s %s'%(p, m, overlap, out_name, input_file)
    elif flag == 3:
        a = 'T'*p
        if flag_discard:
            cmd = 'cutadapt -g %s -m %s -n 2 -e 0.2 --discard-untrimmed --quiet -O %d -o %s %s'%(a, m, overlap, out_name, input_file)
        else:
            cmd = 'cutadapt -g %s -m %s -n 2 -e 0.2 --quiet -O %d -o %s %s'%(a, m, overlap, out_name, input_file)
    print("     ", cmd)
    os.system(cmd)
    return out_name

def mark_polyA_lens(raw_file_name, trim_file_name, window_size, ratio, outdir):
    """
    Func: mark polyA lens of a read.
    Parm raw_file_name: str, <*.fq> file before trim polyA.
    Parm trim_file_name: str, <*.fq> file after trim polyA.
    Parm window_size: int, the window size to judge polyA lens.
    Parm ratio: float, the ratio of A in a read.
    Parm outdir: str, the output path.
    return: the file name after mark the polyA length.
    """

    def calc_polyA_len(sequence, window_size ,ratio):
        if len(sequence) <= window_size:
            pA_len = sequence.rindex('A') + 1
        else:
            pA_len = 0
            for i in range(len(sequence) - window_size + 1):
                if (sequence[i:i+window_size].upper().count('A') / window_size) >= ratio:
                    pA_len += 1
                else:
                    break
            if (i+window_size) == len(sequence):
                pA_len += sequence[i:i+window_size].rindex('A')
            else:
                pA_len += sequence[i:i+window_size].rindex('A') + 1
        return pA_len

    out_name = os.path.join(outdir, "%s.clean.fq"%(os.path.split(trim_file_name)[-1].split('.')[0]))
    out = open(out_name, 'w')
    with open(raw_file_name, 'r') as raw_file, open(trim_file_name, 'r') as trim_file:
        raw_name = raw_file.readline().strip()
        raw_seq = raw_file.readline().strip()
        raw_name2 = raw_file.readline().strip()
        raw_qual = raw_file.readline().strip()
        trim_name = trim_file.readline().strip()
        trim_seq = trim_file.readline().strip()
        trim_name2 = trim_file.readline().strip()
        trim_qual = trim_file.readline().strip()
        while trim_name:
            if trim_name == raw_name:
                sequence = raw_seq[len(trim_seq):]
                if len(sequence) > 0:
                    pA_lens = calc_polyA_len(sequence, window_size, ratio)
                    out.writelines([trim_name + '_pA=%d\n'%(pA_lens), trim_seq+"\n", "+\n", trim_qual+"\n"])
                else:
                    out.writelines([trim_name+"\n", trim_seq+"\n", "+\n", trim_qual+"\n"])

                trim_name = trim_file.readline().strip()
                trim_seq = trim_file.readline().strip()
                trim_name2 = trim_file.readline().strip()
                trim_qual = trim_file.readline().strip()
            else:
                raw_name = raw_file.readline().strip()
                raw_seq = raw_file.readline().strip()
                raw_name2 = raw_file.readline().strip()
                raw_qual = raw_file.readline().strip()
    out.close()
    return out_name

parser = argparse.ArgumentParser(description=HELP, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i", "--input_file", type=str, required=True, nargs="+",
                    help="A list of input file in <*.fastq> or <*.fq.gz> format.")
parser.add_argument("-o", "--output_dir", type=str, default="result",
                    help="The path of output files.")
parser.add_argument("-s", "--seq_direct", type=int, default=5, choices=[5, 3],
                    help="Specify sequencing direction based on mRNA.")
parser.add_argument("-v", "--version", action='version', version='%(prog)s v3.0')
a_step_parser = parser.add_argument_group("#A. MARK 5' PRIME SEQUENCE")
a_step_parser.add_argument("-l", "--prime_length", default=0, type=int,
                    help="the length of expexted artificial sequence, like UMI, barcode, etc.")
b_step_parser = parser.add_argument_group("#B. TRIM 3' ADAPTOR")
b_step_parser.add_argument("-a", "--adaptor_sequence", type=str, required=True,
                           help="The 3' adapter sequence, detail in cutadapt. (illumina:agatcggaagagc)")
b_step_parser.add_argument("-m", "--min_read_length", type=int, default=18,
                           help="The minimum read length, detail in cutadapt.")
c_step_parser = parser.add_argument_group("#C. REMOVE LOW COMPLEXITY READS")
c_step_parser.add_argument("-c", "--cutoff", type=float, default=0.8,
                           help="The cutoff of A/T/C/G ratio in a reads.")
d_step_parser = parser.add_argument_group("#D. TRIM POLY(A) FROM 3' END OR TRIM POLY(T) FROM 5' END")
d_step_parser.add_argument("-p", "--poly_length", type=int, default=15,
                           help="The length of poly(A) or poly(T) sequence used for cutadapt.")
d_step_parser.add_argument("-ol", "--overlap", type=int, default=3,
                           help="The min overlap between reads and adaptor.")
d_step_parser.add_argument("-du", "--discard_untrimmed", type=int, default=0, choices=[0, 1],
                          help="Whether discard reads that without polyT. 1 means YES, 0 means NO.")
e_step_parser = parser.add_argument_group("#E. MARK POLY(A) LENGTH")
e_step_parser.add_argument("-w", "--window_size", type=int, default=10,
                          help="The window size to judge polyA lens.")
e_step_parser.add_argument("-r", "--ratio", type=float, default=0.9,
                          help="The ratio of A in a window size.")
args = parser.parse_args()
check_requirement(args)

base = {5:'A', 3:'T'}
for i in args.input_file:
    outfile_name = i
    print("\n####PROCESSING FILE: %s"%(outfile_name))
    if args.prime_length > 0:
        print("------a. mark 5' prime sequence ...")
        trim_prime = mark_5prime_sequence(outfile_name, args.prime_length, args.output_dir)
        print(outfile_name)
    else:
        trim_prime = outfile_name
    if args.adaptor_sequence:
        print("------b. trim 3' adaptor using cutadapt ...")
        trim_adaptor = cut_adaptor(trim_prime, args.adaptor_sequence, args.min_read_length, args.output_dir, args.seq_direct)
    else:
        trim_adaptor = trim_prime
    if args.cutoff < 1:
        print("------c. remove low complexity reads ...")
        remove_high_density = remove_low_complexity_read(trim_adaptor, args.cutoff, args.output_dir)
    else:
        remove_high_density = trim_adaptor
    if args.poly_length > 0:
        print("------d. trim poly(%s) reads ..."%(base[args.seq_direct]))
        trim_poly = trim_polyA_or_polyT_tail(remove_high_density, args.min_read_length, args.poly_length, args.output_dir, args.seq_direct, base[args.seq_direct], args.overlap, args.discard_untrimmed)
    else:
        trim_poly = remove_high_density
    if args.seq_direct == 5:
        outname = mark_polyA_lens(remove_high_density, trim_poly, args.window_size, args.ratio, args.output_dir)
    elif args.seq_direct == 3:
        prefix = os.path.split(trim_poly)[-1].split('.')[0]
        outname = os.path.join(args.output_dir, prefix+".clean.fq")
        cmd = 'mv %s %s'%(trim_poly, outname)
        print("     ", cmd)
        os.system(cmd)
    cmd = 'gzip -f %s'%(outname)
    print("     ", cmd)
    os.system(cmd)
    cmd = 'rm %s'%(os.path.join(args.output_dir, "*.fq"))
    print("     ", cmd)
    os.system(cmd)

