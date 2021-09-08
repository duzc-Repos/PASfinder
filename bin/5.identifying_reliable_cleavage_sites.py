import argparse
import os
import numpy as np
import scipy as sp
import scipy.optimize as opt
import multiprocessing as mp
import datetime


HELP = """
###############################################################################
-------------------------------------------------------------------------------
This script is to identify reliable cleavage sites.
-------------------------------------------------------------------------------
###############################################################################
"""

def gaussian_kernel(x, u, b):
    """
    f(x) = exp( -(x-u)^2/(2*b^2) )
    """
    return np.exp( ((x-u)**2) / (-2*(b**2)) )

def background_model(args):
    """
    Func: using gaussain kernel to construct background model.
    Parm args = [gene_id, tags, length, bandwitdth, options, threhold, seeds]
         gene_id: str,
         tags: int, the number of reads in 3'UTR.
         length: int, the length of 3'UTR.
         bandwidth: float,
         options: dict, options for SLSQP method.
         seeds: list, seed for random number generation.
    return: list or None, the simulated distribution of the background.
    """

    gene_id, tags, length, bandwidth, options, seeds = args
    if tags == 0:
        return [gene_id, None]
    if tags == 1:
        return [gene_id, [1]*1000]
 
    res = []
    for seed in seeds:
        np.random.seed(seed)
        random_u = np.random.randint(length, size=tags)
        fun = lambda x: -gaussian_kernel(x, random_u, bandwidth).sum()
        dfun = lambda x: (gaussian_kernel(x, random_u, bandwidth) * (x-random_u)/bandwidth**2).sum()
        #res.append( -opt.minimize(fun, x0=random_u[0],  method='SLSQP', bounds=[(0, length)], options=options, jac=dfun).fun  )
        res.append( max( -opt.minimize(fun, x0=u, method='SLSQP', bounds=[(0, length)], options=options, jac=dfun).fun
                        for u in np.random.choice(random_u, size=int(np.log2(tags))//2+1) ) )
    return [gene_id, res ]

def args_check(args):
    if not os.path.exists(args.out_path):
        os.makedirs(args.out_path)
    if args.p_value >= 1 or args.p_value >= 1:
        print('Error: the p value should be >0 and <1.')
        exit()
    if args.core >= mp.cpu_count():
        print('Warning: the setting of CPU core beyond the max CPU. %s are used'%(mp.cpu_count()-1))
    if args.seed != None and (args.seed < 0 or args.seed > 2**32-1):
        print('Error: Seed must be between 0 and 2**32 - 1')
        exit()


parser = argparse.ArgumentParser(description=HELP, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-r", "--ref", type=str, required = True,
                    help = "Annotation file in <*.gtf> format, must include UTR feature.")
parser.add_argument("-i", "--input_file", type=str, nargs="+",required=True,
                    help="a list of <*.bed> file whose 5th column is the height.")
parser.add_argument("-o", "--out_path", type=str, default="./",
                    help="The output path.")
parser.add_argument("-e", "--extend", type=int, default = 1000,
                    help="extend 1000 nt downstream of 3'UTR.")
parser.add_argument("-j", "--core", type=int, default = 1,
                    help="The number of CPU cores to use.")
parser.add_argument("-fiter", "--fun_maxiter", type=int, default = 200,
                    help="The maximum number of iteration for SLSQP method.")
parser.add_argument("-bgiter", "--background_maxiter", type=int, default=1000,
                    help="The maximum number of iteration for background model generation.")
parser.add_argument("-t", "--tolerance", type=float, default = 10e-3,
                    help="Precision goal for the value of f in the stopping criterion.")
parser.add_argument("-b", "--bandwidth", type=float, default=20,
                    help="The bandwidth of gaussian kernal function.")
parser.add_argument("-p", "--p_value", type=float, default=0.05,
                    help="The threshold of significance, should be within (0, 1).")
parser.add_argument("-s", "--seed", type=int, default=None,
                    help="seed for random number generator")
args = parser.parse_args()
args_check(args)

#0. generate 3'UTR region
utr_name = os.path.join(args.out_path, os.path.split(args.ref)[-1] + '.utr_bed')
with open(utr_name, 'w') as out:
    for line in open(args.ref, 'r'):
        s = line.strip().split('\t')
        for i in s[-1].split('|'):
            out.write( '\t'.join([s[0], '\t'.join(i.split(',')), s[3], s[4], s[5] ])+'\n' )

#1. extend downstream of genes
extend_name = os.path.join(args.out_path, os.path.split(args.ref)[-1] + '.extend_bed')
cmd = r"""less -S %s |awk 'BEGIN{OFS="\t"}{if($6=="+") print $1, $3, $3+%d, $4, $5, $6; else if($6=="-" && $2<%d) print $1, 0, $2, $4, $5, $6; else print $1, $2-%d, $2, $4, $5, $6}' |bedtools intersect -s -v -a - -b %s >%s"""%(args.ref, args.extend, args.extend, args.extend, args.ref, extend_name)
os.system(cmd)


for i in args.input_file:
    expr = {}
    cmd = 'bedtools intersect -s -a %s -b %s -wa -wb -loj'%(utr_name, i)
    for line in os.popen(cmd).read().strip().split('\n'):
        s = line.strip().split('\t')
        if not s[3] in expr:
            expr[s[3]] = [0, 0]
        num = 0 if s[6] == '.' else int(s[10])
        length = int(s[2]) - int(s[1])
        expr[s[3]][0] += num
        expr[s[3]][1] += length

    #1. parameter preparation
    gene_id = [ j for j in expr if expr[j][0] > 0 ]
    tags = [ int(expr[j][0]) for j in gene_id]
    length = [ int(expr[j][1]) for j in gene_id]
    bandwidth = [ args.bandwidth ] * len(gene_id)
    options = [ {'ftol':args.tolerance, 'maxiter':args.fun_maxiter} ] * len(gene_id)
    if args.seed != None:
        np.random.seed(args.seed)
    random_seeds = list(np.random.randint(2**32 - 1, size=(len(gene_id), args.background_maxiter) ))
    arg_list = list(zip(gene_id, tags, length, bandwidth, options, random_seeds))

    #2. backgroud generation
    starttime = datetime.datetime.now()
    core = mp.cpu_count()-1 if args.core >= mp.cpu_count() else args.core
    pool = mp.Pool(processes=core)
    result = pool.map(background_model, arg_list)
    pool.close()
    pool.join()
    for j in result:
        expr[j[0]].append(j[-1])
    endtime = datetime.datetime.now()
    print( 'Time for sample %s: '%(i), (endtime - starttime).seconds )

    out_name = os.path.join(args.out_path, os.path.split(i)[-1].split('.')[0]+'.rCS.bed')
    with open(out_name, 'w') as out:
        cmd = r"""cat %s %s |awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $5, $6}' |sortBed |bedtools intersect -s -a %s -b - -wa -wb """%(args.ref, extend_name, i)
        for line in os.popen(cmd).read().strip().split('\n'):
            s = line.strip().split('\t')
            if type(expr[s[-3]][-1]) == list:
                dist = np.array(expr[s[-3]][-1], dtype="float64")
                p = sum( dist >= int(s[4]) ) / args.background_maxiter
                out.write( '\t'.join(['\t'.join(s[:6]), str(p), s[-3], str(expr[s[-3]][0]), str(expr[s[-3]][1]), str(expr[s[-3]][0]/expr[s[-3]][1]) ])+'\n' )
            #else:
            #    out.write( '\t'.join(['\t'.join(s[:6]), 'None', s[-3], str(expr[s[-3]][0]), str(expr[s[-3]][1]), str(expr[s[-3]][0]/expr[s[-3]][1]) ])+'\n' )
    #break

