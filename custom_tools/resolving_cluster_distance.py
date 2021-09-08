import sys
import os

path = sys.argv[1]
files = [ os.path.join(path, i) for i in os.listdir(path) if i.endswith('.rCS.bed') ]

stat = {}
distance = 51
for i in range(distance):
    stat[i] = []

for i in files:
    total_sites = int(os.popen(r'''less -S %s |awk '{if($7<0.05) print $0}' |wc -l'''%(i) ).read().strip())
    for dist in range(distance):
        cmd = r"""less -S %s |awk '{if($7!="None" && $7<0.05) print $0}' |bedtools cluster -s -d %d -i - |awk '{print $NF}' |uniq -c |awk '{if($1==1) print $0}' |wc -l"""%(i, dist)
        non_clustered = int(os.popen(cmd).read().strip())
        stat[dist].append( (total_sites-non_clustered)/total_sites )

for i in stat:
    print(i, sum(stat[i])/len(stat[i]), sep="\t" )








