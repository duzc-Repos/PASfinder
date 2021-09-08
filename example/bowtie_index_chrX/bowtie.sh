wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chrX.fa.gz
gzip -d chrX.fa.gz
samtools faidx chrX.fa
bowtie-build ./chrX.fa chrX.fa

