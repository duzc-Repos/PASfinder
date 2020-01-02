# PASfinder

This pipeline is designed for identification of cleavage sites and detection of alternative polyadenylation from 3' end sequencing.

## 0. Prerequisites:
This pipeline is implemented on __Liunx__ using PYTHON 3.x with SciPy, NumPy, os, multiprocessing, datatime and argparse. 
Some softwares need to be installed and add to environment variables ($PATH):
* cutadapt (version 1.18)
* bowtie (version 1.2.2)
* samtools (Version: 1.6 (using htslib 1.6)), availiable in version 1.2 and above.
* bedtools (v2.26.0). Only v2.26.0 is avaliable.      

some R packages:
* cleanUpdTSeq (v1.16.0, optional): filtering internal primed events by bayes classifer
* DEXSeq (v1.24.4, optional): for alternative polyadenylation events identification.

## 1. Introduction
### 1.1 Main Steps
Generally, there are 5 steps:
* preprocessing include processing artificial sequence (1.processing_artificial_sequence.py) and mapping to the reference genome (2.bowtie_SE_mapping.py).
* collapasing reads with same end (3.collapasing.py) and filter internal primed events (4.filter_internal_primed_events.py, cleanUpdTSeq.r).
* identification of reliable cleavage sites based on dynamic background model (5.identifying_reliable_cleavage_sites.py).
* clustering the close cleavage sites (6.clustering_cleavage_sites.py).
* detecting alternative polyadenylation (7.detecting_alternative_polyadenylation.py, DEXSeq.r)
### 1.2 Directory Structure
The following shows the Directory Structure:  
* __bin/__: The main scripts for this pipeline.  
* __custom_tools/__: Some useful tools to deal with well-organised annotation files (such as .gtf file from ensembl) and determine the distance for clustering.  
* __database/__: The necessary input files in the pipeline. Four species (human(hg19), mouse(mm9), fly(BDGP5) and Arabidopsis_thaliana(TAIR10)) are provided.  
* __example/__: Some demos on usage. see __Usage__ part.  
* __PASfinder.sh__: The overall script for user usage. Some paths need to be set by users.  
```
PASfinder_v1.1
├── bin
│   ├── 1.processing_artificial_sequence.py
│   ├── 2.bowtie_SE_mapping.py
│   ├── 3.collapasing.py
│   ├── 4.filter_internal_primed_events.py
│   ├── 5.identifying_reliable_cleavage_sites.py
│   ├── 6.clustering_cleavage_sites.py
│   ├── 7.detecting_alternative_polyadenylation.py
│   ├── cleanUpdTSeq.r
│   └── DEXSeq.r
├── custom_tools
│   ├── prepare_annotation.py
│   └── resolving_cluster_distance.py
├── database
│   ├── Arabidopsis_thaliana.TAIR10.43.bed
│   ├── BDGP5.chrom.sizes
│   ├── Drosophila_melanogaster.BDGP5.70.bed
│   ├── gencode.v19.annotation.bed
│   ├── gencode.vM1.annotation.bed
│   ├── hg19.chrom.sizes
│   ├── mm9.chrom.sizes
│   └── Tair10.chrom.sizes
├── example
│   ├── 0.rawdata
│   │   ├── testchrX_3seq_5.fq.gz
│   │   └── testchrX_PAS-seq_3.fq.gz
│   ├── bowtie_index_chrX
│   │   └── bowtie.sh
│   ├── PASfinder_quickStart.sh
├── PASfinder.sh
└── README
```

## 2. Usage
### 2.0 Preparation
For usage, users need to know:
* the type of library (sense or antisense) ("-s 5" for sense, "-s 3" for anti-sense)
* how long the artificial sequence existed at the 5' end of reads ("-l 6": the length of artificial sequence)
* set the pathes to needed file, like bowtie index, genome annotation, genemo fasta, geome size, etc.
* input <*.fastq or *.fastq.gz> and output path
### 2.1 Quick start
An example is provided in __./example/__ for quickstart using data from 3'-seq, PAS-seq on chrX. Bowtie index should be build and a script __bowtie.sh__ is provided for this process in __./example/bowtie_index__.
```
cd /path/to/example/bowtie_index_chrX
bash bowtie.sh
cd ../
bash PASfinder_quickStart.sh
```
### 2.2 Step by Step
#### Step 0: preparation
```
PASfinder=<path to "PASfinder_script_v1.1">                    # /path/to/example
sense_library_input_files=<path to raw fastq>                  # ${PASfinder}/example/*5.fq.gz
antisense_library_input_files=<path to raw fastq>              # ${PASfinder}/example/*3.fq.gz
output_path=<paht to output path>                              # ${PASfinder}/example
adaptor=<sequence of adaptor>                                  # agatcggaagagc (illumina)
bowtie_index=<path to bowtie index>                            # ${PASfinder}/example/bowtie_index_chrX/chrX.fa
genome_size=<path to genome size, like "hg19.chrom.sizes">     # ${PASfinder}/database/hg19.chrom.sizes
genome_fa=<path to genome fasta file>                          # ${PASfinder}/example/bowtie_index_chrX/chrX.fa
annotation=<path to processed species annoataion file>         # ${PASfinder}/database/gencode.v19.annotation.bed
end5seq=<artificial sequence at 5 end of reads >               # 6
core=<#CPU>                                                    # 20
extend=<extend __bp downstream of transcripts>                 # 1000
polyA_length=<reads with at least __As are regarded as CS>     # 8 #only use for sense library
p=<P-value>                                                    # 0.05
```
#### Step 1: preprocessing include processing artificial sequence
For sense library
```
python3 ${PASfinder}/bin/1.processing_artificial_sequence.py -l ${end5seq} -i ${sense_library_input_files} -s 5 -a ${adaptor} -o ${output_path}/1.preprocessing
```
For antisense library
```
python3 ${PASfinder}/bin/1.processing_artificial_sequence.py -i ./0.rawdata/*3.fq.gz -a agatcggaagagc -s 3 -o ./1.preprocessing
```
For detail:
```
python3 ${PASfinder}/bin/1.processing_artificial_sequence.py -h
```
#### Step2: mapping to the reference genome
```
python3 ${PASfinder}/2.bowtie_SE_mapping.py -in ${bowtie_index} -i ${output_path}/1.preprocessing/*.clean.fq.gz -o  ${output_path}/2.mapping
```
For detail:
```
python3 ${PASfinder}/2.bowtie_SE_mapping.py -h
```
#### Step 3: collapasing reads with same end
For sense library
```
python3 ${PASfinder}/bin/3.collapasing.py -i ${output_path}/2.mapping/*5.sort.bam -s 5 -l 8 -o ${output_path}/3.filter_internal_primed_events -gs ${genome_size}
```
For antisense library
```
python3 ${PASfinder}/bin/3.collapasing.py -i ${output_path}/2.mapping/*3.sort.bam -s 3 -o ${output_path}/3.filter_internal_primed_events -gs ${genome_size}
```
For detail:
```
python3 ${PASfinder}/bin/3.collapasing.py -h
```
#### Step 4: filter internal primed events
Using heuristic method (default):
```
python3 ${PASfinder}/bin/4.filter_internal_primed_events.py  -i ${output_path}/3.filter_internal_primed_events/*.pCS.bed -g ${genome_fa} -o ${output_path}/3.filter_internal_primed_events/ -m heuristic
```
Using bayes classifer:
```
python3 ${PASfinder}/bin/4.filter_internal_primed_events.py  -i ${output_path}/3.filter_internal_primed_events/*.pCS.bed -g ${genome_fa} -o ${output_path}/3.filter_internal_primed_events/ -m cleanUpdTseq
```
For details:
```
python3 ${PASfinder}/bin/4.filter_internal_primed_events.py -h
```
#### Step 5: identification of reliable cleavage sites based on dynamic background model
```
python3 ${PASfinder}/bin/5.identifying_reliable_cleavage_sites.py -j ${core} -e ${extend} -r ${annotation} -i ${output_path}/3.filter_internal_primed_events/*.heuristic -o ${output_path}/4.reliable_cleavage_sites -p ${p}
```
or
```
python3 ${PASfinder}/bin/5.identifying_reliable_cleavage_sites.py -j ${core} -e ${extend} -r ${annotation} -i ${output_path}/3.filter_internal_primed_events/*.cleanUpdTSeq -o ${output_path}/4.reliable_cleavage_sites -p ${p}
```
#### Step 6: clustering the close cleavage sites
```
python3 ${PASfinder}/bin/6.clustering_cleavage_sites.py -i ${output_path}/4.reliable_cleavage_sites/*.rCS.bed -o ${output_path}/5.cluster -p ${p}
```
### 2.3 Detecting alternative polyadenylation
After identifying the cleavage sites, the results can be subject to detect the alternative polyadenylation. At least 2 samples are required by both control and treatment conditions.
```
python3 7.detecting_alternative_polyadenylation.py -c /path/to/control_1.cluster.bed /path/to/control_2.cluster.bed -t /path/to/treatment_1.cluster.bed /path/to/treatment_2.cluster.bed -p 0.05 -o /6.alternative_polyadenylation
```
```
usage: 7.detecting_alternative_polyadenylation.py [-h] -c CONTROL
                                                  [CONTROL ...] -t TREATMENT
                                                  [TREATMENT ...]
                                                  [-o OUT_PATH] [-d DISTANCE]
                                                  [-p P_VALUE]
                                                  [-padj P_ADJUST]
                                                  [-fc LOG2FC]

optional arguments:
  -h, --help            show this help message and exit
  -c CONTROL [CONTROL ...], --control CONTROL [CONTROL ...]
                        a list of control files generated by
                        <6.clustering_cleavage_sites.py>, al least 2.
  -t TREATMENT [TREATMENT ...], --treatment TREATMENT [TREATMENT ...]
                        a list of treatment files generated by
                        <6.clustering_cleavage_sites.py>, at least 2.
  -o OUT_PATH, --out_path OUT_PATH
                        The output path.
  -d DISTANCE, --distance DISTANCE
                        The distance between clusters in each sample to merge
                        as a single site.
  -p P_VALUE, --p_value P_VALUE
                        The threshold of significance in DEXSeq, should be
                        within (0, 1)
  -padj P_ADJUST, --p_adjust P_ADJUST
                        The threshold of significance in DEXSeq, should be
                        within (0, 1)
  -fc LOG2FC, --log2fc LOG2FC
                        The threshold of fold change (log2 transformed).
```
