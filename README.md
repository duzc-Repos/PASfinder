# PASfinder

This pipeline is designed for identification of cleavage sites from 3' end sequencing.

## Implementation:
This pipeline is implemented on Liunx using PYTHON 3.x with SciPy, NumPy, os, multiprocessing, datatime and argparse. 
Some softwares need to be installed and add to environment variables ($PATH):
* cutadapt (version 1.18)
* bowtie (version 1.2.2)
* samtools (Version: 1.6 (using htslib 1.6)), availiable in version 1.2 and above.
* bedtools (v2.26.0). Only v2.26.0 is avaliable.

## Introduction
The following shows the Directory Structure:
__bin/__: The main scripts for this pipeline.
__custom_tools/__: Some useful tools to deal with well-organised annotation files (such as .gtf file from ensembl) and determine the distance for clustering.
__database/__: The necessary input files in the pipeline. Four species (human(hg19), mouse(mm9), fly(BDGP5) and Arabidopsis_thaliana(TAIR10)) are provided.
__example/__: Some demos on usage. see __Example__ part.
__PASfinder.sh__: The overall script for user usage. Some paths need to be set by users.
```
PASfinder_v1.1
├── bin
│   ├── 1.processing_artificial_sequence.py
│   ├── 2.bowtie_SE_mapping.py
│   ├── 3.collapasing.py
│   ├── 4.filter_internal_primed_events.py
│   ├── 5.identifying_reliable_cleavage_sites.py
│   ├── 6.clustering_cleavage_sites.py
│   └── cleanUpdTSeq.r
├── custom_tools
│   ├── prepare_annotation.py
│   ├── README
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
│   └── PASfinder_stepBystep.sh
├── PASfinder.sh
└── README
```

## Steps
Generally, there are 4 steps:
* preprocessing include processing artificial sequence and mapping to the reference genome
* collapasing reads with same end and filter internal primed events
* identification of reliable cleavage sites based on dynamic background model
* clustering the close cleavage sites

## Usage
For usage, users need to know:
* the type of library (sense or antisense) ("-s 5" for sense, "-s 3" for anti-sense)
* how long the artificial sequence existed at the 5' end of reads ("-l 6": the length of artificial sequence)
* set the pathes to needed file, like bowtie index, genome annotation, genemo fasta, geome size, etc.
* input <*.fastq or *.fastq.gz> and output path

## Example
__./PSItools.sh__ is provided for run with default parameter, some paths need to be set by users. 
For custom usage, user can see scripts in __./bin/__ for each steps and adjust parameters according to the help in scripts. 

An example is provided in __./example/__ for quickstart and step by step using data from 3'-seq, PAS-seq on chrX. Bowtie index should be build and a script __bowtie.sh__ is provided for this process in __./example/bowtie_index__.
