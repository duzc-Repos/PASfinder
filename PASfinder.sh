#config: Absolute path is recommanded. 
# you should change to your path
PASfinder=<path to "PASfinder_script_v1.1">                    # /path/to/PASfinder_script_v1.1/
input_files=<path to raw fastq>                                # ${PASfinder}/example/0.rawdata/*.fq.gz
output_path=<paht to output path>                              # ${PASfinder}/example/result
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

control=<path to control samples>                              # 
treatment=<path to treatment samples>                          #

python3 ${PASfinder}/bin/1.processing_artificial_sequence.py -l ${end5seq} -i ${input_file} -s ${library_type} -a ${adaptor} -o ${output_path}/1.preprocessing
python3 ${PASfidner}/bin/2.bowtie_SE_mapping.py -in ${bowtie_index} -i ${output_path}/1.preprocessing/*.clean.fq.gz -o  ${output_path}/2.mapping
python3 ${PASfinder}/bin/3.collapasing.py -i ${output_path}/2.mapping/*.bam -s ${library_type} -l ${polyA_length} -o ${output_path}/3.filter_internal_primed_events -gs ${genome_size}
python3 ${PASfinder}/bin/4.filter_internal_primed_events.py -i ${output_path}/3.filter_internal_primed_events/*.pCS.bed -o ${output_path}/3.filter_internal_primed_events -g ${genome_fa}
python3 ${PASfinder}/bin/5.identifying_reliable_cleavage_sites.py -j ${core} -e ${extend} -r ${annotation} -i ${output_path}/3.filter_internal_primed_events/*.heuristic -o ${output_path}/4.reliable_cleavage_sites -p ${p}
python3 ${PASfinder}/bin/6.clustering_cleavage_sites.py -i ${output_path}/4.reliable_cleavage_sites/*.rCS.bed -o ${output_path}/5.cluster -p ${p}
python3 ${PASfinder}/bin/7.detecting_alternative_polyadenylation.py -c ${control} -t ${treatment} -p 0.05 -o ${output_path}/6.alternative_polyadenylation

