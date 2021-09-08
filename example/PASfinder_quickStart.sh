PASfinder=${PWD}/../bin
input_file=${PWD}/0.rawdata/*3.fq.gz
output_path=${PWD}/result
library_type=3
adaptor=agatcggaagagc
bowtie_index=${PWD}/bowtie_index_chrX/chrX.fa
genome_size=${PWD}/../database/hg19.chrom.sizes
genome_fa=${PWD}/bowtie_index_chrX/chrX.fa
annotation=${PWD}/../database/gencode.v19.annotation.bed
core=20
extend=1000
polyA_length=8
p=0.05


python3 ${PASfinder}/1.processing_artificial_sequence.py -i ${input_file} -s ${library_type} -a ${adaptor} -o ${output_path}/1.preprocessing
python3 ${PASfinder}/2.bowtie_SE_mapping.py -p ${core} -in ${bowtie_index} -i ${output_path}/1.preprocessing/*.clean.fq.gz -o  ${output_path}/2.mapping
python3 ${PASfinder}/3.collapasing.py -i ${output_path}/2.mapping/*.bam -s ${library_type} -l ${polyA_length} -o ${output_path}/3.filter_internal_primed_events -gs ${genome_size}
python3 ${PASfinder}/4.filter_internal_primed_events.py -i ${output_path}/3.filter_internal_primed_events/*.pCS.bed -o ${output_path}/3.filter_internal_primed_events -g ${genome_fa}
python3 ${PASfinder}/5.identifying_reliable_cleavage_sites.py -j ${core} -e ${extend} -r ${annotation} -i ${output_path}/3.filter_internal_primed_events/*.heuristic.bed -o ${output_path}/4.reliable_cleavage_sites -p ${p}
python3 ${PASfinder}/6.clustering_cleavage_sites.py -i ${output_path}/4.reliable_cleavage_sites/*.rCS.bed -o ${output_path}/5.cluster -p ${p}






