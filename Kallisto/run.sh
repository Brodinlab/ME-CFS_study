#/bin/bash
path='path'
suff_1='R1_001.fastq'
suff_2='R2_001.fastq'
index='path_to_index'
gtf='gft_from_ensembl'
python step1.py $path $suff_1 $suff_2 $index
sh kallisto.sh
python step3.py $path $gtf
