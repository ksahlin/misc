#!/bin/bash

# RUN scripts e.g. as: ./run_all.sh dataset1/

if [ $# -lt 1 ]; then
    # TODO: print usage
    echo "./run_all.sh <dataset_folder>"
    exit 1
fi

folder=$1
genome=$folder/genome.fa
query=$folder/reads.fa
# annot=$folder/annotation.gtf
output=$folder/output
mkdir -p $output

#echo -n -e  "Parameters\t\t\t\t\t#E\t#A\tE_min\tA_min"$'\n'
mm2_params="--eqx -ax sr"
minimap2 $mm2_params $genome $query  1> $output/mm2.sam 2> /dev/null
python evaluate_instance.py $output/mm2.sam $output/mm2.csv minimap2

bwa index -p $output/bwa_mem_idx $genome &> /dev/null
bwa mem -t 1 -o $output/bwa_mem.sam $output/bwa_mem_idx $query 2> /dev/null
python evaluate_instance.py $output/bwa_mem.sam $output/bwa_mem.csv bwa_mem

strobealign -t 1 $genome $query 1> $output/strobealign.sam 2> /dev/null
python evaluate_instance.py $output/strobealign.sam $output/strobealign.csv strobealign


#python evaluate_instance.py --gtf $annot --samfile $mm_out/$name.sam

