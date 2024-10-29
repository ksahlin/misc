#!/bin/bash

# RUN scripts e.g. as: ./run_all.sh dataset1/

if [ $# -lt 1 ]; then
    # TODO: print usage
    echo "./run_all.sh <dataset_folder>"
    exit 1
fi

folder=$1
genome=$folder/genome.fa
query=$folder/query.fa
annot=$folder/annotation.gtf
mm_out=$folder/minimap2/
mkdir -p $mm_out

echo -n -e  "Parameters\t\t\t\t\t#E\t#A\tE_min\tA_min"$'\n'

for params in "--eqx -a" "--eqx -ax splice" "--eqx -ax splice -k 10" "--eqx -ax splice -k 10 -w 1" "--eqx -ax splice -u n" "--eqx -ax splice -u n -k 10" "--eqx -ax splice -u n -k 10 -w 1" "--eqx -ax splice -u n -k 10 -B 10 -O 4,12" "--eqx -ax splice -u n -k 10 -w 1 -B 10 -O 4,12" #"--eqx -ax splice -k 7 -w 1 -u n -B 10 -O 4,12"
do
	# params="--eqx -ax splice"
	name="mm2" #$params"\t\t\t\t\t\t"
	minimap2 $params $genome $query  1> $mm_out/$name.sam 2> /dev/null
	tab="\t"
	le=${#params}
	n=$(((54-$le) / 8))
	result=$(printf "%0.s$tab" $(seq 1 $n))
	echo -n $params
	echo -n -e "$result"
	python evaluate_instance.py --gtf $annot --samfile $mm_out/$name.sam
done


# Ultra
ultra_out=$folder/ultra/
mkdir -p $ultra_out
rm -r $ultra_out/*
echo -n -e  "Parameters\t#E\t#A\tE_min\tA_min"$'\n'

for params in "default"
do
	name="ultra" #$params"\t\t\t\t\t\t" $params
	uLTRA pipeline --prefix $name $genome $annot $query $ultra_out &> /dev/null # $ultra_out/$name.sam 
	tab="\t"
	le=${#params}
	n=2 #$(((54-$le) / 8))
	result=$(printf "%0.s$tab" $(seq 1 $n))
	echo -n $params
	echo -n -e "$result"
	python evaluate_instance.py --gtf $annot --samfile $ultra_out/$name.sam
done
