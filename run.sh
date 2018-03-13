#!/bin/bash
clear

file1=$1
file2=$2
ksize=$3
hash_count=$4
false_positive_rate=$5
example=$6
long_reads=1
len_small_string=$7
len_large_string=$8
long_reads_file=${file1}
reference_file=${file2}
compare=$9

echo "---------------------------------------------------------------"
echo "$(tput setaf 2)[+] Cleaning old output files$(tput sgr0)"
make clean

echo "$(tput setaf 2)[+] Compiling code$(tput sgr0)"
make run

if [ $# -ne 6 ]
then
  rm ${file1} ${file2}
  echo "$(tput setaf 2)[+] Generating Data$(tput sgr0)"
  python generate_data.py ${ksize} ${long_reads} ${len_small_string} ${len_large_string} ${long_reads_file} ${reference_file}
fi

echo "$(tput setaf 2)[+] Long Reads file$(tput sgr0)"
cat ${file1}
echo ""

echo "$(tput setaf 2)[+] Reference Genome file$(tput sgr0)"
cat ${file2}
echo ""

echo "$(tput setaf 2)[+] Estimate Jaccard Index with both approaches$(tput sgr0)"
echo "$(tput setaf 3)EXAMPLE ${example}: using ${file1} and ${file2}$(tput sgr0)"
echo "./main -g <smaller_genome_file> -r <reference_genome_file> -k <kmer_size> -h <number of hash_functions> -f <false_positive_rate>"
echo ""
./main -g ${file1} -r ${file2} -k ${ksize} -h ${hash_count} -f ${false_positive_rate}

if [ $# -ne 6 ]
then
  exit 1
fi
echo ""
echo "---------------------------------------------------------------"
read  -n 1 -p "To compare both approaches: PRESS ANY KEY TO CONTINUE"

echo "$(tput setaf 2)[+] Comparing both approaches$(tput sgr0)"
echo "$(tput setaf 3)EXAMPLE ${example}: using ${file1} and ${file2}$(tput sgr0)"
./main -g ${file1} -r ${file2} -h ${hash_count} -f ${false_positive_rate} -d draw

echo "$(tput setaf 2)[+] Drawing graphs$(tput sgr0)"
echo "$(tput setaf 2)[+] Saving graphs$(tput sgr0)"
python draw_graphs.py ${example}

echo "$(tput setaf 2)[+] Graphs generated are: $(tput sgr0)"
ls results/*.pdf

echo "$(tput setaf 2)[+] Showing graphs$(tput sgr0)"
open results/min${example}.pdf
open results/containment${example}.pdf
