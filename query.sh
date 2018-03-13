#!/bin/bash
clear

file1=$1
file2=$2
threshold=$3
ksize=$4
hash_count=$5
false_positive_rate=$6
long_reads=$7
len_small_string=100
len_large_string=1000
long_reads_file=data/longreads.txt
reference_file=data/reference.txt

echo "---------------------------------------------------------------"
echo "$(tput setaf 2)[+] Cleaning old data files$(tput sgr0)"
rm data/reference.txt data/longreads.txt

echo "$(tput setaf 2)[+] Compiling code$(tput sgr0)"
make mapper

echo "$(tput setaf 2)[+] Generating Data$(tput sgr0)"
python generate_data.py ${ksize} ${long_reads} ${len_small_string} ${len_large_string} ${long_reads_file} ${reference_file}

echo "$(tput setaf 2)[+] Reference Genome file$(tput sgr0)"
cat data/reference.txt
echo ""

echo "$(tput setaf 2)[+] Long Reads file$(tput sgr0)"
cat data/longreads.txt
echo ""

echo "$(tput setaf 2)[+] Running Query$(tput sgr0)"
./mapper -r ${file1} -q ${file2} -t ${threshold} -k ${ksize} -h ${hash_count} -f ${false_positive_rate}
