clone repo and type following:

------------------

- Estimate Jaccard Index by running Min Hash and Containment Min Hash.

./run.sh <file1> <file2> <kmer_size> <hash_functions>
 <false_positive_rate> <number_appended_to_results_files> <order_of_len_A> <order_of_len_B>

./run.sh data/file3.txt data/file4.txt 20 200 0.01 2 150 1000

./run.sh data/file1.txt data/file2.txt 18 1000 0.02 3 15 1000

Sample output: [img1]

- Compare both approaches. Reproduce the graphs in the report.

./run.sh data/filex3.txt data/filex4.txt 20 150 0.04 3

./run.sh data/filex1.txt data/filex2.txt 25 1000 0.01 4

Sample output: [img2] [img3]

------------------

- Generate reference genome and long reads.
- Mapper: maps long reads to reference genome
  if they have similarity higher than threshold

./query.sh <long_reads_file> <reference_genome_file> <threshold>
 <kmer_size> <hash_functions> <false_positive_rate> <number_of_long_reads_to_generate>

./query.sh data/reference.txt data/longreads.txt 0.05 18 100 0.01 5

./query.sh data/reference.txt data/longreads.txt 0.05 20 200 0.02 10

Sample output: [img4] [img5]
