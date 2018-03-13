import numpy as np
import random
import sys, os

def generate_input_files(ksize, len_small_string, len_large_string, long_reads_count, long_reads_file, reference_file):
    # Create data: small set called A, large set called B

    small_string = ""
    len_small = len_small_string
    # small string to form the small set A
    random_numbers = []
    for i in range(long_reads_count):
        random_numbers.append(random.randint(75,150))
    random.shuffle(random_numbers)

    for i in range(long_reads_count):
        #within 75% to 150% of len small string
        len_small = len_small_string + (random_numbers[i] / 10)

        #generate string
        small_string = ''.join(np.random.choice(['A', 'C', 'T', 'G'], len_small))

        #append string to long reads file
        option = 'w'
        if os.path.exists(long_reads_file):
            option = 'a' # append if already exists
        with open(long_reads_file, option) as the_file:
            the_file.write(">long_read" + str(i+1) + "\n")
            the_file.write(small_string)
            the_file.write("\n")

    # size of smaller set, used to convert containment index to Jaccard index
    size_A = len(set([small_string[i:i+ksize] for i in range(len(small_string) - ksize + 1)]))

    # large string to form the larger set B
    large_string = ''.join(np.random.choice(['A', 'C', 'T', 'G'], len_large_string)) + small_string

    with open(reference_file, 'w') as the_file:
        the_file.write(">reference" + "\n")
        the_file.write(large_string)

if __name__ == '__main__':
    ksize = 0
    long_reads_count = 0

    if len(sys.argv) == 7:
        ksize = int(sys.argv[1])  # k-mer length
        long_reads_count = int(sys.argv[2])
        len_small_string = int(sys.argv[3])
        len_large_string = int(sys.argv[4])
        long_reads_file = sys.argv[5]
        reference_file = sys.argv[6]
    else:
        print "\ngenerate_data.py: incorrect commandline args\n"

    generate_input_files(ksize, len_small_string, len_large_string, long_reads_count, long_reads_file, reference_file)
