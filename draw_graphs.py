import matplotlib.pyplot as plt
import sys

def draw_graph(input_file, output_file, true_jaccard_file, title, example=""):
    x = []
    y = []
    with open(input_file) as f:
        y = f.read().split(',')
    y = y[:-1]
    with open(true_jaccard_file) as f:
        x = f.read().split(',')
    x = x[:-1]

    plt.gcf().clear()
    plt.xlabel('True Jaccard')
    plt.ylabel('Estimated Jaccard')
    if "Error" in title:
        plt.xlabel('Kmer Size')
        plt.ylabel('Relative Error')
    if "hash" in input_file:
        plt.xlabel('#Hash functions')
        plt.ylabel('Relative Error')
    plt.title(title)
    plt.plot(x,y)
    plt.savefig("results/" + output_file + str(example) + ".pdf")

if __name__ == '__main__':
    example = 0

    if len(sys.argv) == 2:
        example = int(sys.argv[1])

    input_min = "results/min.txt"
    input_containment = "results/containment.txt"
    true_jaccard = "results/true_jaccard.txt"
    input_min_rel_err = "results/min_relative_error.txt"
    input_containment_rel_err = "results/containment_relative_error.txt"
    input_kmer_sizes = "results/kmer_sizes.txt"
    input_hash_count = "results/hash_count.txt"
    input_hash_min_rel_err = "results/hash_min_relative_error.txt"
    input_hash_containment_rel_err = "results/hash_containment_relative_error.txt"

    output_min = "min"
    output_containment = "containment"
    min_output_relative_err = "min_relative_error"
    containment_output_relative_err = "containment_relative_error"

    draw_graph(input_min, output_min, true_jaccard, "Min Hash", example)
    draw_graph(input_containment, output_containment, true_jaccard, "Containment Hash", example)

    #draw_graph(input_min_rel_err, min_output_relative_err, input_kmer_sizes, "Relative Error : Min Hash", example)
    #draw_graph(input_containment_rel_err, containment_output_relative_err, input_kmer_sizes, "Relative Error : Containment", example)
