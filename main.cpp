#include "hash.h"

using namespace std;
extern vector <int> random_hashes;

int main(int argc, char *argv[]) {

  int option;
  int ksize; //length of kmers
  int hash_count; //number of hashes to use
  float false_positive_rate; //false positive rate for bloom filter
  std::string A;
  std::string B;
  std::string filename1, filename2;
  bool relative_error_hash = false;
  bool draw = false;

  while((option = getopt(argc,argv,"g:r:k:h:f:z:d:")) != -1) {
    switch(option) {
      case 'g':
        filename1 = optarg;
        break;
      case 'r':
        filename2 = optarg;
        break;
      case 'k':
        ksize = atoi(optarg);
        break;
      case 'h':
        hash_count = atoi(optarg);
        break;
      case 'f':
        false_positive_rate = atof(optarg);
        break;
      case 'z':
        relative_error_hash = true;
        break;
      case 'd':
        draw = true;
        break;
      default:
        std::cout << "\nPlease enter commandline args correctly: \nmain -g file1.txt -h file2.txt -k 4 -h 10 -f 0.01 \n";
        return -1;
    }
  }

  A = read_input(filename1);
  B = read_input(filename2);

  std::cout << "\nfile1 = " << filename1 << std::endl;
  //std::cout << "A = " << A << std::endl;
  std::cout << "file2 = " << filename2 << std::endl;
  //std::cout << "B = " << B << std::endl;
  //std::cout << "kmer size = " << ksize << std::endl;
  std::cout << "hash count = " << hash_count << std::endl;
  std::cout << "false positive rate = " << false_positive_rate << std::endl;

  if (A.size() > B.size()) {
    string swap = B;
    B = A;
    A = swap;
  }

  if (!draw) {
    std::cout << "A size = " << A.size() << std::endl;
    std::cout << "B size = " << B.size() << std::endl;

    vector<string> kmers_A, kmers_B;
    vector<int> hashes_A, hashes_B;
    string larger_str, smaller_str;

    //***
    create_random_hashes(random_hashes, hash_count); //
    create_kmers(ksize, A, kmers_A);
    create_kmers(ksize, B, kmers_B); //
    get_hashes(kmers_A, hashes_A, hash_count);
    get_hashes(kmers_B, hashes_B, hash_count); //
    //***

    float min_jaccard = estimate_jaccard_index_min_hash(hashes_A, hashes_B);
    float containment_jaccard = estimate_containment_index(ksize, kmers_A, kmers_B,
      false_positive_rate, hash_count);
    float true_jaccard = get_true_jaccard_index(kmers_A, kmers_B);
    float min_relative_error = std::abs(min_jaccard-true_jaccard) / true_jaccard;
    float containment_relative_error = std::abs(containment_jaccard-true_jaccard) / true_jaccard;

    std::cout << "kmer size = " << ksize <<std::endl;
    std::cout << std::endl;
    std::cout << "Calcuting Similarity: " << std::endl;
    std::cout << "True Jaccard = " << true_jaccard <<std::endl;
    std::cout << "Estimated Jaccard from Min Hash = " << min_jaccard <<std::endl;
    std::cout << "Estimated Jaccard from Containment Min Hash = " << containment_jaccard <<std::endl;
    std::cout << std::endl;
    std::cout << "Calcuting Relative Error: " << std::endl;
    std::cout << "Relative Error with Min Hash = " << min_relative_error <<std::endl;
    std::cout << "Relative Error with Containment Min Hash = " << containment_relative_error <<std::endl;

    return 1;
  }

  //-----------------------------------------------------------
  //compare both approaches:
  ofstream outfile1, outfile2, outfile3, outfile4, outfile5, outfile6;
  if (relative_error_hash){
    goto compare_hash_count;
  }
  outfile1.open("results/min.txt");
  outfile2.open("results/containment.txt");
  outfile3.open("results/true_jaccard.txt");
  outfile4.open("results/min_relative_error.txt");
  outfile5.open("results/containment_relative_error.txt");
  outfile6.open("results/kmer_sizes.txt");

  std::cout << "A size = " << A.size() << std::endl;
  std::cout << "B size = " << B.size() << std::endl;

  std::cout << "kmer size | est jaccard (min) | (containment)  ";
  std::cout << "| true jaccard  | relative error (min) | (containment) " << std::endl;
  for (int i = A.size() / 3; i >= 15; i--) {
    ksize = i;

    vector<string> kmers_A, kmers_B;
    vector<int> hashes_A, hashes_B;
    string larger_str, smaller_str;

    //***
    create_random_hashes(random_hashes, hash_count); //
    create_kmers(ksize, A, kmers_A);
    create_kmers(ksize, B, kmers_B); //
    get_hashes(kmers_A, hashes_A, hash_count);
    get_hashes(kmers_B, hashes_B, hash_count); //
    //***

    float min_jaccard = estimate_jaccard_index_min_hash(hashes_A, hashes_B);
    float containment_jaccard = estimate_containment_index(ksize, kmers_A, kmers_B,
      false_positive_rate, hash_count);
    float true_jaccard = get_true_jaccard_index(kmers_A, kmers_B);
    float min_relative_error = std::abs(min_jaccard-true_jaccard) / true_jaccard;
    float containment_relative_error = std::abs(containment_jaccard-true_jaccard) / true_jaccard;

    std::cout << " ---> " << std::fixed << std::setw(3)<< std::setprecision(6) << ksize;
    std::cout << " | " << std::fixed << std::setw(15)<< std::setprecision(6) << min_jaccard;
    std::cout << " |  " << std::fixed << std::setw(15)<< std::setprecision(6)
    << containment_jaccard;
    std::cout << " | " << std::fixed << std::setw(15)
    << std::setprecision(6) << true_jaccard;
    std::cout << " | " << std::fixed << std::setw(15)
    << std::setprecision(6) << min_relative_error;
    std::cout << " | " << std::fixed << std::setw(15)
    << std::setprecision(6) << containment_relative_error;
    std::cout << std::endl;

    outfile1 << min_jaccard << ",";
    outfile2 << containment_jaccard << ",";
    outfile3 << true_jaccard << ",";
    outfile4 << min_relative_error << ",";
    outfile5 << containment_relative_error << ",";
    outfile6 << ksize << ",";
  }
  outfile1.close();
  outfile2.close();
  outfile3.close();
  outfile4.close();
  outfile5.close();
  outfile6.close();

  //-----------------------------------------------------------

  if (!relative_error_hash) return 1;

  compare_hash_count:

  outfile1.open("results/min.txt");
  outfile2.open("results/containment.txt");
  outfile3.open("results/true_jaccard.txt");
  outfile4.open("results/hash_min_relative_error.txt");
  outfile5.open("results/hash_containment_relative_error.txt");
  outfile6.open("results/hash_count.txt");

  std::cout << "A size = " << A.size() << std::endl;
  std::cout << "B size = " << B.size() << std::endl;

  std::cout << "hashes# | est jaccard (min) | (containment)  ";
  std::cout << "| true jaccard  | relative error (min) | (containment) " << std::endl;

  for (int i = 50; i <= 1000; i = i + 10) {
    hash_count = i;

    vector<string> kmers_A, kmers_B;
    vector<int> hashes_A, hashes_B;
    string larger_str, smaller_str;

    //***
    create_random_hashes(random_hashes, hash_count); //
    create_kmers(ksize, A, kmers_A);
    create_kmers(ksize, B, kmers_B); //
    get_hashes(kmers_A, hashes_A, hash_count);
    get_hashes(kmers_B, hashes_B, hash_count); //
    //***

    float min_jaccard = estimate_jaccard_index_min_hash(hashes_A, hashes_B);
    float containment_jaccard = estimate_containment_index(ksize, kmers_A, kmers_B,
      false_positive_rate, hash_count);
    float true_jaccard = get_true_jaccard_index(kmers_A, kmers_B);
    float min_relative_error = std::abs(min_jaccard-true_jaccard) / true_jaccard;
    float containment_relative_error = std::abs(containment_jaccard-true_jaccard) / true_jaccard;

    /*
    std::cout << " ---> " << std::fixed << std::setw(3)<< std::setprecision(6) << hash_count;
    std::cout << " | " << std::fixed << std::setw(15)<< std::setprecision(6) << min_jaccard;
    std::cout << " |  " << std::fixed << std::setw(15)<< std::setprecision(6)
    << containment_jaccard;
    std::cout << " | " << std::fixed << std::setw(15)
    << std::setprecision(6) << true_jaccard;
    std::cout << " | " << std::fixed << std::setw(15)
    << std::setprecision(6) << min_relative_error;
    std::cout << " | " << std::fixed << std::setw(15)
    << std::setprecision(6) << containment_relative_error;
    std::cout << std::endl;
    */

    outfile1 << min_jaccard << ",";
    outfile2 << containment_jaccard << ",";
    outfile3 << true_jaccard << ",";
    outfile4 << min_relative_error << ",";
    outfile5 << containment_relative_error << ",";
    outfile6 << hash_count << ",";
  }

  outfile1.close();
  outfile2.close();
  outfile3.close();
  outfile4.close();
  outfile5.close();
  outfile6.close();

  return 1;
}
