#include "hash.h"

using namespace std;
extern vector <int> random_hashes;

void read_input_long_reads(string filename, vector<string> &long_reads) {
  ifstream inFile;
  std::string line;

  inFile.open(filename);

  if (!inFile) {
    cerr << "Unable to open file: " << filename;
    exit(1);
  }

  while(std::getline(inFile, line)) {
    if (line.size() > 0 && line[0] == '>') continue;
    //std::cout << "long read: " << line << std::endl;
    long_reads.push_back(line);
  }

  return;
}

std::string read_reference_genome(string filename) {
  std::string str = "";
  ifstream inFile;
  std::string line;

  inFile.open(filename);

  if (!inFile) {
    cerr << "Unable to open file: " << filename;
    exit(1);
  }

  while(std::getline(inFile, line)) {
    if (line.size() > 0 && line[0] == '>') continue;
    str = str + line;
  }
  return str;
}

int get_similarity() {
  return 1;
}

float calculate_threshold(int hash_count) {
  /*
  Pick a number
of bands b and a number of rows r such that br = n

  hash_count # no of bands of sizeof(int) rows each
  b = hash_count, r = sizeof(int)
  threshold = (1/b) to the power (1/r)
  according to threshold calculation in: http://infolab.stanford.edu/~ullman/mmds/ch3.pdf
  */

  int r = sizeof(hash_count);
  int b = hash_count;
  float threshold;

  //std::cout << "b,r = " << b << ", " << r << std::endl;
  threshold = std::pow( (1.0/b), (1.0/r) );
  //std::cout << "threshold = " << threshold << std::endl;

  return threshold;
}

/*
Maps long reads to the reference
only if the Min Hash sketch of long read and genome
have a similarity higher than a threshold
*/
vector<float> map_long_reads(int ksize,
  int hash_count,
  float false_positive_rate,
  vector<string> long_reads,
  string reference_string,
  float threshold) {

  vector<string> kmers_B;
  vector<int> hashes_B;
  vector<float> jaccard_est;

  create_random_hashes(random_hashes, hash_count);
  create_kmers(ksize, reference_string, kmers_B);
  get_hashes(kmers_B, hashes_B, hash_count);

  //populate bloom filter
  bloom_parameters parameters;
  parameters.projected_element_count = kmers_B.size();
  parameters.false_positive_probability = false_positive_rate * 1.0;
  parameters.random_seed = rand();
  if (!parameters) {
    std::cout << "Error: Invalid set of bloom filter parameters!" << std::endl;
    return jaccard_est;
  }
  parameters.compute_optimal_parameters();

  bloom_filter filter(parameters); //instantiate

  //insert into bloom filter
  int filter_size = 0;
  for (int i=0; i < kmers_B.size(); i++) {
    string kmer = kmers_B.at(i);
    if (filter.contains(kmer)) continue;
    filter.insert(kmer);
    filter_size ++;
  }

  vector<float> true_jaccard_indices, min_jaccard_indices;

  for (int i =0; i < long_reads.size(); i++) {
    //find if kmers from set A which give min hash are in bloom filter
    std::string A = long_reads.at(i);
    vector<string> kmers_A;
    vector<int> hashes_A;
    create_kmers(ksize, A, kmers_A);
    get_hashes(kmers_A, hashes_A, hash_count);
    vector<string> min_kmers_A;

    //find min kmers for hashes of A
    min_kmers_A = get_min_kmers_for_hashes (kmers_A, hashes_A, hash_count);
    float intersection = 0.0;
    for (int k=0; k < min_kmers_A.size(); k++) {
      if (! filter.contains(min_kmers_A.at(k))) continue;
      intersection = intersection + 1;
    }

    //similarity - containment approach
    intersection = intersection - (false_positive_rate*hash_count);
    float containment_index = intersection / ((float)hash_count);
    float jaccard_index = std::abs( (kmers_A.size() * containment_index)
    / (kmers_A.size() + filter_size - kmers_A.size() * containment_index) );
    jaccard_est.push_back(jaccard_index);

    float min_jaccard = estimate_jaccard_index_min_hash(hashes_A, hashes_B);
    float true_jaccard = get_true_jaccard_index(kmers_A, kmers_B);

    true_jaccard_indices.push_back(true_jaccard);
    min_jaccard_indices.push_back(min_jaccard);

    int long_read_count = i + 1;
    //print results
    std::cout << std::endl <<"long read " << long_read_count << ":- " << std::endl;
    std::cout << A << std::endl;
    //std::cout << "containment_index = " << containment_index << std::endl;
    //std::cout << "(Min Hash) Est. jaccard_index = " << min_jaccard << std::endl;
    if (jaccard_index >= threshold) {
      /*
      "\033[1;31mbold red text\033[0m\n"
      */
      std::cout << "\033[1;32m(Containment) Est. jaccard_index = " << jaccard_index  << "\033[0m" << std::endl;
      std::cout << "Above Threshold. This Long Read maps to the Reference Genome" << std::endl;
    } else {
      std::cout << "\033[1;31m(Containment) Est. jaccard_index = " << jaccard_index  << "\033[0m" << std::endl;
      std::cout << "Below Threshold. This Long Read \033[1;31mDOES NOT\033[0m map to the Reference Genome" << std::endl;
    }
    //std::cout << "True Jaccard = " << true_jaccard << std::endl;
    std::cout << "Similarity = " << (jaccard_index * 100) << std::endl;
  }

  return jaccard_est;
}

int main(int argc, char * argv[]) {
  int option;
  std::string reference_file, long_reads_file;
  int ksize, hash_count;
  float threshold;
  float false_positive_rate;

  threshold = calculate_threshold(hash_count);

  while ( (option = getopt(argc,argv,"r:q:t:k:h:f:"))  != -1) {
    switch(option) {
      case 'r':
        reference_file = optarg;
        break;
      case 'q':
        long_reads_file = optarg;
        break;
      case 't':
        threshold = atof(optarg);
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
      default:
        std::cout << "Please enter cmdline args correctly:" <<std::endl;
        std::cout << "mapper -r <reference> -q <long reads>" <<
        "-t <threshold> -k <kmer size> -h <hash count>" <<
        "-f <false positive rate>" << std::endl;
        return 1;
    }
  }

  std::cout << "reference_file = " << reference_file << std::endl;
  std::cout << "long_reads_file = " << long_reads_file << std::endl;
  std::cout << "threshold = " << (threshold*100) << " percent"<< std::endl;
  std::cout << "kmer size = " << ksize << std::endl;
  std::cout << "hash count = " << hash_count << std::endl;
  std::cout << "false positive rate = " << false_positive_rate << std::endl;

  vector<string> long_reads;
  std::string reference_string;

  read_input_long_reads(long_reads_file, long_reads);
  reference_string = read_reference_genome(reference_file);
  //std::cout << "reference string = " << reference_string << std::endl;

  vector<float> jaccard_estimates =
  map_long_reads(ksize, hash_count, false_positive_rate,
    long_reads, reference_string, threshold);

  return 1;
}
