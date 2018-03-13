#include "hash.h"

using namespace std;
vector <int> random_hashes;

/*
Read input strings from files
*/
std::string read_input(string filename) {
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

/*
Create kmers from input string
*/
void create_kmers(int ksize, string genome, vector<string> &kmers) {
  for (int i=0; i <= genome.size() - ksize; i++) {
    string kmer = genome.substr(i, ksize);
    kmers.push_back(kmer);
  }
  return;
}

/*
Generate Random hash values
to be used while calculating hashes for kmers
*/
int create_random_hashes(vector<int> &random_hashes, int hash_count) {
  for (int i =0; i < hash_count; i++) {
    random_hashes.push_back(rand());
  }
  return 0;
}

/*
Create hashes for kmers and store kmers that give min hashes
- Using MurmurHash3
*/
vector<string> get_hash_from_kmers(vector<string> &kmers,
  int index_random_hashes, int hash_count) {
  int min_hash = 0;
  string min_kmer = "";
  for (int i =0; i < kmers.size(); i++) {
    string kmer = kmers.at(i);

    uint64_t murmur_hash_output[2];

    MurmurHash3_x64_128(kmer.c_str(),
    kmer.size(),
    random_hashes.at(index_random_hashes),
    murmur_hash_output);

    /*
    suggested hash functions:
      gi(x) =  ( h1(x) + i * h2(x) + i^2 ) mod m,
    page 188 from paper:
      https://www.eecs.harvard.edu/~michaelm/postscripts/rsa2008.pdf
    */
    int hash =  (murmur_hash_output[0] +  hash_count * murmur_hash_output[1]) % PRIME; //according to paper, best way to build hashes from bloom filters

    if (hash < min_hash) {
      min_kmer = kmer;
    }
    min_hash = std::min(hash, min_hash);
  }
  vector<string> min_hash_and_kmer;
  min_hash_and_kmer.push_back(std::to_string(min_hash));
  min_hash_and_kmer.push_back(min_kmer);
  return min_hash_and_kmer;
}

/*
Get hashes from kmers
*/
int get_hashes(vector<string> &kmers,
  vector<int> &hashes, int hash_count) {
  for (int i =0;i < hash_count; i++) {
    vector<string> min_hash_and_kmer = get_hash_from_kmers(kmers, i, hash_count);
    int hash = stoi(min_hash_and_kmer.at(0));
    hashes.push_back(hash);
  }
  return 0;
}

/*
Get kmers that give min hashes
*/
vector<string> get_min_kmers_for_hashes (vector<string> &kmers,
  vector<int> &hashes, int hash_count) {
    vector<string> min_kmers;
    for (int i =0;i < hash_count; i++) {
      vector<string> min_hash_and_kmer = get_hash_from_kmers(kmers, i, hash_count);
      int hash = stoi(min_hash_and_kmer.at(0));
      min_kmers.push_back(min_hash_and_kmer.at(1));
      hashes.push_back(hash);
    }
    return min_kmers;
}

/*
Print hashes of kmers
*/
int print_hashes(vector<int> &hashes) {
  for (int i =0; i < hashes.size(); i++) {
    std::cout << hashes.at(i) << std::endl;
  }
  return 0;
}

/*
Estimate Jaccard Index from Min hash approach

Uses min hashes for kmers of A and B
according to the paper: https://www.biorxiv.org/content/early/2017/09/04/184150

(Calculation similar to finding true jaccard
  EXCEPT that we use min hashes for kmers of A and B with Min Hash approach
  whereas in true jaccard calculation we use kmers of A and B)
*/
float estimate_jaccard_index_min_hash(vector<int> &hashes_A, vector<int> &hashes_B) {
  std::set<int> set_A(hashes_A.begin(), hashes_A.end());
  std::set<int> set_B(hashes_B.begin(), hashes_B.end());

  std::set<int> intersect;
  set_intersection(set_A.begin(),set_A.end(),set_B.begin(),set_B.end(),
                  std::inserter(intersect,intersect.begin()));
  int AUB = (set_A.size() + set_B.size() - intersect.size());
  float intersection = intersect.size() * 1.0;
  float jaccard_index = intersection / AUB;

  return jaccard_index;
}

/*
Calculate true jaccard index

True jaccard calculation uses kmers of A and B
according to author's calculation of true jaccard in: https://github.com/dkoslicki/MinHashMetagenomics
*/
float get_true_jaccard_index(vector<string> &kmers_A, vector<string> &kmers_B) {
  std::set<string> set_A(kmers_A.begin(), kmers_A.end());
  std::set<string> set_B(kmers_B.begin(), kmers_B.end());

  std::set<string> intersect;
  set_intersection(set_A.begin(),set_A.end(),set_B.begin(),set_B.end(),
                  std::inserter(intersect,intersect.begin()));
  int AUB = (set_A.size() + set_B.size() - intersect.size());
  float intersection = intersect.size() * 1.0;
  float jaccard_index = intersection / AUB;

  return jaccard_index;
}

/*
Estimate Jaccard Index from Containment Hash Approach
*/
float estimate_containment_index(int ksize,
  vector<string> &kmers_A,
  vector<string> &kmers_B,
  float false_positive_rate,
  int hash_count) {

  vector<string> min_kmers_A;
  //populate bloom filter with kmers of larger set

  //  - set params
  bloom_parameters parameters;
  parameters.projected_element_count = kmers_B.size();

  parameters.false_positive_probability = false_positive_rate * 1.0;
  parameters.random_seed = rand();
  if (!parameters) {
    std::cout << "Error: Invalid set of bloom filter parameters!" << std::endl;
    return 0.0;
  }
  parameters.compute_optimal_parameters();

  //  - instantiate bloom filter
  bloom_filter filter(parameters);

  //  - insert into bloom filter
  int filter_size = 0;
  for (int i=0; i < kmers_B.size(); i++) {
    string kmer = kmers_B.at(i);
    if (filter.contains(kmer)) continue;
    filter.insert(kmer);
    filter_size ++;
  }

  //find if kmers from set A which give min hash are in bloom filter
  vector<int> hashes_A;
  min_kmers_A = get_min_kmers_for_hashes (kmers_A, hashes_A, hash_count);
  float intersection = 0.0;
  for (int i=0; i < min_kmers_A.size(); i++) {
    if (! filter.contains(min_kmers_A.at(i))) continue;
    intersection = intersection + 1;
  }

  intersection = intersection - (false_positive_rate*hash_count);
  float containment_index = intersection / ((float)hash_count);
  float jaccard_index = (kmers_A.size() * containment_index)
  /(kmers_A.size() + filter_size - kmers_A.size() * containment_index);

  return jaccard_index;
}
