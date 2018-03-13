#include <iostream>
#include <vector>
#define PRIME 9999999999971
#include <set>
#include "include/MurmurHash3.h"
#include "include/bloom_filter.hpp"
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

std::string read_input(string filename);

void create_kmers(int ksize, string genome, vector<string> &kmers);

int create_random_hashes(vector<int> &random_hashes, int hash_count);

vector<string> get_hash_from_kmers(vector<string> &kmers,
  int index_random_hashes, int hash_count);

int get_hashes(vector<string> &kmers,
  vector<int> &hashes, int hash_count);

vector<string> get_min_kmers_for_hashes (vector<string> &kmers,
  vector<int> &hashes, int hash_count);

int print_hashes(vector<int> &hashes);

float estimate_jaccard_index_min_hash(vector<int> &hashes_A, vector<int> &hashes_B);

float get_true_jaccard_index(vector<string> &hashes_A, vector<string> &hashes_B);

float estimate_containment_index(int ksize,
  vector<string> &kmers_A,
  vector<string> &kmers_B,
  float false_positive_rate,
  int hash_count);
