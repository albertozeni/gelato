#ifndef KMERCOUNT_HPP_
#define KMERCOUNT_HPP_

#include <string>
#include <vector>
#include <omp.h>

#include "robin_hood.h"
#include "datatypes.hpp"

void countSeqKmers(int k, std::string sequence, hash_k &localHashTable, bool no_N);

bool searchKmerInAssembly(std::string kmer1, std::string kmer2, std::vector<std::string> &assembly, std::string &patch, int kmer_len);

void kMerCountingParallel(int k, std::vector<std::string> &sequences, int n_threads, hash_k &globalHashTable,
                          graph_k &globalGraph, bool no_N, bool buildGraph, unsigned n, unsigned coverage, float &max, float &min);

void filterMap(hash_k &m, unsigned n, unsigned coverage);

void filterGraph(hash_k &m, graph_k &g, unsigned k);

std::string rev_comp(std::string kmer);

std::uint8_t to_hotbit(char c);

char from_hotbit(std::uint8_t u);

std::uint8_t encode(char c);

char decode(std::uint8_t u);

void encode_kmer(std::string kmer, std::uint64_t *encoded_kmer);

std::uint64_t encode_kmer_minimum(std::string kmer);

std::uint64_t reencode_kmer_minimum(std::uint64_t kmer_en, std::uint32_t k);

std::uint64_t reencode_kmer_minimum_fast(std::uint64_t kmer_en, std::uint32_t k);

std::uint64_t encode_kmer_for(std::string kmer);

std::string decode_kmer(std::uint64_t encode_kmer, unsigned k);

void adaptCoverage(hash_k &m, unsigned c);

#endif //KMERCOUNT_HPP_