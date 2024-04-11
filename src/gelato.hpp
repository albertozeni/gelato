#ifndef GELATO_HPP_
#define GELATO_HPP_

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <zlib.h>
#include <filesystem>
#include <omp.h>
#include <cassert>
#include <zlib.h>
#include <queue>
#include <cmath>

#include "kmercount.hpp"
#include "graphexplore.hpp"
#include "datatypes.hpp"

void printMap(hash_k const &m);

std::string rev_comp(std::string kmer);

void adaptCoverage(hash_k &m, unsigned c);

void parseFile(std::string inputFile, std::vector<std::string> &sequences);

void graphToGFA(hash_k &globalHash, graph_k &globalGraph, std::string file_name, unsigned kmer_len);

void computeSimilarityRefReads(hash_k &readsHash, std::vector<hash_k> &refsHash, std::vector<float> &similarityScores, unsigned k, unsigned nThreads);

void printMessage(std::string s);

#endif // GELATO_HPP_