#ifndef GRAPHEXPLORE_HPP_
#define GRAPHEXPLORE_HPP_

#include <vector>
#include <string>
#include "datatypes.hpp"

void refGuidedGenAssembly(hash_k &readsHash, hash_k &refHash, graph_k &readsGraph, const std::string referenceGenome, std::vector<std::string> &assembly_support,
                    std::string file_name, unsigned kmer_len, bool isfastq, float max_count, float min_count); //, unsigned support)

#endif //GRAPHEXPLORE_HPP_