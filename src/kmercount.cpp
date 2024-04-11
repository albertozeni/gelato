#include "kmercount.hpp"
#include "gelato.hpp"

bool searchKmerInAssembly(std::string kmer1, std::string kmer2, std::vector<std::string> &assembly, std::string &patch, int kmer_len)
{
    for (size_t i = 0; i < assembly.size(); i++)
    {
        std::size_t pos1 = assembly[i].find(kmer1);
        std::size_t pos2;
        if (pos1 != std::string::npos)
        {
            pos2 = assembly[i].find(kmer2, pos1);
            if (pos2 != std::string::npos)
            {
                patch = assembly[i].substr(pos1 + kmer_len, pos2 + kmer_len - (pos1 + kmer_len));
                if (patch.find("N") != std::string::npos)
                    continue;
                else
                    return true;
            }
        }
    }
    return false;
}

// Merge a local hash table into a global hash table
void mergeLocalHashTable(hash_k &localHashTable, hash_k &globalHashTable, float &max, float &min)
{
    for (const auto &[kmer, count] : localHashTable){
        float to_save = globalHashTable[kmer] + count;
        globalHashTable[kmer] = to_save;
        max = to_save > max ? to_save : max;
        min = to_save < min ? to_save : min;
    }
    // Clear the space of the temp hash table
    localHashTable.clear();
}

// Merge a local graph table into a global graph table
void mergeLocalGraph(graph_k &localGraph, graph_k &globalGraph)
{
    for (const auto &[kmer, hotbit] : localGraph)
    {
        if (globalGraph.count(kmer) == 0)
            globalGraph[kmer] = hotbit;
        else
            globalGraph[kmer] |= hotbit;
    }
    // Clear the space of the temp hash table
    localGraph.clear();
}

// Count the kmers for the sequence and also account for building the seq graph
void countSeqKmersAndBuildGraph(int k, std::string sequence, hash_k &localHashTable,
                                graph_k &localGraph, bool no_N)
{
    // Set all char of the sequence to upper case
    std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
    assert(sequence.length() >= k);
    for (unsigned i = 0; i < sequence.length() - k + 1; i++)
    {
        std::string kmer = sequence.substr(i, k);
        // std::cerr << kmer << std::endl;
        std::string kmer_rev = rev_comp(kmer);
        /*
         * If the kmer is valid, then it's added to the hashmap
         */
        if (no_N && kmer.find('N') < kmer.length())
            continue;
        std::uint64_t kmer_coded[2];
        encode_kmer(kmer, kmer_coded);
        localHashTable[kmer_coded[0] < kmer_coded[1] ? kmer_coded[0] : kmer_coded[1]]++;

        /*
         * Also add the kmer to the graph with no outgoing edges
         */
        if (localGraph.count(kmer_coded[0]) == 0)
        {
            localGraph[kmer_coded[0]] = 0;
            localGraph[kmer_coded[1]] = 0;
        }

        if (i < sequence.length() - k)
        {
            std::string kmer_next = sequence.substr(i + 1, k);
            if (no_N && kmer_next.find('N') < kmer_next.length())
                continue;
            /*
             * Add kmer to kmer next connection
             */
            std::uint64_t kmer_next_coded[2];
            encode_kmer(kmer_next, kmer_next_coded);
            localGraph[kmer_coded[0]] |= to_hotbit(kmer_next[k - 1]);
            /*
             * Add the same connection to the reverse complemented kmers
             */
            localGraph[kmer_next_coded[1]] |= to_hotbit(kmer_rev[k - 1]);
        }
    }
}

// Count the kmers for a sequence
void countSeqKmers(int k, std::string sequence, hash_k &localHashTable, bool no_N)
{
    // Set all char of the sequence to upper case
    std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
    assert(sequence.length() >= k);
    for (unsigned i = 0; i < sequence.length() - k + 1; i++)
    {
        std::string kmer = sequence.substr(i, k);
        /*
         * If the kmer is valid, then it's added to the hashmap
         */
        if (no_N && kmer.find('N') < kmer.length())
            continue;
        std::uint64_t kmer_coded = encode_kmer_minimum(kmer);
        localHashTable[kmer_coded]++;
    }
}

// Build a local hash table of k-mers and their adjacent k-mers from a chunk of input sequences
void buildLocalHashTable(int k, std::vector<std::string> const &sequences, hash_k &localHashTable,
                         graph_k &localGraph, bool no_N, bool buildGraph)
{
    // Create a sliding window of size k for each sequence and add each k-mer to the local hash table
    if (buildGraph)
        for (const auto &seq : sequences)
            if(seq.length()>=k)
                countSeqKmersAndBuildGraph(k, seq, localHashTable, localGraph, no_N);
    else
        for (const auto &seq : sequences)
            if(seq.length()>=k)
                countSeqKmers(k, seq, localHashTable, no_N);
}

// Count kmers using multiple threads
void kMerCountingParallel(int k, std::vector<std::string> &sequences, int n_threads, hash_k &globalHashTable,
                          graph_k &globalGraph, bool no_N, bool buildGraph, unsigned n, unsigned coverage, float &max, float &min)
{
    // Divide the input sequences into chunks (1 chunk per thread)
    std::vector<std::vector<std::string>> chunks(n_threads);
    for (int i = 0; i < sequences.size(); i++)
        chunks[i % n_threads].push_back(sequences[i]);

    // Build the local hash tables using multiple threads
    std::vector<hash_k> localHashTables(n_threads);

    // Build the graphs using multiple threads
    std::vector<graph_k> localGraph(n_threads);

#pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < n_threads; i++)
        buildLocalHashTable(k, chunks[i], localHashTables[i], localGraph[i], no_N, buildGraph);

    printMessage("Done counting, merging");

    // Merge the local hash tables into a global hash table
    for (int i = 0; i < n_threads; i++)
        mergeLocalHashTable(localHashTables[i], globalHashTable, max, min);

    // Account for kmer coverage
    adaptCoverage(globalHashTable, coverage);
    // Filter out error kmers
    printMessage("Filter read hash table");
    filterMap(globalHashTable, n, coverage);

    if (buildGraph)
    {
        printMessage("Filter graph");
#pragma omp parallel for num_threads(n_threads)
        for (int i = 0; i < n_threads; i++)
            filterGraph(globalHashTable, localGraph[i], k);
    }

    // Reserve hash table for graph
    if (buildGraph)
    {
        printMessage("Merge graphs");
        for (int i = 0; i < n_threads; i++)
            mergeLocalGraph(localGraph[i], globalGraph);
    }
}

void filterMap(hash_k &m, unsigned n, unsigned coverage)
{
    for (auto it = begin(m); it != end(m);)
    {
        float min = ((float)n) / ((float)coverage);
        if (it->second < min)
            it = m.erase(it);
        else
            it++;
    }
}

void filterGraph(hash_k &m, graph_k &g, unsigned k)
{
    for (auto it = begin(g); it != end(g);)
    {
        /*
         * Kmer 1 is always the same, regardless of link
         * If kmer 1 is not in the kmer set, we delete the node in the graph directly
         */
        std::string kmer1_s = decode_kmer(it->first, k);
        std::uint64_t kmer1_min = encode_kmer_minimum(kmer1_s);
        if (m.count(kmer1_min) == 0)
            it = g.erase(it);
        else
            it++;
    }

    for (auto &node : g)
    {
        /*
         * For all the remaining kmers check all kmer connections
         */
        for (std::uint8_t i = 1; i <= 8; i <<= 1)
        {
            if (node.second & i)
            {
                std::string kmer2_s = decode_kmer(node.first, k);
                kmer2_s = kmer2_s.substr(1, kmer2_s.length() - 1);
                kmer2_s.push_back(from_hotbit(i));
                std::uint64_t kmer2_min = encode_kmer_minimum(kmer2_s);
                if (m.count(kmer2_min) == 0)
                    node.second = node.second - i;
            }
        }
    }
}

std::string rev_comp(std::string kmer)
{
    std::string kmer_rev = kmer;
    std::reverse(kmer_rev.begin(), kmer_rev.end());
    for (unsigned i = 0; i < kmer_rev.length(); i++)
    {
        if (kmer_rev[i] == 'A')
            kmer_rev[i] = 'T';
        else if (kmer_rev[i] == 'T')
            kmer_rev[i] = 'A';
        else if (kmer_rev[i] == 'C')
            kmer_rev[i] = 'G';
        else if (kmer_rev[i] == 'G')
            kmer_rev[i] = 'C';
    }
    return kmer_rev;
}

/*
 * Alphabet for the next kmer letter for the graph: (hotbit)
 * 0001 = A
 * 0010 = T
 * 0100 = C
 * 1000 = G
 */
std::uint8_t to_hotbit(char c)
{
    if (c == 'A')
        return 1;
    else if (c == 'T')
        return 2;
    else if (c == 'C')
        return 4;
    else if (c == 'G')
        return 8;
    else // this state exists only to check we don't encounter erroneus characters
        exit(-1);
}

char from_hotbit(std::uint8_t u)
{
    if (u == 1)
        return 'A';
    else if (u == 2)
        return 'T';
    else if (u == 4)
        return 'C';
    else if (u == 8)
        return 'G';
    else // this state exists only to check we don't encounter erroneus characters
        exit(-1);
}

/*
 * Alphabet for coding and decoding kmers (counting)
 * 0 = A
 * 1 = T
 * 2 = C
 * 3 = G
 */
std::uint8_t encode(char c)
{
    if (c == 'A')
        return 0;
    else if (c == 'T')
        return 1;
    else if (c == 'C')
        return 2;
    else if (c == 'G')
        return 3;
    else // this state exists only to check we don't encounter erroneus characters
        exit(-1);
}

char decode(std::uint8_t u)
{
    if (u == 0)
        return 'A';
    else if (u == 1)
        return 'T';
    else if (u == 2)
        return 'C';
    else if (u == 3)
        return 'G';
    else // this state exists only to check we don't encounter erroneus characters
        exit(-1);
}

void encode_kmer(std::string kmer, std::uint64_t *encoded_kmer)
{
    std::uint64_t encoded_c = 0;
    std::uint64_t mask = (1ULL << (kmer.length() * 2)) - 1, shift = (kmer.length() - 1) * 2;
    for (unsigned i = 0; i < kmer.length(); ++i)
    {
        encoded_c = encode(kmer[i]);
        encoded_kmer[0] = (encoded_kmer[0] << 2 | encoded_c) & mask;                 // forward strand
        encoded_kmer[1] = encoded_kmer[1] >> 2 | (uint64_t)(3 - encoded_c) << shift; // reverse strand
    }
}

std::uint64_t encode_kmer_minimum(std::string kmer)
{
    std::uint64_t encoded_kmer[2];
    std::uint64_t encoded_c = 0;
    std::uint64_t mask = (1ULL << (kmer.length() * 2)) - 1, shift = (kmer.length() - 1) * 2;
    for (unsigned i = 0; i < kmer.length(); ++i)
    {
        encoded_c = encode(kmer[i]);
        encoded_kmer[0] = (encoded_kmer[0] << 2 | encoded_c) & mask;                 // forward strand
        encoded_kmer[1] = encoded_kmer[1] >> 2 | (uint64_t)(3 - encoded_c) << shift; // reverse strand
    }
    return encoded_kmer[0] < encoded_kmer[1] ? encoded_kmer[0] : encoded_kmer[1];
}

std::uint64_t reencode_kmer_minimum(std::uint64_t kmer_en, std::uint32_t k)
{
    std::string kmer = decode_kmer(kmer_en, k);
    std::uint64_t encoded_kmer[2];
    std::uint64_t encoded_c = 0;
    std::uint64_t mask = (1ULL << (kmer.length() * 2)) - 1, shift = (kmer.length() - 1) * 2;
    for (unsigned i = 0; i < kmer.length(); ++i)
    {
        encoded_c = encode(kmer[i]);
        encoded_kmer[0] = (encoded_kmer[0] << 2 | encoded_c) & mask;                 // forward strand
        encoded_kmer[1] = encoded_kmer[1] >> 2 | (uint64_t)(3 - encoded_c) << shift; // reverse strand
    }
    return encoded_kmer[0] < encoded_kmer[1] ? encoded_kmer[0] : encoded_kmer[1];
}

std::uint64_t reencode_kmer_minimum_fast(std::uint64_t kmer_en, std::uint32_t k)
{
    std::uint64_t kmer_rev_en = 0, shift = (k - 1) * 2;
    for (unsigned i = 0; i < k; ++i)
    {
        kmer_rev_en = kmer_rev_en >> 2 | (uint64_t)(3 - ((kmer_en >> (k * 2 - i * 2)) & 3)) << shift;
    }
    return kmer_en < kmer_rev_en ? kmer_en : kmer_rev_en;
}

std::uint64_t encode_kmer_for(std::string kmer)
{
    std::uint64_t encoded_kmer;
    std::uint64_t encoded_c = 0;
    std::uint64_t mask = (1ULL << (kmer.length() * 2)) - 1, shift = (kmer.length() - 1) * 2;
    for (unsigned i = 0; i < kmer.length(); ++i)
    {
        encoded_c = encode(kmer[i]);
        encoded_kmer = (encoded_kmer << 2 | encoded_c) & mask; // forward strand
    }
    return encoded_kmer;
}

std::string decode_kmer(std::uint64_t encode_kmer, unsigned k)
{
    std::string s;
    for (size_t i = 0; i < k; i++)
    {
        std::uint8_t c = encode_kmer & 3;
        s.insert(0, 1, decode(c));
        encode_kmer >>= 2;
    }
    return s;
}

void adaptCoverage(hash_k &m, unsigned c)
{
    for (auto &[kmer, count] : m)
    {
        float new_count = (count) / ((float)c);
        m[kmer] = new_count;
    }
}