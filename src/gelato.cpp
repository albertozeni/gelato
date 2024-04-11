#include "gelato.hpp"
#include <queue>
#define LENS 1000000

void printMessage(std::string s)
{
    std::cerr << "[GELATO 🍦] :: " << s << std::endl;
}

void printMap(hash_k const &m)
{
    for (auto const &pair : m)
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
}

void parseFile(std::string inputFile, std::vector<std::string> &sequences)
{
    std::filesystem::path filePath = inputFile;
    bool data = false;

    if (filePath.extension() == ".gz"){
        std::string tmp = inputFile;
        for(int i = 0; i < 3; i++) tmp.pop_back();
        filePath = tmp;
    }

    if (filePath.extension() == ".fastq" || filePath.extension() == ".fq")
    {
        char *tmp_s = new char[LENS];
        gzFile inFile = gzopen(inputFile.data(), "r");
        while (0!=gzgets(inFile,tmp_s,LENS))
        {
            std::string s(tmp_s);
            s.pop_back();
            if (s[0] == '@')
            {
                data = true;
                continue;
            }
            char *tmp_add_s = new char[LENS];
            while (data && 0!=gzgets(inFile, tmp_add_s, LENS) && tmp_add_s[0] != '+'){
                std::string add_s(tmp_add_s);
                add_s.pop_back();
                s += add_s;
            }
            delete[] tmp_add_s;
            if (data)
                sequences.push_back(s);
            data = false;
        }
        delete[] tmp_s;
        gzclose(inFile);
    }
    else if (filePath.extension() == ".fasta" || filePath.extension() == ".fa")
    {
        char *tmp_s = new char[LENS];
        gzFile inFile = gzopen(inputFile.data(), "r" );
        while (0!=gzgets(inFile,tmp_s,LENS))
        {
            std::string s(tmp_s);
            s.pop_back();
            if (s[0] == '>')
                continue;
            std::string add_s;
            char *tmp_add_s = new char[LENS];
            while (0!=gzgets(inFile,tmp_add_s,LENS) && add_s[0]!='>')
            {
                std::string add_s(tmp_add_s);
                add_s.pop_back();
                s += add_s;
            }
            delete[] tmp_add_s;
            sequences.push_back(s);
        }
        delete[] tmp_s;
        gzclose(inFile);
    }
}

void graphToGFA(hash_k &globalHash, graph_k &globalGraph, std::string file_name, unsigned kmer_len)
{
    std::ofstream outFile;
    outFile.open(file_name);
    outFile << "H\tVN:Z:1.2\n";
    // Save the segments
    for (auto const &pair : globalHash)
        outFile << "S\t" << pair.first << "\t" << decode_kmer(pair.first, kmer_len) << "\n";
    // // Save the links
    for (auto const &node : globalGraph)
    {
        // kmer 1 is always the same, regardless of link
        std::string tmp = decode_kmer(node.first, kmer_len);
        std::uint64_t kmer1[2];
        std::string kmer1_s = tmp;
        encode_kmer(kmer1_s, kmer1);
        char kmer1_orientation = kmer1[0] < kmer1[1] ? '+' : '-';
        kmer1_s = kmer1[0] < kmer1[1] ? kmer1_s : rev_comp(kmer1_s);
        std::uint64_t min_k1 = kmer1[0] < kmer1[1] ? kmer1[0] : kmer1[1];
        // iterate and write all kmer 1 connections
        for (std::uint8_t i = 1; i <= 8; i <<= 1)
        {
            if (i & node.second)
            {
                std::string kmer2_s = tmp.substr(1, kmer_len - 1);
                kmer2_s.push_back(from_hotbit(i));
                // assert(kmer2.length() == kmer_len);
                std::uint64_t kmer2[2];
                encode_kmer(kmer2_s, kmer2);
                char kmer2_orientation = kmer2[0] < kmer2[1] ? '+' : '-';
                kmer2_s = kmer2[0] < kmer2[1] ? kmer2_s : rev_comp(kmer2_s);
                std::uint64_t min_k2 = kmer2[0] < kmer2[1] ? kmer2[0] : kmer2[1];
                outFile << "L\t" << min_k1 << "\t" << kmer1_orientation << "\t" << min_k2 << "\t" << kmer2_orientation << "\t" << kmer_len - 1 << "M"
                        << "\n";
            }
        }
    }
    outFile.close();
}

void computeSimilarityRefReads(hash_k &readsHash, std::vector<hash_k> &refsHash, std::vector<float> &similarityScores, unsigned k, unsigned nThreads)
{
#pragma omp parallel for num_threads(nThreads)
    for (int i = 0; i < refsHash.size(); i++)
    {
        float m = 0; // number of ref kmers in common with the reads
        for (auto const &pair : refsHash[i])
            if (readsHash.count(pair.first))
                m++;

        float n = readsHash.size(); // number of reads kmers

        similarityScores[i] = pow(m / n, (1 / (double(k))));
    }
}