#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <zlib.h>
#include <filesystem>
#include <omp.h>
#include <cassert>
#include <zlib.h>
#include "gelato.hpp"

// Example usage
int main(int argc, char *argv[])
{
    // Boolean vars
    bool help = false, no_N = true, filter_referece = false, verbose = false;
    // Default kmer size
    unsigned k = 31;
    // Default minimum kmer appearence
    unsigned n = 2;
    // Default coverage value
    unsigned coverage = 1;
    // Input path string
    std::string ref_in, ref_list, seq_in, seq_list, assm_in, assm_list;
    // List of input path string
    std::string list_seq_in;
    // Output file name (graph & genome)
    std::string output, outputGraph;
    // Flag to set if we want to build the sequence graph or the genome
    bool buildGraph = false;
    bool buildGenome = false;
    bool isfastq = false;
    // Strings for outputfiles for the kmer counting
    std::string k_refs, k_reads;
    // Set the num of threads
    unsigned n_threads = 32; // std::thread::hardware_concurrency();
    // Program options
    char c;
    while ((c = getopt(argc, argv, "k:i:I:r:hn:R:p:Nfc:o:g:vs:S:a:A:O:")) >= 0)
    {
        if (c == 'h')
            help = true;
        else if (c == 'k'){
            k = atoi(optarg);
            printMessage("Maximum countable kmer length is 32");
            return 0;
        }
        else if (c == 's')
            k_reads = optarg;
        else if (c == 'S')
            k_refs = optarg;
        else if (c == 'i')
            seq_in = optarg;
        else if (c == 'I')
            seq_list = optarg;
        else if (c == 'r')
            ref_in = optarg;
        else if (c == 'n')
            n = atoi(optarg);
        else if (c == 'p')
            n_threads = atoi(optarg);
        else if (c == 'R')
            ref_list = optarg;
        else if (c == 'N')
            no_N = false;
        else if (c == 'f')
            filter_referece = true;
        else if (c == 'c')
            coverage = atoi(optarg);
        else if (c == 'a')
            assm_in = optarg;
        else if (c == 'A')
            assm_list = optarg;
        else if (c == 'v')
            verbose = true;
        else if (c == 'g')
        {
            std::string format(".gfa");
            outputGraph = optarg + format;
            buildGraph = true;
        }
        else if (c == 'o')
        {
            std::string format(".fasta");
            output = optarg + format;
            buildGraph = true;
            buildGenome = true;
        }
        else if (c == 'O')
        {
            std::string format(".fastq");
            output = optarg + format;
            buildGraph = true;
            buildGenome = true;
            isfastq = true;
        }
    }
    // Print out manual
    if (help)
    {
        std::cerr << "Usage: gelato [options]" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  -h               help" << std::endl;
        std::cerr << "  -v               verbose prints" << std::endl;
        std::cerr << "  -N               count kmers containing the N character (default do not count)" << std::endl;
        std::cerr << "  -s STR           save reads kmer counting info" << std::endl;
        std::cerr << "  -S STR           save references kmer counting info" << std::endl;
        std::cerr << "  -k INT           kmer len used for counting and graph building (default 31, max 32)" << std::endl;
        std::cerr << "  -r STR           single reference input file" << std::endl;
        std::cerr << "  -R STR           list of reference input files (on per line)" << std::endl;
        std::cerr << "  -i STR           single reads input file" << std::endl;
        std::cerr << "  -I STR           list of reads input files (on per line)" << std::endl;
        std::cerr << "  -a STR           additional assembly input file" << std::endl;
        std::cerr << "  -A STR           list of additional assembly input files (on per line)" << std::endl;
        std::cerr << "  -n INT           num of times a kmer has to appear to not be considered error" << std::endl;
        std::cerr << "  -p INT           set n. of threads to run gelato with (default 32)" << std::endl;
        std::cerr << "  -c INT           set read coverage (default 1)" << std::endl;
        std::cerr << "  -g STR           output filename for graph (gfa format and extension are automatically applied)" << std::endl;
        std::cerr << "  -o STR           output filename for genome (fasta format and extension are automatically applied)" << std::endl;
        std::cerr << "  -O STR           output filename for genome (fastq format and extension are automatically applied)" << std::endl;
        std::cerr << "Note that gelato always counts for canonical k-mers (the count is the combined number of occurrences of both a k-mer and it reverse complement)" << std::endl;
        return 0;
    }
    // Define the input sequences to compute the kmer counting from
    std::vector<std::string> reads, refs, supp_assm;

    if (seq_in.empty() && seq_list.empty())
    {
        printMessage("Missing read input file");
        return 0;
    }
    if (ref_in.empty() && ref_list.empty())
    {
        printMessage("Missing reference input file");
        return 0;
    }
    // if (assm_in.empty() && assm_list.empty())
    // {
    //     printMessage("Missing support assembly file");
    //     return 0;
    // }
    /*
     * Declare a filestream
     */
    std::ifstream inFile;

    /*
     * Parse the sequence file(s)
     */
    if (!seq_list.empty())
    {
        inFile.open(seq_list);
        while (std::getline(inFile, seq_in))
            parseFile(seq_in, reads);
    }
    else
        parseFile(seq_in, reads);

    if (!seq_list.empty())
        inFile.close();

    printMessage("Finished reading sequence inputs, got " + std::to_string(reads.size()) + " sequences");
    printMessage("Now counting reads kmers");
    if (buildGraph)
        printMessage("And constructing the kmer graph");

    /*
     * Declare hash tables and graph, count read kmers first and then proceed to parse refs
     */
    hash_k readsHashtable;
    graph_k readsGraph;
    std::vector<hash_k> refHashTables;
    std::vector<std::string> refPaths;
    float max_count = 0, min_count = 0;

    kMerCountingParallel(k, reads, n_threads, readsHashtable, readsGraph, no_N, buildGraph, n, coverage, max_count, min_count);

    printMessage(std::to_string(max_count));
    printMessage(std::to_string(min_count));

    unsigned id = 0;
    // Parse the reference file(s)
    if (!ref_list.empty())
    {
        inFile.open(ref_list);
        while (std::getline(inFile, ref_in))
        {
            refPaths.push_back(ref_in);
            parseFile(ref_in, refs);
            id++;
        }
    }
    else
    {
        refPaths.push_back(ref_in);
        parseFile(ref_in, refs);
    }
    // Resize the has table vector to the number of actual reference genomes
    refHashTables.resize(refs.size());

    printMessage("Finished reading reference inputs, got " + std::to_string(refs.size()) + " references");
    printMessage("Now counting reference kmers");

    // Parallelize kmer counting using multiple threads (1 per genome)
#pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < refs.size(); i++)
    {
        countSeqKmers(k, refs[i], refHashTables[i], no_N);
        if (filter_referece)
        {
            // filter reference hashmap
            if (i == 0) // so that is printed only once and not by all threads
                printMessage("Filtering ref");
            filterMap(refHashTables[i], n, coverage);
        }
    }

    if (!ref_list.empty())
        inFile.close();

    printMessage("Compute similarity between reads and reference datasets");
    // Compute similarity
    std::vector<float> similarity_score(refs.size());
    computeSimilarityRefReads(readsHashtable, refHashTables, similarity_score, k, n_threads);

    // Compute the most similar dataset to the set of reads
    std::vector<std::pair<unsigned, float>> ids_similar;
    for (unsigned i = 0; i < similarity_score.size(); i++)
        ids_similar.push_back(std::make_pair(i, similarity_score[i]));

    std::sort(ids_similar.begin(), ids_similar.end(), [](auto &left, auto &right)
              { return left.second > right.second; });

#ifndef RANDOM_REF
    unsigned closest_ref_id = 0;
    std::cout << refPaths[ids_similar[closest_ref_id].first] << std::endl;
#else
    srand(time(NULL));
    unsigned closest_ref_id = rand() % refs.size();
    std::cout << refPaths[ids_similar[closest_ref_id].first] << std::endl;
#endif

    if (verbose)
    {
        printMessage("-------------------------------------------------------------------");
        for (unsigned i = 0; i < similarity_score.size(); i++)
            printMessage("Path/Reference: " + refPaths[ids_similar[i].first] + " score: " + std::to_string(ids_similar[i].second));
        printMessage("-------------------------------------------------------------------");
    }

    // TODO: delete unecessary refs hash and plain variables

    // Parse support assembly file if present
    printMessage("Parsing support assembly file");

    if (!assm_list.empty())
    {
        inFile.open(assm_list);
        while (std::getline(inFile, assm_in))
            parseFile(assm_in, supp_assm);
    }
    else if(!assm_in.empty())
        parseFile(assm_in, supp_assm);

    if (buildGraph && !outputGraph.empty())
    {
        printMessage("Saving the sequence k-mer DBG graph");
        // TODO: print out waring in RED that tells the user the fixed kmer length and other info of the saved GFA
        // TODO: speed up GFA save
        graphToGFA(readsHashtable, readsGraph, outputGraph, k);
    }
    if (buildGenome)
    {
        printMessage("Generating the genome");
        refGuidedGenAssembly(readsHashtable, refHashTables[ids_similar[closest_ref_id].first], readsGraph, refs[ids_similar[closest_ref_id].first], supp_assm, output, k, isfastq, max_count, min_count);
        std::cerr << std::endl;
        printMessage("Finished genome generation");
    }
    return 0;
}
