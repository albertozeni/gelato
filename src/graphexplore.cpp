#include <cassert>
#include <queue>

#include "graphexplore.hpp"
#include "kmercount.hpp"
#include "gelato.hpp"

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
#define QUALITY_CHAR_NO '!' - '~'
#define EMULATED_QUALITY 'I'

void printProgress(double i, double tot)
{
    double percentage = i / tot;
    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    fprintf(stderr, "\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

std::uint64_t hash_func(const std::vector<std::uint64_t> &v)
{
    std::uint64_t answer = v.size();
    for (auto x : v)
    {
        x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
        x = x ^ (x >> 31);
        answer ^= x + 0x9e3779b9 + (answer << 6) + (answer >> 2); // this last line might not be ideal since it is thought for 32 bits vectors
    }
    return answer;
}

int get_length_similarity_score(std::vector<std::vector<std::uint64_t>> &valid_paths, std::uint32_t max_len)
{

    /*
     * Look into all path candidates and choose the one whose length is closer to the missing one in the reference assembly
     */

    int similarity_score = (valid_paths[0].size() - 1) - max_len >= 0 ? (valid_paths[0].size() - 1) - max_len : max_len - (valid_paths[0].size() - 1);
    for (int i = 0; i < valid_paths.size(); i++)
    {
        std::vector<std::uint64_t> path_tmp = valid_paths[i];
        // smaller is better
        int sim_score_tmp = (valid_paths[i].size() - 1) - max_len >= 0 ? (valid_paths[i].size() - 1) - max_len : max_len - (valid_paths[i].size() - 1);
        if (sim_score_tmp <= similarity_score)
        {
            if (sim_score_tmp < similarity_score)
            {
                similarity_score = sim_score_tmp;
            }
        }
    }

    return similarity_score;

    /*
     * End length method
     */
}

void evaluate_paths(hash_k &h, std::vector<std::vector<std::uint64_t>> &valid_paths, std::uint32_t max_len, std::uint32_t kmer_len, std::uint64_t *selected_end, std::string &patch)
{
    // TODO: it might make sense to have the length relative to the selected end, as different anchors are at different legnths with respect to the start kmer
    std::vector<std::vector<std::uint64_t>> optimal_paths;
    /*
     * Look into all path candidates and choose the one whose length is closer to the missing one in the reference assembly
     */
    int similarity_score = (valid_paths[0].size() - 1) - max_len >= 0 ? (valid_paths[0].size() - 1) - max_len : max_len - (valid_paths[0].size() - 1);
    for (int i = 0; i < valid_paths.size(); i++)
    {
        std::vector<std::uint64_t> path_tmp = valid_paths[i];
        // the smaller the similarity score is, the closer is the path to the one searched
        int sim_score_tmp = (valid_paths[i].size() - 1) - max_len >= 0 ? (valid_paths[i].size() - 1) - max_len : max_len - (valid_paths[i].size() - 1);
        if (sim_score_tmp <= similarity_score)
        {
            if (sim_score_tmp < similarity_score)
            {
                similarity_score = sim_score_tmp;
                optimal_paths.clear();
            }
            optimal_paths.push_back(valid_paths[i]);
        }
    }
    /*
     * End length method
     */

    /*
     * Look into all path candidates and choose the coverage is higher
     */
    similarity_score = 0;
    for (int i = 0; i < optimal_paths.size(); i++)
    {
        std::vector<std::uint64_t> path_tmp = optimal_paths[i];
        std::string patch_tmp;
        float sim_score_tmp = 0;
        /*
         * For every path we ignore the first char,
         * as it comes from a previously found kmer and has been already accounted for
         */
        for (int j = 1; j < path_tmp.size(); j++)
        {
            sim_score_tmp += h[reencode_kmer_minimum(path_tmp[j], kmer_len)];
            // assert(reencode_kmer_minimum(path_tmp[j], kmer_len) == reencode_kmer_minimum_fast(path_tmp[j], kmer_len));
        }

        if (sim_score_tmp > similarity_score)
        {
            similarity_score = sim_score_tmp;
            for (int j = 1; j < path_tmp.size(); j++)
            {
                std::string kmer_tmp = decode_kmer(path_tmp[j], kmer_len);
                patch_tmp.push_back(kmer_tmp[kmer_len - 1]);
            }
            patch = patch_tmp;
            *selected_end = path_tmp[path_tmp.size() - 1];
        }
    }
}

bool compatiblePath(std::vector<std::uint64_t> &pathToCheck, std::vector<std::uint64_t> &alternativeBranch)
{
    std::uint32_t i, j;
    for (i = 0; i < pathToCheck.size() && pathToCheck[i] == alternativeBranch[0]; i++)
        ;
    for (j = i + 1; j < pathToCheck.size() && pathToCheck[j] == alternativeBranch[alternativeBranch.size() - 1]; j++)
        ;
    return (i != pathToCheck.size() && j != pathToCheck.size());
}

std::vector<std::uint64_t> replaceBranch(std::vector<std::uint64_t> &pathToCheck, std::vector<std::uint64_t> &alternativeBranch)
{
    std::uint32_t i, j;
    std::vector<std::uint64_t> new_path = pathToCheck;
    for (i = 0; i < new_path.size() && new_path[i] == alternativeBranch[0]; i++)
        ;
    for (j = i + 1; j < new_path.size() && new_path[j] == alternativeBranch[alternativeBranch.size() - 1]; j++)
        ;
    new_path.erase(new_path.begin() + i, new_path.begin() + j);
    new_path.insert(new_path.begin() + i, alternativeBranch.begin(), alternativeBranch.end());
    return new_path;
}

void generate_alternative_paths(std::vector<std::vector<std::uint64_t>> &possible_paths, std::vector<std::vector<std::uint64_t>> &alternative_branches, std::uint32_t max_len)
{
    assert(possible_paths.size() > 0);
    float best_similarity = (float)get_length_similarity_score(possible_paths, max_len);
    float eval = best_similarity / 10;
    eval += best_similarity;

    std::vector<std::vector<std::uint64_t>> alternative_path_chosen(possible_paths.size()); // contains all the ids of the various branch chosen per particular path
    set_k alternative_paths_hash_set;

    for (std::uint32_t i = 0; i < possible_paths.size(); i++)
        alternative_paths_hash_set.insert(hash_func(possible_paths[i]));

    for (std::uint32_t i = 0; i < possible_paths.size(); i++)
    {
        // std::cerr << possible_paths.size() << " " << i << std::endl;
        std::vector<std::uint64_t> curr_path = possible_paths[i];
        for (std::uint32_t j = 0; j < alternative_branches.size(); j++)
        {
            if (!std::count(alternative_path_chosen[i].begin(), alternative_path_chosen[i].end(), j) && compatiblePath(curr_path, alternative_branches[j]))
            {
                std::vector<std::uint64_t> new_path = replaceBranch(curr_path, alternative_branches[j]);
                float similarity_score = (new_path.size() - 1) - max_len >= 0 ? (new_path.size() - 1) - max_len : max_len - (new_path.size() - 1);

                if ((similarity_score < eval) && !alternative_paths_hash_set.count(hash_func(new_path))) // check if this path has been already created (somehow)
                {
                    possible_paths.push_back(new_path);
                    alternative_path_chosen.push_back(alternative_path_chosen[i]);
                    alternative_path_chosen[alternative_path_chosen.size() - 1].push_back(j);
                }
            }
        }
    }
}

bool searchPaths(graph_k &g, hash_k &h, std::uint64_t start, set_k &ends, std::string &patch, std::uint32_t kmer_len, std::uint32_t max_len, std::uint64_t *selected_end)
{
    // std::deque<std::uint64_t> queue;
    std::vector<std::vector<std::uint64_t>> valid_paths;
    std::vector<std::vector<std::uint64_t>> alternative_branches;
    std::deque<std::vector<std::uint64_t>> queue;

    // set_v possible_paths;
    // set_v alternative_branches;
    set_k visited_nodes;

    // hash_k visited_nodes;

    size_t curr_len = 0, curr_lev, next_lev;
    bool found = false;

    std::vector<std::uint64_t> path;
    path.push_back(start);
    queue.push_back(path);
    curr_lev = 1;
    next_lev = 0;

    while ((queue.size() > 0) && (curr_len < max_len * 2))
    {
        if (curr_lev == 0)
        {
            curr_len++;
            curr_lev = next_lev;
            next_lev = 0;
        }
        /*
         * Take the first path in the queue and continue to explore it,
         * the last node contains the last explored kmer.
         * From that we take the connected kmers, and continue to add new paths to the queue.
         */
        std::vector<std::uint64_t> path_to_explore = queue.front();
        std::uint64_t kmer_to_explore = path_to_explore[path_to_explore.size() - 1];
        std::uint8_t to_explore = g[kmer_to_explore];
        std::string kmer_to_explore_s = decode_kmer(kmer_to_explore, kmer_len);
        queue.pop_front();

        for (std::uint8_t i = 1; i <= 8; i <<= 1)
        {
            if (to_explore & i)
            {
                /*
                 * Observe the connected kmer, if it is an anchor than mark the path as valid,
                 * otherwise continue to explore
                 */
                std::string new_expl = kmer_to_explore_s.substr(1, kmer_len - 1);
                new_expl.push_back(from_hotbit(i));
                std::uint64_t new_expl_encoded = encode_kmer_for(new_expl);

                std::vector<std::uint64_t> new_path = path_to_explore;
                new_path.push_back(new_expl_encoded);
                if (visited_nodes.count(new_expl_encoded))
                {
                    if(alternative_branches.size() < MAX_START)
                        alternative_branches.push_back(new_path);
                    continue;
                }
                else
                    visited_nodes.insert(new_expl_encoded);
                if (ends.count(new_expl_encoded))
                {
                    /*
                     * Found at least one path that connects the prev found kmer and the anchor kmer
                     */
                    found = true;
                    valid_paths.push_back(new_path);
                }
                /*
                 * Put it back in queue as it might lead to another destination
                 */
                queue.push_back(new_path);
                next_lev++;
            }
        }
        curr_lev--;
    }

    if (found)
    {
        generate_alternative_paths(valid_paths, alternative_branches, max_len);
        evaluate_paths(h, valid_paths, max_len, kmer_len, selected_end, patch);
    }

    valid_paths.clear();
    alternative_branches.clear();

    return found;
}

std::uint64_t findAssemblyStart(hash_k &readsHash, hash_k &refHash, const std::string referenceGenome)
{
    std::uint32_t currentMax = 0;
    std::uint32_t arg_max = 0;
    while (true)
    {
        for (auto it = readsHash.cbegin(); it != readsHash.cend(); ++it)
        {
            if (it->second > currentMax && refHash.count(it->first) != 0)
            {
                arg_max = it->first;
                currentMax = it->second;
            }
        }
    }
    return arg_max;
}

inline char computeQuality(float max, float min, float k_count)
{
    // if(k_count > 0)
        return 'I';
    // else
    //     return 'J';
    // unsigned interval = (max - min) / QUALITY_CHAR_NO;
    // unsigned id = k_count;//(k_count - min);// interval;
    // // std::cerr << id << " ";
    // // std::cerr << k_count << std::endl;
    // return (char)('!' + id);
}

void refGuidedGenAssembly(hash_k &readsHash, hash_k &refHash, graph_k &readsGraph, const std::string referenceGenome, std::vector<std::string> &assembly_support,
                          std::string file_name, std::uint32_t kmer_len, bool isfastq, float max_count, float min_count) //, std::uint32_t support)
{
    std::deque<std::string> prev_found_kmers;
    std::deque<std::uint64_t> start_positions;

    /*
     * The out genome will always start and end with the first and last kmers of the chosen reference respectively
     */
    std::string outputGenome = referenceGenome.substr(0, kmer_len);
    std::string outputQual;

    /*
     * check if previous kmer is in the graph, this might not be the case for the very first kmer
     */
    bool first_kmer_in_hash = readsHash.count(encode_kmer_minimum(outputGenome));

    /*
     *  compute quality for the first kmer
     */
    char qual_first_kmer = '!';
    if (first_kmer_in_hash)
        qual_first_kmer = computeQuality(max_count, min_count, readsHash[encode_kmer_minimum(outputGenome)]);
    for (int i = 0; i < kmer_len; i++)
        outputQual.push_back(qual_first_kmer);

    /*
     * First kmer and first i value
     */
    prev_found_kmers.push_front(outputGenome);
    start_positions.push_front(0);

    /*
     *  Navigate the reference genome in kmer-size batches
     */
    for (std::uint32_t i = 1; i <= referenceGenome.length() - kmer_len * 2;)
    {
        // std::cerr << "Progression: " << i << std::endl;
        printProgress(i, referenceGenome.length() - kmer_len * 2);
        // float perc = (((float)i/referenceGenome.length() - kmer_len * 2);
        std::string curr_kmer = referenceGenome.substr(i, kmer_len);
        std::uint64_t saved_encoded_kmer = encode_kmer_minimum(curr_kmer);
        /*
         * If the kmer is also in the reads we mark its support on the output genome
         */
        if (readsHash.count(saved_encoded_kmer))
        {
            outputGenome.push_back(curr_kmer[kmer_len - 1]);
            outputQual.push_back(computeQuality(max_count, min_count, readsHash[encode_kmer_minimum(curr_kmer)]));
            prev_found_kmers.push_front(curr_kmer);
            start_positions.push_front(i);

            if (prev_found_kmers.size() > MAX_START)
            {
                prev_found_kmers.pop_back();
                start_positions.pop_back();
            }
            i++;
        }
        /*
         * If there is no support from the reads, we look into the supporting assembly
         * Otherwise we have to explore the reads graph and reattach to the next common kmer
         */
        else
        {
            std::uint32_t j = i;

            std::string patch;                     /* will contain the final patch to be set in the genome, regardless of the method used */
            std::vector<std::string> anchor_kmers; /* contains all the anchor kmers */
            set_k anchor_kmers_set;                /* hash map containing the possible ends in the graph, e.g. all the anchors */
            std::vector<std::uint32_t> end_positions;
            /*
             * Look in the reference and in the reads for the N next
             * common kmers to look for possible paths
             */
            for (; j <= referenceGenome.length() - kmer_len * 2 && anchor_kmers.size() < MAX_DEST; j++)
            {
                curr_kmer = referenceGenome.substr(j, kmer_len);
                saved_encoded_kmer = encode_kmer_minimum(curr_kmer);
                if (readsHash.count(saved_encoded_kmer))
                {
                    anchor_kmers.push_back(curr_kmer);
                    anchor_kmers_set.insert(encode_kmer_for(curr_kmer));
                    end_positions.push_back(j);
                }
            }

            /*
             *  First we search for the anchor kmer and the previous kmer inside the supporting assembly
             *  if we do not find them we move on searching the path in the reads graph
             *  If we found an anchor kmer, we search for a path in the reads graph starting from a previously found kmer,
             *  we search for the path with the most similar length (i.e. basepairs) between the two found kmers.
             */
            // TODO: remove the first_kmer_in_hash stuff, otherwise we never search in some rare instances
            if (anchor_kmers.size() > 0 && first_kmer_in_hash)
            {
                std::uint32_t selected_start_idx;

                std::vector<std::string> patch_candidates(prev_found_kmers.size());
                std::vector<bool> contig_found_start(prev_found_kmers.size());
                std::vector<std::uint32_t> selected_dest_per_start(prev_found_kmers.size());
                bool path_found_in_contigs = false;

                if (assembly_support.size() > 0) // check if there is any assembly to use
                {
#pragma omp parallel for
                    for (std::uint32_t m = 0; m < prev_found_kmers.size(); m++)
                    {
                        for (std::uint32_t l = 0; l < anchor_kmers.size(); l++)
                        {
                            contig_found_start[m] = searchKmerInAssembly(prev_found_kmers[m], anchor_kmers[l], std::ref(assembly_support), std::ref(patch_candidates[m]), kmer_len);
                            if (contig_found_start[m])
                            {
                                path_found_in_contigs |= contig_found_start[m];
                                selected_dest_per_start[m] = l;
                                break;
                            }
                        }
                    }
                }

                std::string graph_patch;
                std::vector<std::string> graph_candidate_patches(prev_found_kmers.size());
                std::vector<std::uint64_t> selected_ends(prev_found_kmers.size());
                std::vector<std::uint32_t> selected_anchor_idx(prev_found_kmers.size()); // value that indicates which end kmer has been taken for a particular start kmer
                std::vector<bool> path_found_start(prev_found_kmers.size());
                bool path_found_in_read_graph = false;
                
                assert(prev_found_kmers.size()<=8);

                if (!path_found_in_contigs)
                {
// TODO: check if this is correct as it might be that we have multiple identical kmer anchors
#pragma omp parallel for
                    for (std::uint32_t m = 0; m < prev_found_kmers.size(); m++)
                    {
                        path_found_start[m] = searchPaths(std::ref(readsGraph), std::ref(readsHash), encode_kmer_for(prev_found_kmers[m]),
                                                          anchor_kmers_set, std::ref(graph_candidate_patches[m]),
                                                          kmer_len, end_positions[anchor_kmers.size() - 1] - start_positions[m] + 1,
                                                          &selected_ends[m]);
                        if (path_found_start[m])
                            for (std::uint32_t l = 0; l < anchor_kmers.size(); l++)
                                selected_anchor_idx[m] = (encode_kmer_for(anchor_kmers[l]) == selected_ends[m]) ? l : selected_anchor_idx[m];
                    }

                    for (std::uint32_t m = 0; m < prev_found_kmers.size(); m++)
                    {
                        if (path_found_start[m])
                        {
                            path_found_in_read_graph = true;
                            graph_patch = graph_candidate_patches[m];
                            selected_start_idx = m;
                            break;
                        }
                    }
                }

                for (std::uint32_t m = 0; m < prev_found_kmers.size(); m++)
                {
                    if (path_found_in_contigs && contig_found_start[m])
                    {
                        patch = patch_candidates[m];
                        outputQual.resize(outputGenome.size() - (i - start_positions[m] - 1));
                        outputGenome.resize(outputGenome.size() - (i - start_positions[m] - 1));
                        for (int j = 0; j < patch.length(); j++)
                            outputQual.push_back(computeQuality(max_count, min_count, readsHash[encode_kmer_minimum(curr_kmer)]));
                        i = end_positions[selected_dest_per_start[m]] + 1;
                        break;
                    }
                    else if (path_found_in_read_graph)
                    {
                        patch = graph_patch;
                        outputQual.resize(outputGenome.size() - (i - start_positions[selected_start_idx] - 1));
                        outputGenome.resize(outputGenome.size() - (i - start_positions[selected_start_idx] - 1));
                        for (int j = 0; j < patch.length(); j++)
                            outputQual.push_back(computeQuality(max_count, min_count, readsHash[encode_kmer_minimum(curr_kmer)]));
                        i = end_positions[selected_anchor_idx[selected_start_idx]] + 1;
                        break;
                    }
                    else if (m == prev_found_kmers.size() - 1)
                    {
                        /*
                         *  If there is no way to support the genome we take part of the assembly from the closest reference and
                         *  we also mark this path as untrustworthy, since there was no trace of it in our reads,
                         *  the only thing we keep as trustworthy is the part supported by the anchor kmer
                         *  we use the first anchor kmer to reduce the lost info
                         */
                        /* We use the last first found destination kmer and the last found prev kmer to patch up when copying only */
                        patch = referenceGenome.substr(i + kmer_len - 1, end_positions[0] - i + 1);
                        if (end_positions[0] - i + 1 > kmer_len)
                            std::transform(patch.begin(), patch.end() - (kmer_len), patch.begin(), ::tolower);
                        for (int j = 0; j < patch.length(); j++)
                            outputQual.push_back(EMULATED_QUALITY);
                        prev_found_kmers.push_front(anchor_kmers[0]);
                        start_positions.push_front(end_positions[0]);
                        if (prev_found_kmers.size() > MAX_START)
                        {
                            prev_found_kmers.pop_back();
                            start_positions.pop_back();
                        }
                        i = end_positions[0] + 1;
                    }
                }
            }
            else
            {
                /*
                 * If no anchor kmer is found we copy the reference genome,
                 * again marking this as untrustworthy
                 */
                patch = referenceGenome.substr(i + kmer_len - 1, referenceGenome.length() - kmer_len - (i + kmer_len - 1));
                std::transform(patch.begin(), patch.end(), patch.begin(), ::tolower);
                i += (referenceGenome.length() - kmer_len - (kmer_len - 1));
                for (int j = 0; j < patch.length(); j++)
                    outputQual.push_back(EMULATED_QUALITY);
            }
            outputGenome.append(patch);
        }
    }
    /*
     * We copy the last kmer of the reference too
     */
    outputGenome.append(referenceGenome.substr(referenceGenome.length() - kmer_len, kmer_len));
    for (int i = 0; i < kmer_len; i++)
        outputQual.push_back(EMULATED_QUALITY);
    /*
     * Save the output assembly
     */
    assert(outputGenome.size() == outputQual.size());
    std::ofstream outFile;
    outFile.open(file_name);
    outFile << ">ASSEMBLY\n";
    outFile << outputGenome << "\n";
    if (isfastq)
    {
        outFile << "+\n";
        outFile << outputQual << "\n";
    }
    outFile.close();
}