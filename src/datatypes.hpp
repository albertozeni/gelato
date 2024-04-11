#ifndef DATATYPES_HPP_
#define DATATYPES_HPP_

#include "robin_hood.h"
#include <bits/stdc++.h>

// Hash function
struct hash
{
    std::uint64_t operator()(const std::vector<std::uint64_t>
                                 &v) const
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
};

typedef robin_hood::unordered_flat_map<std::uint64_t, float> hash_k;
typedef robin_hood::unordered_flat_map<std::uint64_t, std::uint8_t> graph_k;
typedef robin_hood::unordered_flat_set<std::uint64_t> set_k;
typedef robin_hood::unordered_flat_set<std::vector<std::uint64_t>, hash> set_v;

#define N_BASES 4
#define MAX_DEST 16
#define MAX_START 8

#endif // DATATYPES_HPP_