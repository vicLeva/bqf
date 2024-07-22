#pragma once
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

template<class T = uint64_t>
inline void vec_random(T* dst, const int length, const int rand_max = 65536)
{
    for(int x = 0; x < length; x += 1)
        dst[x] = rand() % rand_max;
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

template<class T = uint64_t>
inline void vec_random(std::vector<T>& dst, const int rand_max = 65536)
{
    vec_random(dst.data(), dst.size(), rand_max);
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
