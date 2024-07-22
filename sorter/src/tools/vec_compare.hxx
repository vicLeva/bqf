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

inline bool vec_compare(std::vector<uint64_t>& dst, std::vector<uint64_t>& src)
{
    if(dst.size() != src.size())
        return false;

    const int length = src.size();
    for(int x = 0; x < length; x += 1)
        if( dst[x] != src[x] )
            return false;
    return true;
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
