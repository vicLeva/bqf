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
inline void vec_copy(T* dst, const T* src, const int length)
{
    for(int x = 0; x < length; x += 1)
        dst[x] = src[x];
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

template<class T = uint64_t>
inline void vec_copy(std::vector<T>& dst, uint64_t* src, const int length)
{
    if(dst.size() != length)
        dst.resize( length );

    vec_copy(dst.data(), src, length);
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

template<class T = uint64_t>
inline void vec_copy(std::vector<T>& dst, const std::vector<T>& src)
{
    if(dst.size() != src.size())
        dst.resize( src.size() );

    vec_copy(dst.data(), src.data(), dst.size());
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

