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

inline void vec_print(std::vector<uint64_t>& src)
{
    printf("(II) Printed values :");
    for (int i = 0; i < src.size(); i += 1) {
        if ((i % 32) == 0) printf("\n%3d | ", i);
        printf("%3llu ", src[i]);
    }
    printf("\n");
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

inline void vec_print(uint64_t* src, int length)
{
    printf("(II) Printed values :");
    for (int i = 0; i < length; i += 1)
    {
        if ((i % 32) == 0) printf("\n%3d | ", i);
        printf("%3llu ", src[i]);
    }
    printf("\n");
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
