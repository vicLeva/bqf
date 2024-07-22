#pragma once
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <omp.h>

#include "../tools/vec_copy.hxx"
#include "../tools/vec_print.hxx"

#include <thread>
#include <execution>

static void local_sort(
        std::vector<uint64_t>* dst,
        std::vector<uint64_t>* src,
        uint64_t start, int length)
{
    vec_copy(*dst, src->data() + start, length );
    std::sort( dst->begin(), dst->end() );
}

void std_2cores(std::vector<uint64_t>& test)
{
    const int size = test.size();
    const int half = size / 2;

    std::vector<uint64_t> v1( half );
    std::vector<uint64_t> v2( half );

    std::thread t1 (&local_sort, &v1, &test, 0, half);
    std::thread t2 (&local_sort, &v2, &test, half, half);

    t1.join();
    t2.join();

    int ptr_1 = 0;
    int ptr_2 = 0;

    int i;
    for(i = 0; i < size; i += 1)
    {
        if( v1[ptr_1] < v2[ptr_2] ) {
            test[i] = v1[ptr_1++];
        } else {
            test[i] = v2[ptr_2++];
        }
        if( (ptr_1 == half) || (ptr_2 == half)  )
            break;
    }
    i += 1;

    if( ptr_1 == half )
    {
        for(; i < size; i += 1)
            test[i] = v2[ptr_2++];
    }else{
        for(; i < size; i += 1)
            test[i] = v1[ptr_1++];
    }

    v1.clear();
    v2.clear();
}
