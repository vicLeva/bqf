#pragma once
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <thread>

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

template <class T = uint64_t>
void seq_qsort(int p, int r, T* data) {
    if (p < r) {
        int q = partition<T>(p, r, data);
        seq_qsort(p, q - 1, data);
        seq_qsort(q + 1, r, data);
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

template <class T = uint64_t>
void q_sort_tasks(int p, int r, T* data, int low_limit) {
    if (p < r) {
        if (r - p < low_limit)
        {
            // small list => sequential.
            return seq_qsort<T>(p, r, data);
        }else{

            int q = partition<T>(p, r, data);
            // create two tasks
#pragma omp task shared(data)
            q_sort_tasks<T>(p, q - 1, data, low_limit);
#pragma omp task shared(data)
            q_sort_tasks<T>(q + 1, r, data, low_limit);
        }
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

template <class T = uint64_t>
void par_q_sort_tasks(int p, int r, T* data){
    const auto processor_count = std::thread::hardware_concurrency();
#pragma omp parallel
    {
#pragma omp single
        q_sort_tasks<T>(p, r, data, processor_count - 1);
#pragma omp taskwait
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
