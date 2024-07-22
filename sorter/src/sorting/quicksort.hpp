#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <thread>

#define TASK_SIZE 256

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

int partition(int * a, int p, int r)
{
    int lt[r-p];
    int gt[r-p];
    int i;
    int j;
    int key = a[r];
    int lt_n = 0;
    int gt_n = 0;

    for(i = p; i < r; i++){
        if(a[i] < a[r]){
            lt[lt_n++] = a[i];
        }else{
            gt[gt_n++] = a[i];
        }
    }

    for(i = 0; i < lt_n; i++){
        a[p + i] = lt[i];
    }

    a[p + lt_n] = key;

    for(j = 0; j < gt_n; j++){
        a[p + lt_n + j + 1] = gt[j];
    }

    return p + lt_n;
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

void quicksort(int * a, int p, int r)
{
    int div;
    if(p < r)
    {
        div = partition(a, p, r);
#pragma omp task shared(a) if(r - p > TASK_SIZE)
        quicksort(a, p, div - 1);
#pragma omp task shared(a) if(r - p > TASK_SIZE)
        quicksort(a, div + 1, r);
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

template<class T>
void QuickSort_omp(std::vector<T> vec)
{
    const auto numThreads = std::thread::hardware_concurrency();
    omp_set_dynamic(0);              /** Explicitly disable dynamic teams **/
    omp_set_num_threads(numThreads); /** Use N threads for all parallel regions **/
#pragma omp parallel
    {
#pragma omp single
        quicksort(vec, 0, vec.size());
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

