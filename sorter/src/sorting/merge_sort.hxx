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


template<class Iter>
void v2_merge_sort(Iter first, Iter last)
{
    if (last - first > 1)
    {
        Iter middle = first + (last - first) / 2;
        v2_merge_sort(first, middle);
        v2_merge_sort(middle, last);
        std::inplace_merge(first, middle, last);
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

template<class T>
void v2_merge_sort(std::vector<T>& v)
{
    v2_merge_sort(v.begin(), v.end());
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

void v2_mergeSortSerial(int aux[], int left, int right){
    if (left < right){
        int middle = (left + right)/2;
        v2_mergeSortSerial(aux,left,middle); //call 1
        v2_mergeSortSerial(aux,middle+1,right); //call 2
        merge(aux, left, middle, right);
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

void v2_mergeSort (int aux[], int left, int right){
    if (left < right){
        if ((right-left) > 1000)
        {
            int middle = (left + right)/2;
#pragma omp task firstprivate (aux, left, middle)
            v2_mergeSort(aux,left,middle); //call 1
#pragma omp task firstprivate (aux, middle, right)
            v2_mergeSort(aux,middle+1,right); //call 2
#pragma omp taskwait
            std::merge(aux,left,middle,right);
        } else{
            v2_mergeSortSerial(aux, left, right);
        }
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

template<class Iter>
void v4_merge_sort(Iter first, Iter last)
{
    if (last - first > 1)
    {
        Iter middle = first + (last - first) / 2;
        v4_merge_sort(first, middle);
        v4_merge_sort(middle, last);
        std::inplace_merge(first, middle, last);
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
