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
void merge_sort(Iter first, Iter last)
{
    if ((last - first) > 1)
    {
        Iter middle = first + (last - first) / 2;
        merge_sort(first, middle);
        merge_sort(middle, last);
        std::inplace_merge(first, middle, last);
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
