#pragma once
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <thread>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

void mergeSortSerial(int aux[], int left, int right){
    if (left < right){
        int middle = (left + right)/2;
        mergeSortSerial(aux,left,middle); //call 1
        mergeSortSerial(aux,middle+1,right); //call 2
        std::merge(aux,left,middle,right);
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

void mergeSort (int aux[], int left, int right){
    if (left < right){
        if ((right-left) > 1000){
            int middle = (left + right)/2;
#pragma omp task firstprivate (aux, left, middle)
            mergeSort(aux,left,middle); //call 1
#pragma omp task firstprivate (aux, middle, right)
            mergeSort(aux,middle+1,right); //call 2
#pragma omp taskwait
            merge(aux,left,middle,right);
        } else{mergeSortSerial(aux, left, right);}
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
