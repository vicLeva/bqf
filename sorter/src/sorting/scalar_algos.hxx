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

template<class T>
void Bubble_Sort(std::vector<T> arr)
{
    const int n = arr.size();
    for(int i = 0;i < n; i++)
    {
        for(int j=0;j < n - 1; j++)
        {
            if(arr[j] > arr[j+1])
            {
                int temp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = temp;
            }
        }
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

template<class T>
void Revised_Bubble_Sort(std::vector<T> arr)
{
    const int n = arr.size();
    for(int i = 0;i < n; i++){
        bool isSwapped = false;
        for(int j=0;j < n - 1; j++){
            if(arr[j] > arr[j+1]){
                int temp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = temp;
                isSwapped = true;
            }
            if(!isSwapped){
                break;
            }
        }
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

template<class T>
void Selection_Sort(std::vector<T> arr)
{
    const int n = arr.size();
    for(int i = 0; i < n; i++) {
        int minIndex = i;
        for(int j = i+1; j < n; j++)
        {
            if( arr[j] < arr[minIndex] )
            {
                minIndex = j;
            }
        }
        int temp = arr[i];
        arr[i] = arr[minIndex];
        arr[minIndex] = temp;
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

template<class T>
void Insertion_Sort(std::vector<T> arr)
{
    const int n = arr.size();
    for(int i = 1; i < n; i++) {
        int j = i;
        while(j > 0 && arr[j] < arr[j-1]) {
            int temp = arr[j];
            arr[j]   = arr[j-1];
            arr[j-1] = temp;
            j--;
        }
    }
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

int partition(int* arr, const int low, const int high)
{
    int pivot = arr[high];
    int left = low, right = high-1;
    while(left < right) {
        while(arr[left]<pivot) {
            left++;
        }
        while(arr[right]>pivot) {
            right--;
        }
        if(left >= right) {
            break;
        }
        int temp = arr[left];
        arr[left] = arr[right];
        arr[right] = temp;
    }
    int temp = arr[left];
    arr[left] = arr[high];
    arr[high] = temp;
    return left;
}

 void quicksort(int* arr, int low, int high)
 {
    if(low >= high) return;
    int pivotPosition = partition(arr, low, high);
    quicksort(arr,low, pivotPosition-1);
    quicksort(arr, pivotPosition+1, high);
}

//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
#if 0
void merge(int* arr, int l, int m, int r) {
    int n1 = m-l+1;
    int n2 = r-m;
    int[] L = new int[n1];
    int[] R = new int[n2];
    for(int i = 0;i < n1; i++) {
        L[i] = arr[l+i];
    }
    for(int i = 0;i < n2; i++) {
        R[i] = arr[m+1+i];
    }
    int i = 0, j = 0, k =l;
    while(i < n1 && j < n2) {
        if(L[i] <= R[j]) {
            arr[k++] = L[i++];
        }
        else {
            arr[k++] = R[j++];
        }
    }
    while(i < n1) {
        arr[k++] = L[i++];
    }
    while(j < n2) {
        arr[k++] = R[j++];
    }
}

void mergeSort(int* arr, int l, int r)
{
    if (l < r) {
        int m = l + (r-l)/2;
        sort(arr, l, m);
        sort(arr , m+1, r);
        merge(arr, l, m, r);
    }
}
#endif