#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <omp.h>

#include "sorting/std_2cores.hpp"
#include "sorting/std_4cores.hpp"
#include "tools/tools.hpp"

#include "../../src/headers/additional_methods.hpp"


std::chrono::steady_clock::time_point begin;

int main(int argc, char* argv[])
{
    int64_t minv = 256;
    int64_t maxv = 1073741824;

    if( argc >= 2 )
        minv = atoi( argv[1] );

    if( argc >= 3 )
        maxv = atoi( argv[1] );


    for (int64_t v = minv; v <= maxv; v *= 2)
    {
        std::vector<uint64_t> list_std(v);
        std::vector<uint64_t> list_custom(v);

        for (int64_t i = 0; i < v; i += 1)
        {
            list_std[i]    = rand() % 999;
            list_custom[i] = list_std[i];
        }

        if (v == 256)
        {
            printf("(II) Unsorted values :");
            for (int i = 0; i < 256; i += 1) {
                if ((i % 32) == 0) printf("\n%3d | ", i);
                printf("%3llu ", list_std[i]);
            }
            printf("\n\n");
        }

        double start_time_std = omp_get_wtime();
        std::sort(list_std.begin(), list_std.end());
        double end_time_std = omp_get_wtime();

        printf("%10lld : It took %g seconds\n", v, end_time_std - start_time_std);

        double start_time_cstm = omp_get_wtime();
        p_sort_4x(list_custom);
        double end_time_cstm = omp_get_wtime();

        printf("       :         %g seconds\n",    end_time_cstm - start_time_cstm);

        //boost::sort::spreadsort::integer_sort(v.begin(), v.end());
        if (v == 256)
        {
            printf("(II) Sorted values :");
            for (int i = 0; i < 256; i += 1) {
                if ((i % 32) == 0) printf("\n%3d | ", i);
                printf("%3llu ", list_std[i]);
            }
            printf("\n");

            printf("(II) Sorted values :");
            for (int i = 0; i < 256; i += 1) {
                if ((i % 32) == 0) printf("\n%3d | ", i);
                printf("%3llu ", list_custom[i]);
            }
            printf("\n");


        }

        if( vec_compare(list_std, list_custom) == false )
        {
            printf("(EE) Error, the comparison failed\n");
            exit( EXIT_FAILURE );
        }

    }

    return 0;
}