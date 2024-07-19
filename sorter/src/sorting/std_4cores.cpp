#include "std_2cores.hpp"
#include <thread>

static void local_sort_thread_4x(
        uint64_t* dst,
        uint64_t* src,
        uint64_t start, int length)
{
    vec_copy (dst + start,  src + start, length );
    std::sort(dst + start, dst + start + length );
}

void local_merge(uint64_t* dst, uint64_t* src_1, uint64_t* src_2, int size)
{
    int half  = size / 2;
    int ptr_1 = 0;
    int ptr_2 = 0;

    int i;
    for(i = 0; i < size; i += 1)
    {
        if( src_1[ptr_1] < src_2[ptr_2] ) {
            dst[i] = src_1[ptr_1++];
        } else {
            dst[i] = src_2[ptr_2++];
        }
        if( (ptr_1 == half) || (ptr_2 == half)  )
            break;
    }
    i += 1;

    if( ptr_1 == half )
    {
        for(; i < size; i += 1)
            dst[i] = src_2[ptr_2++];
    }else{
        for(; i < size; i += 1)
            dst[i] = src_1[ptr_1++];
    }
}

void p_sort_4x(std::vector<uint64_t>& test)
{
    double A = omp_get_wtime();
    const int size = test.size();
    const int half = size / 2;
    const int quar = half / 2;
    const int qua3 = half + quar;

    uint64_t* tmp = new uint64_t[size];

    std::thread t1 (&local_sort_thread_4x, tmp, test.data(), 0,    quar);
    std::thread t2 (&local_sort_thread_4x, tmp, test.data(), quar, quar);
    std::thread t3 (&local_sort_thread_4x, tmp, test.data(), half, quar);
    std::thread t4 (&local_sort_thread_4x, tmp, test.data(), qua3, quar);

    t1.join();
    t2.join();
    t3.join();
    t4.join();
    double B = omp_get_wtime();

    std::thread z1(local_merge, test.data(),        tmp,        tmp + quar, half);
    std::thread z2(local_merge, test.data() + half, tmp + half, tmp + qua3, half);
    z1.join();
    z2.join();

    //
    double C = omp_get_wtime();
    local_merge(tmp,                test.data(),test.data() + half, size);

    double D = omp_get_wtime();

    for(int x = 0; x < size; x += 1)
        test[x] = tmp[x];

    double E = omp_get_wtime();

    printf(" (1) %g seconds\n",    B - A);
    printf(" (2) %g seconds\n",    C - B);
    printf(" (3) %g seconds\n",    D - C);
    printf(" (4) %g seconds\n",    E - D);

    delete[] tmp;
}
