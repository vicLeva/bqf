#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <omp.h>
#include <sstream>

#include "../../src/headers/additional_methods.hpp"

#include "progress/progressbar.h"

#include "./sorting/std_2cores.hxx"
#include "./sorting/std_4cores.hxx"
#include "./sorting/crumsort_2cores.hxx"


void vec_copy(std::vector<uint64_t>& dst, std::vector<uint64_t>& src)
{
    if(dst.size() != src.size())
        dst.resize( src.size() );

    const int length = src.size();
    for(int x = 0; x < length; x += 1)
        dst[x] = src[x];
}

void vec_copy(std::vector<uint64_t>& dst, uint64_t* src, const int length)
{
    if((int)dst.size() != length)
        dst.resize( length );

    for(int x = 0; x < length; x += 1)
        dst[x] = src[x];
}

void my_sort(std::vector<uint64_t>& test)
{
    const int size = test.size();
    const int half = size / 2;

    std::vector<uint64_t> v1( half );
    std::vector<uint64_t> v2( half );

    vec_copy(v1, test.data(),        half );
    vec_copy(v2, test.data() + half, half );

    std::sort( v1.begin(), v1.end() );
    std::sort( v2.begin(), v2.end() );

    int ptr_1 = 0;
    int ptr_2 = 0;

    int i;
    for(i = 0; i < size; i += 1)
    {
        if( v1[ptr_1] < v2[ptr_2] )
        {
            test[i] = v1[ptr_1++];
        } else {
            test[i] = v2[ptr_2++];
        }
    }
    if( ptr_1 == half )
    {
        for(; i < size; i += 1)
            test[i] = v2[ptr_2++];
    }else{
        for(; i < size; i += 1)
            test[i] = v2[ptr_1++];
    }
    v1.clear();
    v2.clear();
}



//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

std::chrono::steady_clock::time_point begin;

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

int count_lines(std::string filename)
{
    printf("(II) Counting the number of sequences (original)\n");
    std::ifstream ifile( filename );

    if( ifile.is_open() == false )
    {
        printf("(EE) File does not exist !\n");
        return EXIT_FAILURE;
    }
    double start_time = omp_get_wtime();
    int n_sequences = 0;
    std::vector<uint64_t> list_hash;
    std::string line, seq;
    while( ifile.eof() == false ) {
        getline(ifile, line);
        n_sequences += 1;
    }
    double end_time = omp_get_wtime();
    ifile.close();
    printf("(II) - %d lines in file\n", n_sequences);
    printf("(II) - It took %g seconds\n", end_time - start_time);
    return n_sequences;
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

int read_k_value(std::string filename)
{
    std::ifstream ifile( filename );
    if( ifile.is_open() == false )
    {
        printf("(EE) File does not exist !\n");
        return EXIT_FAILURE;
    }
    std::string line, seq;

    getline(ifile, line);

    ifile.close();

    std::stringstream ss(line);
    getline(ss, seq, '\t');

    int k_length = seq.length();
    printf("(II) Sequence line = %s\n", seq.c_str());
    printf("(II) k_length      = %d\n", k_length);

    return k_length;
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

int count_lines_cpp(std::string filename)
{
    printf("(II) Counting the number of sequences (cpp)\n");
    std::ifstream ifile( filename );
    if( ifile.is_open() == false )
    {
        printf("(EE) File does not exist !\n");
        return EXIT_FAILURE;
    }

    int n_sequences = 0;
    std::vector<char> buffer;
    double start_time = omp_get_wtime();
    buffer.resize(1024 * 1024); // buffer of 1MB size
    while (ifile)
    {
        ifile.read( buffer.data(), buffer.size() );
        n_sequences += std::count( buffer.begin(),
                                 buffer.begin() + ifile.gcount(), '\n' );
    }
    ifile.close();
    double end_time = omp_get_wtime();
    printf("(II) - %d lines in file\n", n_sequences);
    printf("(II) - It took %g seconds\n", end_time - start_time);
    return n_sequences;
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
int fast_atoi(const char *a)
{
    int res = (*a++) - '0';
    while( *a != '\0' )
    {
        res = res * 10 + (*a - 48);
        a  += 1;
    }
    return res;
}

class fast_kmer_file
{
private:
    const int kmer_size;
    const int buff_size = 64 * 1024 * 1024;
    char* buffer;
    int n_data   = 0;
    int c_ptr    = 0;
    int n_lines        = 0;
    bool no_more_load  = false;
    bool file_ended    = false;
    FILE* f;

public:
    fast_kmer_file(const std::string filen, const int km_size) : kmer_size(km_size)
    {
        printf("(II) Creating fast_kmer_file object\n");

        buffer = new char[buff_size];

        f = fopen( filen.c_str(), "r" );
        if( f == NULL )
        {
            printf("(EE) File does not exist !\n");
            exit( EXIT_FAILURE );
        }
        n_data = fread(buffer, sizeof(char), buff_size, f);
        printf("(II) Creation done\n");
    }

    //////////////////////////////////////////////////////////////////////////////////////////

    ~fast_kmer_file()
    {
        delete[] buffer;
        fclose( f );
    }

    //////////////////////////////////////////////////////////////////////////////////////////

    bool next_kmer(char* n_kmer, char* n_value)
    {
        if( file_ended == true )
        {
            return false;
        }

        //
        // On regarde si l'on a encore une entrée dans le buffer
        //
        int pos_nline = -1;
        for(int i = c_ptr; i < n_data; i += 1)
        {
            if(buffer[i] == '\n')
            {
                pos_nline = i;
                break;
            }
        }

        if( pos_nline == -1 )
        {
            if( is_eof() == true )
                return false;

            if( reload() == true )
            {
                return next_kmer(n_kmer, n_value);
            }else{
                return false;
            }
        }

        //
        // On peut recopier le k-mer dans la structure
        //

        for(int i = 0; i < kmer_size; i += 1)
        {
            n_kmer[i] = buffer[c_ptr + i];
        }
        c_ptr += kmer_size;
        n_kmer[kmer_size] = 0;

        //
        // On verifie que l'on a bien un espace / tabulation
        //

        c_ptr += 1;

        //
        // On recopie la valeur numérique de l'abondance
        //

        int cnt = 0;
        while( buffer[c_ptr] != '\n' )
        {
            n_value[cnt++] = buffer[c_ptr++];
        }
        n_value[cnt] = '\0';
        c_ptr += 1;

        n_lines += 1;

        file_ended = (c_ptr == n_data) && (n_data != buff_size);

        return true;
    }

    //////////////////////////////////////////////////////////////////////////////////////////

    bool reload()
    {
//        printf("[%6d] START OF DATA RELOADING !\n", n_lines);
        int reste = n_data - c_ptr;
        for(int i = c_ptr; i < n_data; i += 1)
        {
            buffer[i - c_ptr] = buffer[i];
        }
        int nread = fread(buffer + reste, sizeof(char), buff_size - reste, f);
//        printf("[%6d] ASKED        : %d\n", n_lines, buff_size - reste);
//        printf("[%6d] READ         : %d\n", n_lines, nread);
        no_more_load = ( n_data != buff_size ); // a t'on atteint la fin du fichier ?
//        printf("[%6d] no_more_load : %d\n", n_lines, no_more_load);
        c_ptr        = 0;                       // on remet a zero le pointeur de lecture
        n_data       = nread + reste;           // on met a jour le nombre de données dans le buffer
//        printf("[%6d] n_data       : %d\n", n_lines, n_data);
//        printf("[%6d] END OF DATA RELOADING !\n", n_lines);
        return true;
    }

    //////////////////////////////////////////////////////////////////////////////////////////

    bool is_eof()
    {
        if( (no_more_load == true) && (c_ptr == n_data) )
            return true;
        return false;
    }

    //////////////////////////////////////////////////////////////////////////////////////////
};

int count_lines_c(std::string filename)
{
    printf("(II) Counting the number of sequences (c)\n");
    FILE* f = fopen( filename.c_str(), "r" );
    if( f == NULL )
    {
        printf("(EE) File does not exist !\n");
        return EXIT_FAILURE;
    }

    int n_sequences = 0;
    char buffer[1024 * 1024];
    double start_time = omp_get_wtime();
    while ( !feof(f) )
    {
        int n = fread(buffer, sizeof(char), 1024 * 1024, f);
        for(int x = 0; x < n; x += 1)
            n_sequences += (buffer[x] == '\n');
    }
    fclose( f );
    double end_time = omp_get_wtime();
    printf("(II) - %d lines in file\n", n_sequences);
    printf("(II) - It took %g seconds\n", end_time - start_time);
    return n_sequences;
}

bool vec_compare(std::vector<uint64_t> dst, std::vector<uint64_t> src)
{
    if(dst.size() != src.size())
        return false;

    const int length = src.size();
    for(int x = 0; x < length; x += 1)
        if( dst[x] != src[x] )
            return false;
    return true;
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

void fast_decode(uint64_t revhash, int k, char* ptr)
{
    char lut[4] = {'A', 'C', 'T', 'G'};
    for(int i = 0; i < k; i += 1)
    {
        ptr[ k-i-1 ] = lut[ revhash & 0x03 ];
        revhash >>= 2;
    }
}


int int_to_char(int value, char* ptr)
{
    if( value >= 100 ){
        int unit = value%10;
        int diza = (value/10)%10;
        int cent = (value/100);
        ptr[0] = '0' + cent;
        ptr[1] = '0' + diza;
        ptr[2] = '0' + unit;
        return 3;

    }else if( value >= 10){
        int unit = value%10;
        int diza = (value/10)%10;
        ptr[0] = '0' + diza;
        ptr[1] = '0' + unit;
        return 2;

    }else{
        ptr[0] = '0' + value;
        return 1;
    }
}

bool SaveToFile(const std::string filename, const std::vector<uint64_t> list_hash, const int smer_size, int q, int r)
{
    printf("(II)\n");
    printf("(II) Saving data set in [%s]\n", filename.c_str());
    progressbar *progress = progressbar_new("Loading k-mers",100);
    double start_time = omp_get_wtime();

    std::string n_file = filename;
    FILE* f = fopen( n_file.c_str(), "w" );
    if( f == NULL )
    {
        printf("(EE) File does not exist !\n");
        exit( EXIT_FAILURE );
    }

    const int buffer_size = 64 * 1024 * 1024;
    char* buffer = new char[ buffer_size ];

    int  n_data = 0;
    const int n_lines = list_hash.size();
    const int prog_step = n_lines / 100;
    for(int y = 0; y < n_lines; y += 1)
    {
        const uint64_t element = list_hash[ y ];
        const uint64_t s_hash  = (element   >> 8);
        const int      r_abon  = (element & 0xFF);

        const uint64_t mask   = {(~(0ULL)) >> (64 - q - r)};
        const uint64_t r_hash = (s_hash >> r | s_hash << q) & mask;

//      if( y < 4 )  printf( "(2) %16.16lX\n", s_hash);
//      if( y < 4 )  printf( "(1) %16.16lX\n", r_hash);

        //std::string toto = hash_to_kmer(r_hash, smer_size);
        uint64_t smer = bfc_hash_64_inv(r_hash, mask_right(smer_size * 2));

        fast_decode(smer, smer_size, buffer + n_data);
        n_data += smer_size;

        buffer[n_data] = '\t';
        n_data += 1;

        n_data += int_to_char(r_abon, buffer + n_data);

        buffer[n_data] = '\n';
        n_data += 1;

        if( n_data >= (buffer_size - 256))
        {
            fwrite(buffer, sizeof(char), n_data, f);
            n_data = 0;
        }

        if( y%prog_step == 0)
            progressbar_inc(progress);
    }
    // on flush le dernier lot de data aant de quitter
    fwrite(buffer, sizeof(char), n_data, f);

    fclose( f );
    double end_time = omp_get_wtime();

    progressbar_finish(progress);

    printf("(II) - Execution time    = %f\n", end_time - start_time);
    return true;
}

bool SaveToFileRAW(const std::string filename, const std::vector<uint64_t> list_hash, const int smer_size, int q, int r)
{
    printf("(II)\n");
    printf("(II) Saving data set in [%s]\n", filename.c_str());
    progressbar *progress = progressbar_new("Saving RAW k-mer data",100);
    double start_time = omp_get_wtime();

    std::string n_file = filename;
    FILE* f = fopen( n_file.c_str(), "w" );
    if( f == NULL )
    {
        printf("(EE) File does not exist !\n");
        exit( EXIT_FAILURE );
    }

    char buff_smer[128];

    const int n_lines = list_hash.size();
    const int prog_step = n_lines / 100;
    for(int y = 0; y < n_lines; y += 1)
    {
        const uint64_t element = list_hash[ y ];
              uint64_t s_hash  = (element   >> 8);
        const int      r_abon  = (element & 0xFF);

        const uint64_t mask   = {(~(0ULL)) >> (64 - q - r)};
        const uint64_t r_hash = (s_hash >> r | s_hash << q) & mask;

        //std::string toto = hash_to_kmer(r_hash, smer_size);
        uint64_t smer = bfc_hash_64_inv(r_hash, mask_right(smer_size * 2));

        fast_decode(smer, smer_size, buff_smer);
        buff_smer[smer_size] = '\0';

        fprintf(f, "%16.16llX | %s | %d\n", list_hash[y], buff_smer, r_abon);

        if( y%prog_step == 0)
            progressbar_inc(progress);
    }
    // on flush le dernier lot de data aant de quitter
    fclose( f );
    double end_time = omp_get_wtime();

    progressbar_finish(progress);

    printf("(II) - Execution time    = %f\n", end_time - start_time);
    return true;
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

int main(int argc, char* argv[])
{

    std::string filename = "../AHX_ACXIOSF_6_1_23_all.txt";

    if(argc == 2)
        filename = argv[1];

    int n_lines = count_lines_c( filename );

    const int k = read_k_value(filename);
    const int s = k;
    const int z = k - s;         //

    const int q = 30;
    const int r = 2 * s - q;

    std::vector<uint64_t> list_hash( n_lines );
    const int MB = (list_hash.size() / 1024 / 1024 * sizeof(uint64_t));

    printf("(II)\n");
    printf("(II)\n");
    printf("(II) # of sequences : %d\n", n_lines);
    printf("(II) k-mer length   : %d\n", k);
    printf("(II) s-mer length   : %d\n", s);
    printf("(II) # s-mer/k-mer  : %d\n", z);
    printf("(II) q value        : %d\n", q);
    printf("(II) r value        : %d\n", r);
    printf("(II)\n");
    printf("(II) memory occupancy = %6d MB (%2.3f%%)\n", MB, 100.f * (float)n_lines / (float)n_lines );
    printf("(II)\n");

    char kmer_value[256]; bzero(kmer_value, 256);
    char kmer_abond[256]; bzero(kmer_abond, 256);

    fast_kmer_file fast_ifile(filename, k);

    const uint64_t s_mask = mask_right(2 * s);

    progressbar *progress = progressbar_new("Loading k-mers",100);
    const int prog_step = n_lines / 100;
    for(int l_number = 0; l_number < n_lines; l_number += 1)
    {
        bool not_oef = fast_ifile.next_kmer(kmer_value, kmer_abond);
        if( not_oef == false )
        {
            break;
        }

//      printf("encoded > kmer_value = %s\n", kmer_value);
        uint64_t current_smer = 0;
        for (auto i = 0; i < k; i++)
        {
            current_smer <<= 2;
            current_smer |= nucl_encode(kmer_value[i]);
        }

        const uint64_t canon  = canonical(current_smer, 2 * s);
//      const uint64_t vflip  = flip(canon, 2*s);
              uint64_t s_hash = bfc_hash_64(canon, s_mask);

//        if( l_number < 4 )  printf("(canon ) %16.16lX\n", canon);
//        if( l_number < 4 )  printf("(s_hash) %16.16lX\n", s_hash);

//        if( l_number < 4 )  printf("(1) %16.16lX\n", s_hash);
        const uint64_t mask {(~(0ULL)) >> (64 - q - r)};
        s_hash = (s_hash << r | s_hash >> q) & mask;
//        if( l_number < 4 )  printf( "(2) %16.16lX\n", s_hash);

        const uint64_t abondan   = fast_atoi( kmer_abond );
        const uint64_t element   = (s_hash << 8 ) | abondan;
        list_hash[ l_number ] = element;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        //
#if 0
        const uint64_t v_hash  = (element   >> 8);
        const int      v_abon  = (element & 0xFF);
        const uint64_t q_hash  = (v_hash >> r | v_hash << q) & mask;

        if( l_number < 4 )  printf( "(2) %16.16lX\n", v_hash);
        if( l_number < 4 )  printf( "(1) %16.16lX\n", q_hash);

        const uint64_t x_smer = bfc_hash_64_inv(q_hash, mask_right(s * 2));

        if( l_number < 4 )  printf( "(0) %16.16lX\n", x_smer);

        fast_decode(x_smer, k, kmer_value);  // inv. nucleotide

        printf("kmer_value = %s\n", kmer_value);  //
        exit( 0 );
#endif
        //
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if( l_number%prog_step == 0)
            progressbar_inc(progress);
    }

    progressbar_finish(progress);

    printf("(II) Number of ADN sequences : %d\n",  n_lines);
    printf("(II) Number of s-mer         : %zu\n", list_hash.size());

//    SaveToFile(filename + ".loaded", list_hash, s, q, r);

    printf("(II)\n");
    printf("(II) Launching the sorting step\n");
    printf("(II) - Number of samples = %ld\n", list_hash.size());
    double start_time = omp_get_wtime();

//    std::sort( list_hash.begin(), list_hash.end() );
//    std_2cores( list_hash );
//    std_4cores( list_hash );
//    crumsort_2cores( list_hash );
    crumsort_prim( list_hash.data(), list_hash.size(), 9 /*uint64*/ );

    double end_time = omp_get_wtime();
    printf("(II) - Execution time    = %f\n", end_time - start_time);

    SaveToFile(filename + ".ordered", list_hash, s, q, r);

//    SaveToFileRAW(filename + ".hash", list_hash, s, q, r);

    return 0;
}