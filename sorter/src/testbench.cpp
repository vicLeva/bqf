#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>

#include "./tools/vec_compare.hxx"
#include "./tools/vec_copy.hxx"
#include "./tools/vec_print.hxx"
#include "./tools/vec_random.hxx"

#include "./sorting/scalar_algos.hxx"
#include "./sorting/std_2cores.hxx"
#include "./sorting/std_4cores.hxx"
#include "./sorting/crumsort_2cores.hxx"

//
//
///////////////////////////////////////////////////////////////////////////////////////////
//
//

TEST_CASE( "Testing ", "[std::2cores]" )
{
    const int v_min = 256;
    const int v_max = 262144;

    for(int i = v_min; i < v_max; i *= 2)
    {
        std::vector<uint64_t> v1( i );
        std::vector<uint64_t> v2( i );
        vec_random(v1, i);
        vec_copy(v2, v1);
        std::sort(v1.begin(), v1.end());
        std_2cores(v2);
        const bool is_ok = vec_compare(v1, v2);
        if( is_ok == false )
        {
            vec_print(v1);
            vec_print(v2);
        }
        REQUIRE( is_ok );
    }
}

//
//
///////////////////////////////////////////////////////////////////////////////////////////
//
//

TEST_CASE( "Testing ", "[std::4cores]" )
{
    const int v_min = 256;
    const int v_max = 262144;

    for(int i = v_min; i < v_max; i *= 2)
    {
        std::vector<uint64_t> v1( i );
        std::vector<uint64_t> v2( i );
        vec_random(v1, i);
        vec_copy(v2, v1);
        std::sort(v1.begin(), v1.end());
        std_4cores(v2);
        const bool is_ok = vec_compare(v1, v2);
        if( is_ok == false )
        {
            vec_print(v1);
            vec_print(v2);
        }
        REQUIRE( is_ok );
    }
}

//
//
///////////////////////////////////////////////////////////////////////////////////////////
//
//

TEST_CASE( "Testing ", "[quadsort]" )
{
    const int v_min = 256;
    const int v_max = 262144;

    for(int i = v_min; i < v_max; i *= 2)
    {
        std::vector<uint64_t> v1( i );
        std::vector<uint64_t> v2( i );
        vec_random(v1, i);
        vec_copy(v2, v1);
        std::sort(v1.begin(), v1.end());
        quadsort_prim(v2.data(), v2.size(), 8); // code => int64_t
        const bool is_ok = vec_compare(v1, v2);
        if( is_ok == false )
        {
            vec_print(v1);
            vec_print(v2);
        }
        REQUIRE( is_ok );
    }
}

//
//
///////////////////////////////////////////////////////////////////////////////////////////
//
//

TEST_CASE( "Testing ", "[crumsort]" )
{
    const int v_min = 256;
    const int v_max = 262144;

    for(int i = v_min; i < v_max; i *= 2)
    {
        std::vector<uint64_t> v1( i );
        std::vector<uint64_t> v2( i );
        vec_random(v1, i);
        vec_copy(v2, v1);
        std::sort(v1.begin(), v1.end());
        crumsort_prim(v2.data(), v2.size(), 8); // code => int64_t
        const bool is_ok = vec_compare(v1, v2);
        if( is_ok == false )
        {
            vec_print(v1);
            vec_print(v2);
        }
        REQUIRE( is_ok );
    }
}

//
//
///////////////////////////////////////////////////////////////////////////////////////////
//
//

TEST_CASE( "Testing ", "[crumsort_2cores]" )
{
    const int v_min = 256;
    const int v_max = 262144;

    for(int i = v_min; i < v_max; i *= 2)
    {
        std::vector<uint64_t> v1( i );
        std::vector<uint64_t> v2( i );
        vec_random(v1, i);
        vec_copy(v2, v1);
        std::sort(v1.begin(), v1.end());
        crumsort_2cores(v2 ); // code => uint64_t
        const bool is_ok = vec_compare(v1, v2);
        if( is_ok == false )
        {
            vec_print(v1);
            vec_print(v2);
        }
        REQUIRE( is_ok );
    }
}
