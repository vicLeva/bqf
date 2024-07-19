#include "tools.hpp"

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

void vec_copy(std::vector<uint64_t>& dst, std::vector<uint64_t>& src)
{
    if(dst.size() != src.size())
        dst.resize( src.size() );

    const int length = src.size();
    for(int x = 0; x < length; x += 1)
        dst[x] = src[x];
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

void vec_copy(std::vector<uint64_t>& dst, uint64_t* src, int length)
{
    if(dst.size() != length)
        dst.resize( length );

    for(int x = 0; x < length; x += 1)
        dst[x] = src[x];
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

void vec_copy(uint64_t* dst, uint64_t* src, int length)
{
    for(int x = 0; x < length; x += 1)
        dst[x] = src[x];
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

bool vec_compare(std::vector<uint64_t>& dst, std::vector<uint64_t>& src)
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

void vec_print(std::vector<uint64_t>& src)
{
    printf("(II) Printed values :");
    for (int i = 0; i < src.size(); i += 1) {
        if ((i % 32) == 0) printf("\n%3d | ", i);
        printf("%3llu ", src[i]);
    }
    printf("\n");
}

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

void vec_print(uint64_t* src, int length)
{
    printf("(II) Printed values :");
    for (int i = 0; i < length; i += 1)
    {
        if ((i % 32) == 0) printf("\n%3d | ", i);
        printf("%3llu ", src[i]);
    }
    printf("\n");
}

