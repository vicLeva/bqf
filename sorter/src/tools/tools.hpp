#pragma once
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>

//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//

extern void vec_copy(std::vector<uint64_t>& dst, std::vector<uint64_t>& src            );
extern void vec_copy(std::vector<uint64_t>& dst, uint64_t* src,              int length);
extern void vec_copy(            uint64_t*  dst, uint64_t* src,              int length);

extern bool vec_compare(std::vector<uint64_t>& dst, std::vector<uint64_t>& src         );
extern void vec_print  (std::vector<uint64_t>& src);
extern void vec_print  (uint64_t* src, int length);
