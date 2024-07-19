#ifndef __bitselectasm__
#define __bitselectasm__
#include "pdep.hpp"
#include "tzcnt.hpp"
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
inline uint64_t bitselectasm_u64_c(uint64_t num, uint64_t rank){
    const uint64_t i = 1ULL << (rank - 1); // i = 2^rank
    const uint64_t j = pdep_u64_c(i, num);
    const uint64_t k = tzcnt_u64_c(j);
	return k;
}
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
#if defined(__aarch64__) || defined(_M_ARM64)
inline uint64_t bitselectasm_u64_arm(uint64_t num, uint64_t rank){
    const uint64_t i = 1ULL << (rank - 1); // i = 2^rank
    const uint64_t j = pdep_u64_arm(i, num);
    const uint64_t k = tzcnt_u64_arm(j);
	return k;
}
#endif
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
#if defined(__SSE4_2__)
inline uint64_t bitselectasm_u64_builtin(uint64_t num, uint64_t rank){
#if 0
    const uint64_t i = 1ULL << (rank - 1); // i = 2^rank
	const uint64_t j = _pdep_u64(i, num);
    const uint64_t k = __builtin_ctzll(j);
	return k;
#else
    const uint64_t i = 1ULL << (rank - 1); // i = 2^rank
	const uint64_t j = pdep_u64_builtin(i, num);
    const uint64_t k = tzcnt_u64_builtin(j);
	return k;
#endif
}
#endif
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
#if defined(__SSE4_2__)
inline uint64_t bitselectasm_u64_x86(uint64_t num, uint64_t rank){
#if 0
    uint64_t i = 1ULL << (rank - 1); // i = 2^rank
	asm("pdep %[num], %[mask], %[num]"
			: [num] "+r" (num)
			: [mask] "r" (i));
	asm("tzcnt %[bit], %[index]" : [index] "=r" (i) : [bit] "g" (num) : "cc");
	return i;
#else
    const uint64_t i = 1ULL << (rank - 1); // i = 2^rank
    const uint64_t j = pdep_u64_x86(i, num);
    const uint64_t k = tzcnt_u64_x86(j);
    return k;
#endif
}
#endif
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
#endif
