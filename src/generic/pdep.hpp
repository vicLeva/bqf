#ifndef __pdep__
#define __pdep__
#include <bit>
#include <bitset>
#include <cstdint>
#include <iostream>
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
#if defined(__aarch64__) || defined(_M_ARM64)
    #include <arm_neon.h>
#else
    #include <immintrin.h>
#endif
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
inline int pdep_u32_c(uint32_t val, uint32_t mask)
{
  unsigned int res = 0;
  for (uint32_t bb = 1; mask; bb += bb) {
    if (val & bb)
      res |= mask & -mask;
    mask &= mask - 1;
  }
  return res;
}
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
inline int64_t pdep_u64_c(const uint64_t val, uint64_t mask)
{
  int64_t res = 0;
  for (uint64_t bb = 1; mask; bb += bb) {
    if (val & bb)
      res |= mask & -mask;
    mask &= mask - 1;
  }
  return res;
}
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
#if defined(__SSE4_2__)
inline int pdep_u32_builtin(uint32_t val, uint32_t mask)
{
  return _pdep_u32(val, mask);
}
#endif
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
#if defined(__SSE4_2__)
inline int64_t pdep_u64_builtin(const uint64_t val, const uint64_t mask)
{
  return _pdep_u64(val, mask);
}
#endif
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
#if defined(__SSE4_2__)
int pdep_u32_x86(uint32_t num, uint32_t rank)
{
    int rsu;
	  asm("pdep %[mask], %[num], %[rsu]"
		    : [rsu]  "=r" (rsu)
        : [num]  "r" (val),
          [mask] "r" (mask)
    );
    return rsu;
}
#endif
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
#if defined(__SSE4_2__)
inline int64_t pdep_u64_x86 (const uint64_t val, const uint64_t mask)
{
    int64_t rsu;
	asm("pdep %[mask], %[num], %[rsu]"
			: [rsu]  "=r" (rsu)
            : [num]  "r" (val),
              [mask] "r" (mask)
        );
    return rsu;
}
#endif
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
inline int64_t pdep_u64_arm (const uint64_t val, const uint64_t _mask)
{
  uint64_t mask = _mask;
  int64_t res = 0;
  for (uint64_t bb = 1; mask; bb += bb) {
    if (val & bb)
      res |= mask & -mask;
    mask &= mask - 1;
  }
  return res;
}
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
#endif