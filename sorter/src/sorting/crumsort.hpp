// crumsort 1.2.1.3 - Igor van den Hoven ivdhoven@gmail.com

#ifndef CRUMSORT_H
#define CRUMSORT_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include <stdalign.h>
#include <float.h>
#include <string.h>

typedef int CMPFUNC (const void *a, const void *b);

//#define cmp(a,b) (*(a) > *(b))

        #ifndef QUADSORT_H
  #include "quadsort.hpp"
#endif

// When sorting an array of pointers, like a string array, the QUAD_CACHE needs
// to be set for proper performance when sorting large arrays.
// crumsort_prim() can be used to sort arrays of 32 and 64 bit integers
// without a comparison function or cache restrictions.

// With a 6 MB L3 cache a value of 262144 works well.

#ifdef cmp
  #define QUAD_CACHE 4294967295
#else
//#define QUAD_CACHE 131072
#define QUAD_CACHE 262144
//#define QUAD_CACHE 524288
//#define QUAD_CACHE 4294967295
#endif

//////////////////////////////////////////////////////////
// ┌───────────────────────────────────────────────────┐//
// │       ██████┐ ██████┐    ██████┐ ██████┐████████┐ │//
// │       └────██┐└────██┐   ██┌──██┐└─██┌─┘└──██┌──┘ │//
// │        █████┌┘ █████┌┘   ██████┌┘  ██│     ██│    │//
// │        └───██┐██┌───┘    ██┌──██┐  ██│     ██│    │//
// │       ██████┌┘███████┐   ██████┌┘██████┐   ██│    │//
// │       └─────┘ └──────┘   └─────┘ └─────┘   └─┘    │//
// └───────────────────────────────────────────────────┘//
//////////////////////////////////////////////////////////

#define VAR int
#define FUNC(NAME) NAME##32

#include "crumsort.hxx"

#undef VAR
#undef FUNC

// crumsort_prim

#define VAR int
#define FUNC(NAME) NAME##_int32
#ifndef cmp
  #define cmp(a,b) (*(a) > *(b))
  #include "crumsort.hxx"
  #undef cmp
#else
  #include "crumsort.c"
#endif
#undef VAR
#undef FUNC

#define VAR unsigned int
#define FUNC(NAME) NAME##_uint32
#ifndef cmp
  #define cmp(a,b) (*(a) > *(b))
  #include "crumsort.hxx"
  #undef cmp
#else
  #include "crumsort.c"
#endif
#undef VAR
#undef FUNC

//////////////////////////////////////////////////////////
// ┌───────────────────────────────────────────────────┐//
// │        █████┐ ██┐  ██┐   ██████┐ ██████┐████████┐ │//
// │       ██┌───┘ ██│  ██│   ██┌──██┐└─██┌─┘└──██┌──┘ │//
// │       ██████┐ ███████│   ██████┌┘  ██│     ██│    │//
// │       ██┌──██┐└────██│   ██┌──██┐  ██│     ██│    │//
// │       └█████┌┘     ██│   ██████┌┘██████┐   ██│    │//
// │        └────┘      └─┘   └─────┘ └─────┘   └─┘    │//
// └───────────────────────────────────────────────────┘//
//////////////////////////////////////////////////////////

#define VAR long long
#define FUNC(NAME) NAME##64

#include "crumsort.hxx"

#undef VAR
#undef FUNC

// crumsort_prim

#define VAR long long
#define FUNC(NAME) NAME##_int64
#ifndef cmp
  #define cmp(a,b) (*(a) > *(b))
  #include "crumsort.hxx"
  #undef cmp
#else
  #include "crumsort.c"
#endif
#undef VAR
#undef FUNC

#define VAR unsigned long long
#define FUNC(NAME) NAME##_uint64
#ifndef cmp
  #define cmp(a,b) (*(a) > *(b))
  #include "crumsort.hxx"
  #undef cmp
#else
  #include "crumsort.c"
#endif
#undef VAR
#undef FUNC

// This section is outside of 32/64 bit pointer territory, so no cache checks
// necessary, unless sorting 32+ byte structures.

#undef QUAD_CACHE
#define QUAD_CACHE 4294967295

//////////////////////////////////////////////////////////
//┌────────────────────────────────────────────────────┐//
//│                █████┐    ██████┐ ██████┐████████┐  │//
//│               ██┌──██┐   ██┌──██┐└─██┌─┘└──██┌──┘  │//
//│               └█████┌┘   ██████┌┘  ██│     ██│     │//
//│               ██┌──██┐   ██┌──██┐  ██│     ██│     │//
//│               └█████┌┘   ██████┌┘██████┐   ██│     │//
//│                └────┘    └─────┘ └─────┘   └─┘     │//
//└────────────────────────────────────────────────────┘//
//////////////////////////////////////////////////////////

#define VAR char
#define FUNC(NAME) NAME##8

#include "crumsort.hxx"

#undef VAR
#undef FUNC

//////////////////////////////////////////////////////////
//┌────────────────────────────────────────────────────┐//
//│           ▄██┐   █████┐    ██████┐ ██████┐████████┐│//
//│          ████│  ██┌───┘    ██┌──██┐└─██┌─┘└──██┌──┘│//
//│          └─██│  ██████┐    ██████┌┘  ██│     ██│   │//
//│            ██│  ██┌──██┐   ██┌──██┐  ██│     ██│   │//
//│          ██████┐└█████┌┘   ██████┌┘██████┐   ██│   │//
//│          └─────┘ └────┘    └─────┘ └─────┘   └─┘   │//
//└────────────────────────────────────────────────────┘//
//////////////////////////////////////////////////////////

#define VAR short
#define FUNC(NAME) NAME##16

#include "crumsort.hxx"

#undef VAR
#undef FUNC

//////////////////////////////////////////////////////////
//┌────────────────────────────────────────────────────┐//
//│  ▄██┐  ██████┐  █████┐    ██████┐ ██████┐████████┐ │//
//│ ████│  └────██┐██┌──██┐   ██┌──██┐└─██┌─┘└──██┌──┘ │//
//│ └─██│   █████┌┘└█████┌┘   ██████┌┘  ██│     ██│    │//
//│   ██│  ██┌───┘ ██┌──██┐   ██┌──██┐  ██│     ██│    │//
//│ ██████┐███████┐└█████┌┘   ██████┌┘██████┐   ██│    │//
//│ └─────┘└──────┘ └────┘    └─────┘ └─────┘   └─┘    │//
//└────────────────────────────────────────────────────┘//
//////////////////////////////////////////////////////////

// 128 reflects the name, though the actual size of a long double is 64, 80,
// 96, or 128 bits, depending on platform.

#if (DBL_MANT_DIG < LDBL_MANT_DIG)
  #define VAR long double
  #define FUNC(NAME) NAME##128
    #include "crumsort.c"
  #undef VAR
  #undef FUNC
#endif

///////////////////////////////////////////////////////////
//┌─────────────────────────────────────────────────────┐//
//│ ██████┐██┐   ██┐███████┐████████┐ ██████┐ ███┐  ███┐│//
//│██┌────┘██│   ██│██┌────┘└──██┌──┘██┌───██┐████┐████││//
//│██│     ██│   ██│███████┐   ██│   ██│   ██│██┌███┌██││//
//│██│     ██│   ██│└────██│   ██│   ██│   ██│██│└█┌┘██││//
//│└██████┐└██████┌┘███████│   ██│   └██████┌┘██│ └┘ ██││//
//│ └─────┘ └─────┘ └──────┘   └─┘    └─────┘ └─┘    └─┘│//
//└─────────────────────────────────────────────────────┘//
///////////////////////////////////////////////////////////

/*
typedef struct {char bytes[32];} struct256;
#define VAR struct256
#define FUNC(NAME) NAME##256

#include "crumsort.c"

#undef VAR
#undef FUNC
*/

 //////////////////////////////////////////////////////////////////////////
//┌─────────────────────────────────────────────────────────────────────┐//
//│ ██████┐██████┐ ██┐   ██┐███┐  ███┐███████┐ ██████┐ ██████┐ ████████┐│//
//│██┌────┘██┌──██┐██│   ██│████┐████│██┌────┘██┌───██┐██┌──██┐└──██┌──┘│//
//│██│     ██████┌┘██│   ██│██┌███┌██│███████┐██│   ██│██████┌┘   ██│   │//
//│██│     ██┌──██┐██│   ██│██│└█┌┘██│└────██│██│   ██│██┌──██┐   ██│   │//
//│└██████┐██│  ██│└██████┌┘██│ └┘ ██│███████│└██████┌┘██│  ██│   ██│   │//
//│ └─────┘└─┘  └─┘ └─────┘ └─┘    └─┘└──────┘ └─────┘ └─┘  └─┘   └─┘   │//
//└─────────────────────────────────────────────────────────────────────┘//
//////////////////////////////////////////////////////////////////////////

void crumsort(void *array, size_t nmemb, size_t size, CMPFUNC *cmp)
{
	if (nmemb < 2)
	{
		return;
	}

	switch (size)
	{
		case sizeof(char):
			crumsort8(array, nmemb, cmp);
			return;

		case sizeof(short):
			crumsort16(array, nmemb, cmp);
			return;

		case sizeof(int):
			crumsort32(array, nmemb, cmp);
			return;

		case sizeof(long long):
			crumsort64(array, nmemb, cmp);
			return;
#if (DBL_MANT_DIG < LDBL_MANT_DIG)
		case sizeof(long double):
			crumsort128(array, nmemb, cmp);
			return;
#endif
//		case sizeof(struct256):
//			crumsort256(array, nmemb, cmp);
			return;

		default:
#if (DBL_MANT_DIG < LDBL_MANT_DIG)
			assert(size == sizeof(char) || size == sizeof(short) || size == sizeof(int) || size == sizeof(long long) || size == sizeof(long double));
#else
			assert(size == sizeof(char) || size == sizeof(short) || size == sizeof(int) || size == sizeof(long long));
#endif
//			qsort(array, nmemb, size, cmp);
	}
}

// suggested size values for primitives:

//		case  0: unsigned char
//		case  1: signed char
//		case  2: signed short
//		case  3: unsigned short
//		case  4: signed int
//		case  5: unsigned int
//		case  6: float
//		case  7: double
//		case  8: signed long long
//		case  9: unsigned long long
//		case  ?: long double, use sizeof(long double):

void crumsort_prim(void *array, size_t nmemb, size_t size)
{
	if (nmemb < 2)
	{
		return;
	}

	switch (size)
	{
		case 4:
			crumsort_int32(array, nmemb, NULL);
			return;
		case 5:
			crumsort_uint32(array, nmemb, NULL);
			return;
		case 8:
			crumsort_int64(array, nmemb, NULL);
			return;
		case 9:
			crumsort_uint64(array, nmemb, NULL);
			return;
		default:
			assert(size == sizeof(int) || size == sizeof(int) + 1 || size == sizeof(long long) || size == sizeof(long long) + 1);
			return;
	}
}

#undef QUAD_CACHE

#endif
