#include <iostream>
#include <bitset>
#include <stdint.h> 
#include <cassert>
#include <vector>

#include "additional_methods.hpp"


#define MEM_UNIT 64ULL
#define MET_UNIT 3ULL
#define OFF_POS 0ULL
#define OCC_POS 1ULL
#define RUN_POS 2ULL


void print_bits(uint64_t x) {
  std::bitset<MEM_UNIT> bits(x);
  std::cout << bits << std::endl;
}


uint64_t mask_right(uint64_t numbits){
    uint64_t mask = -(numbits >= MEM_UNIT) | ((1ULL << numbits) - 1ULL);
    return mask;
}

uint64_t mask_left(uint64_t numbits){
    return ~(mask_right(MEM_UNIT-numbits));
}

uint64_t shift_left(uint64_t value, uint64_t shift){
    if (shift == MEM_UNIT) return 0;
    else return (value << shift);
}

uint64_t shift_right(uint64_t value, uint64_t shift){
    if (shift == MEM_UNIT) return 0;
    else return (value >> shift);
}

uint64_t rebuild_number(uint64_t quot, uint64_t rem, uint64_t shift){
    rem = shift_left(rem, shift);
    return (rem | quot);
}

uint64_t get_block_id(uint64_t position){
    return position / MEM_UNIT;
}

uint64_t get_shift_in_block(uint64_t position){
    return position % MEM_UNIT;
}

uint64_t get_quot_from_block_shift(uint64_t block, uint64_t shift){
    return block*MEM_UNIT + shift;
}

uint64_t bitselectasm(uint64_t num, uint64_t rank){
    if( !((rank != 0) && (rank <= MEM_UNIT)) ) { 
        std::cout << "rank " << rank << "\n";
    }
    assert((rank != 0) && (rank <= MEM_UNIT));
    uint64_t i = 1ULL << (rank - 1); // i = 2^rank

    // SELECT(v,i) = TZCNT(PDEP(2^rank,num))     
	asm("pdep %[num], %[mask], %[num]"
			: [num] "+r" (num)
			: [mask] "r" (i));

	asm("tzcnt %[bit], %[index]"
			: [index] "=r" (i)
			: [bit] "g" (num)
			: "cc");

	return i;
}

uint64_t bitrankasm(uint64_t val, uint64_t pos) {
    assert(pos < MEM_UNIT);
	val = val & ((2ULL << pos) - 1);

    // POPCOUNT(v & (2^i − 1)
	asm("popcnt %[val], %[val]"
			: [val] "+r" (val)
			:
			: "cc");
	return val;
}

uint64_t get_bit_from_word(uint64_t word, uint64_t pos_bit){
    return ((word >> pos_bit) & 0b1);
}

uint64_t get_bits(std::vector<uint64_t>& vec, uint64_t pos, uint64_t len){

    if (!len) return 0;

    uint64_t word = get_block_id(pos);
    uint64_t shift = get_shift_in_block(pos);
    uint64_t mask = mask_right(len);

    if (shift + len <= MEM_UNIT) return (vec[word] >> shift) & mask;

    return (vec[word] >> shift) | ((vec[word+1] << (MEM_UNIT - shift)) & mask);
    //return (cqf[block] >> shift) | ((cqf[block+1] & mask_right(len-(MEM_UNIT-shift))) << (MEM_UNIT - shift));
}

using namespace std;
void set_bits(std::vector<uint64_t>& vec, uint64_t pos, uint64_t value, uint64_t len) {
    assert(pos + len <= vec.size() * 64);
    assert(len == MEM_UNIT or (value >> len) == 0);
    if (len == 0) return;

    uint64_t mask = mask_right(len);
    uint64_t word = get_block_id(pos); //not rly block id, more like word id (pos is exact bit pos)
    uint64_t shift = get_shift_in_block(pos);

    value &= mask;

    vec[word] &= ~(mask << shift); //set to 0 all bits in the range of modified ones
    vec[word] |= (value << shift); //OR op between the previous 0s and bits of value

    uint64_t stored = MEM_UNIT - shift;

    if (len > stored){
        vec[word+1] &= ~(mask_right(len-stored));
        vec[word+1] |= (value >> stored);
    }
}