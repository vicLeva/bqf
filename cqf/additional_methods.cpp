#include <iostream>
#include <bitset>
#include <stdint.h> 
#include <cassert>
#include <vector>
#include <string.h>

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
    /*if( !((rank != 0) && (rank <= MEM_UNIT)) ) { 
        std::cout << "rank " << rank << "\n";
    }
    assert((rank != 0) && (rank <= MEM_UNIT));*/
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

uint64_t encode(string kmer){
    uint64_t encoded = 0;
    for(char& c : kmer) {
        if (c=='A'){
            encoded <<= 2;
            encoded |= 3;
        }
        else if (c=='C'){
            encoded <<= 2;
            encoded |= 2;
        }
        else if (c=='G'){
            encoded <<= 2;
            encoded |= 1;
        }
        else{ //T is 00 so with reverse complementarity we won't get 0000000000000 as input for xorshift
            encoded <<= 2;
        }
    }

    return encoded;
}

string decode(uint64_t revhash, uint64_t size){
    string kmer;

    for (size_t i=0; i<size; i++){
        switch(revhash & mask_right(2)){
            case 3:
                kmer = 'A' + kmer;
                break;
            case 2:
                kmer = 'C' + kmer;
                break;
            case 1:
                kmer = 'G' + kmer;
                break;
            default:
                kmer = 'T' + kmer;
        }
        revhash >>= 2;
    }

    return kmer;
}



uint64_t bfc_hash_64(uint64_t key, uint64_t mask) {
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}


uint64_t bfc_hash_64_inv(uint64_t key, uint64_t mask){
	uint64_t tmp;
 
	// Invert key = key + (key << 31)
	tmp = (key - (key << 31));
	key = (key - (tmp << 31)) & mask;
 
	// Invert key = key ^ (key >> 28)
	tmp = key ^ key >> 28;
	key = key ^ tmp >> 28;
 
	// Invert key *= 21
	key = (key * 14933078535860113213ull) & mask;
 
	// Invert key = key ^ (key >> 14)
	tmp = key ^ key >> 14;
	tmp = key ^ tmp >> 14;
	tmp = key ^ tmp >> 14;
	key = key ^ tmp >> 14;
 
	// Invert key *= 265
	key = (key * 15244667743933553977ull) & mask;
 
	// Invert key = key ^ (key >> 24)
	tmp = key ^ key >> 24;
	key = key ^ tmp >> 24;
 
	// Invert key = (~key) + (key << 21)
	tmp = ~key;
	tmp = ~(key - (tmp << 21));
	tmp = ~(key - (tmp << 21));
	key = ~(key - (tmp << 21)) & mask;
 
	return key;
}


uint64_t kmer_to_hash(string kmer, uint64_t k){
    return bfc_hash_64(encode(kmer), mask_right(k*2));
}

string hash_to_kmer(uint64_t hash, uint64_t k){
    return decode(bfc_hash_64_inv(hash, mask_right(k*2)), k);
}


char complement(char nucl){
  //Returns the complement of a nucleotide
  switch (nucl){
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'C':
      return 'G';
    default:
      return 'C';
  }
}

std::string revcomp(const std::string& kmer) {
    std::string revComp(kmer.rbegin(), kmer.rend());
    for (auto& c : revComp) {
        c = complement(c);
    }
    return revComp;
}

// Returns an array of canonical k-mers in a given sequence s
std::vector<std::string> canonical_kmers(const std::string& s, int len, int k) {
    std::vector<std::string> kmers{};
    string kmer = s.substr(0, k);
    string rc = revcomp(kmer);

    if (kmer < rc) {
        kmers.push_back(kmer);
    } else {
        kmers.push_back(rc);
    }

    for (int i = k; i < len; i++) {
        kmer = kmer.substr(1) + s[i];
        rc = complement(s[i]) + rc.substr(0, k-1);
        if (kmer < rc) {
            kmers.push_back(kmer);
        } else {
            kmers.push_back(rc);
        }
    }

    return kmers;
}