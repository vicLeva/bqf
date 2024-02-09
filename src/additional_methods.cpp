#include "additional_methods.hpp"

#include "./generic/bitrankasm.hpp"
#include "./generic/bitselectasm.hpp"

//BITS & KMER MANIPULATION

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

uint64_t bitselectasm(uint64_t num, uint64_t rank)
{
#if defined(__aarch64__) || defined(_M_ARM64)
    return bitselectasm_u64_arm(num, rank);
#else
    return bitselectasm_u64_builtin(num, rank);
#endif
//    uint64_t i = 1ULL << (rank - 1); // i = 2^rank
//    asm("pdep %[num], %[mask], %[num]"
//            : [num] "+r" (num)
//            : [mask] "r" (i));
//
//    asm("tzcnt %[bit], %[index]"
//            : [index] "=r" (i)
//            : [bit] "g" (num)
//            : "cc");
//
//    return i;
}

uint64_t bitrankasm(uint64_t val, uint64_t pos)
{
#if defined(__aarch64__) || defined(_M_ARM64)
    return bitrankasm_arm(val, pos);
#else
    return bitrankasm_u64_builtin(val, pos);
#endif
//    assert(pos < MEM_UNIT);
//    val = val & ((2ULL << pos) - 1);
//
//    // POPCOUNT(v & (2^i − 1)
//    asm("popcnt %[val], %[val]"
//            : [val] "+r" (val)
//            :
//            : "cc");
//    return val;
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

string decode(uint64_t revhash, uint64_t k){
    string kmer;

    for (size_t i=0; i<k; i++){
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


// Thomas Wang's integer hash functions. See <https://gist.github.com/lh3/59882d6b96166dfc3d8d> for a snapshot.
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


uint64_t nucl_encode(char nucl){
  //Returns the binary encoding of a nucleotide
  //different from encode() function because this is for query, so we have to check for canonical smer
  //and this encoding allows fast comparison (<) of lexico order 
  switch (nucl){
    case 'A':
      return 0;
    case 'T':
      return 3;
    case 'C':
      return 1;
    case 'G':
      return 2;
    default :
        cout << "non nucl : " << nucl << endl;
        throw std::invalid_argument( "received non nucleotidic value" );
  }
}

uint64_t flip(uint64_t encoding, size_t bitsize){
  //Used so Ts are 0s so we dont get 0000.. in xorshift (through revcompl)
    return ~encoding << (64-bitsize) >> (64-bitsize);
}

uint64_t revcomp64 (const uint64_t v, size_t bitsize){
  return (((uint64_t)rev_table[v & 0xff] << 56) | 
    ((uint64_t)rev_table[(v >> 8) & 0xff] << 48) | 
    ((uint64_t)rev_table[(v >> 16) & 0xff] << 40) | 
    ((uint64_t)rev_table[(v >> 24) & 0xff] << 32) | 
    ((uint64_t)rev_table[(v >> 32) & 0xff] << 24) | 
    ((uint64_t)rev_table[(v >> 40) & 0xff] << 16) |
    ((uint64_t)rev_table[(v >> 48) & 0xff] << 8) |
    ((uint64_t)rev_table[(v >> 56) & 0xff])) >> (64-bitsize);
}

uint64_t canonical(uint64_t smer, size_t size){
    uint64_t revcomp = revcomp64(smer, size);
    if (revcomp < smer) { return revcomp; }
    else { return smer; }
}

std::string canonical(const std::string& smer, size_t s){
  //debug purpose only
  return decode(flip(canonical(flip(encode(smer), 2*s), 2*s), 2*s), s);
}

std::ostream& operator<<(std::ostream& os, result_query const& res) {
    return os << "(min:" << res.minimum << ", max:" << res.maximum << ", average:" << res.average << ", presence ratio:" << res.kmer_present_ratio << ")" << endl;
}
