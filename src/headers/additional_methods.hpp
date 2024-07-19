#ifndef ADDITIONAL_METHODS_HPP
#define ADDITIONAL_METHODS_HPP

#include <stdint.h> 
#include <string>
#include <iostream>
#include <bitset>
#include <cassert>
#include <vector>


// STATIC VARIABLES 
#define MEM_UNIT 64ULL
#define BLOCK_SIZE 64ULL
#define MET_UNIT 3ULL
#define OFF_POS 0ULL
#define OCC_POS 1ULL
#define RUN_POS 2ULL
#define SCALE_INPUT 8388608ULL

/** 
 * \brief Print the bits of a given uint64 word
 * It prints on the standard output the bits of the uint64 word
 * \param x the uint64 of which you want to see the bits
 */
void print_bits(uint64_t x);

/** 
 * \brief returns a mask of 1s of the leftmost Nth bits
 * e.g. mask_left(5) gives 11111000000000000...
 * \param num_bits number of 1s from the left in the mask
 * \return the uint64 corresponding to the mask
 */
uint64_t mask_left(uint64_t num_bits);

/** 
 * \brief returns a mask of 1s of the rightmost Nth bits
 * e.g. mask_right(5) gives 000...00011111
 * \param num_bits number of 1s from the right in the mask
 * \return the uint64 corresponding to the mask
 */
uint64_t mask_right(uint64_t num_bits);

/** 
 * \brief returns the value shifted left of shift bits
 * if shift is 64 returns a 0
 * \param value to shift
 * \param shift applied to value
 * \return the uint64 corresponding to the shifted value
 */
uint64_t shift_left(uint64_t value, uint64_t shift);

/** 
 * \brief returns the value shifted right of shift bits
 * if shift is 64 returns a 0
 * \param value to shift
 * \param shift applied to value
 * \return the uint64 corresponding to the shifted value
 */
uint64_t shift_right(uint64_t value, uint64_t shift);

/** 
 * \brief gets a hash back from quotient & remainder
 * Used for enumerating and resizing the QF, because remainder and quotient were
 * computed from splitting the hash, cancatenating them together gives us the original hash
 * \param quotient the quotient whose run including the remainder belongs to
 * \param remainder the remainder of the original hash
 * \param shift the quotient size (in the QF)
 * \return the uint64 corresponding to original hash
 */
uint64_t rebuild_number(uint64_t quotient, uint64_t remainder, uint64_t shift);

/** 
 * \brief returns the id of the block associated to the requested quotient
 * e.g. if blocks have 64 quotients inside each, 
 * get_block_id(155) = 2 (there is block 0)
 * \param position the quotient of wich you want to know the block
 * \return the block id
 */
uint64_t get_block_id(uint64_t position);

/** 
 * \brief returns the position in the (corresponding) block of the quotient
 * e..g. if blocks have 64 quotients inside,
 * get_shit_in_block(155) = 28
 * \param position the quotient of which you want to know the position in its block
 * \return the position in which it is in the block
 */
uint64_t get_shift_in_block(uint64_t position);

/** 
 * \brief gives back the quotient corresponding to a particular block and position in the block
 * e.g. get_quot_from_block_shift(2,28) = 155
 * \param block the particular block 
 * \param shift the shitft in the block (eg the position in the block)
 * \return a uint64 corresponding to the quotient associated with them
 */
uint64_t get_quot_from_block_shift(uint64_t block, uint64_t shift);

/** 
 * \brief bitwise select operation
 * \param num word of the bitvector of runends
 * \param rank rank position
 * \return the select operation
 */
uint64_t bitselectasm(uint64_t num, uint64_t rank);

/** 
 * \brief bitwise rank operation
 * \param val   word of the bitvector of occupieds
 * \param pos   position of the bit in the word
 * \return the rank operation
 */
uint64_t bitrankasm(uint64_t val, uint64_t pos);

/** 
 * \brief get a specific bit from a uint64 word given its position in it
 * \param word  word from which you want a bit
 * \param pos_bit   position of the bit in the word
 * \return uint64 with the bit (in the rightmost pos) corresponding to the extracted bit
 */
uint64_t get_bit_from_word(uint64_t word, uint64_t pos_bit);

/** 
 * \brief get 'len' bits of the vector from the bitposition pos 
 * \param vec   vector from which you want to extract bits
 * \param pos   starting position of the bits you want to extract (in bits)
 * \param len   how many bits you want to extract (up to 64)
 * \return a uint64 with the bits you extracted
 */
uint64_t get_bits(std::vector<uint64_t>& vec, uint64_t pos, uint64_t len);

/** 
 * \brief set 'len' bits of the vector from the bitposition pos. Up to 64 bits at a time (value is a uint64)
 * \param vec vector you want to modify
 * \param pos   position of the bit where you want to start editing the vector
 * \param value uint64 with bits that you want to cast in the vector
 * \param len   number of bits that are gonna be effectively casted in the vector
 */
void set_bits(std::vector<uint64_t>& vec, uint64_t pos, uint64_t value, uint64_t len);

/** 
 * \brief encode a kmer into a number by encoding each nucl to 2bits
 * \param data the kmer you want to encode
 * A->11 C->10 G->01 T->00, T is 00 because we'll never have a kmer full of T 
 * as its reverse complement will be smaller, so we never have 0 as encoded value,
 * which wd be problematic for hashing
 * \return a uint64, encoding the kmer
 */
uint64_t encode(std::string data);
uint64_t encode(const char* kmer, const int kmer_size);

/** 
 * \brief decode a number into a kmer
 * \param hash the number to turn into a kmer
 * \param size the size of the kmer, to know where to stop in the 64bits number
 * \return the kmer string decoded
 */
std::string decode(uint64_t hash, uint64_t size);

/** 
 * \brief modular and reversible xorshift64
 * \param key the number to hash (usually an encoded kmer)
 * \param mask the hash output size (in bits) as a mask (see mask_right())
 * from https://github.com/lh3/bfc/blob/master/kmer.h
 * \return uint64_t, the hash
 */
uint64_t bfc_hash_64(uint64_t key, uint64_t mask);

/** 
 * \brief the reverse of xorshift64 bfc_hash_64()
 * \param key the number to dehash
 * \param mask the hash output size (in bits) as a mask (see mask_right())
 * from https://github.com/lh3/bfc/blob/master/kmer.h
 * \return uint64_t, the originally hashed number
 */
uint64_t bfc_hash_64_inv(uint64_t key, uint64_t mask);

/** 
 * \brief encode the smallest between the kmer and its revcompl, then hash it
 * \param data the kmer to hash
 * \param size the kmer size (in base)
 * \return uint64_t, the hash
 */
uint64_t kmer_to_hash(std::string data, uint64_t size);

/** 
 * \brief reverse the hash, then decode the number into a kmer (canonical)
 * \param data the hash to reverse
 * \param size the kmer size (in base)
 * \return the kmer (string)
 */
std::string hash_to_kmer(uint64_t hash, uint64_t size);

uint64_t nucl_encode(char nucl);

//uint64_t flip(uint64_t encoding, size_t bitsize);

const uint64_t rev_table [256] = {170, 234, 42, 106, 186, 250, 58, 122, 138, 202, 10, 74, 154, 218, 26, 90, 174, 238, 46, 110, 190, 254, 62, 126, 142, 206, 14, 78, 158, 222, 30, 94, 162, 226, 34, 98, 178, 242, 50, 114, 130, 194, 2, 66, 146, 210, 18, 82, 166, 230, 38, 102, 182, 246, 54, 118, 134, 198, 6, 70, 150, 214, 22, 86, 171, 235, 43, 107, 187, 251, 59, 123, 139, 203, 11, 75, 155, 219, 27, 91, 175, 239, 47, 111, 191, 255, 63, 127, 143, 207, 15, 79, 159, 223, 31, 95, 163, 227, 35, 99, 179, 243, 51, 115, 131, 195, 3, 67, 147, 211, 19, 83, 167, 231, 39, 103, 183, 247, 55, 119, 135, 199, 7, 71, 151, 215, 23, 87, 168, 232, 40, 104, 184, 248, 56, 120, 136, 200, 8, 72, 152, 216, 24, 88, 172, 236, 44, 108, 188, 252, 60, 124, 140, 204, 12, 76, 156, 220, 28, 92, 160, 224, 32, 96, 176, 240, 48, 112, 128, 192, 0, 64, 144, 208, 16, 80, 164, 228, 36, 100, 180, 244, 52, 116, 132, 196, 4, 68, 148, 212, 20, 84, 169, 233, 41, 105, 185, 249, 57, 121, 137, 201, 9, 73, 153, 217, 25, 89, 173, 237, 45, 109, 189, 253, 61, 125, 141, 205, 13, 77, 157, 221, 29, 93, 161, 225, 33, 97, 177, 241, 49, 113, 129, 193, 1, 65, 145, 209, 17, 81, 165, 229, 37, 101, 181, 245, 53, 117, 133, 197, 5, 69, 149, 213, 21, 85};

uint64_t revcomp64 (const uint64_t v, size_t bitsize);

uint64_t canonical(uint64_t smer, size_t size);

std::string canonical(const std::string& smer, size_t s);

/** 
 * \struct result_query
 * \brief Result of a sequence query : statistics over all kmers queries of this sequence
 */
struct result_query {
  int minimum; ///< minimum of all kmer abundances queried
  int maximum; ///< maximum of all kmer abundances queried
  float average; ///< average of all kmer abundances queried
  float kmer_present_ratio; ///< Ratio of present kmers (abundance > 0) over all kmers of queried sequence
};
std::ostream& operator<<(std::ostream& os, result_query const& res);


#endif
