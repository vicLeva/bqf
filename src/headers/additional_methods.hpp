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
 * \brief hashes the encoded kmer
 * \param data the kmer to hash
 * \param size the kmer size (in base)
 * \return uint64_t, the hash
 */
uint64_t kmer_to_hash(uint64_t coded_kmer, uint64_t k);

/** 
 * \brief reverse the hash, then decode the number into a kmer (canonical)
 * \param data the hash to reverse
 * \param size the kmer size (in base)
 * \return the kmer (string)
 */
std::string hash_to_kmer(uint64_t hash, uint64_t size);

/**
 * \brief tests if a character is a valid encoding of a nucleotide, i.e. 'A', 'C', 'G' or 'T'
 * \param c the character to encode
 * \returns if it is valid
 */
bool is_valid(char c);

uint64_t nucl_encode(char nucl);

uint64_t flip(uint64_t encoding, size_t bitsize);

const uint8_t rev_table[256] = {
  0xff, 0xbf, 0x7f, 0x3f, 0xef, 0xaf, 0x6f, 0x2f, 0xdf, 0x9f, 0x5f, 0x1f, 0xcf, 0x8f, 0x4f, 0xf, 
  0xfb, 0xbb, 0x7b, 0x3b, 0xeb, 0xab, 0x6b, 0x2b, 0xdb, 0x9b, 0x5b, 0x1b, 0xcb, 0x8b, 0x4b, 0xb, 
  0xf7, 0xb7, 0x77, 0x37, 0xe7, 0xa7, 0x67, 0x27, 0xd7, 0x97, 0x57, 0x17, 0xc7, 0x87, 0x47, 0x7, 
  0xf3, 0xb3, 0x73, 0x33, 0xe3, 0xa3, 0x63, 0x23, 0xd3, 0x93, 0x53, 0x13, 0xc3, 0x83, 0x43, 0x3, 
  0xfe, 0xbe, 0x7e, 0x3e, 0xee, 0xae, 0x6e, 0x2e, 0xde, 0x9e, 0x5e, 0x1e, 0xce, 0x8e, 0x4e, 0xe, 
  0xfa, 0xba, 0x7a, 0x3a, 0xea, 0xaa, 0x6a, 0x2a, 0xda, 0x9a, 0x5a, 0x1a, 0xca, 0x8a, 0x4a, 0xa, 
  0xf6, 0xb6, 0x76, 0x36, 0xe6, 0xa6, 0x66, 0x26, 0xd6, 0x96, 0x56, 0x16, 0xc6, 0x86, 0x46, 0x6, 
  0xf2, 0xb2, 0x72, 0x32, 0xe2, 0xa2, 0x62, 0x22, 0xd2, 0x92, 0x52, 0x12, 0xc2, 0x82, 0x42, 0x2, 
  0xfd, 0xbd, 0x7d, 0x3d, 0xed, 0xad, 0x6d, 0x2d, 0xdd, 0x9d, 0x5d, 0x1d, 0xcd, 0x8d, 0x4d, 0xd, 
  0xf9, 0xb9, 0x79, 0x39, 0xe9, 0xa9, 0x69, 0x29, 0xd9, 0x99, 0x59, 0x19, 0xc9, 0x89, 0x49, 0x9, 
  0xf5, 0xb5, 0x75, 0x35, 0xe5, 0xa5, 0x65, 0x25, 0xd5, 0x95, 0x55, 0x15, 0xc5, 0x85, 0x45, 0x5, 
  0xf1, 0xb1, 0x71, 0x31, 0xe1, 0xa1, 0x61, 0x21, 0xd1, 0x91, 0x51, 0x11, 0xc1, 0x81, 0x41, 0x1, 
  0xfc, 0xbc, 0x7c, 0x3c, 0xec, 0xac, 0x6c, 0x2c, 0xdc, 0x9c, 0x5c, 0x1c, 0xcc, 0x8c, 0x4c, 0xc, 
  0xf8, 0xb8, 0x78, 0x38, 0xe8, 0xa8, 0x68, 0x28, 0xd8, 0x98, 0x58, 0x18, 0xc8, 0x88, 0x48, 0x8, 
  0xf4, 0xb4, 0x74, 0x34, 0xe4, 0xa4, 0x64, 0x24, 0xd4, 0x94, 0x54, 0x14, 0xc4, 0x84, 0x44, 0x4, 
  0xf0, 0xb0, 0x70, 0x30, 0xe0, 0xa0, 0x60, 0x20, 0xd0, 0x90, 0x50, 0x10, 0xc0, 0x80, 0x40, 0x0};

uint64_t revcomp64 (const uint64_t v, size_t bitsize);

uint64_t canonical(uint64_t smer, size_t size);

/** 
 * \brief computes the reverse kmer of a string no matter its encoding
 * 
 * Is rather inefficient so should only be used for debug purposes
 */
std::string revcomp(std::string kmer);

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