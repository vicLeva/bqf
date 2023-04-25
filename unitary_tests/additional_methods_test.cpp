#include <gtest/gtest.h>
#include "additional_methods.hpp"


TEST(AddMethods, mask_right) {
  EXPECT_EQ (mask_right(0),  0);
  EXPECT_EQ (mask_right(1),  1);
  EXPECT_EQ (mask_right(10),  1023);
  EXPECT_EQ (mask_right(400),  18446744073709551615ULL);
  EXPECT_EQ (mask_right(64),  18446744073709551615ULL);
}

TEST(AddMethods, mask_left) {
  EXPECT_EQ (mask_left(0),  0);
  EXPECT_EQ (mask_left(1),  9223372036854775808ULL);
  EXPECT_EQ (mask_left(10),  18428729675200069632ULL);
  EXPECT_EQ (mask_left(40),  18446744073692774400ULL);
}

TEST(AddMethods, shift_left_right) {
  EXPECT_EQ (shift_left(0, 2),  0);
  EXPECT_EQ (shift_left(1, 10),  1024);
  EXPECT_EQ (shift_left(122563, 64),  0);
  EXPECT_EQ (shift_right(0, 2),  0);
  EXPECT_EQ (shift_right(1ULL << 35, 10),  (1ULL << 35)>>10);
  EXPECT_EQ (shift_right(122563, 64),  0);
}

TEST(AddMethods, get_block_id) {
  EXPECT_EQ (get_block_id(0),  0);
  EXPECT_EQ (get_block_id(1),  0);
  EXPECT_EQ (get_block_id(63),  0);
  EXPECT_EQ (get_block_id(64),  1);
  EXPECT_EQ (get_block_id(63598423),  993725);
}

TEST(AddMethods, get_shift_in_block) {
  EXPECT_EQ (get_shift_in_block(0),  0);
  EXPECT_EQ (get_shift_in_block(1),  1);
  EXPECT_EQ (get_shift_in_block(63),  63);
  EXPECT_EQ (get_shift_in_block(64),  0);
  EXPECT_EQ (get_shift_in_block(63598423),  23);
}

TEST(AddMethods, get_quot_from_block_shift) {
  EXPECT_EQ (get_quot_from_block_shift(0,0),  0);
  EXPECT_EQ (get_quot_from_block_shift(0, 32),  32);
  EXPECT_EQ (get_quot_from_block_shift(1, 10),  74);
  EXPECT_EQ (get_quot_from_block_shift(365, 24),  23384);
}

TEST(AddMethods, bitselectasm) {
  EXPECT_EQ (bitselectasm(0, 1),  64);
  EXPECT_EQ (bitselectasm(0, 5),  64);
  EXPECT_EQ (bitselectasm(1, 1), 0);
  EXPECT_EQ (bitselectasm(1, 2),  64);
  EXPECT_EQ (bitselectasm(2, 1),  1);
  EXPECT_EQ (bitselectasm(2, 2),  64);
}

TEST(AddMethods, bitrankasm) {
  EXPECT_EQ (bitrankasm(20, 0),  0);
  EXPECT_EQ (bitrankasm(20, 2),  1);
  EXPECT_EQ (bitrankasm(20, 4),  2);
  EXPECT_EQ (bitrankasm(20, 60),  2);
  EXPECT_EQ (bitrankasm(365<<10, 14),  3);
  EXPECT_EQ (bitrankasm(365<<10, 20),  6);
}

TEST(AddMethods, get_bits) {
  std::vector<uint64_t> v = {((31<<25) + 365), 1023};
  EXPECT_EQ (get_bits(v, 0, 0),  0);
  EXPECT_EQ (get_bits(v, 25, 30),  31);
  EXPECT_EQ (get_bits(v, 24, 30),  31<<1);
  EXPECT_EQ (get_bits(v, 0, 64),  ((31<<25) + 365));
  EXPECT_EQ (get_bits(v, 50, 64),  16760832);
}

TEST(AddMethods, set_bits) {
  std::vector<uint64_t> v = {((31<<25) + 365), 1023};
  set_bits(v, 11, 15, 9);
  EXPECT_EQ (v[0], ((31<<25) + (15 << 11) + 365));
  EXPECT_EQ (v[1], 1023);

  v = {((31<<25) + 365), 1023};
  set_bits(v, 35, ((1ULL<<63ULL)-1), 64);
  EXPECT_EQ (v[0], 18446744040390001005ULL);
  EXPECT_EQ (v[1], 17179869183ULL);
}





