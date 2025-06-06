#include "rsqf.hpp"
#include "abstract_bqf.hpp" 
#include "bqf_ec.hpp"

#include <chrono>
#include <iomanip>
#include <functional>
#include <random>
#include <iostream>
#include <set>

std::string red(std::string s, bool bold = false){
  if (bold) {
    return "\033[1;31m" + s + "\033[0m";
  } else {
    return "\033[0;31m" + s + "\033[0m";
  }
}

std::string green(std::string s, bool bold = false){
  if (bold) {
    return "\033[1;32m" + s + "\033[0m";
  } else {
    return "\033[0;32m" + s + "\033[0m";
  }
}

std::string to_bits(uint64_t num, uint size){
  std::bitset<64ULL> b(num);
  return b.to_string().substr(64-size, 64);
}

std::string debug_filter(Bqf* bqf){
  std::stringstream ss;

  uint64_t curr_occ;
  uint64_t quotient;
  std::pair<uint64_t, uint64_t> bounds;
  uint64_t cursor;
  uint64_t remainder;
  uint64_t count;

  for(uint block = 0; block < bqf->number_blocks; ++block){
    curr_occ = bqf->get_occupied_word(block);
    if (curr_occ == 0) continue;

    for (uint64_t i=0; i<BLOCK_SIZE; i++){
      if (curr_occ & 1ULL){
        quotient = block*BLOCK_SIZE + i;
        ss << "Run at q = " << to_bits(bounds.first, bqf->quotient_size);

        bounds = bqf->get_run_boundaries(quotient);
        ss << " [" << to_bits(bounds.first, bqf->quotient_size) << ", " << to_bits(bounds.second, bqf->quotient_size) << "] :" << std::endl << "(";

        cursor = bounds.first;
        while (cursor != (bounds.second)){
          remainder = bqf->get_remainder(cursor, true);
          count = remainder & mask_right(bqf->count_size);
          remainder >>= bqf->count_size;
          ss << "q : " << to_bits(cursor, bqf->quotient_size) << " = " << to_bits(remainder, bqf->remainder_size) << " " << to_bits(count, bqf->count_size) << ", ";
          cursor = bqf->get_next_quot(cursor);
        }
        remainder = bqf->get_remainder(cursor, true);
        count = remainder & mask_right(bqf->count_size);
        remainder >>= bqf->count_size;
        ss << "q : " << to_bits(cursor, bqf->quotient_size) << " = " << to_bits(remainder, bqf->remainder_size) << " " << to_bits(count, bqf->count_size) << ")" << std::endl;
      }
      curr_occ >>= 1ULL;
    }
  }
  return ss.str();
}

std::pair<bool, std::string> debug_filters_and_compare(Bqf* bqf1, Bqf* bqf2){
  bool same = true;
  bool first;
  std::stringstream ss;
  std::set<uint64_t> word_set1;
  std::set<uint64_t> word_set2;

  uint64_t curr_occ1;
  uint64_t curr_occ2;
  uint64_t quotient;
  std::pair<uint64_t, uint64_t> bounds1;
  std::pair<uint64_t, uint64_t> bounds2;
  uint64_t cursor;
  uint64_t end_cursor;
  uint64_t remainder;
  uint64_t count;
  uint64_t word;

  for(uint block = 0; block < bqf1->number_blocks; ++block){
    curr_occ1 = bqf1->get_occupied_word(block);
    curr_occ2 = bqf2->get_occupied_word(block);
    if (curr_occ1 == 0 && curr_occ2 == 0) continue;

    for (uint64_t i=0; i<BLOCK_SIZE; i++){
      if (curr_occ1 & 1ULL && curr_occ2 & 1ULL){
        quotient = block*BLOCK_SIZE + i;
        ss << "Run at q = " << quotient;

        bounds1 = bqf1->get_run_boundaries(quotient);
        bounds2 = bqf2->get_run_boundaries(quotient);
        ss << " expected [" << bounds1.first << ", " << bounds1.second << "], resize [" << bounds2.first << ", " << bounds2.second << "] :" << std::endl << "(expected : ";

        cursor = bounds1.first;
        end_cursor = bqf1->get_next_quot(bounds1.second);
        while (cursor != end_cursor){
          remainder = bqf1->get_remainder(cursor, true);
          count = remainder & mask_right(bqf1->count_size);
          remainder >>= bqf1->count_size;
          word = (((remainder << bqf1->count_size) | count) << bqf1->quotient_size) | cursor;
          word_set1.insert(word);
          cursor = bqf1->get_next_quot(cursor);
        }

        cursor = bounds2.first;
        end_cursor = bqf2->get_next_quot(bounds2.second);
        while (cursor != end_cursor){
          remainder = bqf2->get_remainder(cursor, true);
          count = remainder & mask_right(bqf2->count_size);
          remainder >>= bqf2->count_size;
          word = (((remainder << bqf2->count_size) | count) << bqf2->quotient_size) | cursor;
          word_set2.insert(word);
          cursor = bqf2->get_next_quot(cursor);
        }

        first = true;
        for (auto word : word_set1){
          quotient = word & mask_right(bqf2->quotient_size);
          count = (word >> bqf2->quotient_size) & mask_right(bqf2->count_size);
          remainder = (word >> (bqf2->quotient_size + bqf2->count_size));
          if (!first){
           ss << ", "; 
          }
          first = false;
          if(word_set2.count(word))
            ss << "q : " << to_bits(quotient, bqf1->quotient_size) << " = " << to_bits(remainder, bqf1->remainder_size) << " " << to_bits(count, bqf1->count_size);
          else {
            ss << "q : " << green(to_bits(quotient, bqf1->quotient_size)) << " = " << green(to_bits(remainder, bqf1->remainder_size)) << " " << green(to_bits(count, bqf1->count_size));
            same = false;
          }
        }
        ss << ")" << std::endl << "(resize   : ";
        first = true;
        for (auto word : word_set2){
          quotient = word & mask_right(bqf2->quotient_size);
          count = (word >> bqf2->quotient_size) & mask_right(bqf2->count_size);
          remainder = (word >> (bqf2->quotient_size + bqf2->count_size));
          if (!first){
           ss << ", "; 
          }
          first = false;
          if(word_set1.count(word))
            ss << "q : " << to_bits(quotient, bqf2->quotient_size) << " = " << to_bits(remainder, bqf2->remainder_size) << " " << to_bits(count, bqf2->count_size);
          else {
            ss << "q : " << red(to_bits(quotient, bqf2->quotient_size)) << " = " << red(to_bits(remainder, bqf2->remainder_size)) << " " << red(to_bits(count, bqf2->count_size));
            same = false;
          }
        }
        ss << ")" << std::endl;
        word_set1.clear();
        word_set2.clear();
      } else if (curr_occ1 & 1ULL){
        quotient = block*BLOCK_SIZE + i;
        ss << "Run at q = " << quotient;

        bounds1 = bqf1->get_run_boundaries(quotient);
        ss << " expected [" << bounds1.first << ", " << bounds1.second << "] :" << std::endl << "(";

        cursor = bounds1.first;
        end_cursor = bqf1->get_next_quot(bounds1.second);
        while (cursor != end_cursor){
          remainder = bqf1->get_remainder(cursor, true);
          count = remainder & mask_right(bqf1->count_size);
          remainder >>= bqf1->count_size;
          ss << "q : " << green(to_bits(quotient, bqf1->quotient_size)) << " = " << green(to_bits(remainder, bqf1->remainder_size)) << " " << green(to_bits(count, bqf1->count_size));
          cursor = bqf1->get_next_quot(cursor);
        }
        ss << ")" << std::endl;
        same = false;
      } else if (curr_occ2 & 1ULL){
        quotient = block*BLOCK_SIZE + i;
        ss << "Run at q = " << quotient;

        bounds2 = bqf2->get_run_boundaries(quotient);
        ss << " expected [" << bounds2.first << ", " << bounds2.second << "] :" << std::endl << "(";

        cursor = bounds2.first;
        end_cursor = bqf2->get_next_quot(bounds2.second);
        while (cursor != end_cursor){
          remainder = bqf2->get_remainder(cursor, true);
          count = remainder & mask_right(bqf2->count_size);
          remainder >>= bqf2->count_size;
          ss << "q : " << red(to_bits(quotient, bqf2->quotient_size)) << " = " << red(to_bits(remainder, bqf2->remainder_size)) << " " << red(to_bits(count, bqf2->count_size));
          cursor = bqf2->get_next_quot(cursor);
        }
        ss << ")" << std::endl;
        same = false;
      }
      curr_occ1 >>= 1ULL;
      curr_occ2 >>= 1ULL;
    }
  }

  return std::make_pair(same, ss.str());
}

void mock_resize(Bqf* bqf, int n){
    std::map<uint64_t, uint64_t> inserted_elements = bqf->enumerate();

    bqf->quotient_size += n;
    bqf->remainder_size -= n;
    
    uint64_t num_quots = 1ULL << bqf->quotient_size; 
    uint64_t num_of_words = num_quots * (MET_UNIT + bqf->remainder_size) / MEM_UNIT; 

    bqf->size_limit = num_quots * 0.95;

    // In machine words
    bqf->number_blocks = std::ceil(num_quots / BLOCK_SIZE);

    bqf->filter = std::vector<uint64_t>(num_of_words);

    bqf->elements_inside = 0;

    for (auto const& elem : inserted_elements){
        bqf->insert(elem.first, elem.second);
    }
}

std::pair<bool, std::string> compare(Bqf* mock, Bqf* resize){
  bool result = true;
  std::stringstream ss;

  if (mock->quotient_size != resize->quotient_size) {
    result = false;
    ss << red("Error") << " : different quotient_size : expected " << green(std::to_string(mock->quotient_size)) << " but received " << red(std::to_string(resize->quotient_size)) << std::endl;
  }
  if (mock->remainder_size != resize->remainder_size) {
    result = false;
    ss << red("Error") << " : different remainder_size : expected " << green(std::to_string(mock->remainder_size)) << " but received " << red(std::to_string(resize->remainder_size)) << std::endl;
  }
  if (mock->size_limit != resize->size_limit) {
    result = false;
    ss << red("Error") << " : different size_limit : expected " << green(std::to_string(mock->size_limit)) << " but received " << red(std::to_string(resize->size_limit)) << std::endl;
  }
  if (mock->number_blocks != resize->number_blocks) {
    result = false;
    ss << red("Error") << " : different number_blocks : expected " << green(std::to_string(mock->number_blocks)) << " but received " << red(std::to_string(resize->number_blocks)) << std::endl;
  }
  if (mock->elements_inside != resize->elements_inside) {
    result = false;
    ss << red("Error") << " : different elements_inside : expected " << green(std::to_string(mock->elements_inside)) << " but received " << red(std::to_string(resize->elements_inside)) << std::endl;
  }
  if (mock->filter.size() != resize->filter.size()) {
    result = false;
    ss << red("Error") << " : wrong filter size : expected " << green(std::to_string(mock->filter.size())) << " but received " << red(std::to_string(resize->filter.size())) << std::endl;
  }
  
  if(result){
    std::pair<bool, std::string> compared = debug_filters_and_compare(mock, resize);
    if(!compared.first){
      result = false;
      ss << compared.second;
    }
  }

  return std::make_pair(result, ss.str());
}

// rajouter n en param et faire diff implementation des insterts

template <typename F>
void test(bool printExceptations, F* insert, std::string test_name, const uint64_t q_size, uint64_t c_size, uint64_t k, uint64_t z, uint n, uint* success, uint* total, double* avg){
  Bqf_ec mock = Bqf_ec(q_size, c_size, k, z, false);
  Bqf_ec resize = Bqf_ec(q_size, c_size, k, z, false);

  (*insert)(&mock);
  (*insert)(&resize);

  if (printExceptations)
    std::cout << test_name << " before resize :" << std::endl << debug_filter(&mock) << std::endl;

  auto mockStart = std::chrono::high_resolution_clock::now();
  mock_resize(&mock, 1);
  double mockTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - mockStart).count();
  
  auto resizeStart = std::chrono::high_resolution_clock::now();
  resize.resize(1);
  double resizeTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - resizeStart).count();

  if (printExceptations)
    std::cout << test_name << " expected resize :" << std::endl << debug_filter(&mock) << std::endl;

  std::pair<bool, std::string> compared = compare(&mock, &resize);
  if (!compared.first) {
    std::cout << "Params pre resize : q_size = " << q_size << ", c_size = " << c_size << ", k = " << k << ", z = " << z << std::endl;
    std::cout << "Missing elements in "<< green("green") << " and elements that shouldn't exists in " << red("red") << std::endl;
    std::cout << compared.second << std::endl;
  } else {
    (*success)++;
  }

  std::string result = compared.first? green("PASSED", true) : red("FAILED", true);
  std::cout << std::setw(50) << std::left << test_name << " : [" << result << "] (mock time : " << std::to_string(mockTime) << "ms, resize time : " << std::to_string(resizeTime) << "ms    x" << std::to_string(mockTime/resizeTime) << ")" << std::endl;
  
  (*total)++;
  (*avg)+=mockTime/resizeTime;
}

std::string make_kmer(uint q, uint r, uint shift, uint k){
  return hash_to_kmer(rebuild_number(q, r, shift), k);
}

template <typename F>
void test_multiple(bool printExceptations, F* insert, std::string test_name, uint* success, uint* total, double* avg){
  uint i = 0;
  test(printExceptations, insert, test_name + std::to_string(++i), 8, 2, 8, 2, 1, success, total, avg);
  test(printExceptations, insert, test_name + std::to_string(++i), 10, 2, 12, 4, 1, success, total, avg);
  test(printExceptations, insert, test_name + std::to_string(++i), 12, 2, 17, 5, 1, success, total, avg);
  test(printExceptations, insert, test_name + std::to_string(++i), 8, 2, 8, 2, 2, success, total, avg);
  test(printExceptations, insert, test_name + std::to_string(++i), 8, 2, 8, 2, 3, success, total, avg);
  test(printExceptations, insert, test_name + std::to_string(++i), 8, 2, 8, 2, 0, success, total, avg);
}

void testEmpty(bool printExceptations, uint* success, uint* total, double* avg){
  auto insert = [](Bqf* bqf) {};
  test_multiple(printExceptations, &insert, "Test Empty ", success, total, avg);
}

void testOneInsert(bool printExceptations, uint* success, uint* total, double* avg){
  auto insert1 = [](Bqf* bqf) {
    bqf->insert(make_kmer(1, 0b1000, bqf->quotient_size, bqf->smer_size), 1);
  };
  auto insert2 = [](Bqf* bqf) {
    bqf->insert(make_kmer(1, 0b1001, bqf->quotient_size, bqf->smer_size), 1);
  };
  test_multiple(printExceptations, &insert1, "Test One Insert 1-", success, total, avg);
  test_multiple(printExceptations, &insert2, "Test One Insert 2-", success, total, avg);
}

void testOneRun(bool printExceptations, uint* success, uint* total, double* avg){
  auto insert1 = [](Bqf* bqf) {
    const uint64_t q = 2;
    bqf->insert(make_kmer(q, 0b0000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q, 0b1000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q, 0b0100, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q, 0b0010, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q, 0b0001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q, 0b1001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q, 0b0101, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q, 0b1110, bqf->quotient_size, bqf->smer_size), 4);
  };
  auto insert2 = [](Bqf* bqf) {
    const uint64_t q = 62;
    bqf->insert(make_kmer(q, 0b0000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q, 0b1000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q, 0b0100, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q, 0b0010, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q, 0b0001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q, 0b1001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q, 0b0101, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q, 0b1110, bqf->quotient_size, bqf->smer_size), 4);
  };
  test_multiple(printExceptations, &insert1, "Test One Run 1-", success, total, avg);
  test_multiple(printExceptations, &insert2, "Test One Run 2-", success, total, avg);
}


void testOneRunBackToZero(bool printExceptations, uint* success, uint* total, double* avg){
  auto insert = [](Bqf* bqf) {
    const uint64_t q = (1 << bqf->quotient_size) - 2;
    bqf->insert(make_kmer(q, 0b0000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q, 0b1000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q, 0b0100, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q, 0b0010, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q, 0b1001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q, 0b0101, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q, 0b0011, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q, 0b1111, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q, 0b1011, bqf->quotient_size, bqf->smer_size), 1);
  };
  test_multiple(printExceptations, &insert, "Test One Run Back To Zero ", success, total, avg);
}

void testSimpleCluster(bool printExceptations, uint* success, uint* total, double* avg){
  auto insert1 = [](Bqf* bqf) {
    const uint64_t q1 = 2;
    const uint64_t q2 = 4;
    bqf->insert(make_kmer(q1, 0b0000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1110, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0100, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0010, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0101, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0011, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q2, 0b0000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0100, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0010, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0101, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0011, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1101, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1111, bqf->quotient_size, bqf->smer_size), 2);
  };
  auto insert2 = [](Bqf* bqf) {
    const uint64_t q1 = 63;
    const uint64_t q2 = 65;
    bqf->insert(make_kmer(q1, 0b0000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1110, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0100, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0010, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0101, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0011, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q2, 0b0000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0100, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0010, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0101, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0011, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1101, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1111, bqf->quotient_size, bqf->smer_size), 2);
  };
  test_multiple(printExceptations, &insert1, "Test Simple Cluster 1-", success, total, avg);
  test_multiple(printExceptations, &insert2, "Test Simple Cluster 2-", success, total, avg);
}

void testComplexCluster(bool printExceptations, uint* success, uint* total, double* avg){
  auto insert1 = [](Bqf* bqf) {
    const uint64_t q1 = 2;
    const uint64_t q2 = 4;
    const uint64_t q3 = 6;
    const uint64_t q4 = 8;
    bqf->insert(make_kmer(q1, 0b0000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0100, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0010, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0101, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0011, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1111, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q2, 0b0000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0100, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0010, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0101, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q3, 0b0000, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b1000, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0100, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0010, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0001, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b1001, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0101, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q4, 0b0000, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b1000, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0100, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0010, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0001, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b1001, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0101, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0011, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b1111, bqf->quotient_size, bqf->smer_size), 4);
  };
  auto insert2 = [](Bqf* bqf) {
    const uint64_t q1 = 61;
    const uint64_t q2 = 63;
    const uint64_t q3 = 65;
    const uint64_t q4 = 67;
    bqf->insert(make_kmer(q1, 0b0000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0100, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0010, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0101, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0011, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1111, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q2, 0b0000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0100, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0010, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0101, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q3, 0b0000, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b1000, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0100, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0010, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0001, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b1001, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0101, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q4, 0b0000, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b1000, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0100, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0010, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0001, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b1001, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0101, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0011, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b1111, bqf->quotient_size, bqf->smer_size), 4);
  };
  test_multiple(printExceptations, &insert1, "Test Complex Cluster 1-", success, total, avg);
  test_multiple(printExceptations, &insert2, "Test Complex Cluster 2-", success, total, avg);
}

void testSimpleClusterBackToZero(bool printExceptations, uint* success, uint* total, double* avg){
  auto insert1 = [](Bqf* bqf) {
    const uint64_t q1 = (1 << bqf->quotient_size) - 4;
    const uint64_t q2 = (1 << bqf->quotient_size) - 2;
    bqf->insert(make_kmer(q1, 0b0000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0100, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0010, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0101, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0011, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q2, 0b0000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0100, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0010, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0101, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0011, bqf->quotient_size, bqf->smer_size), 2);
  };
  auto insert2 = [](Bqf* bqf) {
    const uint64_t q1 = (1 << bqf->quotient_size) - 2;
    const uint64_t q2 = 2;
    bqf->insert(make_kmer(q1, 0b0000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0100, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0010, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0101, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0011, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q2, 0b0000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0100, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0010, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0101, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0011, bqf->quotient_size, bqf->smer_size), 2);
  };
  auto insert3 = [](Bqf* bqf) {
    const uint64_t q1 = (1 << bqf->quotient_size) - 8;
    const uint64_t q2 = (1 << bqf->quotient_size) - 4;
    bqf->insert(make_kmer(q1, 0b0000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0100, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0010, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0101, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0011, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1111, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q2, 0b0000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0100, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0010, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0101, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0011, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1111, bqf->quotient_size, bqf->smer_size), 2);
  };
  auto insert4 = [](Bqf* bqf) {
    const uint64_t q1 = (1 << bqf->quotient_size) - 4;
    const uint64_t q2 = 0;
    bqf->insert(make_kmer(q1, 0b0000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0100, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0010, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0101, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0011, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q2, 0b0000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0100, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0010, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0101, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0011, bqf->quotient_size, bqf->smer_size), 2);
  };

  test_multiple(printExceptations, &insert1, "Test Simple Cluster Back To Zero 1-", success, total, avg);
  test_multiple(printExceptations, &insert2, "Test Simple Cluster Back To Zero 2-", success, total, avg);
  test_multiple(printExceptations, &insert3, "Test Simple Cluster Back To Zero 3-", success, total, avg);
  test_multiple(printExceptations, &insert4, "Test Simple Cluster Back To Zero 4-", success, total, avg);
}

void testComplexClusterBackToZero(bool printExceptations, uint* success, uint* total, double* avg){
  auto insert = [](Bqf* bqf) {
    const uint64_t q1 = (1 << bqf->quotient_size) - 4;
    const uint64_t q2 = (1 << bqf->quotient_size) - 2;
    const uint64_t q3 = 0;
    const uint64_t q4 = 2;
    bqf->insert(make_kmer(q1, 0b0000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1000, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0100, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0010, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b1001, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0101, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q1, 0b0011, bqf->quotient_size, bqf->smer_size), 1);
    bqf->insert(make_kmer(q2, 0b0000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1000, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0100, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0010, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b1001, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0101, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q2, 0b0011, bqf->quotient_size, bqf->smer_size), 2);
    bqf->insert(make_kmer(q3, 0b0000, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b1000, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0100, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0010, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0001, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b1001, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0101, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q3, 0b0011, bqf->quotient_size, bqf->smer_size), 3);
    bqf->insert(make_kmer(q4, 0b0000, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b1000, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0100, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0010, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0001, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b1001, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0101, bqf->quotient_size, bqf->smer_size), 4);
    bqf->insert(make_kmer(q4, 0b0011, bqf->quotient_size, bqf->smer_size), 4);
  };
  test_multiple(printExceptations, &insert, "Test Complex Cluster Back To Zero ", success, total, avg);
}

std::string generateRandomKMer(int k, std::mt19937* gen) {
  const char alphabet[] = "ACGT";
  const int alphabetSize = sizeof(alphabet) - 1;
  std::uniform_int_distribution<> dis(0, alphabetSize - 1);

  std::string randomKMer;
  for (int i = 0; i < k; ++i) {
      randomKMer += alphabet[dis(*gen)];
  }

  return randomKMer;
}

void testRandomInserts(bool printExceptations, uint* success, uint* total, double* avg, uint n){
  std::random_device rd;
  uint_fast32_t seed = 0;
  const uint64_t q_sizes[] = {10, 12, 14};
  const uint64_t c_sizes[] = {2, 4, 6};
  const uint64_t ks[] = {20, 25, 31};
  const uint64_t zs[] = {2, 6, 11};
  const uint64_t ns[] = {1, 2, 3};

  auto insert = [seed](Bqf* bqf) {
    std::mt19937 gen;
    gen.seed(seed);
    while(bqf->elements_inside < bqf->size_limit - 1) {
      bqf->insert(generateRandomKMer(bqf->kmer_size, &gen), 1);
    }
  };
  
  for (uint i = 0 ; i < n ; ++i){
    seed = rd();
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> dis(0, 2);
    const uint64_t q_size = q_sizes[dis(gen)];
    const uint64_t c_size = c_sizes[dis(gen)];
    const uint64_t k = ks[dis(gen)];
    const uint64_t z = zs[dis(gen)];
    const uint64_t n = ns[dis(gen)];
    test(printExceptations, &insert, "Test Random Inserts " + std::to_string(i) + " seed : " + std::to_string(seed), q_size, c_size, k, z, n, success, total, avg);
  }
}

int main() {
  std::cout << std::endl << "\033[1;37mRESIZING TESTS\033[0m : " << std::endl;
  
  const bool printExceptations = false;
  uint success = 0;
  uint total = 0;
  double avg = 0;
  testEmpty(printExceptations, &success, &total, &avg);
  testOneInsert(printExceptations, &success, &total, &avg);
  testOneRun(printExceptations, &success, &total, &avg);
  testOneRunBackToZero(printExceptations, &success, &total, &avg);
  testSimpleCluster(printExceptations, &success, &total, &avg);
  testComplexCluster(printExceptations, &success, &total, &avg);
  testSimpleClusterBackToZero(printExceptations, &success, &total, &avg);
  testComplexClusterBackToZero(printExceptations, &success, &total, &avg);

  std::cout << std::endl << "\033[1;37mRESIZING RANDOM TESTS\033[0m : " << std::endl;

  uint success_rdm = 0;
  uint total_rdm = 0;
  double avg_rdm = 0;
  testRandomInserts(false, &success_rdm, &total_rdm, &avg_rdm, 50);

  std::cout << std::endl << "\033[1;37mResults\033[0m : " << std::endl;

  std::cout << "Tests :        " << success << "/" << total << " [" << green("PASSED", true) << "]" << " (avg : x" << std::to_string(avg/total) << ")" << std::endl;
  std::cout << "Random tests : " << success_rdm << "/" << total_rdm << " [" << green("PASSED", true) << "]" << " (avg : x" << std::to_string(avg_rdm/total_rdm) << ")" << std::endl;
}
