#include "rsqf.hpp"
#include "additional_methods.hpp" 
#include "abstract_bqf.hpp" 
#include "bqf_ec.hpp"
#include "bqf_ec_resize.hpp"

#include <chrono>
#include <regex>
#include <iomanip>
#include <functional>

std::string prettyFilter(Bqf* bqf){
  std::stringstream ss;
  ss << "OFF   : ";
  for (uint block = 0; block < bqf->number_blocks; ++block){
    ss << std::setw(64) << (bqf->filter[block *(MET_UNIT + bqf->remainder_size)]) << " ";
  }
  ss << std::endl;
  for (uint position = 1; position < MET_UNIT; ++position){
    switch (position) {
      case OCC_POS:
        ss << "OCC   : ";
        break;
      case RUN_POS:
        ss << "RUN   : ";
        break;
    }
    for (uint block = 0; block < bqf->number_blocks; ++block){
      std::bitset<64ULL> bytes(bqf->filter[(block *(MET_UNIT + bqf->remainder_size)) + position]);
      auto str = bytes.to_string();
      std::reverse(str.begin(), str.end());
      ss << str << " ";
    }
    ss << std::endl;
  }
  for (uint position = 0; position < bqf->remainder_size; ++position){
    if (position == bqf->remainder_size - bqf->count_size){
        ss << "CNT   : ";
    }else if (position == 0){
        ss << "REM   : ";
    }else{
        ss << "        ";
    }
    for (uint q = 0; q < (1ULL << bqf->quotient_size); ++q){
      if (q % 64 == 0 && q != 0){
        ss << " ";
      }
      ss << char('0' + ((bqf->get_remainder(q, true) & (1 << (bqf->remainder_size - 1 - position))) >> (bqf->remainder_size - 1 - position)));
    }
    ss << std::endl;
  }

  return ss.str();
}

std::pair<std::string, std::string> showDifferences(std::string s1, std::string s2){
  std::stringstream ss1;
  std::stringstream ss2;
  for(std::string::size_type i = 0; i < s1.size(); ++i) {
    if (s1[i] != s2[i]){
      ss1 << "\033[1;32m" << s1[i] << "\033[0m";
      ss2 << "\033[1;31m" << s2[i] << "\033[0m";
    } else if (s1[i] == '1'){
      ss1 << "\033[1;33m" << s1[i] << "\033[0m";
      ss2 << "\033[1;33m" << s2[i] << "\033[0m";
    } else {
      ss1 << s1[i];
      ss2 << s2[i];
    }
  }
  return std::make_pair(ss1.str(), ss2.str());
}

void prettyPrint(Bqf* bqf){
  std::string str = prettyFilter(bqf);
  str = std::regex_replace(str, std::regex("1"), "\033[1;33m1\033[0m");
  std::cout << str;
}

bool compare(Bqf* old, Bqf* revised){
  bool ok = true;
  if (old->quotient_size != revised->quotient_size) {
    ok = false;
    std::cout << "different quotient_size : " << old->quotient_size << " and " << revised->quotient_size << std::endl;
  }
  if (old->remainder_size != revised->remainder_size) {
    ok = false;
    std::cout << "different remainder_size : " << old->remainder_size << " and " << revised->remainder_size << std::endl;
  }
  if (old->size_limit != revised->size_limit) {
    ok = false;
    std::cout << "different size_limit : " << old->size_limit << " and " << revised->size_limit << std::endl;
  }
  if (old->number_blocks != revised->number_blocks) {
    ok = false;
    std::cout << "different number_blocks : " << old->number_blocks << " and " << revised->number_blocks << std::endl;
  }
  if (old->elements_inside != revised->elements_inside) {
    ok = false;
    std::cout << "different elements_inside : " << old->elements_inside << " and " << revised->elements_inside << std::endl;
  }
  if (old->filter.size() != revised->filter.size()) {
    ok = false;
    std::cout << "wrong filter size : " << old->filter.size() << " and " << revised->filter.size() << std::endl;
  }
  
  if(ok){
    std::string strOld = prettyFilter(old);
    std::string strRevised = prettyFilter(revised);
    if (strOld != strRevised){
      std::pair<std::string, std::string> strs = showDifferences(strOld, strRevised);
      std::cout << "filters have different values : " << std::endl;
      std::cout << "expected filter : " << std::endl;
      std::cout << strs.first;
      std::cout << "result filter : " << std::endl;
      std::cout << strs.second;
      ok = false;
    }
  }

  return ok;
}

std::string makeKmer(uint q, uint r, uint shift, uint k){
  return hash_to_kmer(rebuild_number(q, r, shift), k);
}

template <typename F>
void test(bool debug, F* insert, std::string name){
  // Setting the parameters
  uint64_t q_size = 8;
  uint64_t c_size = 2;
  uint64_t k = 8;
  uint64_t z = 2;
  Bqf_ec old = Bqf_ec(q_size, c_size, k, z, false);
  Bqf_ec_resize revised = Bqf_ec_resize(q_size, c_size, k, z, false);

  // Inserting kmers
  (*insert)(&old, &revised, q_size, k, z);

  if (debug){
    std::cout << "Before resize :" << std::endl;
    prettyPrint(&revised);
  }
  
  // Calculating times
  auto oldStart = std::chrono::high_resolution_clock::now();
  old.resize(1);
  double oldTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - oldStart).count();
  
  auto revisedStart = std::chrono::high_resolution_clock::now();
  revised.resize(1);
  double revisedTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - revisedStart).count();

  if (debug){
    std::cout << "After resize :" << std::endl;
    prettyPrint(&revised);
  }

  // Verification
  std::string result = compare(&old, &revised)? "Passed" : "Failed";
  std::cout << name << " " << result << " (old time : " << std::to_string(oldTime) << "ms, revised time : " << revisedTime << "ms)" << std::endl;
}

void testEmpty(bool debug){
  auto insert = [](Bqf* old, Bqf* revised, uint64_t q_size, uint64_t k, uint64_t z) {};
  test(debug, &insert, "Test Empty");
}

void testOneInsert(bool debug){
  auto insert = [](Bqf* old, Bqf* revised, uint64_t q_size, uint64_t k, uint64_t z) {
    revised->insert(makeKmer(0b0, 0b1000, q_size, k - z), 1);
    old->insert(makeKmer(0b0, 0b1000, q_size, k - z), 1);
  };
  test(debug, &insert, "Test One Insert");
}

void testOneRun(bool debug){
  auto insert = [](Bqf* old, Bqf* revised, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q = 0b0;
    revised->insert(makeKmer(q, 0b0000, q_size, k - z), 1);
    revised->insert(makeKmer(q, 0b1000, q_size, k - z), 2);
    revised->insert(makeKmer(q, 0b0100, q_size, k - z), 3);
    revised->insert(makeKmer(q, 0b0010, q_size, k - z), 4);
    revised->insert(makeKmer(q, 0b0001, q_size, k - z), 1);
    revised->insert(makeKmer(q, 0b1001, q_size, k - z), 2);
    revised->insert(makeKmer(q, 0b0101, q_size, k - z), 3);
    revised->insert(makeKmer(q, 0b0011, q_size, k - z), 4);
    old->insert(makeKmer(q, 0b0000, q_size, k - z), 1);
    old->insert(makeKmer(q, 0b1000, q_size, k - z), 2);
    old->insert(makeKmer(q, 0b0100, q_size, k - z), 3);
    old->insert(makeKmer(q, 0b0010, q_size, k - z), 4);
    old->insert(makeKmer(q, 0b0001, q_size, k - z), 1);
    old->insert(makeKmer(q, 0b1001, q_size, k - z), 2);
    old->insert(makeKmer(q, 0b0101, q_size, k - z), 3);
    old->insert(makeKmer(q, 0b0011, q_size, k - z), 4);
  };
  test(debug, &insert, "Test One Run");
}


void testOneRunBackToZero(bool debug){
  auto insert = [](Bqf* old, Bqf* revised, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q = 0b11111110;
    revised->insert(makeKmer(q, 0b0000, q_size, k - z), 1);
    revised->insert(makeKmer(q, 0b1000, q_size, k - z), 1);
    revised->insert(makeKmer(q, 0b0100, q_size, k - z), 1);
    revised->insert(makeKmer(q, 0b0010, q_size, k - z), 1);
    revised->insert(makeKmer(q, 0b0001, q_size, k - z), 1);
    revised->insert(makeKmer(q, 0b1001, q_size, k - z), 1);
    revised->insert(makeKmer(q, 0b0101, q_size, k - z), 1);
    revised->insert(makeKmer(q, 0b0011, q_size, k - z), 1);
    old->insert(makeKmer(q, 0b0000, q_size, k - z), 1);
    old->insert(makeKmer(q, 0b1000, q_size, k - z), 1);
    old->insert(makeKmer(q, 0b0100, q_size, k - z), 1);
    old->insert(makeKmer(q, 0b0010, q_size, k - z), 1);
    old->insert(makeKmer(q, 0b0001, q_size, k - z), 1);
    old->insert(makeKmer(q, 0b1001, q_size, k - z), 1);
    old->insert(makeKmer(q, 0b0101, q_size, k - z), 1);
    old->insert(makeKmer(q, 0b0011, q_size, k - z), 1);
  };
  test(debug, &insert, "Test One Run Back To Zero");
}

void testTwoRunsWithOverlap(bool debug){
  auto insert = [](Bqf* old, Bqf* revised, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q1 = 0b00000000;
    uint64_t q2 = 0b00000010;
    revised->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    revised->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
    old->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    old->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
  };
  test(debug, &insert, "Test Two Run With Overlap");
}

void testManyRunsWithOverlap(bool debug){
  auto insert = [](Bqf* old, Bqf* revised, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q1 = 0b00000000;
    uint64_t q2 = 0b00000010;
    uint64_t q3 = 0b00000100;
    uint64_t q4 = 0b00001000;
    revised->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    revised->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
    revised->insert(makeKmer(q3, 0b0000, q_size, k - z), 3);
    revised->insert(makeKmer(q3, 0b1000, q_size, k - z), 3);
    revised->insert(makeKmer(q3, 0b0100, q_size, k - z), 3);
    revised->insert(makeKmer(q3, 0b0010, q_size, k - z), 3);
    revised->insert(makeKmer(q3, 0b0001, q_size, k - z), 3);
    revised->insert(makeKmer(q3, 0b1001, q_size, k - z), 3);
    revised->insert(makeKmer(q3, 0b0101, q_size, k - z), 3);
    revised->insert(makeKmer(q3, 0b0011, q_size, k - z), 3);
    revised->insert(makeKmer(q4, 0b0000, q_size, k - z), 4);
    revised->insert(makeKmer(q4, 0b1000, q_size, k - z), 4);
    revised->insert(makeKmer(q4, 0b0100, q_size, k - z), 4);
    revised->insert(makeKmer(q4, 0b0010, q_size, k - z), 4);
    revised->insert(makeKmer(q4, 0b0001, q_size, k - z), 4);
    revised->insert(makeKmer(q4, 0b1001, q_size, k - z), 4);
    revised->insert(makeKmer(q4, 0b0101, q_size, k - z), 4);
    revised->insert(makeKmer(q4, 0b0011, q_size, k - z), 4);
    old->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    old->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
    old->insert(makeKmer(q3, 0b0000, q_size, k - z), 3);
    old->insert(makeKmer(q3, 0b1000, q_size, k - z), 3);
    old->insert(makeKmer(q3, 0b0100, q_size, k - z), 3);
    old->insert(makeKmer(q3, 0b0010, q_size, k - z), 3);
    old->insert(makeKmer(q3, 0b0001, q_size, k - z), 3);
    old->insert(makeKmer(q3, 0b1001, q_size, k - z), 3);
    old->insert(makeKmer(q3, 0b0101, q_size, k - z), 3);
    old->insert(makeKmer(q3, 0b0011, q_size, k - z), 3);
    old->insert(makeKmer(q4, 0b0000, q_size, k - z), 4);
    old->insert(makeKmer(q4, 0b1000, q_size, k - z), 4);
    old->insert(makeKmer(q4, 0b0100, q_size, k - z), 4);
    old->insert(makeKmer(q4, 0b0010, q_size, k - z), 4);
    old->insert(makeKmer(q4, 0b0001, q_size, k - z), 4);
    old->insert(makeKmer(q4, 0b1001, q_size, k - z), 4);
    old->insert(makeKmer(q4, 0b0101, q_size, k - z), 4);
    old->insert(makeKmer(q4, 0b0011, q_size, k - z), 4);
  };
  test(debug, &insert, "Test Many Run With Overlap");
}

void testTwoRunsWithOverlapBackToZero1(bool debug){
  auto insert = [](Bqf* old, Bqf* revised, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q1 = 0b11111100;
    uint64_t q2 = 0b11111110;
    revised->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    revised->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
    old->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    old->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
  };
  test(debug, &insert, "Test Two Run With Overlap Back To Zero 1");
}

void testTwoRunsWithOverlapBackToZero2(bool debug){
  auto insert = [](Bqf* old, Bqf* revised, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q1 = 0b11111110;
    uint64_t q2 = 0b00000010;
    revised->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    revised->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    revised->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    revised->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
    old->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    old->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    old->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    old->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
  };
  test(debug, &insert, "Test Two Run With Overlap Back To Zero 2");
}

int main() {
  std::cout << "RESIZING TESTS :" << std::endl;
  bool debug = false;
  testEmpty(debug);
  testOneInsert(debug);
  testOneRun(debug);
  testOneRunBackToZero(debug);
  testTwoRunsWithOverlap(debug);
  testManyRunsWithOverlap(debug);
  testTwoRunsWithOverlapBackToZero1(debug);
  testTwoRunsWithOverlapBackToZero2(debug);
}
