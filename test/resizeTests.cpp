#include "rsqf.hpp"
#include "additional_methods.hpp" 
#include "abstract_bqf.hpp" 
#include "bqf_ec.hpp"

#include <chrono>
#include <regex>
#include <iomanip>
#include <functional>
#include <random>

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
      ss1 << "\033[;33m" << s1[i] << "\033[0m";
      ss2 << "\033[;33m" << s2[i] << "\033[0m";
    } else {
      ss1 << s1[i];
      ss2 << s2[i];
    }
  }
  return std::make_pair(ss1.str(), ss2.str());
}

void prettyPrint(Bqf* bqf){
  std::string str = prettyFilter(bqf);
  str = std::regex_replace(str, std::regex("1"), "\033[;33m1\033[0m");
  std::cout << str;
}

bool compare(Bqf* mock, Bqf* resize){
  bool ok = true;
  if (mock->quotient_size != resize->quotient_size) {
    ok = false;
    std::cout << "different quotient_size : " << mock->quotient_size << " and " << resize->quotient_size << std::endl;
  }
  if (mock->remainder_size != resize->remainder_size) {
    ok = false;
    std::cout << "different remainder_size : " << mock->remainder_size << " and " << resize->remainder_size << std::endl;
  }
  if (mock->size_limit != resize->size_limit) {
    ok = false;
    std::cout << "different size_limit : " << mock->size_limit << " and " << resize->size_limit << std::endl;
  }
  if (mock->number_blocks != resize->number_blocks) {
    ok = false;
    std::cout << "different number_blocks : " << mock->number_blocks << " and " << resize->number_blocks << std::endl;
  }
  if (mock->elements_inside != resize->elements_inside) {
    ok = false;
    std::cout << "different elements_inside : " << mock->elements_inside << " and " << resize->elements_inside << std::endl;
  }
  if (mock->filter.size() != resize->filter.size()) {
    ok = false;
    std::cout << "wrong filter size : " << mock->filter.size() << " and " << resize->filter.size() << std::endl;
  }
  
  if(ok){
    std::string strmock = prettyFilter(mock);
    std::string strresize = prettyFilter(resize);
    if (strmock != strresize){
      std::pair<std::string, std::string> strs = showDifferences(strmock, strresize);
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
void test(bool printExceptations, F* insert, std::string name, uint64_t q_size, uint64_t c_size, uint64_t k, uint64_t z){
  // Setting the parameters
  Bqf_ec mock = Bqf_ec(q_size, c_size, k, z, false);
  Bqf_ec resize = Bqf_ec(q_size, c_size, k, z, false);

  // Inserting kmers
  (*insert)(&mock, q_size, k, z);
  (*insert)(&resize, q_size, k, z);

  if (printExceptations){
    std::cout << "Before resize :" << std::endl;
    prettyPrint(&resize);
  }
  
  // Calculating times
  auto mockStart = std::chrono::high_resolution_clock::now();
  mock_resize(&mock, 1);
  double mockTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - mockStart).count();
  
  auto resizeStart = std::chrono::high_resolution_clock::now();
  resize.resize(1);
  double resizeTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - resizeStart).count();

  if (printExceptations){
    std::cout << "Expected resize :" << std::endl;
    prettyPrint(&mock);
  }

  // Verification
  std::string result = compare(&mock, &resize)? "\033[1;32mPassed\033[0m" : "\033[1;31mFailed\033[0m";
  std::cout << std::setw(50) << std::left <<name << " : " << result << " (mock time : " << std::to_string(mockTime) << "ms, resize time : " << std::to_string(resizeTime) << "ms x" << std::to_string(mockTime/resizeTime) << ")" << std::endl;
}

void testEmpty(bool printExceptations){
  auto insert = [](Bqf* bqf, uint64_t q_size, uint64_t k, uint64_t z) {};
  test(printExceptations, &insert, "Test Empty", 8, 2, 8, 2);
}

void testOneInsert(bool printExceptations){
  auto insert = [](Bqf* bqf, uint64_t q_size, uint64_t k, uint64_t z) {
    bqf->insert(makeKmer(1, 0b1000, q_size, k - z), 1);
  };
  test(printExceptations, &insert, "Test One Insert 1", 8, 2, 8, 2);
  test(false, &insert, "Test One Insert 2", 9, 2, 6, 1);
  test(false, &insert, "Test One Insert 3", 10, 2, 10, 2);
}

void testOneRun(bool printExceptations){
  auto insert = [](Bqf* bqf, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q = 2;
    bqf->insert(makeKmer(q, 0b0000, q_size, k - z), 1);
    bqf->insert(makeKmer(q, 0b1000, q_size, k - z), 2);
    bqf->insert(makeKmer(q, 0b0100, q_size, k - z), 3);
    bqf->insert(makeKmer(q, 0b0010, q_size, k - z), 4);
    bqf->insert(makeKmer(q, 0b0001, q_size, k - z), 1);
    bqf->insert(makeKmer(q, 0b1001, q_size, k - z), 2);
    bqf->insert(makeKmer(q, 0b0101, q_size, k - z), 3);
    bqf->insert(makeKmer(q, 0b0011, q_size, k - z), 4);
  };
  test(printExceptations, &insert, "Test One Run 1", 8, 2, 8, 2);
  test(false, &insert, "Test One Run 2", 9, 2, 6, 1);
  test(false, &insert, "Test One Run 3", 11, 1, 9, 2);
}


void testOneRunBackToZero(bool printExceptations){
  auto insert = [](Bqf* bqf, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q = (1 << q_size) - 2;
    bqf->insert(makeKmer(q, 0b0000, q_size, k - z), 1);
    bqf->insert(makeKmer(q, 0b1000, q_size, k - z), 1);
    bqf->insert(makeKmer(q, 0b0100, q_size, k - z), 1);
    bqf->insert(makeKmer(q, 0b0010, q_size, k - z), 1);
    bqf->insert(makeKmer(q, 0b0001, q_size, k - z), 1);
    bqf->insert(makeKmer(q, 0b1001, q_size, k - z), 1);
    bqf->insert(makeKmer(q, 0b0101, q_size, k - z), 1);
    bqf->insert(makeKmer(q, 0b0011, q_size, k - z), 1);
  };
  test(printExceptations, &insert, "Test One Run Back To Zero 1", 8, 2, 8, 2);
  test(false, &insert, "Test One Run Back To Zero 2", 9, 2, 6, 1);
  test(false, &insert, "Test One Run Back To Zero 3", 10, 1, 9, 2);
}

void testTwoRunsWithOverlap(bool printExceptations){
  auto insert = [](Bqf* bqf, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q1 = 2;
    uint64_t q2 = 4;
    bqf->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    bqf->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
  };
  test(printExceptations, &insert, "Test Two Run With Overlap 1", 8, 2, 8, 2);
  test(false, &insert, "Test Two Run With Overlap 2", 9, 2, 6, 1);
  test(false, &insert, "Test Two Run With Overlap 3", 10, 2, 10, 1);
}

void testManyRunsWithOverlap(bool printExceptations){
  auto insert = [](Bqf* bqf, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q1 = 2;
    uint64_t q2 = 4;
    uint64_t q3 = 6;
    uint64_t q4 = 8;
    bqf->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    bqf->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
    bqf->insert(makeKmer(q3, 0b0000, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b1000, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b0100, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b0010, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b0001, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b1001, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b0101, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b0011, q_size, k - z), 3);
    bqf->insert(makeKmer(q4, 0b0000, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b1000, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b0100, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b0010, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b0001, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b1001, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b0101, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b0011, q_size, k - z), 4);
  };
  test(printExceptations, &insert, "Test Many Run With Overlap 1", 8, 3, 8, 2);
  test(false, &insert, "Test Many Run With Overlap 2", 9, 2, 6, 1);
  test(false, &insert, "Test Many Run With Overlap 3", 10, 3, 10, 3);
}

void testTwoRunsWithOverlapBackToZero(bool printExceptations){
  auto insert1 = [](Bqf* bqf, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q1 = (1 << q_size) - 4;
    uint64_t q2 = (1 << q_size) - 2;
    bqf->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    bqf->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
  };
  auto insert2 = [](Bqf* bqf, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q1 = (1 << q_size) - 2;
    uint64_t q2 = 2;
    bqf->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    bqf->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
  };
  auto insert3 = [](Bqf* bqf, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q1 = (1 << q_size) - 8;
    uint64_t q2 = (1 << q_size) - 4;
    bqf->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1111, q_size, k - z), 1);
    bqf->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1111, q_size, k - z), 2);
  };
  test(printExceptations, &insert1, "Test Two Run With Overlap Back To Zero 1-1", 8, 2, 8, 2);
  test(false, &insert1, "Test Two Run With Overlap Back To Zero 1-2", 9, 2, 11, 2);
  test(false, &insert1, "Test Two Run With Overlap Back To Zero 1-3", 10, 3, 20, 1);
  test(printExceptations, &insert2, "Test Two Run With Overlap Back To Zero 2-1", 8, 2, 8, 2);
  test(false, &insert2, "Test Two Run With Overlap Back To Zero 2-2", 9, 2, 11, 2);
  test(false, &insert2, "Test Two Run With Overlap Back To Zero 2-3", 10, 3, 20, 1);
  test(printExceptations, &insert3, "Test Two Run With Overlap Back To Zero 3-1", 8, 2, 8, 2);
  test(false, &insert3, "Test Two Run With Overlap Back To Zero 3-2", 9, 2, 11, 2);
  test(false, &insert3, "Test Two Run With Overlap Back To Zero 3-3", 10, 3, 20, 1);
}

void testManyRunsWithOverlapBackToZero(bool printExceptations){
  auto insert = [](Bqf* bqf, uint64_t q_size, uint64_t k, uint64_t z) {
    uint64_t q1 = (1 << q_size) - 4;
    uint64_t q2 = (1 << q_size) - 2;
    uint64_t q3 = 0;
    uint64_t q4 = 2;
    bqf->insert(makeKmer(q1, 0b0000, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1000, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0100, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0010, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0001, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b1001, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0101, q_size, k - z), 1);
    bqf->insert(makeKmer(q1, 0b0011, q_size, k - z), 1);
    bqf->insert(makeKmer(q2, 0b0000, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1000, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0100, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0010, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0001, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b1001, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0101, q_size, k - z), 2);
    bqf->insert(makeKmer(q2, 0b0011, q_size, k - z), 2);
    bqf->insert(makeKmer(q3, 0b0000, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b1000, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b0100, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b0010, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b0001, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b1001, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b0101, q_size, k - z), 3);
    bqf->insert(makeKmer(q3, 0b0011, q_size, k - z), 3);
    bqf->insert(makeKmer(q4, 0b0000, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b1000, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b0100, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b0010, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b0001, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b1001, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b0101, q_size, k - z), 4);
    bqf->insert(makeKmer(q4, 0b0011, q_size, k - z), 4);
  };
  test(printExceptations, &insert, "Test Many Run With Overlap Back To Zero 1", 8, 3, 8, 2);
  test(false, &insert, "Test Many Run With Overlap Back To Zero 2", 9, 2, 6, 1);
  test(false, &insert, "Test Many Run With Overlap Back To Zero 3", 10, 3, 10, 3);
}

std::string generateRandomKMer(int k, std::mt19937* gen) {
  const char alphabet[] = "ACGT";
  const int alphabetSize = sizeof(alphabet) - 1;
  //static std::random_device rd;
  //static std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, alphabetSize - 1);

  std::string randomKMer;
  for (int i = 0; i < k; ++i) {
      randomKMer += alphabet[dis(*gen)];
  }

  return randomKMer;
}

void testRandomInserts(bool printExceptations){
  std::random_device rd;
  const uint_fast32_t seed = rd();

  auto insert = [seed](Bqf* bqf, uint64_t q_size, uint64_t k, uint64_t z) {
    std::mt19937 gen;
    gen.seed(seed);
    while(bqf->elements_inside < bqf->size_limit - 1) {
      bqf->insert(generateRandomKMer(k, &gen), 1);
    }
  };
  test(printExceptations, &insert, "Test Random Inserts 1 seed : " + std::to_string(seed), 8, 2, 8, 2);
  test(false, &insert, "Test Random Inserts 2 seed : " + std::to_string(seed), 9, 2, 6, 1);
  test(false, &insert, "Test Random Inserts 3 seed : " + std::to_string(seed), 10, 3, 7, 2);
}

int main() {
  std::cout << "RESIZING TESTS :" << std::endl;
  const bool printExceptations = false;
  testEmpty(printExceptations);
  testOneInsert(printExceptations);
  testOneRun(printExceptations);
  testOneRunBackToZero(printExceptations);
  testTwoRunsWithOverlap(printExceptations);
  testManyRunsWithOverlap(printExceptations);
  testTwoRunsWithOverlapBackToZero(printExceptations);
  testManyRunsWithOverlapBackToZero(printExceptations);
  testRandomInserts(printExceptations);
}
