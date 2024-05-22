#include "rsqf.hpp"
#include "additional_methods.hpp" 
#include "abstract_bqf.hpp" 
#include "bqf_ec.hpp"
#include "bqf_ec_resize.hpp"

#include <chrono>


#define MEM_UNIT 64ULL
#define MET_UNIT 3ULL
#define OFF_POS 0ULL
#define OCC_POS 1ULL
#define RUN_POS 2ULL

void printFilter(Bqf* bqf){
  std::cout << "{";
  for (uint64_t i : bqf->filter){
    std::cout << i << ", ";
  }
  std::cout << "}" << std::endl;
}

void prettyPrintFilter(Bqf* bqf){
  for (uint position = 0; position < MET_UNIT; ++position){
    switch (position) {
      case OFF_POS:
        std::cout << "OFF   : ";
        break;
      case OCC_POS:
        std::cout << "OCC   : ";
        break;
      case RUN_POS:
        std::cout << "RUN   : ";
        break;
    }
    for (uint block = 0; block < bqf->number_blocks; ++block){
      std::bitset<64ULL> bytes(bqf->filter[(block *(MET_UNIT + bqf->remainder_size)) + position]);
      std::cout << bytes << " ";
    }
    std::cout << std::endl;
  }
  for (uint position = 0; position < bqf->remainder_size; ++position){
    std::cout << "        ";
    for (uint q = 0; q < (1ULL << bqf->quotient_size); ++q){
      if (q % 64 == 0 && q != 0){
        std::cout << " ";
      }
      std::cout << ((bqf->get_remainder(q, true) & (1 << position)) == 0? "0" : "1");
    }
    std::cout << std::endl;
  }
}

bool compare(Bqf* bqf1, Bqf* bqf2){
  bool ok = true;
  if (bqf1->quotient_size != bqf2->quotient_size) {
    ok = false;
    std::cout << "different quotient_size : " << bqf1->quotient_size << " and " << bqf2->quotient_size << std::endl;
  }
  if (bqf1->remainder_size != bqf2->remainder_size) {
    ok = false;
    std::cout << "different remainder_size : " << bqf1->remainder_size << " and " << bqf2->remainder_size << std::endl;
  }
  if (bqf1->size_limit != bqf2->size_limit) {
    ok = false;
    std::cout << "different size_limit : " << bqf1->size_limit << " and " << bqf2->size_limit << std::endl;
  }
  if (bqf1->number_blocks != bqf2->number_blocks) {
    ok = false;
    std::cout << "different number_blocks : " << bqf1->number_blocks << " and " << bqf2->number_blocks << std::endl;
  }
  if (bqf1->elements_inside != bqf2->elements_inside) {
    ok = false;
    std::cout << "different elements_inside : " << bqf1->elements_inside << " and " << bqf2->elements_inside << std::endl;
  }
  if (bqf1->filter.size() != bqf2->filter.size()) {
    ok = false;
    std::cout << "wrong filter size : " << bqf1->filter.size() << " and " << bqf2->filter.size() << std::endl;
  }
  
  if(ok){
    for (uint i = 0; i < bqf1->filter.size(); ++i){
      if (bqf1->filter[i] != bqf1->filter[i]){
        ok = false;
        std::cout << "wrong bit in filter : " << bqf1->filter[i] << " and " << bqf2->filter[i] << " at " << i << std::endl;
        prettyPrintFilter(bqf1);
        prettyPrintFilter(bqf2);
      }
    }
  }

  return ok;
}

void testEmpty(){
  // Setting the parameters
  uint64_t q_size = 8;
  uint64_t c_size = 2;
  uint64_t k = 6;
  uint64_t z = 0;
  Bqf_ec original = Bqf_ec(q_size, c_size, k, z, false);
  Bqf_ec_resize revised = Bqf_ec_resize(q_size, c_size, k, z, true);
  revised.insert("ACTGAG", 1);
  prettyPrintFilter(&revised);
  printFilter(&revised);
  prettyPrintFilter(&original);
  printFilter(&original);
  
  // Calculating times
  auto originalStart = std::chrono::high_resolution_clock::now();
  original.resize(1);
  double originalTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - originalStart).count();
  
  auto revisedStart = std::chrono::high_resolution_clock::now();
  revised.resize(1);
  double revisedTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - revisedStart).count();

  // prettyPrintFilter(&revised);
  // prettyPrintFilter(&original);

  // Verification
  std::string result = compare(&original, &revised)? "Passed" : "Failed";
  std::cout << "Test Empty " << result << " (original time : " << std::to_string(originalTime) << "ms, revised time : " << revisedTime << "ms)" << std::endl;
}

int main() {
  std::cout << "RESIZING TESTS :" << std::endl;
  testEmpty();
}
