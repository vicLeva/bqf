/*
PRINTING, DEBUGGING AND TESTING
*/
#include <iostream>
#include <bitset>
#include <stdint.h> 
#include <limits>
#include <random>
#include <ctime>
#include <getopt.h>
#include <vector>
#include <algorithm>
#include <chrono>
#include "sys/types.h"
#include "sys/sysinfo.h"

#include "filter.hpp" 

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include <iostream>
#include <fstream>
using namespace std;


#define MEM_UNIT 64ULL
#define MET_UNIT 3ULL
#define OFF_POS 0ULL
#define OCC_POS 1ULL
#define RUN_POS 2ULL

using namespace std;

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int getValue(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

void print_pair(std::pair<uint64_t,uint64_t> bound){
  std::cout << bound.first << "-" << bound.second << " ";
}

void show(uint64_t value, std::string name){
  std::cout << name << std::endl;
  print_bits(value);
}



void test_lots_of_full_cqf(){
  uint64_t seed = time(NULL);
  default_random_engine generator;
  generator.seed(seed);
  uniform_int_distribution<uint64_t> distribution;


  int qsize = 19; 
  Rsqf small_qf(qsize, 64-qsize, false);
  Rsqf usual_qf(4); //524288, qsize=19

  for (size_t j=0 ; j<1000 ; j++) {
    uint64_t seedTMP = distribution(generator);
    cout << "\nseed " << seedTMP << endl;

    
    default_random_engine generatorTMP;
    generatorTMP.seed(seedTMP);
    uniform_int_distribution<uint64_t> distributionTMP;

    cout << "j " << j << endl;
    Rsqf small_qf(qsize, 64-qsize, false);

    for (size_t i=0 ; i<(1ULL<<qsize) ; i++) { //fill to 2^qsize elements (100%-1)
      //cout << "llooooopppp" << endl;
      //cout << i << " ";
      uint64_t val = distributionTMP(generatorTMP);      
      small_qf.insert(val);
    }     
  }
}


void test_lots_of_full_cqf_enumerate() {
  uint64_t seed = time(NULL);
  default_random_engine generator;
  generator.seed(seed);
  uniform_int_distribution<uint64_t> distribution;

  int qsize = 24; 

  for (size_t j=0 ; j<1000 ; j++) {
    
    uint64_t seedTMP = distribution(generator);
    cout << "\nseed " << seedTMP << endl;

    default_random_engine generatorTMP;
    generatorTMP.seed(seedTMP);
    uniform_int_distribution<uint64_t> distributionTMP;

    cout << "j " << j << endl;
    Rsqf small_qf(qsize, 64-qsize, false);
    std::unordered_set<uint64_t> verif;

    for (size_t i=0 ; i<(1ULL<<qsize)-1 ; i++) { //fill to 2^qsize elements (100%-1)
      uint64_t val = distributionTMP(generatorTMP);      
      small_qf.insert(val);
      verif.insert(val);
    }   

    if (verif != small_qf.enumerate()) {
      cout << "error verif != enum" << endl;
      exit(0);
    }
  }  
}



void test_lots_of_full_cqf_remove() {
  uint64_t seed = time(NULL);
  default_random_engine generator;
  generator.seed(seed);
  uniform_int_distribution<uint64_t> distribution;

  int qsize = 19; 
  uint64_t val;
  std::unordered_set<uint64_t> enu;

  for (size_t j=0 ; j<10000 ; j++) {
    
    uint64_t seedTMP = distribution(generator);
    cout << "\nseed " << seedTMP << endl;

    default_random_engine generatorTMP;
    generatorTMP.seed(seedTMP);
    uniform_int_distribution<uint64_t> distributionTMP;

    cout << "j " << j << endl;
    Rsqf small_qf(qsize, 64-qsize, false);
    std::vector<uint64_t> verif;

    //INSERT
    for (size_t i=0 ; i<(1ULL<<qsize)-1 ; i++) { //fill to 2^qsize elements (100%-1) (1ULL<<qsize)-1
      val = distributionTMP(generatorTMP);      
      small_qf.insert(val);
      verif.push_back(val);
    }   

    //REMOVE ELEMS
    for (size_t i=0 ; i<(1ULL<<qsize)-1 ; i++) { //(1ULL<<qsize)/2
      val = verif.back();
      verif.pop_back();
      small_qf.remove(val);
    } 

    //CHECK ENUMERATE
    enu = small_qf.enumerate();
    if (verif.size() != enu.size()) {
      cout << "error verif != enum" << endl;
      exit(0);
    }
    for ( auto it = verif.begin(); it != verif.end(); ++it ){
      if (enu.find(*it) == enu.end()){
        cout << "error verif != enum" << endl;
        exit(0);
      }
    }
  }  
}




void test_one_cqf(){
  uint64_t seed = 9915022754138594861ULL; 
  default_random_engine generator;
  generator.seed(seed);
  uniform_int_distribution<uint64_t> distribution;


  int qsize = 7;
  uint64_t val;
  std::unordered_set<uint64_t> enu;
  Rsqf small_qf(qsize, 64-qsize, true);

  std::vector<uint64_t> verif;

  for (size_t i=0 ; i<127 ; i++) {
    //std::cout << "\ni " << i << endl;
    uint64_t val = distribution(generator);
    /* val &= mask_right(qsize);
    if (val == 0) { val += (1ULL << 45); }
    else { val += (val<<qsize); } */
    
    small_qf.insert(val);
    if (find(verif.begin(), verif.end(), val) == verif.end()) { verif.push_back(val); } 
  }


  //REMOVE ELEMS

  
  //CHECK ENUMERATE
  enu = small_qf.enumerate();
  cout << "done inserting, verif size " << verif.size() << " enum size " << enu.size() << endl;

  std::cout << small_qf.block2string(0) << "\n" << small_qf.block2string(1);

  if (verif.size() != enu.size()) {
    cout << "error verif != enum" << endl;
    exit(0);
  }
  for ( auto it = verif.begin(); it != verif.end(); ++it ){
    if (enu.find(*it) == enu.end()){
      cout << "error verif != enum (diff)" << endl;
      exit(0);
    }
  }

  

   //REMOVE ELEMS
  for (size_t i=0 ; i<127 ; i++) { //(1ULL<<qsize)/2
    val = verif.back();
    verif.pop_back();
    small_qf.remove(val);
  }


  std::cout << small_qf.block2string(0) << "\n" << small_qf.block2string(1);

  //CHECK ENUMERATE
  enu = small_qf.enumerate();
  if (verif.size() != enu.size()) {
    cout << "error verif != enum (post remove, size) verif:" << verif.size() << "  " << enu.size() << endl;
    exit(0);
  }
  for ( auto it = verif.begin(); it != verif.end(); ++it ){
    if (enu.find(*it) == enu.end()){
      cout << "error verif != enum (post remove, diff)" << endl;
      exit(0);
    }
  }

  
}


void test_time_fill_cqf(int q, int n){
  uint64_t seed = time(NULL);
  default_random_engine generator;
  generator.seed(seed);
  uniform_int_distribution<uint64_t> distribution;

  ofstream myfile;
  myfile.open ("/tmp/tmp");
  
  auto ttot = std::chrono::high_resolution_clock::now();

  

  for (int j=0 ; j<n ; j++) {

    uint64_t seedTMP = distribution(generator);
    default_random_engine generatorTMP;
    generatorTMP.seed(seedTMP);
    uniform_int_distribution<uint64_t> distributionTMP;

    Rsqf small_qf(q, 64-q, false);

    for (size_t i=0 ; i<(1ULL<<q)-1 ; i++) { //fill to 2^qsize elements (100%-1)

      auto t1 = std::chrono::high_resolution_clock::now();

      uint64_t val = distributionTMP(generatorTMP);      
      small_qf.insert(val);

      myfile << to_string( std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - t1
      ).count()) << "\n";
    }     
  }

  myfile.close();

  cout << to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms\n";
}




int main(int argc, char** argv) {
    //test_one_cqf();

    //test_lots_of_full_cqf();

    //test_lots_of_full_cqf_enumerate();

    //test_lots_of_full_cqf_remove();

    test_time_fill_cqf(22, 1);



    
    

    
    
  
  return 0;
}