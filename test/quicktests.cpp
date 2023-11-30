/*
PRINTING, DEBUGGING AND TESTING
*/
#include "rsqf.hpp" 
#include "additional_methods.hpp" 
#include "bqf_ec.hpp" 
#include "bqf_oom.hpp" 

#include <random>
#include <ctime>
#include <getopt.h>
#include <chrono>



#define MEM_UNIT 64ULL
#define MET_UNIT 3ULL
#define OFF_POS 0ULL
#define OCC_POS 1ULL
#define RUN_POS 2ULL

using namespace std;

void print_pair(std::pair<uint64_t,uint64_t> bound){
  std::cout << bound.first << "-" << bound.second << " ";
}

void show(uint64_t value, std::string name){
  std::cout << name << std::endl;
  print_bits(value);
}


std::string generateRandomKMer(int k) {
    static const char alphabet[] = "ACGT";
    static const int alphabetSize = sizeof(alphabet) - 1;
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<> dis(0, alphabetSize - 1);

    std::string randomKMer;
    for (int i = 0; i < k; ++i) {
        randomKMer += alphabet[dis(gen)];
    }

    return randomKMer;
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
    std::cout << "\nseed " << seedTMP << endl;

    
    default_random_engine generatorTMP;
    generatorTMP.seed(seedTMP);
    uniform_int_distribution<uint64_t> distributionTMP;

    std::cout << "j " << j << endl;
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
    std::cout << "\nseed " << seedTMP << endl;

    default_random_engine generatorTMP;
    generatorTMP.seed(seedTMP);
    uniform_int_distribution<uint64_t> distributionTMP;

    std::cout << "j " << j << endl;
    Rsqf small_qf(qsize, 64-qsize, false);
    std::unordered_set<uint64_t> verif;

    for (size_t i=0 ; i<(1ULL<<qsize)-1 ; i++) { //fill to 2^qsize elements (100%-1)
      uint64_t val = distributionTMP(generatorTMP);      
      small_qf.insert(val);
      verif.insert(val);
    }   

    if (verif != small_qf.enumerate()) {
      std::cout << "error verif != enum" << endl;
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
    std::cout << "\nseed " << seedTMP << endl;

    default_random_engine generatorTMP;
    generatorTMP.seed(seedTMP);
    uniform_int_distribution<uint64_t> distributionTMP;

    std::cout << "j " << j << endl;
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
      std::cout << "error verif != enum" << endl;
      exit(0);
    }
    for ( auto it = verif.begin(); it != verif.end(); ++it ){
      if (enu.find(*it) == enu.end()){
        std::cout << "error verif != enum" << endl;
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

  std::map<uint64_t, uint64_t> enu;
  std::map<uint64_t, uint64_t> verif;

  int qsize = 8;
  int hashsize = 56;
  Bqf_ec cqf(qsize, 5, 32, 32-(hashsize/2), true);

  

  //INSERT ELEMS
  for (size_t i=0 ; i<3 ; i++) {
    //std::cout << "\ni " << i << endl;
    uint64_t val = distribution(generator);
    //val &= mask_right(qsize);
    //if (val == 0) { val += (1ULL << 45); }
    //else { val += (val<<qsize); }

    val &= mask_right(hashsize); 
    val &= mask_left(64-qsize);
    val |= 64; 

    print_bits(val);
    std::cout << "inserting " << val << " => " << (val%63) << endl; 
    cqf.insert(val, val%63);
    verif.insert({ val, val%63 });
  }

  

  //REMOVE ELEMS

  
  //CHECK ENUMERATE
  enu = cqf.enumerate();
  std::cout << "done inserting, verif size " << verif.size() << " enum size " << enu.size() << endl;

  std::cout << cqf.block2string(0) << "\n" << cqf.block2string(1) << "\n" << cqf.block2string(2) << "\n" << cqf.block2string(3);

  if (verif.size() != enu.size()) {
    std::cout << "error verif != enum" << endl;
    exit(0);
  }
  for ( auto it = verif.begin(); it != verif.end(); ++it ){
    std::cout << "checking if " << (*it).first << " => " << verif[(*it).first] << " is in enu " << endl;
    if (enu.find((*it).first) == enu.end()){
      std::cout << "error verif != enum (diff)" << endl;
      /* std::cout << "enu atm:" << endl;
      for ( auto ite = enu.begin(); ite != enu.end(); ++ite ){
        std::cout << (*ite).first << " => " << (*ite).second << endl;
      } */
      exit(0);
    }
  }


   //REMOVE ELEMS
  for (std::map<uint64_t,uint64_t>::iterator it = verif.begin(); it != verif.end(); it++){
    std::cout << "removing " << (*it).first << " => " << (*it).second << endl; 
    cqf.remove((*it).first, (*it).second);
  }
  verif.clear();


  //std::cout << cqf.block2string(0) << "\n" << cqf.block2string(1);

  //CHECK ENUMERATE
  enu = cqf.enumerate();
  if (verif.size() != enu.size()) {
    std::cout << "error verif != enum (post remove, size) verif:" << verif.size() << "  " << enu.size() << endl;
    exit(0);
  }
  for ( auto it = verif.begin(); it != verif.end(); ++it ){
    if (enu.find((*it).first) == enu.end()){
      std::cout << "error verif != enum (post remove, diff)" << endl;
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

  std::cout << to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms\n";
}


void experiments(){
  std::cout << "START EXPERIMENT\n";

  std::vector<std::string> random32MerList; //pour requêtes positives 
  std::vector<std::string> randomQueryList; //pour requetes situation réelles
  std::vector<std::string> positive32MerList; //requêtes supposément négatives

  auto ttot = std::chrono::high_resolution_clock::now(); //timer build structure + inserts

  //Bqf_ec cqf(q, c, k, z, debug_print); 
  //on choisit q tel que : 2^(q-1) < #nb_kmers_uniques_insérés < 2^q
  //on choisit r tel que : taille_hash - q (avec taille_hash = 2s = 2*taille_s-mers)
  //on choisit c selon la précision voulue sur les compteurs et l'espace qu'on souhaite économiser

  Bqf_ec cqf(31, 5, 32, 5, false); 
  //Insertion des s-mers comptés avec KMC
  //cqf.insert("/scratch/vlevallois/data/AHX_ACXIOSF_6_1_28_all.txt");
  cqf.insert("/scratch/vlevallois/data/AHX_ACXIOSF_6_1_27_all.txt");
	
  std::cout << to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms (inserts)\n"; 
    


  //Verification (pour l'expé) des sur-estimations / sous-estimations avec les valeurs de 32-mers commptés,
  //qu'on essaie de retrouver via les 27-mers insérés dans le BQF
  /*
  std::ifstream infile("/scratch/vlevallois/data/AHX_ACXIOSF_6_1_32_all.txt");

  std::cout << "start verif" << std::endl;
  std::string a;
  uint64_t b;

  int surestim = 0;
  int totsurestim = 0;
  int bug = 0;
  
  int i = 0;
  uint64_t query;
  int count_limit = 0;

  while (infile >> a >> b) {
      if (i < 100000){
          positive32MerList.push_back(a);
      }
      i++;

      query = cqf.query(a);
      if (query > b){
          surestim ++;
          totsurestim += query - b;
      } 
      if (query < b && b < 31 && query !=-1){
          std::cout << a << " a--b " << b << "  vs query : " << query << std::endl;
          bug ++;
      } 
      if (b > 31){
          count_limit ++;
      } 
  }
  
  std::cout << "nb elems inserted : " << i << std::endl;
  std::cout << "nb count limited : " << count_limit << std::endl;
  std::cout << "nb surestim : " << surestim << std::endl;
  std::cout << "avg surestim : " << totsurestim / (double)surestim << std::endl;
  std::cout << "nb bug : " << bug << std::endl;
  std::cout << "end verif" << std::endl;

  ttot = std::chrono::high_resolution_clock::now();
  for (const auto& mer : positive32MerList) {
      cqf.query(mer);
  }
  std::cout << to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms (100k 32-mers positive query)\n";
  */

  ttot = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < 10000; ++i) {
      std::string random32Mer = generateRandomKMer(32);
      random32MerList.push_back(random32Mer);
  }
  for (const auto& mer : random32MerList) {
	    cqf.query(mer);
  }
  std::cout << to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms (100k 32-mers negative (random) query)\n";

  ttot = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < 10; ++i) {
      std::string randomQuery = generateRandomKMer(100);
      randomQueryList.push_back(randomQuery);
  }
  for (const auto& mer : randomQueryList) {
	    cout << cqf.query(mer) << endl;
  }
  std::cout << to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms (100k 32-mers negative (random) query)\n";
}


int main(int argc, char** argv) {
    //test_one_cqf();

    //test_lots_of_full_cqf();

    //test_lots_of_full_cqf_enumerate();

    ////test_lots_of_full_cqf_remove();

    //test_time_fill_cqf(22, 1);

    experiments();

}
