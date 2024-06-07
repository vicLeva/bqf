#include "rsqf.hpp"
#include "abstract_bqf.hpp" 
#include "bqf_ec.hpp"

#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>

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

std::string generateRandomKMer(int kmer_size) {
    const char alphabet[] = "ACGT";
    const int alphabetSize = sizeof(alphabet) - 1;
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<> dis(0, alphabetSize - 1);

    std::stringstream randomKmer;
    for (int i = 0; i < kmer_size; ++i) {
        randomKmer << alphabet[dis(gen)];
    }

    return randomKmer.str();
}

struct Perf
{
    uint64_t q_size;
    uint64_t inserted_elements;
    double mock_time;
    double resize_time;
};

Perf time_for_n_insertions(Bqf_ec* mock, Bqf_ec* resize){
    // almost full
    while(resize->elements_inside < resize->size_limit - 1) {
        std::string kmer = generateRandomKMer(resize->kmer_size);
        resize->insert(kmer, 1);
        mock->insert(kmer, 1);
    }

    // calculating times
    auto mockStart = std::chrono::high_resolution_clock::now();
    mock_resize(mock, 1);
    double mockTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - mockStart).count();
    
    auto resizeStart = std::chrono::high_resolution_clock::now();
    resize->resize(1);
    double resizeTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - resizeStart).count();

    return Perf { resize->quotient_size, resize->elements_inside, mockTime, resizeTime};
}

int main() {
    const uint64_t q = 20;
    const uint64_t c = 5;
    const uint64_t k = 32;
    const uint64_t z = 11;
    Bqf_ec mock = Bqf_ec(q, c, k, z, false);
    Bqf_ec resize = Bqf_ec(q, c, k, z, false);
    
    std::cout << "Quotient size, Number of inserted elements, Old time, new time" << std::endl;
    
    while(resize.remainder_size != 0){
        Perf p = time_for_n_insertions(&mock, &resize);

        //std::cout << ' ' << p.q_size;
        //std::cout.flush();

        std::cout << (p.q_size - 1) << ',' << (p.inserted_elements + 1) << ',' << p.mock_time << ',' << p.resize_time << std::endl;
       
    	//checking the bqfs
    	bool same = true;
    	for (uint i = 0; i < mock.filter.size(); ++i){
        	if (mock.filter[i] != resize.filter[i]) {
            	same = false;
            	break;
        	}
    	}
    	std::cout << "same stuff : " << (same? "yes!" : "no...") << std::endl; 
    
    	if(resize.quotient_size == 25) return 0;
    }

    // myfile.close();
    // std::cout << std::endl;

}
