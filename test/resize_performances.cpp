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

    bqf->number_blocks = std::ceil(num_quots / BLOCK_SIZE);

    bqf->filter = std::vector<uint64_t>(num_of_words);

    bqf->elements_inside = 0;

    for (auto const& elem : inserted_elements){
        bqf->insert(elem.first, elem.second);
    }
}

int main() {
    const uint64_t q = 8;
    const uint64_t c = 5;
    const uint64_t k = 32;
    const uint64_t z = 11;
    const uint64_t kmer_size = (k - z);
    Bqf_ec mock = Bqf_ec(q, c, k, z, false);
    Bqf_ec resize = Bqf_ec(q, c, k, z, false);
    
    std::random_device rd;
    std::mt19937_64 gen(rd()); // random seed
    std::uniform_int_distribution<uint64_t> dis(0ULL, (1ULL << (2 * kmer_size)) - 1ULL);
    uint64_t hashed_kmer;

    std::cout << "Quotient size, Number of inserted elements, Old time, new time" << std::endl;

    while (resize.remainder_size != resize.count_size) {
        // fill up
        while(resize.elements_inside < resize.size_limit - 1) {
            // inserting random hashes
            hashed_kmer = dis(gen);
            resize.insert(hashed_kmer, 1);
            mock.insert(hashed_kmer, 1);
        }

        // calculating times
        auto mockStart = std::chrono::high_resolution_clock::now();
        mock_resize(&mock, 1);
        double mockTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - mockStart).count();
        
        auto resizeStart = std::chrono::high_resolution_clock::now();
        resize.resize(1);
        double resizeTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - resizeStart).count();

        std::cout << (resize.quotient_size - 1) << ',' << (resize.elements_inside + 1) << ',' << mockTime << ',' << resizeTime << std::endl;

        if(resize.quotient_size == 30) break;

        // resize.save_on_disk("./bqfs/rdm_resize_q" + std::to_string(resize.quotient_size));
        // mock.save_on_disk("./bqfs/rdm_mock_q" + std::to_string(mock.quotient_size));
    }
    // resize.save_on_disk("./bqfs/rdm_resize_done" + std::to_string(resize.quotient_size));
    // mock.save_on_disk("./bqfs/rdm_mock_done" + std::to_string(mock.quotient_size));
    std::cout << "All Done !" << std::endl;
}