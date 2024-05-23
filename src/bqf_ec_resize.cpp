#include "bqf_ec_resize.hpp"

Bqf_ec_resize::Bqf_ec_resize(uint64_t q_size, uint64_t c_size, uint64_t k, uint64_t z, bool verb): Bqf_ec(q_size, c_size, k, z, verb) {}

void Bqf_ec_resize::resize(int n){
    assert(n > 0);

    std::map<uint64_t, uint64_t> inserted_elements = this->enumerate();

    this->quotient_size += n;
    this->remainder_size -= n;
    
    uint64_t num_quots = 1ULL << this->quotient_size; 
    uint64_t num_of_words = num_quots * (MET_UNIT + this->remainder_size) / MEM_UNIT; 

    this->size_limit = num_quots * 0.95;

    // In machine words
    this->number_blocks = std::ceil(num_quots / BLOCK_SIZE);

    this->filter = std::vector<uint64_t>(num_of_words);

    this->elements_inside = 0;

    for (auto const& elem : inserted_elements){
        this->insert(elem.first, elem.second);
    }

    // std::cout << "testing" << std::endl;
}