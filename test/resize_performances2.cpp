#include "rsqf.hpp"
#include "abstract_bqf.hpp" 
#include "bqf_ec.hpp"

#include <chrono>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>

int main() {
    const uint64_t q = 8;
    const uint64_t c = 5;
    const uint64_t k = 32;
    const uint64_t z = 11;
    const uint64_t kmer_size = (k - z);
    Bqf_ec resize = Bqf_ec(q, c, k, z, false);
    
    std::random_device rd;
    std::mt19937_64 gen(rd()); // random seed
    std::uniform_int_distribution<uint64_t> dis(0ULL, (1ULL << (2 * kmer_size)) - 1ULL);
    uint64_t hashed_kmer;

    std::cout << "Quotient size, Number of inserted elements, new time" << std::endl;

    while (resize.remainder_size != resize.count_size) {
        // fill up
        while(resize.elements_inside < resize.size_limit - 1) {
            // inserting random hashes
            hashed_kmer = dis(gen);
            resize.insert(hashed_kmer, 1);
        }

        if(resize.quotient_size == 25) break;

        auto resizeStart = std::chrono::high_resolution_clock::now();
        resize.resize(1);
        double resizeTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - resizeStart).count();

        std::cout << (resize.quotient_size - 1) << ',' << (resize.elements_inside + 1) << ',' << resizeTime << std::endl;
    }
    resize.save_on_disk("./new2");
    std::cout << "All Done !" << std::endl;
}